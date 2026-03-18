#!/usr/bin/env python3
"""
ks9.py  -  PacBio subread BAM pre-conditioner for xz compression
=================================================================
Usage (compress):
    samtools view -h subreads.bam | python3 ks9.py -c | xz -9e > out.sq9.xz

Usage (decompress):
    xz -d < out.sq9.xz | python3 ks9.py -u | samtools view -b - > restored.bam

Changes from ks8 (derived from analysis of real Revio/Sequel II SAM data)
--------------------------------------------------------------------------
1. QUAL ELIMINATION
   PacBio subread QUAL is universally all '!' (Phred 0) -- a placeholder.
   Stored as a 1-byte mode flag instead of seq_len identical bytes:
     0x00 = all '!'     -> restored as '!' * seq_len
     0x01 = constant V  -> 1 extra byte; restored as chr(V) * seq_len
     0x02 = literal     -> seq_len bytes of real per-base quality
   Frees ~30% of the xz input that was previously all-identical bytes,
   giving the LZMA dictionary more room for ip/pw cross-read matching.

2. REDUNDANT TAG ELIMINATION
   zm:i / qs:i / qe:i are byte-for-byte redundant with QNAME field
   (movie/zmw/qs_qe). Dropped on compression, reconstructed on decompress.
   qt:Z: is always all '!' (barcode quality placeholder). Dropped; length
   reconstructed from the paired bt:Z: tag on decompress.
   Tag insertion order is preserved via a compact per-read order record.

3. BLOCK-STRIPED OUTPUT  (O(block) RAM, no temp files)
   Accumulate BLOCK_SIZE reads into per-channel buffers, flush channel-
   by-channel so xz sees homogeneous byte streams:
     [all meta] [all dna] [all qual-flags] [all tags]
   RAM bounded at ~300 MB for default BLOCK_SIZE=5000 on Revio reads.
"""

import sys
import struct

MAGIC      = b'KS9D'
BLOCK_SIZE = 5000

K4_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
K4_REV = ['A', 'C', 'G', 'T']

QUAL_ALL_BANG = 0x00
QUAL_CONSTANT = 0x01
QUAL_LITERAL  = 0x02

# These are dropped during compression and rebuilt from context on decompress
_DROP_TAGS = {'zm', 'qs', 'qe'}   # redundant with QNAME
# qt:Z: is dropped only when all '!' -- handled inline


# ── DNA ───────────────────────────────────────────────────────────────────────

def pack_dna_k4(seq):
    seq = seq.upper()
    packed = bytearray()
    exceptions = []
    for i in range(0, len(seq), 4):
        byte_val = 0
        block = seq[i:i+4]
        for j, base in enumerate(block):
            if base in K4_MAP:
                val = K4_MAP[base]
            else:
                val = 0
                exceptions.append((i + j, base))
            byte_val |= (val << (6 - j * 2))
        packed.append(byte_val)
    exc_bin = struct.pack(">H", len(exceptions))
    for pos, char in exceptions:
        exc_bin += struct.pack(">IB", pos, ord(char))
    return bytes(packed), exc_bin


def unpack_dna_k4(packed, length, exc_bin):
    seq = []
    for i in range(length):
        byte_idx = i // 4
        sub_idx  = i % 4
        if byte_idx >= len(packed):
            break
        val = (packed[byte_idx] >> (6 - sub_idx * 2)) & 0x03
        seq.append(K4_REV[val])
    if len(exc_bin) >= 2:
        exc_count = struct.unpack(">H", exc_bin[:2])[0]
        offset = 2
        for _ in range(exc_count):
            if offset + 5 > len(exc_bin):
                break
            pos, char_code = struct.unpack(">IB", exc_bin[offset:offset+5])
            if pos < len(seq):
                seq[pos] = chr(char_code)
            offset += 5
    return "".join(seq)


# ── QUAL ──────────────────────────────────────────────────────────────────────

def encode_qual(qual_bytes):
    unique = set(qual_bytes)
    if len(unique) == 1:
        v = next(iter(unique))
        if v == 33:
            return bytes([QUAL_ALL_BANG])
        return bytes([QUAL_CONSTANT, v])
    return bytes([QUAL_LITERAL]) + qual_bytes


def decode_qual(qual_enc, seq_len):
    mode = qual_enc[0]
    if mode == QUAL_ALL_BANG:
        return '!' * seq_len
    if mode == QUAL_CONSTANT:
        return chr(qual_enc[1]) * seq_len
    return qual_enc[1:seq_len + 1].decode('latin-1')


# ── TAGS ──────────────────────────────────────────────────────────────────────

def encode_tags(extra_cols):
    """
    Returns (tag_payload_bytes, order_bytes).
    order_bytes: 2 bytes per original tag, encoding what was stored and how.
      High nibble: 0=stored-normal, 1=stored-ip/pw, 2=dropped-redundant, 3=dropped-qt
      Low nibble + next byte: the 2-char tag name as 2 ASCII bytes
    Packed as: [n_tags: 1 byte] [per-tag: 1 flag byte + 2 name bytes]
    """
    payload = bytearray()
    order   = bytearray()

    n_tags = len(extra_cols)
    order.append(n_tags)

    qt_dropped = False

    for t in extra_cols:
        name = t[:2]

        if name in _DROP_TAGS:
            order.append(0x20); order.extend(name.encode())
            continue

        if t.startswith('qt:Z:') and all(c == '!' for c in t[5:]):
            order.append(0x30); order.extend(name.encode())
            qt_dropped = True
            continue

        if t.startswith(("ip:B:C,", "pw:B:C,")):
            vals = bytes(int(x) for x in t[7:].split(','))
            payload += struct.pack(">B",  2)
            payload += name.encode()
            payload += struct.pack(">I",  len(vals))
            payload += vals
            order.append(0x10); order.extend(name.encode())
        else:
            t_b = t.encode('utf-8')
            payload += struct.pack(">BH", 3, len(t_b))
            payload += t_b
            order.append(0x00); order.extend(name.encode())

    return bytes(payload), bytes(order)


def decode_tags(tp, order_bytes, qname):
    """Reconstruct tags in original order."""
    # Parse the order record
    n_tags = order_bytes[0]
    order  = []
    for i in range(n_tags):
        base = 1 + i * 3
        flag = order_bytes[base]
        name = order_bytes[base+1:base+3].decode()
        order.append((flag, name))

    # Decode the payload into a dict by name
    stored = {}   # name -> tag_string (for ip/pw) or tag_string (for normal)
    off = 0
    while off < len(tp):
        kind = tp[off]; off += 1
        if kind == 2:
            name  = tp[off:off+2].decode(); off += 2
            v_len = struct.unpack(">I", tp[off:off+4])[0]; off += 4
            vals  = tp[off:off+v_len]; off += v_len
            stored[name] = f"{name}:B:C," + ",".join(map(str, list(vals)))
        elif kind == 3:
            g_len = struct.unpack(">H", tp[off:off+2])[0]; off += 2
            tag_str = tp[off:off+g_len].decode('utf-8'); off += g_len
            stored[tag_str[:2]] = tag_str

    # Reconstruct zm/qs/qe from QNAME
    parts = qname.split('/')
    reconstructed = {}
    if len(parts) == 3:
        reconstructed['zm'] = f"zm:i:{parts[1]}"
        coords = parts[2].split('_')
        if len(coords) == 2:
            reconstructed['qs'] = f"qs:i:{coords[0]}"
            reconstructed['qe'] = f"qe:i:{coords[1]}"

    # Find bt tag to reconstruct qt length
    bt_tag = stored.get('bt')
    if bt_tag and bt_tag.startswith('bt:Z:'):
        reconstructed['qt'] = f"qt:Z:{'!' * len(bt_tag[5:])}"

    # Emit in original order
    tags = []
    for flag, name in order:
        if flag == 0x20:   # dropped redundant
            t = reconstructed.get(name)
            if t:
                tags.append(t)
        elif flag == 0x30:  # dropped qt
            t = reconstructed.get(name)
            if t:
                tags.append(t)
        else:
            t = stored.get(name)
            if t:
                tags.append(t)

    return tags


# ── streaming helpers ─────────────────────────────────────────────────────────

def sread(stream, n):
    buf = bytearray()
    while len(buf) < n:
        chunk = stream.read(n - len(buf))
        if not chunk:
            raise EOFError(f"Expected {n} bytes, got {len(buf)}")
        buf.extend(chunk)
    return bytes(buf)


def flush_block(out, ch_meta, ch_dna, ch_qual, ch_tags, n):
    out.write(struct.pack(">BI", 0x10, n))
    out.write(ch_meta)
    out.write(ch_dna)
    out.write(ch_qual)
    out.write(ch_tags)


# ── compress ──────────────────────────────────────────────────────────────────

def compress():
    out = sys.stdout.buffer
    out.write(MAGIC)

    ch_meta = bytearray(); ch_dna  = bytearray()
    ch_qual = bytearray(); ch_tags = bytearray()
    n_in_block = 0

    for line in sys.stdin:
        if line.startswith('@'):
            b = line.encode('utf-8')
            out.write(struct.pack(">BI", 0, len(b)))
            out.write(b)
            continue

        cols = line.rstrip('\r\n').split('\t')
        if len(cols) < 11:
            continue

        meta    = "\t".join(cols[0:9]).encode('utf-8')
        seq_str = cols[9]
        qual_b  = cols[10].encode('latin-1')
        sl      = len(seq_str)

        dna_bin, exc_bin = pack_dna_k4(seq_str)
        qual_enc         = encode_qual(qual_b)
        tag_payload, tag_order = encode_tags(cols[11:])

        ch_meta += struct.pack(">BII", 1, len(meta), sl) + meta
        ch_dna  += (struct.pack(">BI", 2, len(dna_bin)) + dna_bin
                  + struct.pack(">I",  len(exc_bin))    + exc_bin)
        # qual: type + 2-byte enc_len + enc_bytes
        ch_qual += struct.pack(">BH", 3, len(qual_enc)) + qual_enc
        # tags: type + 4-byte payload_len + payload + 2-byte order_len + order
        ch_tags += (struct.pack(">BI", 4, len(tag_payload)) + tag_payload
                  + struct.pack(">H",    len(tag_order))    + tag_order)
        n_in_block += 1

        if n_in_block == BLOCK_SIZE:
            flush_block(out, ch_meta, ch_dna, ch_qual, ch_tags, n_in_block)
            ch_meta = bytearray(); ch_dna  = bytearray()
            ch_qual = bytearray(); ch_tags = bytearray()
            n_in_block = 0

    if n_in_block:
        flush_block(out, ch_meta, ch_dna, ch_qual, ch_tags, n_in_block)

    out.write(struct.pack(">BI", 0, 0))
    out.write(struct.pack(">B",  0xFF))


# ── decompress ────────────────────────────────────────────────────────────────

def decompress():
    inp = sys.stdin.buffer
    out = sys.stdout

    if sread(inp, 4) != MAGIC:
        sys.stderr.write("Bad magic -- is this a ks9 (KS9D) file?\n")
        return

    while True:
        type_b = inp.read(1)
        if not type_b:
            break
        rec_type = type_b[0]

        if rec_type == 0x00:
            l = struct.unpack(">I", sread(inp, 4))[0]
            if l == 0:
                continue
            out.write(sread(inp, l).decode('utf-8'))

        elif rec_type == 0x10:
            n = struct.unpack(">I", sread(inp, 4))[0]
            meta_list  = []
            dna_list   = []
            qual_list  = []
            tag_list   = []

            for _ in range(n):
                assert sread(inp, 1)[0] == 1
                ml, sl = struct.unpack(">II", sread(inp, 8))
                meta_list.append((sread(inp, ml).decode('utf-8'), sl))

            for _ in range(n):
                assert sread(inp, 1)[0] == 2
                dl  = struct.unpack(">I", sread(inp, 4))[0]
                dna = sread(inp, dl)
                el  = struct.unpack(">I", sread(inp, 4))[0]
                exc = sread(inp, el)
                dna_list.append((dna, exc))

            for _ in range(n):
                assert sread(inp, 1)[0] == 3
                qe_len = struct.unpack(">H", sread(inp, 2))[0]
                qual_list.append(sread(inp, qe_len))

            for _ in range(n):
                assert sread(inp, 1)[0] == 4
                tl         = struct.unpack(">I", sread(inp, 4))[0]
                tag_payload = sread(inp, tl)
                ol          = struct.unpack(">H", sread(inp, 2))[0]
                tag_order   = sread(inp, ol)
                tag_list.append((tag_payload, tag_order))

            for i in range(n):
                meta, sl = meta_list[i]
                qname    = meta.split('\t')[0]
                dna      = unpack_dna_k4(dna_list[i][0], sl, dna_list[i][1])
                qual     = decode_qual(qual_list[i], sl)
                tags     = decode_tags(tag_list[i][0], tag_list[i][1], qname)

                line = f"{meta}\t{dna}\t{qual}"
                if tags:
                    line += "\t" + "\t".join(tags)
                out.write(line + "\n")

        elif rec_type == 0xFF:
            break

        else:
            sys.stderr.write(f"Unknown record type 0x{rec_type:02x}\n")
            break


# ── entry point ───────────────────────────────────────────────────────────────

def process():
    if len(sys.argv) < 2 or sys.argv[1] not in ('-c', '-u'):
        sys.stderr.write(__doc__ + "\n")
        sys.exit(1)
    if sys.argv[1] == '-c':
        compress()
    else:
        decompress()

if __name__ == '__main__':
    process()
