# kineticsquash — PacBio Subread BAM Pre-conditioner for xz Compression

A lossless, reversible pre-compression filter for PacBio subread BAM files. It operates as a pipe stage between `samtools view` and `xz`, transforming the SAM byte stream into a binary format that xz can compress significantly more efficiently than either raw SAM or raw BAM.

Python version and RUST version included

---

## Quick Start - python

```bash
# Compress
samtools view -h subreads.bam | python3 kineticsquash.py -c | xz -9e > subreads.sq.xz

# Decompress
xz -d < subreads.sq.xz | python3 kineticsquash.py -u | samtools view -b - > restored.bam
```

**Requirements:** Python 3.6+, samtools, xz. No third-party Python packages.

---

## Quick Start - RUST

```bash
# Compress
./kineticsquash subreads.bam
Successfully compressed to subreads.sq.xz
[output is subreads.sq.xz]

# Decompress
kineticsquash subreads.sq.xz
Successfully decompressed to subreads.bam
```

**Requirements:** Linux OS

---

## Why This Exists

A PacBio Revio run produces an unaligned BAM file of 300–600 GB. BAM is already BGZF-compressed (block gzip), but for long-term archival the goal is maximum compression, typically with `xz -9e` which uses the LZMA2 algorithm. The problem is that naively applying `xz` to the SAM representation of the data achieves poor ratios because:

1. **SAM is ASCII-verbose.** Numeric values like kinetic tags are stored as comma-separated decimal strings. `ip:B:C,5,12,3,0,8,...` uses 3–4 bytes per value as ASCII, vs 1 byte as a raw integer.
2. **The data streams are heterogeneous.** In standard SAM/BAM, each read record contains all its fields interleaved: sequence, quality, kinetics, metadata. When this enters xz, the LZMA2 dictionary is shared across all data types simultaneously. The ip/pw kinetic arrays from one read compete for dictionary space with the sequence from the next, preventing the compressor from finding the long-range cross-read patterns that matter most for kinetic data.
3. **Redundant and zero-entropy data occupies dictionary space.** PacBio subread QUAL is universally all `!` (ASCII 33, Phred 0) — a mandatory SAM field filled with a placeholder because subreads have no per-base quality score. Several tags (`zm:i`, `qs:i`, `qe:i`, `qt:Z:`) duplicate information already present in the QNAME.

kineticsquash addresses all three issues with a binary pre-encoding pass.

---

## Usage

### Compression

```bash
samtools view -h subreads.bam | python3 kineticsquash.py -c | xz -9e > subreads.sq.xz
```

The `-h` flag to `samtools view` includes the SAM header lines (beginning with `@`), which are required for round-trip reconstruction of a valid BAM. Without `-h`, the header is lost and the restored BAM will be invalid.

`xz -9e` uses LZMA2 at maximum compression (preset 9, extreme mode). This is slow but achieves the best ratio for archival. For faster compression at modest cost to ratio, `xz -6` is a reasonable alternative. The kineticsquash pre-encoding is independent of the xz preset — it improves the input quality regardless of which preset is used.

### Decompression

```bash
xz -d < subreads.sq.xz | python3 kineticsquash.py -u | samtools view -b - > restored.bam
```

The decompressed output is valid SAM (text) piped into `samtools view -b` to produce binary BAM. The `-` argument tells samtools to read from stdin.

### Parallelism

xz is single-threaded by default. For faster compression on multi-core systems:

```bash
samtools view -h subreads.bam | python3 kineticsquash.py -c | xz -9e -T 0 > subreads.sq.xz
```

`-T 0` uses all available cores. Note that multithreaded xz uses independent LZMA streams per thread, which slightly reduces the compression ratio compared to single-threaded mode because cross-thread dictionary references are not possible.

---

## Architecture

### The Pipeline

```
subreads.bam
     │
     ▼  samtools view -h
SAM text stream (ASCII, ~3–4× BAM size when decompressed)
     │
     ▼  kineticsquash.py -c
Binary kineticsquashD stream:
  - 2-bit packed DNA
  - 1-byte qual mode flags (not seq_len bytes)
  - raw binary ip/pw arrays
  - deduplicated tags
  - block-striped channel layout
     │
     ▼  xz -9e
subreads.sq9.xz
```

The decompression pipeline runs in exact reverse.

### Memory Model

Both the compressor and decompressor are streaming and bounded in RAM. Records are accumulated into four channel buffers up to `BLOCK_SIZE` reads (default 5000), then flushed. At no point is the entire file held in memory.

At default settings on Revio data (avg ~15,000 bases/read):

| Channel | Size per block |
|---------|----------------|
| META    | ~10 MB         |
| DNA     | ~19 MB         |
| QUAL    | ~0.04 MB (1 byte/read after encoding) |
| TAGS    | ~150 MB        |
| **Total peak RAM** | **~300 MB** |

This is constant regardless of whether the input file is 1 GB or 600 GB.

---

## Encoding Elements in Detail

### 1. File Header and Magic Bytes

Every kineticsquash file begins with the 4-byte magic sequence `kineticsquashD` followed immediately by SAM `@`-header records, each encoded as:

```
[0x00: 1 byte] [length: uint32 big-endian] [SAM header line: length bytes]
```

The header section is terminated by a zero-length sentinel:

```
[0x00: 1 byte] [0x00000000: 4 bytes]
```

This design allows the decompressor to emit SAM header lines immediately as they are read, before any read data is parsed.

### 2. Block Structure

After the header section, the file consists of a sequence of blocks. Each block begins with:

```
[0x10: 1 byte] [n_reads: uint32 big-endian]
```

Followed by four consecutive channel sections, each containing exactly `n_reads` records in the same read order:

```
[Channel META:  n_reads records]
[Channel DNA:   n_reads records]
[Channel QUAL:  n_reads records]
[Channel TAGS:  n_reads records]
```

The file ends with a single `0xFF` byte after the last block.

**Why blocks rather than one fully-striped file?**

Striping the entire file's channels would require buffering every read in RAM — infeasible for a 500 GB Revio BAM. The block approach limits RAM to one block at a time in both the compressor and decompressor. The block size of 5000 reads is chosen so that each channel block (particularly ip/pw at ~75 MB per block for 15 kb reads) exceeds xz's 64 MB LZMA2 dictionary. Once a channel block exceeds the dictionary size, xz's cross-read pattern matching is fully saturated — larger blocks give diminishing returns. Smaller blocks (e.g. 1000 reads) reduce RAM at the cost of shorter xz back-reference distances.

**Why channel separation helps xz**

LZMA2 finds compression by identifying repeated byte sequences (LZ77 matches) and entropy-coding the residuals (range coding). Long matches across reads — where the ip array of read N resembles the ip array of read N+1 — are the most valuable because kinetic profiles on the same flow cell share similar distributions. In the interleaved layout of standard SAM, these ip arrays are separated by ~10 KB of sequence, quality, and tag data from the intervening reads. In the block-striped layout, all ip data within a block is contiguous, so xz's dictionary is maximally populated with comparable data when encoding each ip array.

### 3. DNA Encoding (2-bit Pack with Exception Map)

**Channel record format:**
```
[0x02: 1 byte] [dna_len: uint32] [packed_bytes: dna_len bytes]
[exc_len: uint32] [exception_bytes: exc_len bytes]
```

**2-bit encoding:**

The four canonical bases are mapped to 2-bit values:

| Base | Bits |
|------|------|
| A    | 00   |
| C    | 01   |
| G    | 10   |
| T    | 11   |

Four bases are packed into one byte, most-significant bits first. A sequence of length `L` requires `ceil(L/4)` bytes, a 4× reduction over ASCII. This is the theoretical minimum for a 4-symbol alphabet with uniform distribution (2 bits/base = Shannon entropy of a uniform 4-symbol source).

**Exception map:**

Non-ACGT bases (IUPAC ambiguity codes: `N`, `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`) are encoded as `A` (00) in the 2-bit stream, and their true values are stored separately in the exception map:

```
[count: uint16 big-endian]
[per-exception: position(uint32) char(uint8)]  ×count
```

This handles ambiguous bases without changing the storage cost for the common case. For typical PacBio subread data, non-ACGT bases are uncommon (adapter trimming removes most), so the exception map is usually a few dozen bytes at most per read.

**Why not a more sophisticated DNA encoding?**

After 2-bit packing, the DNA channel looks like uniformly random bytes — xz achieves essentially no further compression on it (it is already at the Shannon limit of 2 bits/base). Alternatives such as k-mer tables, Markov models, or reference-based encoding were considered but ruled out: kineticsquash targets *unaligned* subreads, which have no reference, and the additional complexity is not justified for a channel that is already at its theoretical minimum.

### 4. QUAL Encoding — Critical Optimisation

**Channel record format:**
```
[0x03: 1 byte] [enc_len: uint16] [encoded_qual: enc_len bytes]
```

**The core observation:**

PacBio subread BAM files do not contain per-base quality scores. The QUAL field in SAM (column 11) is a mandatory field, so the PacBio toolchain fills it with the ASCII character `!` (decimal 33, Phred score 0) for every base of every read. A read of 15,000 bases produces 15,000 identical bytes of QUAL data with zero information content.

In ks8, this placeholder was stored verbatim: `seq_len` bytes per read. xz compresses a run of identical bytes efficiently (down to ~27 bytes regardless of length), but those 27 bytes of compressed output are not the problem. The problem is that xz's LZMA dictionary is finite (64 MB at `-9e`) and the 15,000 *input* bytes of `!!!!...!!!!` still occupy that dictionary while being processed, displacing ip and pw data that would otherwise find long-range matches across reads.

**Encoding modes:**

| Mode byte | Meaning | Stored bytes |
|-----------|---------|--------------|
| `0x00`    | All `!` (ASCII 33) | 1 byte total |
| `0x01`    | Constant value V | 2 bytes total (mode + value) |
| `0x02`    | Literal per-base quality | 1 + seq_len bytes |

Mode `0x00` handles 100% of real PacBio subread data observed. Mode `0x01` handles edge cases such as reads where all bases happen to have the same non-`!` quality. Mode `0x02` provides a lossless fallback for HiFi reads or any data where real per-base quality scores are present — the format degrades gracefully to storing the full quality string, and the round-trip is still bit-perfect.

**Impact at scale:** For a 500 GB Revio file at ~30 M reads of 15,000 bases, the QUAL channel would have contributed approximately 450 GB to the xz input in ks8. In kineticsquash it contributes 30 MB (1 byte per read). This is the single largest absolute saving.

### 5. Redundant Tag Elimination — Critical Optimisation

**Tags dropped during compression, reconstructed during decompression:**

| Tag | Reason for dropping | Reconstruction source |
|-----|--------------------|-----------------------|
| `zm:i:N` | ZMW number is QNAME field 2 | `qname.split('/')[1]` |
| `qs:i:N` | Query start is QNAME field 3, before `_` | `qname.split('/')[2].split('_')[0]` |
| `qe:i:N` | Query end is QNAME field 3, after `_` | `qname.split('/')[2].split('_')[1]` |
| `qt:Z:!!!!` | Barcode quality, always all `!` | `'!' * len(bt_tag[5:])` |

The PacBio QNAME format is `{movie}/{zmw}/{qs}_{qe}`, for example `m64018_201129_132425/61867919/8400_11446`. All four values are completely derivable from the QNAME, which must be stored anyway as it is the read identifier. Storing them additionally as tags wastes approximately 60 bytes per read in ASCII.

**Tag order preservation:**

The SAM specification states that tags can appear in any order, so reconstructing them in a different order would produce a valid file. However, to guarantee byte-identical output and protect against any downstream tools sensitive to tag ordering, the original column order is preserved via a compact order record appended to each read's tag payload:

```
[n_tags: 1 byte]
[flag: 1 byte][name: 2 bytes]  ×n_tags
```

The flag byte encodes how the tag was handled:

| Flag | Meaning |
|------|---------|
| `0x00` | Stored normally in payload |
| `0x10` | Stored as binary ip/pw array |
| `0x20` | Dropped (zm/qs/qe); reconstruct from QNAME |
| `0x30` | Dropped (qt); reconstruct from bt length |

On decompression, tags are emitted in their original column order by replaying this record. For a read with 19 tags (the count in the analysed Sequel II data), the order record is 58 bytes — negligible overhead relative to the savings.

### 6. ip and pw Kinetic Tag Encoding

**Tag payload record format (within the tags channel):**
```
[0x02: 1 byte] [name: 2 bytes, 'ip' or 'pw'] [length: uint32] [values: length bytes]
```

The `ip` (inter-pulse duration) and `pw` (pulse width) tags are stored as binary `uint8` arrays. Each value occupies one raw byte. In SAM, these same values are stored as ASCII decimal strings separated by commas — `ip:B:C,5,12,3,0,8,...` — which uses 3–4 bytes per value. The binary encoding is a 3–4× reduction in storage before compression, and presents xz with a byte stream of small integers rather than a mix of digits and punctuation.

**The PacBio CodecV1 distribution:**

PacBio encodes IPD and pulse width using a non-linear quantisation scheme (CodecV1). Codepoints 0–63 represent single-frame resolution (one codepoint per frame), codepoints 64–127 represent two-frame resolution, and so on in doublings. The codepoint distribution is therefore heavily right-skewed: in real Sequel II data, approximately 85% of ip values fall below 64 and 95% of pw values fall below 64, with both distributions peaking in the range 3–15.

**Encodings that were tried and rejected:**

*Delta encoding (first differences):* The intuition is that smooth kinetic traces would produce small deltas that compress better. In practice, delta encoding of ip/pw *increases* compressed size by approximately 13% on real data. The geometric distribution of raw values already compresses well because xz can model the skewed histogram; their first differences lose the skew and approach a uniform distribution, which compresses worse. This was confirmed empirically by measuring zlib compression of raw vs delta-encoded ip arrays extracted from real Sequel II subread data.

*6-bit packing (values 0–63 packed at 6 bits each):* Since 85% of ip values fit in 6 bits, packing the main stream at 6 bits/value and storing overflow values separately in a sidecar stream seems attractive. In practice, the overflow sidecar requires storing both the position (as a 2-byte index) and value (1 byte) for each of the 15% of values exceeding 63. The framing overhead of this sidecar exceeds the saving from 6-bit packing: total compressed size is larger than raw bytes + xz. Measured on real data, the 6-bit approach is worse by approximately 38%.

*Context-based grouping (values grouped by preceding base A/C/G/T):* Bases in the same nucleotide context have similar kinetic profiles in the polymerase model. Reorganising ip values by the base they follow before compression theoretically improves homogeneity. Empirically, this produced no measurable improvement — the compressor found the same patterns regardless of grouping, presumably because the per-context distributions are similar enough that separation does not help.

The conclusion is that raw binary storage of ip/pw values, with block-level channel separation to maximise xz dictionary utilisation, is the optimal approach for a general-purpose lossless scheme without a bespoke entropy coder.

### 7. Other Tags (Pass-through)

All tags not specifically handled above are stored as their full ASCII SAM representation, length-prefixed:

```
[0x03: 1 byte] [string_len: uint16] [tag_string: string_len bytes]
```

This includes: `bx:B:i`, `np:i`, `rq:f`, `sn:B:f`, `we:i`, `ws:i`, `bc:B:S`, `bq:i`, `cx:i`, `bl:Z:`, `bt:Z:`, `ql:Z:`, `RG:Z:`.

The `sn:B:f` tag (signal-to-noise ratio, four IEEE 754 float32 values) is per-ZMW and constant across all subreads of the same ZMW. It was considered as a candidate for per-ZMW deduplication (store once, reference multiple times), but since xz naturally compresses the repeated float byte patterns across reads from the same ZMW via LZ77 back-references, the implementation complexity is not justified for a tag that represents only ~38 bytes per read.

---

## What Helps and What Does Not

### Critical — significant impact

| Technique | Impact |
|-----------|--------|
| **QUAL elimination** | Removes ~15 KB/read of zero-entropy placeholder bytes from the xz input. On a 500 GB Revio file (~30 M reads), this removes ~450 GB from the xz input stream, directly freeing LZMA dictionary capacity for ip/pw. This is the single largest optimisation. |
| **Redundant tag elimination** | Drops `zm:i`, `qs:i`, `qe:i`, `qt:Z:` which are byte-for-byte duplicates of QNAME-derived or zero-entropy data. Small per-read saving (~60 bytes) but requires no approximation or information loss. |
| **2-bit DNA packing** | 4× reduction in sequence storage (8 bits/base ASCII → 2 bits/base binary). DNA is already at Shannon entropy after packing; xz cannot compress it further. |
| **ip/pw as raw binary** | Converts `ip:B:C,5,12,...` (3–4 bytes/value ASCII) to raw uint8 arrays (1 byte/value). 3–4× reduction on the two largest data channels. |

### Marginal — small or file-size-dependent impact

| Technique | Why it is marginal |
|-----------|-------------------|
| **Block-striped channel layout** | In theory: separating channels so xz sees homogeneous data improves dictionary utilisation. On a small test file (< 64 MB total): negligible improvement, because xz's 64 MB dictionary already spans the entire file regardless of layout. On a full Revio file (500 GB), striping is expected to provide a real benefit as each ~75 MB ip channel block fully saturates the dictionary, but the gain is difficult to measure without a full-scale test file. |

### Actively harmful — do not use

| Technique | Why it hurts |
|-----------|-------------|
| **Delta encoding of ip/pw** | ip/pw values follow a geometric (right-skewed) distribution. Their first differences are more uniformly distributed (higher entropy). Measured at approximately 13% worse compressed size than raw bytes. |
| **Delta encoding of QUAL** | Subread QUAL is all identical bytes — there is nothing to delta. For HiFi data with real quality scores, delta encoding ASCII Phred values does not consistently help because the ASCII range (33–126) is already a restricted alphabet that xz handles efficiently. |
| **QNAME delta encoding** | ZMW numbers in a subread BAM are not necessarily monotonically increasing or sorted. Zigzag-encoded deltas can be larger than the original ASCII digits. Complexity overhead is not recovered in practice. |
| **6-bit packing of ip/pw** | The sidecar stream required for the 15% of values exceeding 63 adds more overhead than the 6-bit main stream saves. Measured at approximately 38% worse than raw bytes after compression. |

---

## On-Disk Format Reference

All multi-byte integers are big-endian.

```
[Magic: 4 bytes 'kineticsquashD']

--- Header section ---
[0x00][len: uint32][SAM @header line: len bytes]   (one per header line)
[0x00][0x00000000]                                  (sentinel)

--- Block section (repeated until 0xFF) ---
[0x10][n_reads: uint32]

  Channel META (n_reads records):
    [0x01][meta_len: uint32][seq_len: uint32][meta_bytes: meta_len]
    meta_bytes = TAB-joined SAM columns 0–8 (QNAME through TLEN)

  Channel DNA (n_reads records):
    [0x02][dna_len: uint32][packed_dna: dna_len bytes]
           [exc_len: uint32][exception_bytes: exc_len bytes]
    exception_bytes = [count: uint16][pos: uint32, char: uint8] ×count

  Channel QUAL (n_reads records):
    [0x03][enc_len: uint16][encoded_qual: enc_len bytes]
    encoded_qual[0] = mode:
      0x00 = all '!'  → enc_len=1
      0x01 = constant → enc_len=2, encoded_qual[1]=value byte
      0x02 = literal  → enc_len=1+seq_len, encoded_qual[1:]=raw qual bytes

  Channel TAGS (n_reads records):
    [0x04][payload_len: uint32][tag_payload: payload_len bytes]
          [order_len: uint16][order_bytes: order_len bytes]

    tag_payload records (concatenated, type-dispatched):
      ip/pw:  [0x02][name: 2 bytes][count: uint32][values: count bytes]
      other:  [0x03][str_len: uint16][tag_string: str_len bytes]

    order_bytes:
      [n_tags: 1 byte]
      [flag: 1 byte][name: 2 bytes]  ×n_tags
      flag values:
        0x00 = stored as normal string tag
        0x10 = stored as binary ip/pw array
        0x20 = dropped (zm/qs/qe); reconstruct from QNAME on decompress
        0x30 = dropped (qt all-bang); reconstruct from bt tag length on decompress

[0xFF]   (end marker)
```

---

## Compatibility Notes

**Not compatible with ks8.** The binary format is different; ks8 files must be decompressed with the matching ks8 script.

**SAM spec compliance.** The restored SAM is fully valid per the SAM v1.6 specification. Tag order is preserved exactly; all dropped tags are reconstructed with their original values.

**CCS / HiFi reads.** HiFi reads have real per-base quality scores (QUAL ≠ all `!`). The QUAL encoder handles this via mode `0x02` (literal), storing the full quality string with no loss. The compression benefit for this channel is smaller for HiFi files, but all other optimisations (2-bit DNA, binary ip/pw, redundant tag elimination, channel striping) still apply fully.

**Non-PacBio data.** The `zm`/`qs`/`qe`/`qt` tag elimination assumes PacBio QNAME format (`movie/zmw/qs_qe`). For reads with a different QNAME format, the tag parsing returns the QNAME unchanged and the tags remain in the payload. The script will not corrupt non-PacBio data but will not apply the tag deduplication optimisations.

---

## Tuning `BLOCK_SIZE`

The constant `BLOCK_SIZE = 5000` near the top of the script controls the number of reads accumulated before each channel flush. It trades RAM usage against compression quality:

| BLOCK_SIZE | Approx RAM (15 kb reads) | ip channel block size | Notes |
|------------|--------------------------|----------------------|-------|
| 1000       | ~60 MB                   | ~15 MB               | ip block < 64 MB xz dict; reduced cross-read matching |
| 5000       | ~300 MB                  | ~75 MB               | Default; ip block > dict, fully saturated |
| 20000      | ~1.2 GB                  | ~300 MB              | Diminishing returns beyond 5000 |

For shorter reads (e.g. amplicon sequencing, avg 500 bp), a larger `BLOCK_SIZE` is needed to fill the xz dictionary: at 500 bp/read, 5000 reads yields only 2.5 MB per ip channel block, well below the 64 MB dictionary. In that case, increase `BLOCK_SIZE` to approximately 130,000 to achieve the same dictionary saturation — though at ~4 GB RAM. For short-read data the benefit of striping is modest regardless, as the reads are short enough that xz finds cross-read patterns naturally at any block size.
