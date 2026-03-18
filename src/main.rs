use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use std::collections::HashMap;
use std::env;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::process::{Command, Stdio};
use xz2::read::XzDecoder;
use xz2::write::XzEncoder;

const MAGIC: &[u8; 4] = b"KS9D";
const BLOCK_SIZE: u32 = 5000;

// ── DNA ───────────────────────────────────────────────────────────────────────
fn pack_dna_k4(seq: &str) -> (Vec<u8>, Vec<u8>) {
    let mut packed = Vec::new();
    let mut exceptions = Vec::new();
    let seq_bytes = seq.as_bytes();

    let mut i = 0;
    while i < seq_bytes.len() {
        let mut byte_val = 0u8;
        for j in 0..4 {
            if i + j >= seq_bytes.len() { break; }
            let base = seq_bytes[i + j].to_ascii_uppercase();
            let val = match base {
                b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
                _ => { exceptions.push(((i + j) as u32, base)); 0 }
            };
            byte_val |= val << (6 - j * 2);
        }
        packed.push(byte_val);
        i += 4;
    }

    let mut exc_bin = Vec::new();
    exc_bin.write_u16::<BigEndian>(exceptions.len() as u16).unwrap();
    for (pos, char_code) in exceptions {
        exc_bin.write_u32::<BigEndian>(pos).unwrap();
        exc_bin.write_u8(char_code).unwrap();
    }
    (packed, exc_bin)
}

fn unpack_dna_k4(packed: &[u8], length: usize, exc_bin: &[u8]) -> String {
    let mut seq = Vec::with_capacity(length);
    let rev_map = [b'A', b'C', b'G', b'T'];

    for i in 0..length {
        let byte_idx = i / 4;
        let sub_idx = i % 4;
        if byte_idx >= packed.len() { break; }
        let val = (packed[byte_idx] >> (6 - sub_idx * 2)) & 0x03;
        seq.push(rev_map[val as usize]);
    }

    if exc_bin.len() >= 2 {
        let mut cursor = io::Cursor::new(exc_bin);
        let exc_count = cursor.read_u16::<BigEndian>().unwrap_or(0);
        for _ in 0..exc_count {
            if let (Ok(pos), Ok(char_code)) = (cursor.read_u32::<BigEndian>(), cursor.read_u8()) {
                if pos < seq.len() as u32 { seq[pos as usize] = char_code; }
            }
        }
    }
    String::from_utf8(seq).unwrap()
}

// ── QUAL ──────────────────────────────────────────────────────────────────────
fn encode_qual(qual: &[u8]) -> Vec<u8> {
    if qual.is_empty() { return vec![0x02]; }
    let first = qual[0];
    let all_same = qual.iter().all(|&b| b == first);

    if all_same {
        if first == b'!' { vec![0x00] } else { vec![0x01, first] }
    } else {
        let mut out = vec![0x02];
        out.extend_from_slice(qual);
        out
    }
}

fn decode_qual(qual_enc: &[u8], seq_len: usize) -> String {
    if qual_enc.is_empty() { return String::new(); }
    let mode = qual_enc[0];
    if mode == 0x00 {
        "!".repeat(seq_len)
    } else if mode == 0x01 {
        let c = qual_enc.get(1).copied().unwrap_or(b'!') as char;
        c.to_string().repeat(seq_len)
    } else {
        qual_enc[1..].iter().map(|&b| b as char).collect()
    }
}

// ── TAGS ──────────────────────────────────────────────────────────────────────
fn encode_tags(tags: &[&str]) -> (Vec<u8>, Vec<u8>) {
    let mut payload = Vec::new();
    let mut order = Vec::new();
    order.push(tags.len() as u8);

    for t in tags {
        if t.len() < 2 { continue; }
        let name = &t[0..2];

        if matches!(name, "zm" | "qs" | "qe") {
            order.push(0x20); order.extend_from_slice(name.as_bytes());
            continue;
        }

        if t.starts_with("qt:Z:") && t[5..].chars().all(|c| c == '!') {
            order.push(0x30); order.extend_from_slice(name.as_bytes());
            continue;
        }

        if t.starts_with("ip:B:C,") || t.starts_with("pw:B:C,") {
            let vals: Vec<u8> = t[7..].split(',').filter_map(|x| x.parse().ok()).collect();
            payload.write_u8(2).unwrap();
            payload.extend_from_slice(name.as_bytes());
            payload.write_u32::<BigEndian>(vals.len() as u32).unwrap();
            payload.extend_from_slice(&vals);
            order.push(0x10); order.extend_from_slice(name.as_bytes());
        } else {
            let t_b = t.as_bytes();
            payload.write_u8(3).unwrap();
            payload.write_u16::<BigEndian>(t_b.len() as u16).unwrap();
            payload.extend_from_slice(t_b);
            order.push(0x00); order.extend_from_slice(name.as_bytes());
        }
    }
    (payload, order)
}

fn decode_tags(tp: &[u8], order_bytes: &[u8], qname: &str) -> Vec<String> {
    if order_bytes.is_empty() { return Vec::new(); }
    let n_tags = order_bytes[0] as usize;
    let mut order = Vec::new();
    for i in 0..n_tags {
        let base = 1 + i * 3;
        if base + 3 > order_bytes.len() { break; }
        let flag = order_bytes[base];
        let name = std::str::from_utf8(&order_bytes[base+1..base+3]).unwrap_or("").to_string();
        order.push((flag, name));
    }

    let mut stored = HashMap::new();
    let mut off = 0;
    while off < tp.len() {
        let kind = tp[off]; off += 1;
        if kind == 2 {
            let name = std::str::from_utf8(&tp[off..off+2]).unwrap_or(""); off += 2;
            let mut cur = io::Cursor::new(&tp[off..off+4]);
            let v_len = cur.read_u32::<BigEndian>().unwrap() as usize; off += 4;
            let vals = &tp[off..off+v_len]; off += v_len;
            let vals_str: Vec<String> = vals.iter().map(|b| b.to_string()).collect();
            stored.insert(name.to_string(), format!("{}:B:C,{}", name, vals_str.join(",")));
        } else if kind == 3 {
            let mut cur = io::Cursor::new(&tp[off..off+2]);
            let g_len = cur.read_u16::<BigEndian>().unwrap() as usize; off += 2;
            let tag_str = std::str::from_utf8(&tp[off..off+g_len]).unwrap_or(""); off += g_len;
            stored.insert(tag_str[..2].to_string(), tag_str.to_string());
        }
    }

    let mut reconstructed = HashMap::new();
    let parts: Vec<&str> = qname.split('/').collect();
    if parts.len() == 3 {
        reconstructed.insert("zm".to_string(), format!("zm:i:{}", parts[1]));
        let coords: Vec<&str> = parts[2].split('_').collect();
        if coords.len() == 2 {
            reconstructed.insert("qs".to_string(), format!("qs:i:{}", coords[0]));
            reconstructed.insert("qe".to_string(), format!("qe:i:{}", coords[1]));
        }
    }
    if let Some(bt_tag) = stored.get("bt") {
        if bt_tag.starts_with("bt:Z:") {
            reconstructed.insert("qt".to_string(), format!("qt:Z:{}", "!".repeat(bt_tag.len() - 5)));
        }
    }

    let mut tags = Vec::new();
    for (flag, name) in order {
        if flag == 0x20 || flag == 0x30 {
            if let Some(t) = reconstructed.get(&name) { tags.push(t.clone()); }
        } else if let Some(t) = stored.get(&name) {
            tags.push(t.clone());
        }
    }
    tags
}

fn flush_block(out: &mut impl Write, ch_meta: &[u8], ch_dna: &[u8], ch_qual: &[u8], ch_tags: &[u8], n: u32) {
    out.write_u8(0x10).unwrap();
    out.write_u32::<BigEndian>(n).unwrap();
    out.write_all(ch_meta).unwrap();
    out.write_all(ch_dna).unwrap();
    out.write_all(ch_qual).unwrap();
    out.write_all(ch_tags).unwrap();
}

// ── COMPRESS ──────────────────────────────────────────────────────────────────
fn compress(input_bam: &str) {
    let output_xz = input_bam.replace(".bam", ".sq.xz");
    
    let mut samtools = Command::new("samtools")
        .arg("view").arg("-h").arg(input_bam)
        .stdout(Stdio::piped()).spawn().expect("Failed to start samtools.");
        
    let stdout = samtools.stdout.take().unwrap();
    let reader = BufReader::new(stdout);
    
    let file = std::fs::File::create(&output_xz).unwrap();
    let mut out = XzEncoder::new(file, 9);
    out.write_all(MAGIC).unwrap();

    let mut ch_meta = Vec::new();
    let mut ch_dna = Vec::new();
    let mut ch_qual = Vec::new();
    let mut ch_tags = Vec::new();
    let mut n_in_block = 0u32;

    for line_result in reader.lines() {
        let line = line_result.unwrap();
        if line.starts_with('@') {
            let b = line.as_bytes();
            out.write_u8(0).unwrap();
            out.write_u32::<BigEndian>(b.len() as u32).unwrap();
            out.write_all(b).unwrap();
            continue;
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 11 { continue; }

        let meta = cols[0..9].join("\t");
        let meta_bytes = meta.as_bytes();
        let seq_str = cols[9];
        let qual_b = cols[10].as_bytes();
        let sl = seq_str.len() as u32;

        let (dna_bin, exc_bin) = pack_dna_k4(seq_str);
        let qual_enc = encode_qual(qual_b);
        let (tag_payload, tag_order) = encode_tags(&cols[11..]);

        ch_meta.write_u8(1).unwrap();
        ch_meta.write_u32::<BigEndian>(meta_bytes.len() as u32).unwrap();
        ch_meta.write_u32::<BigEndian>(sl).unwrap();
        ch_meta.write_all(meta_bytes).unwrap();

        ch_dna.write_u8(2).unwrap();
        ch_dna.write_u32::<BigEndian>(dna_bin.len() as u32).unwrap();
        ch_dna.write_all(&dna_bin).unwrap();
        ch_dna.write_u32::<BigEndian>(exc_bin.len() as u32).unwrap();
        ch_dna.write_all(&exc_bin).unwrap();

        ch_qual.write_u8(3).unwrap();
        ch_qual.write_u16::<BigEndian>(qual_enc.len() as u16).unwrap();
        ch_qual.write_all(&qual_enc).unwrap();

        ch_tags.write_u8(4).unwrap();
        ch_tags.write_u32::<BigEndian>(tag_payload.len() as u32).unwrap();
        ch_tags.write_all(&tag_payload).unwrap();
        ch_tags.write_u16::<BigEndian>(tag_order.len() as u16).unwrap();
        ch_tags.write_all(&tag_order).unwrap();

        n_in_block += 1;
        if n_in_block == BLOCK_SIZE {
            flush_block(&mut out, &ch_meta, &ch_dna, &ch_qual, &ch_tags, n_in_block);
            ch_meta.clear(); ch_dna.clear(); ch_qual.clear(); ch_tags.clear();
            n_in_block = 0;
        }
    }

    if n_in_block > 0 {
        flush_block(&mut out, &ch_meta, &ch_dna, &ch_qual, &ch_tags, n_in_block);
    }
    out.write_u8(0).unwrap(); out.write_u32::<BigEndian>(0).unwrap();
    out.write_u8(0xFF).unwrap();
    
    let _ = samtools.wait();
    println!("Successfully compressed to {}", output_xz);
}

// ── DECOMPRESS ────────────────────────────────────────────────────────────────
fn decompress(input_xz: &str) {
    let output_bam = input_xz.replace(".sq.xz", ".bam");
    
    let file = std::fs::File::open(input_xz).unwrap();
    let decoder = XzDecoder::new(file);
    let mut reader = BufReader::new(decoder);
    
    let mut magic_buf = [0u8; 4];
    if reader.read_exact(&mut magic_buf).is_err() || &magic_buf != MAGIC {
        eprintln!("Bad magic -- is this a KS9D file?");
        return;
    }

    let mut samtools = Command::new("samtools")
        .arg("view").arg("-b").arg("-")
        .stdin(Stdio::piped())
        .stdout(std::fs::File::create(&output_bam).unwrap())
        .spawn().expect("Failed to start samtools.");
        
    let mut out = samtools.stdin.take().unwrap();
    
    loop {
        let rec_type = match reader.read_u8() {
            Ok(b) => b,
            Err(_) => break,
        };

        if rec_type == 0x00 {
            let l = reader.read_u32::<BigEndian>().unwrap() as usize;
            if l == 0 { continue; }
            let mut buf = vec![0u8; l];
            reader.read_exact(&mut buf).unwrap();
            out.write_all(&buf).unwrap();
            out.write_all(b"\n").unwrap();
        } else if rec_type == 0x10 {
            let n = reader.read_u32::<BigEndian>().unwrap() as usize;

            let mut meta_list = Vec::with_capacity(n);
            let mut dna_list = Vec::with_capacity(n);
            let mut qual_list = Vec::with_capacity(n);
            let mut tag_list = Vec::with_capacity(n);

            for _ in 0..n {
                let _ = reader.read_u8().unwrap(); // Assert 1
                let ml = reader.read_u32::<BigEndian>().unwrap() as usize;
                let sl = reader.read_u32::<BigEndian>().unwrap() as usize;
                let mut m_buf = vec![0u8; ml];
                reader.read_exact(&mut m_buf).unwrap();
                meta_list.push((String::from_utf8(m_buf).unwrap(), sl));
            }

            for _ in 0..n {
                let _ = reader.read_u8().unwrap(); // Assert 2
                let dl = reader.read_u32::<BigEndian>().unwrap() as usize;
                let mut d_buf = vec![0u8; dl];
                reader.read_exact(&mut d_buf).unwrap();
                let el = reader.read_u32::<BigEndian>().unwrap() as usize;
                let mut e_buf = vec![0u8; el];
                reader.read_exact(&mut e_buf).unwrap();
                dna_list.push((d_buf, e_buf));
            }

            for _ in 0..n {
                let _ = reader.read_u8().unwrap(); // Assert 3
                let ql = reader.read_u16::<BigEndian>().unwrap() as usize;
                let mut q_buf = vec![0u8; ql];
                reader.read_exact(&mut q_buf).unwrap();
                qual_list.push(q_buf);
            }

            for _ in 0..n {
                let _ = reader.read_u8().unwrap(); // Assert 4
                let tl = reader.read_u32::<BigEndian>().unwrap() as usize;
                let mut t_buf = vec![0u8; tl];
                reader.read_exact(&mut t_buf).unwrap();
                let ol = reader.read_u16::<BigEndian>().unwrap() as usize;
                let mut o_buf = vec![0u8; ol];
                reader.read_exact(&mut o_buf).unwrap();
                tag_list.push((t_buf, o_buf));
            }

            for i in 0..n {
                let (meta, sl) = &meta_list[i];
                let qname = meta.split('\t').next().unwrap();
                let dna = unpack_dna_k4(&dna_list[i].0, *sl, &dna_list[i].1);
                let qual = decode_qual(&qual_list[i], *sl);
                let tags = decode_tags(&tag_list[i].0, &tag_list[i].1, qname);

                let mut line = format!("{}\t{}\t{}", meta, dna, qual);
                if !tags.is_empty() {
                    line.push('\t');
                    line.push_str(&tags.join("\t"));
                }
                line.push('\n');
                out.write_all(line.as_bytes()).unwrap();
            }
        } else if rec_type == 0xFF {
            break;
        } else {
            eprintln!("Unknown record type 0x{:02x}", rec_type);
            break;
        }
    }
    drop(out); // Close stdin to signal samtools we are done
    let _ = samtools.wait();
    println!("Successfully decompressed to {}", output_bam);
}

// ── ENTRY POINT ───────────────────────────────────────────────────────────────
fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: ./kineticsquash <subreads.bam | subreads.sq.xz>");
        std::process::exit(1);
    }
    
    let input = &args[1];
    if input.ends_with(".bam") {
        compress(input);
    } else if input.ends_with(".sq.xz") {
        decompress(input);
    } else {
        eprintln!("Error: File must end in .bam (to compress) or .sq.xz (to decompress).");
    }
}
