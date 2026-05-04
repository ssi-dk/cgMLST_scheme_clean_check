#!/usr/bin/env python
import os
import subprocess
import argparse
import sys
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
from Bio.Seq import Seq
from pyfaidx import Fasta
import multiprocessing

import os

LOCAL_CWD = "/tmp/chewbbaca_runs"  # local xfs
os.makedirs(LOCAL_CWD, exist_ok=True)


# Define valid bases and codon sets (uppercase)
VALID_BASES = set(('A', 'T', 'C', 'G'))
START_CODONS = {"ATG", "TTG", "GTG", "CTG", "ATA", "ATT"}
STOP_CODONS  = {"TAA", "TAG", "TGA"}

# Configure logging to file (INFO and above) and preserve warnings/info
logging.basicConfig(
    filename='extraction.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s'
)
console = logging.StreamHandler(sys.stderr)
console.setLevel(logging.INFO)
console.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
logging.getLogger().addHandler(console)



def is_valid_cds(seq: str) -> bool:
    """
    Returns True if the sequence starts with a valid start codon,
    ends with a valid stop codon, and is a multiple of 3 in length.
    """
    return (
        seq[:3] in START_CODONS and
        seq[-3:] in STOP_CODONS and
        len(seq) % 3 == 0
    )


def has_internal_stop(seq: str) -> bool:
    """
    Returns True if there is any internal in-frame stop codon (excluding first and last codons).
    """
    for i in range(3, len(seq) - 3, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            return True
    return False

def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    return str(Seq(seq).reverse_complement())



def orient_and_frame_fix(seq: str, qstart: int, qend: int, genome: Fasta, seq_id: str, header_raw: str, strand: int, max_stop_extend: int):
    """
    Orient (reverse-complement if needed) and ensure seq length % 3 == 0
    by codon-aware padding or trimming. Returns fixed_seq, new_start, new_end, stats.
    """
    stats = Counter()

    # 1) Orientation from BLAST strand
    if strand == -1:
        seq = reverse_complement(seq)
        stats['reoriented'] += 1

    # Single uppercase pass
    seq = seq.upper()

    # 2) Frame correction
    length = len(seq)
    mod = length % 3
    # Frame correction (only the index moves change)

    if mod != 0:
        needed = 3 - mod
        adjusted = False
        contig_len = len(genome[seq_id])

        # a) Try extending 3' end (codon boundary)
        if strand == 1:
            if qend + needed <= contig_len:
                cand = genome[seq_id][qstart:qend + needed].seq.upper()
                if len(cand) % 3 == 0:
                    seq = cand
                    qend += needed
                    stats['padded_stop'] += 1
                    adjusted = True
        else:
            if qstart - needed >= 0:
                cand = genome[seq_id][qstart - needed:qend].seq.upper()
                if len(cand) % 3 == 0:
                    # sequence stays RC relative to genome coords
                    seq = reverse_complement(cand)
                    qstart -= needed
                    stats['padded_stop'] += 1
                    adjusted = True

        # b) Try extending 5' end if not adjusted
        if not adjusted:
            if strand == 1 and qstart - needed >= 0:
                cand = genome[seq_id][qstart - needed:qend].seq.upper()
                if len(cand) % 3 == 0:
                    seq = cand
                    qstart -= needed
                    stats['padded_start'] += 1
                    adjusted = True
            elif strand == -1 and qend + needed <= contig_len:
                cand = genome[seq_id][qstart:qend + needed].seq.upper()
                if len(cand) % 3 == 0:
                    seq = reverse_complement(cand)
                    qend += needed
                    stats['padded_start'] += 1
                    adjusted = True

    # After padding attempts, require frame to be restored
    if len(seq) % 3 != 0:
        logging.warning(f"{header_raw}: out of frame after padding; skipping.")
        return "", qstart, qend, stats  # signal caller to discard


    # After frame fix, if no terminal stop, walk in codon steps up to 90 bp
    if seq[-3:] not in STOP_CODONS:
        max_extend = max_stop_extend
        contig_len = len(genome[seq_id])
        extended = False
        for i in range(3, max_extend + 1, 3):
            if strand == 1:
                if qend + i > contig_len: break
                cand = genome[seq_id][qstart : qend + i].seq.upper()
            else:
                if qstart - i < 0: break
                raw = genome[seq_id][qstart - i : qend].seq.upper()
                cand = reverse_complement(raw)
            if cand[-3:] in STOP_CODONS and len(cand) % 3 == 0:
                seq = cand
                if strand == 1: qend += i
                else:           qstart -= i
                stats['stop_extended'] = i
                extended = True
                break
        if not extended:
            msg = f"{header_raw}: no terminal stop within {max_extend} bp 3' window; skipping."
            print(f"[WARNING] {msg}", file=sys.stderr); logging.warning(msg)
            return "", qstart, qend, stats  # cause caller to skip

    return seq, qstart, qend, stats


def extract_subsequences(genome: Fasta, alleles: list, max_stop_extend: int) -> list:
    """
    Iterate over allele regions, perform validation, orientation, framing,
    and return valid (header, seq) pairs.
    """
    extracted = []
    for allele in alleles:
        seq_id, qstart, qend, saccver, strand = allele
        header_raw = f">{seq_id}_{saccver}_{qstart+1}_{qend}"

        # 0) Extract raw subsequence
        subseq = genome[seq_id][qstart:qend].seq

        # After extracting subseq
        if len(subseq) != (qend - qstart):
            logging.error(f"{header_raw}: slice length {len(subseq)} != qend-qstart ({qend - qstart})")
            continue

        # 1) Check for ambiguous bases and report details
        invalid_positions = [(i, b) for i, b in enumerate(subseq) if b.upper() not in VALID_BASES]
        if invalid_positions:
            for pos, base in invalid_positions[:20]:
                msg = f"{header_raw}: invalid base '{base}' at position {pos+1}"
                print(f"[WARNING] {msg}", file=sys.stderr); logging.warning(msg)
            if len(invalid_positions) > 20:
                logging.warning(f"{header_raw}: ... and {len(invalid_positions)-20} more invalid base(s)")
            continue


        # 2) Orient and frame-fix
        fixed, new_start, new_end, stats = orient_and_frame_fix(
            subseq, qstart, qend, genome, seq_id, header_raw, strand, max_stop_extend)

        if not fixed:
            # already logged inside orient_and_frame_fix
            continue
        
        # Check for internal in-frame stop codons
        if has_internal_stop(fixed):
            msg = f"{header_raw}: contains internal in-frame stop codon, likely truncated or pseudogene."
            print(f"[WARNING] {msg}", file=sys.stderr)
            logging.warning(msg)
            continue  # Skip this sequence

        # 4) Final validation and append if valid CDS
        if is_valid_cds(fixed):
            final_hdr = f">{seq_id}_{saccver}_{new_start+1}_{new_end}"
            extracted.append((final_hdr, fixed))
        else:
            msg = f"{header_raw}: skipped due to invalid CDS (start/stop missing or length not divisible by 3)."
            print(f"[WARNING] {msg}", file=sys.stderr)
            logging.warning(msg)
    return extracted


def parse_blast_output(blast_output: str, genome: Fasta, min_cov_ratio: float) -> list:
    """
    Parse BLAST tabular output and select perfect or best fallback hits.
    Returns alleles as list of tuples: (qaccver, qstart, qend, saccver).
    """
    perfect_hits = []
    fallback_strict = None
    fallback_strict_key = None
    fallback_any = None
    fallback_any_key = None
    with open(blast_output) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 13:
                continue
            rec = {
                'qaccver': cols[0],
                'saccver': cols[1],
                'slen': int(cols[2]),
                'pident': float(cols[3]),
                'length': int(cols[4]),
                'mismatch': int(cols[5]),
                'gapopen': int(cols[6]),
                'qstart': int(cols[7]),
                'qend': int(cols[8]),
                'sstart': int(cols[9]),
                'send': int(cols[10]),
                'evalue': float(cols[11]),
                'bitscore': float(cols[12])
            }
            # if rec['pident'] == 100.0 and rec['length'] == rec['slen']:
            #     perfect_hits.append(rec)
            # else:
            #     if (fallback is None or
            #         rec['length'] > fallback['length'] or
            #         (rec['length'] == fallback['length'] and rec['pident'] > fallback['pident']) or
            #         (rec['length'] == fallback['length'] and rec['pident'] == fallback['pident'] and rec['bitscore'] > fallback['bitscore'])):
            #         fallback = rec
    

            qstart, qend = rec['qstart'], rec['qend']
            query_span = abs(qend - qstart) + 1
            no_gaps = (rec['gapopen'] == 0)

            # PERFECT: clean, full-length, no indels
            if rec['pident'] == 100.0 and no_gaps and query_span == rec['slen']:
                perfect_hits.append(rec)
            else:
                # your policy: consider both "no_gaps" and "length" first, then length, then %id, then bitscore
                rank_key = (no_gaps, query_span, rec['pident'], rec['bitscore'])

                # Track best overall (even if short)
                if fallback_any_key is None or rank_key > fallback_any_key:
                    fallback_any = rec
                    fallback_any_key = rank_key

                # Track best that meets coverage floor
                if query_span >= min_cov_ratio * rec['slen']:
                    if fallback_strict_key is None or rank_key > fallback_strict_key:
                        fallback_strict = rec
                        fallback_strict_key = rank_key

    # If multiple perfect hits share the same qstart or the same qend on the same contig, keep only the longest (then highest bitscore).
    if len(perfect_hits) > 1:
        def span(r): return abs(r['qend'] - r['qstart']) + 1
        # longest first; bitscore as tiebreak (pident is 100% already for perfect hits)
        perfect_hits.sort(key=lambda r: (span(r), r['bitscore']), reverse=True)

        kept = []
        seen_start = set()   # (qaccver, qstart)
        seen_end   = set()   # (qaccver, qend)
        for r in perfect_hits:
            kS = (r['qaccver'], r['saccver'], r['qstart'])
            kE = (r['qaccver'], r['saccver'], r['qend'])
            if kS in seen_start or kE in seen_end:
                continue  # a longer, same-anchored hit already kept
            kept.append(r)
            seen_start.add(kS)
            seen_end.add(kE)
        perfect_hits = kept

    #hits = perfect_hits if perfect_hits else ([fallback] if fallback else [])
    hits = perfect_hits if perfect_hits else (
    [fallback_strict] if fallback_strict else ([fallback_any] if fallback_any else [])
)

    alleles = []
    for hit in hits:
        sstart, send = hit['sstart'], hit['send']
        qstart, qend = hit['qstart'], hit['qend']
        slen = hit['slen']
        # Normalize coordinates
        start, end = min(sstart, send), max(sstart, send)
        strand = 1 if sstart <= send else -1

        # Expand query coords so the subject spans 1..slen
        if start != 1:
            qstart -= (start - 1)   # expand left (was + before; that shrank!)
        if end != slen:
            qend += (slen - end)    # expand right

        # Clamp to contig and convert to 0-based start / 1-based inclusive end
        qstart0 = max(0, qstart - 1)
        qend_incl = min(qend, len(genome[hit['qaccver']]))

        alleles.append((hit['qaccver'], qstart0, qend_incl, hit['saccver'], strand))
    return alleles


def run_blast_locus(assembly_path: str, locus_path: str, genome: Fasta, assembly_name: str, min_cov_ratio: float) -> list:
    """
    Run BLAST for a single locus against an assembly and return allele regions.
    """
    locus_name = os.path.basename(locus_path).replace('.fasta', '')
    #blast_output = f"blast_{assembly_name}_{locus_name}.txt" 
    blast_output = os.path.join(LOCAL_CWD, f"blast_{assembly_name}_{locus_name}.txt")

    cmd = [
        'blastn', '-query', assembly_path,
        '-subject', locus_path,
        '-out', blast_output,
        '-outfmt', '6 qaccver saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        '-max_target_seqs', '1000',
        '-num_threads', '1'
    ]
    subprocess.run(cmd, check=True)
    alleles = []
    if os.path.exists(blast_output):
        alleles = parse_blast_output(blast_output, genome, min_cov_ratio)
        os.remove(blast_output)
    return alleles


def process_single_assembly(assembly_file: str, scheme_dir: str, output_dir: str, assemblies_dir: str, max_workers: int, max_stop_extend: int, min_cov_ratio: float):
    assembly_path = os.path.join(assemblies_dir, assembly_file)
    if not os.path.exists(assembly_path):
        raise FileNotFoundError(f"Assembly file not found: {assembly_path}")

    output_fasta = os.path.join(output_dir, os.path.splitext(os.path.basename(assembly_file))[0] + ".fa")

    genome = Fasta(assembly_path)

    loci = [f for f in os.listdir(scheme_dir) if f.endswith('.fasta')]

    all_alleles = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                run_blast_locus,
                assembly_path,
                os.path.join(scheme_dir, locus),
                genome,
                os.path.splitext(os.path.basename(assembly_file))[0],
                min_cov_ratio
            ): locus for locus in loci
        }
        for future in as_completed(futures):
            all_alleles.extend(future.result())

    extracted = extract_subsequences(genome, all_alleles, max_stop_extend)

    with open(output_fasta, 'w') as out_f:
        for hdr, seq in extracted:
            out_f.write(f"{hdr}\n{seq}\n")

    print(f"Wrote {len(extracted)} alleles to {output_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description="Run per-locus BLAST in parallel and extract valid CDS subsequences."
    )
    parser.add_argument("--scheme", required=True, help="Directory containing locus FASTA files")
    parser.add_argument("--fa", required=True, help="Directory or file of assembly FASTA sequences")
    parser.add_argument("--out", required=True, help="Directory to store the output files")
    parser.add_argument("--max_workers", type=int, default=multiprocessing.cpu_count(), help="Max concurrent BLASTs")
    parser.add_argument("--samples_file", help="Path to a text file listing one SAMPLE ID per line")
    parser.add_argument("--max_stop_extend", type=int, default=90,
                    help="Max 3' in-frame extension (bp) to find a stop [default: 90]")
    parser.add_argument("--min_cov", type=float, default=0.70,
                    help="Min ungapped query coverage vs slen to consider 'strict' fallback [default: 0.70]")


    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)


    args.fa = os.path.abspath(args.fa)
    args.out = os.path.abspath(args.out)
    if args.samples_file:
        args.samples_file = os.path.abspath(args.samples_file)

    # Determine assemblies list
    valid_ext = ('.fa', '.fna', '.fas', '.fasta')

    if os.path.isdir(args.fa):
        assembly_paths = [
            os.path.join(args.fa, f)
            for f in os.listdir(args.fa)
            if f.lower().endswith(valid_ext)
        ]
    elif os.path.isfile(args.fa):
        assembly_paths = [args.fa]
    else:
        raise ValueError(f"{args.fa!r} is not a file or directory")

    for assembly_path in assembly_paths:
        assembly_name = os.path.splitext(os.path.basename(assembly_path))[0]
        print(f"Processing assembly: {assembly_name}")
        process_single_assembly(
            os.path.basename(assembly_path),
            args.scheme,
            args.out,
            os.path.dirname(assembly_path),
            args.max_workers,
            args.max_stop_extend,
            args.min_cov
        )

if __name__ == "__main__":
    main()
