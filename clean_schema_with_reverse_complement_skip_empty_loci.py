#!/usr/bin/env python3


"""
Filter cgMLST scheme FASTA files to keep only "proper" coding-sequence alleles.

Given a directory with thousands of locus FASTA files (each containing many allele
sequences), this script validates each allele using strict CDS checks and writes
only the valid ones to an output directory (one cleaned FASTA per input FASTA).

What counts as a proper allele here (strict mode):
  1) Only unambiguous DNA bases {A,C,G,T}
  2) Length is a multiple of 3 (codon complete)
  3) Starts with a valid bacterial start codon (ATG, TTG, GTG, CTG, ATA, ATT)
  4) Ends with a valid stop codon (TAA, TAG, TGA)
  5) No internal in-frame stop codons

Example:
  python clean_schema.py \
      --scheme path/to/cgmlst_scheme \
      --out path/to/clean_scheme \

"""

import sys
import os
import logging
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
from typing import Tuple, Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# === Global Constants ===
VALID_BASES = set("ACGT")
START_CODONS = {"ATG", "TTG", "GTG", "CTG", "ATA", "ATT"}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def setup_logging(log_name: str):
    """Configure logging to both stderr and a file named after the scheme."""
    log_file = Path.cwd() / f"{log_name}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler(log_file, mode="w")
        ]
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


def open_maybe_gzip(path: Path, mode: str = "rt"):
    """Open a file in normal or gzip mode depending on the file extension."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def is_unambiguous_dna(seq: str) -> bool:
    """Check if a DNA sequence contains only A, C, G, or T."""
    return all(ch in VALID_BASES for ch in seq.upper())


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement()).upper()


def has_internal_stop(seq: str) -> bool:
    """Check whether a DNA sequence contains internal stop codons."""
    for i in range(3, len(seq) - 3, 3):
        if seq[i:i+3] in STOP_CODONS:
            return True
    return False


def check_cds(seq: str) -> Tuple[bool, str]:
    """Check CDS criteria on a single orientation."""
    s = seq.upper()
    if len(s) % 3 != 0:
        return False, "length_not_multiple_of_3"
    if s[:3] not in START_CODONS:
        return False, "missing_start_codon"
    if s[-3:] not in STOP_CODONS:
        return False, "missing_stop_codon"
    if has_internal_stop(s):
        return False, "internal_stop"
    return True, ""


def validate_allele(seq: str) -> Tuple[bool, str, str, bool]:
    """
    Validate an allele sequence.
    Strategy:
      1) Require unambiguous DNA.
      2) Try full CDS validation on forward strand.
      3) If that fails, try full CDS validation on reverse complement.
      4) Accept whichever passes; otherwise report forward failure reason.
    Returns (is_valid, reason, corrected_seq, reoriented)
    """
    s = seq.upper()
    reoriented = False

    if not s:
        return False, "empty", s, reoriented
    if not is_unambiguous_dna(s):
        return False, "ambiguous_base", s, reoriented

    ok_fwd, reason_fwd = check_cds(s)
    if ok_fwd:
        return True, "", s, reoriented

    rc = reverse_complement(s)
    ok_rc, reason_rc = check_cds(rc)
    if ok_rc:
        return True, "", rc, True

    return False, reason_fwd, s, reoriented


def process_fasta_file(in_path: Path, out_dir: Path) -> Tuple[str, Counter, Dict[str, int]]:
    """
    Validate all records in a single FASTA file (per locus).
    Writes only valid alleles to the output FASTA file.
    If no valid alleles are found, no output file is written.
    """
    stats = Counter()
    invalid_reasons = defaultdict(int)

    out_name = in_path.name
    if out_name.endswith(".gz"):
        out_name = out_name[:-3]
    out_path = out_dir / out_name

    valid_records = []

    with open_maybe_gzip(in_path, "rt") as fh_in:
        for rec in SeqIO.parse(fh_in, "fasta"):
            seq_str = str(rec.seq).upper()
            ok, reason, fixed, reoriented = validate_allele(seq_str)
            if ok:
                stats["valid"] += 1
                if reoriented:
                    stats["reoriented"] += 1
                rec.seq = Seq(fixed)
                valid_records.append(rec)
            else:
                stats["invalid"] += 1
                stats[f"invalid_{reason}"] += 1
                invalid_reasons[reason] += 1

    if stats["valid"] == 0:
        stats["empty_locus"] = 1
        return in_path.name, stats, dict(invalid_reasons)

    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh_out:
        SeqIO.write(valid_records, fh_out, "fasta")

    return in_path.name, stats, dict(invalid_reasons)


def main():
    """Main entry point."""
    import argparse
    ap = argparse.ArgumentParser(description="Filter and orient cgMLST allele FASTAs")
    ap.add_argument("--scheme", required=True, help="Directory with locus FASTA files")
    ap.add_argument("--out", required=True, help="Directory to write cleaned FASTA files")
    ap.add_argument("--max_workers", type=int, default=None, help="Threads for parallel processing")
    args = ap.parse_args()

    scheme_dir = Path(args.scheme).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    scheme_name = scheme_dir.name
    setup_logging(scheme_name)

    fasta_paths = sorted(
        p for p in scheme_dir.iterdir()
        if p.suffix in (".fasta", ".fa") or p.name.endswith((".fasta.gz", ".fa.gz"))
    )
    if not fasta_paths:
        logging.error(f"No FASTA files found in {scheme_dir}")
        sys.exit(1)

    workers = args.max_workers or int(
        (os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1)
    )

    global_stats = Counter()
    empty_loci = []
    per_locus_counts = {}
    all_invalid_reasons = set()

    with ThreadPoolExecutor(max_workers=workers) as ex:
        fut2path = {ex.submit(process_fasta_file, p, out_dir): p for p in fasta_paths}
        for fut in as_completed(fut2path):
            p = fut2path[fut]
            try:
                locus_name, stats, invalid_reason_counts = fut.result()
                global_stats.update(stats)
                per_locus_counts[locus_name] = {
                    "valid": stats.get("valid", 0),
                    "invalid": stats.get("invalid", 0),
                    "reoriented": stats.get("reoriented", 0),
                    **{f"invalid_{k}": v for k, v in invalid_reason_counts.items()}
                }
                all_invalid_reasons.update(invalid_reason_counts.keys())
                if stats.get("empty_locus"):
                    empty_loci.append(locus_name)
            except Exception as e:
                logging.exception(f"Error processing {p}: {e}")
                global_stats["errors"] += 1

    logging.info(f"Processed {len(fasta_paths)} FASTA files")
    logging.info(f"Valid alleles total: {global_stats['valid']}")
    logging.info(f"Invalid alleles total: {global_stats['invalid']}")
    if global_stats.get("reoriented"):
        logging.info(f"Alleles auto-oriented (reverse-complemented): {global_stats['reoriented']}")
    if global_stats.get("errors"):
        logging.info(f"Files with errors: {global_stats['errors']}")

    for k, v in sorted(global_stats.items()):
        if k.startswith("invalid_"):
            logging.info(f"  {k}: {v}")

    if empty_loci:
        empty_file = Path.cwd() / f"empty_loci_{scheme_name}.txt"
        logging.warning(f"{len(empty_loci)} loci had NO valid alleles and were skipped.")
        with open(empty_file, "w") as f:
            for locus in sorted(empty_loci):
                f.write(f"{locus}\n")
        logging.warning(f"Empty loci list saved to {empty_file}")
    else:
        logging.info("No empty loci detected 🎉")

    invalid_columns = [f"invalid_{reason}" for reason in sorted(all_invalid_reasons)]
    counts_file = Path.cwd() / f"locus_counts_{scheme_name}.tsv"
    with open(counts_file, "w") as f:
        f.write("locus\tvalid\tinvalid\treoriented")
        for col in invalid_columns:
            f.write(f"\t{col}")
        f.write("\n")

        for locus in sorted(per_locus_counts.keys()):
            row = per_locus_counts[locus]
            f.write(f"{locus}\t{row['valid']}\t{row['invalid']}\t{row['reoriented']}")
            for col in invalid_columns:
                f.write(f"\t{row.get(col, 0)}")
            f.write("\n")
    logging.info(f"Per-locus counts (with invalid breakdown) written to {counts_file}")

    print(
        f"VALID={global_stats['valid']} INVALID={global_stats['invalid']} "
        f"REORIENTED={global_stats.get('reoriented',0)} EMPTY_LOCI={len(empty_loci)} "
        f"ERRORS={global_stats.get('errors',0)}",
        file=sys.stdout
    )


if __name__ == "__main__":
    main()
