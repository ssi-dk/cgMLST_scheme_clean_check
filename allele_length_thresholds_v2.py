#!/usr/bin/env python3
"""
allele_length_thresholds_v2.py

Compute per-locus allele length distributions and suggest length ratio thresholds for cgMLST/chewBBACA QC,
and also produce a per-locus length frequency table.

CHANGES vs v1
- Writes TWO CSVs:
    1) <prefix>_thresholds.csv      (summary statistics and suggested bounds per locus)
    2) <prefix>_length_counts.csv   (rows = locus, length, count)
- New option --out-prefix to control filenames (default: "allele_lengths")

USAGE
    python allele_length_thresholds_v2.py --schema /path/to/dir --out-prefix results/allele_lengths

EXAMPLE OUTPUT FILES
    results/allele_lengths_thresholds.csv
    results/allele_lengths_length_counts.csv
"""
from __future__ import annotations
import argparse
import os
import sys
import math
from collections import Counter, defaultdict
from statistics import median, mean

def parse_args():
    p = argparse.ArgumentParser(description="Compute per-locus allele length thresholds and length counts.")
    p.add_argument("--schema", required=True, help="Directory containing FASTA files (one per locus).")
    p.add_argument("--epsilon", type=float, default=0.03, help="Cushion added to [P1, P99] ratio bounds (default: 0.03)")
    p.add_argument("--min-alleles", type=int, default=5, help="Minimum alleles required to compute thresholds (default: 5)")
    p.add_argument("--global-hard-bound", type=float, default=None, help="Hard bound on |ratio-1|; e.g. 0.6 => [0.4, 1.6]")
    p.add_argument("--extension", action="append", default=[".fasta"], help="FASTA file extension(s) to include (default: .fasta). Can be used multiple times.")
    p.add_argument("--strict-extension", action="store_true", help="If set, only include files ending exactly with provided extensions.")
    p.add_argument("--out-prefix", default="allele_lengths", help="Prefix for output CSV files (default: allele_lengths)")
    return p.parse_args()

def iter_fasta_lengths(path):
    """Return list of lengths for sequences in FASTA file at path."""
    lens = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        seq_len = 0
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if seq_len > 0:
                    lens.append(seq_len)
                    seq_len = 0
            else:
                s = line.strip()
                if s:
                    seq_len += len(s)
        if seq_len > 0:
            lens.append(seq_len)
    return lens

def percentile(xs, p):
    """Compute the p-th percentile (0-100) via linear interpolation (numpy-like)."""
    if not xs:
        return math.nan
    xs_sorted = sorted(xs)
    k = (len(xs_sorted)-1) * (p/100.0)
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return float(xs_sorted[int(k)])
    d0 = xs_sorted[f] * (c - k)
    d1 = xs_sorted[c] * (k - f)
    return float(d0 + d1)

def mad(xs, med=None):
    """Median Absolute Deviation (raw, not scaled)."""
    if not xs:
        return math.nan
    if med is None:
        med = median(xs)
    return median([abs(x - med) for x in xs])

def suggest_bounds(ratio_p1, ratio_p99, epsilon, hard=None):
    lo = max(ratio_p1 - epsilon, 0.0)
    hi = ratio_p99 + epsilon
    if hard is not None:
        lo = max(lo, 1.0 - hard)
        hi = min(hi, 1.0 + hard)
        if lo > hi:
            lo, hi = 1.0 - hard, 1.0 + hard
    return lo, hi

def main():
    args = parse_args()
    schema = args.schema
    if not os.path.isdir(schema):
        print(f"ERROR: schema path is not a directory: {schema}", file=sys.stderr)
        sys.exit(2)

    # Collect relevant FASTA files
    fasta_files = []
    for root, _, files in os.walk(schema):
        for fn in files:
            fpath = os.path.join(root, fn)
            if args.strict_extension:
                if any(fn.endswith(ext) for ext in args.extension):
                    fasta_files.append(fpath)
            else:
                if any(fn.lower().endswith(ext.lower()) for ext in args.extension) or fn.lower().endswith((".fa",".fna",".fas",".fasta")):
                    fasta_files.append(fpath)

    # Prepare output paths
    out_thresh = f"{args.out_prefix}_thresholds.csv"
    out_counts = f"{args.out_prefix}_length_counts.csv"

    # Ensure output directories exist
    out_dir = os.path.dirname(os.path.abspath(out_thresh))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Open writers
    th = open(out_thresh, "w", encoding="utf-8")
    lc = open(out_counts, "w", encoding="utf-8")

    # Write headers
    th_header = [
        "locus","n_alleles",
        "min_bp","p1_bp","p5_bp","q1_bp","median_bp","mean_bp","q3_bp","p95_bp","p99_bp","max_bp",
        "iqr_bp","mad_bp","inframe_fraction", "p1p99diff", "p5p95diff",
        "p1mediandiff", "p5mediandiff", "p95mediandiff", "p99mediandiff", 
        "ratio_p1","ratio_p5","ratio_q1","ratio_q3","ratio_p95","ratio_p99",
        "suggest_ratio_lower","suggest_ratio_upper"
    ]
    th.write(",".join(th_header) + "\n")
    lc.write("locus,length_bp,count\n")

    # Process each file
    for fpath in sorted(fasta_files):
        base = os.path.basename(fpath)
        locus = base
        for ext in [".fasta",".fa",".fna",".fas"] + args.extension:
            if locus.lower().endswith(ext.lower()):
                locus = locus[: -len(ext)]
                break

        lens = iter_fasta_lengths(fpath)
        n = len(lens)
        if n == 0:
            continue

        # --- Length counts ---
        counts = Counter(lens)
        for L, cnt in sorted(counts.items()):
            lc.write(f"{locus},{int(L)},{cnt}\n")

        # --- Stats ---
        min_bp = min(lens)
        max_bp = max(lens)
        p1 = percentile(lens, 1)
        p5 = percentile(lens, 5)
        q1 = percentile(lens, 25)
        med = float(median(lens))
        mu = float(mean(lens))
        q3 = percentile(lens, 75)
        p95 = percentile(lens, 95)
        p99 = percentile(lens, 99)
        iqr = q3 - q1
        mad_val = mad(lens, med)
        inframe_frac = sum(1 for L in lens if (L % 3) == 0) / n
        p1p99diff=100-((p1/p99)*100)
        p5p95diff=100-((p5/p95)*100)
        p1mediandiff=100-((p1/med)*100)
        p5mediandiff=100-((p5/med)*100)
        p95mediandiff=100-((p95/med)*100)
        p99mediandiff=100-((p99/med)*100)

        # Ratios vs median
        def r(x):
            return (x / med) if med > 0 else float("nan")
        ratio_p1  = r(p1)
        ratio_p5  = r(p5)
        ratio_q1  = r(q1)
        ratio_q3  = r(q3)
        ratio_p95 = r(p95)
        ratio_p99 = r(p99)

        if n >= args.min_alleles:
            lo, hi = suggest_bounds(ratio_p1, ratio_p99, args.epsilon, args.global_hard_bound)
        else:
            lo, hi = (1.0 - 0.15, 1.0 + 0.15)
            if args.global_hard_bound is not None:
                lo = max(lo, 1.0 - args.global_hard_bound)
                hi = min(hi, 1.0 + args.global_hard_bound)

        th_row = [
            locus, str(n),
            f"{min_bp:.0f}", f"{p1:.2f}", f"{p5:.2f}", f"{q1:.2f}", f"{med:.2f}", f"{mu:.2f}", f"{q3:.2f}", f"{p95:.2f}", f"{p99:.2f}", f"{max_bp:.0f}",
            f"{iqr:.2f}", f"{mad_val:.2f}", f"{inframe_frac:.3f}", f"{p1p99diff:.2f}", f"{p5p95diff:.2f}",
            f"{p1mediandiff:.2f}", f"{p5mediandiff:.2f}", f"{p95mediandiff:.2f}", f"{p99mediandiff:.2f}",
            f"{ratio_p1:.3f}", f"{ratio_p5:.3f}", f"{ratio_q1:.3f}", f"{ratio_q3:.3f}", f"{ratio_p95:.3f}", f"{ratio_p99:.3f}",
            f"{lo:.3f}", f"{hi:.3f}"
        ]
        th.write(",".join(th_row) + "\n")

    th.close()
    lc.close()

    print(f"Wrote: {out_thresh}")
    print(f"Wrote: {out_counts}")

if __name__ == "__main__":
    main()
