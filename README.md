### Filter and orient cgMLST allele FASTAs
```
usage: clean_schema_with_reverse_complement_skip_empty_loci.py [-h] --scheme SCHEME --out OUT [--max_workers MAX_WORKERS]

options:

  -h, --help            show this help message and exit
  --scheme SCHEME       Directory with locus FASTA files
  --out OUT             Directory to write cleaned FASTA files
  --max_workers MAX_WORKERS
                        Threads for parallel processing

```



### Compute per-locus allele length thresholds and length counts.
```
usage: allele_length_thresholds_v2.py [-h] --schema SCHEMA [--epsilon EPSILON] [--min-alleles MIN_ALLELES] [--global-hard-bound GLOBAL_HARD_BOUND]
                                      [--extension EXTENSION] [--strict-extension] [--out-prefix OUT_PREFIX]

options:

  -h, --help            show this help message and exit  
  --schema SCHEMA       Directory containing FASTA files (one per locus).
  --epsilon EPSILON     Cushion added to [P1, P99] ratio bounds (default: 0.03)
  --min-alleles MIN_ALLELES
                        Minimum alleles required to compute thresholds (default: 5)                        
  --global-hard-bound GLOBAL_HARD_BOUND
                        Hard bound on |ratio-1|; e.g. 0.6 => [0.4, 1.6]                        
  --extension EXTENSION
                        FASTA file extension(s) to include (default: .fasta). Can be used multiple times.
  --strict-extension    If set, only include files ending exactly with provided extensions.
  --out-prefix OUT_PREFIX
                        Prefix for output CSV files (default: allele_lengths)
```
