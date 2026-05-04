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


### Run chewBBACA with blast genecalls and plot dendrogram with R

run chewbbaca_auto.py, 
- this script runs chewbbaca_incremental.py to check if any new assemblies have been added to the assembly folder, 
- then the blast gene call with snakefile_genecall.smk, which runs the actual gene calling script blast_gene_call_per_locus_moreChecks_optimized2_skipNonCDS.py,
- then runs the actual allelecall with snakefile_chewbbaca.smk,
- lastly plots the dendrogram with Snakefile_chewBBACA_dist_plot.smk which uses the R script automatic_chewbbaca_dendrograms.R

If you need to change the plot after the pipeline is done you can use the second R script change_automatic_chewbbaca_dendrograms.R


Remember to change all paths in the scripts. 

