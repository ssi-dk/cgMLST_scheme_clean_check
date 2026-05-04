import os

# Load config variables
sample_map = config["samples"] 
sample_ids = list(sample_map.keys())

results_path = config["results_path"]
fasta_path = config["fasta"]
scheme_path = config["scheme"]
#blastdb_path = config["blastdb"]
blast_out = os.path.join(results_path, "blast_gene_call")

# Handle empty sample list gracefully
if not sample_ids:
    print("No new samples for BLAST gene calling. Skipping this step.")
    sample_ids = []

# Rule: final targets
rule all:
    input:
        expand(os.path.join(blast_out, "{sample}.fa"), sample=sample_ids)

# Rule: run gene calling per sample
rule blast_gene_call:
    input:
        fasta=lambda wc: os.path.join(fasta_path, sample_map[wc.sample])
    output:
        result=os.path.join(blast_out, "{sample}.fa")
    params:
        scheme=scheme_path,
        out=blast_out
    threads: 
        6
    resources:
        mem_mb=6000,
        slurm_partition=config["partition"],
        slurm_mem_mb=6000,
        slurm_threads=6
    shell:
        """
        python /users/data/Projects/FBI_SOFI/proj/chewBBACA/Tools/scripts/new_per_locus/new/blast_gene_call_per_locus_moreChecks_optimized2_skipNonCDS.py --scheme {params.scheme} --fa {input.fasta} --out {params.out}
        """
