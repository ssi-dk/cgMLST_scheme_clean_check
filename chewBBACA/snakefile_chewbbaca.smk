import os
import glob

FOLDER=config["results_path"]
scheme_dir = config["scheme"]

results_dir = os.path.join(FOLDER, "chewBBACA")
blast_out_dir = os.path.join(FOLDER, "blast_gene_call")

rule all:
    input:
        results_dir


rule run_chewbbaca:
    input:
        fasta_dir=blast_out_dir
    output:
        directory(results_dir)
    params:
        scheme=scheme_dir,
        part=config["partition"]
    threads:
        12
    resources:
        mem_mb=10000,
        slurm_partition=config["partition"],
        slurm_mem_mb=10000,
        slurm_threads=12
    conda:
        'chewie_test_raah'
    shell:
        """
        chewBBACA.py AlleleCall -i {input} -g {params.scheme} -o {output} --cpu {threads} --cds
        """
