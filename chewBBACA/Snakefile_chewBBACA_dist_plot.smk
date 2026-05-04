##env FBI_cgMLST

import os
import glob

## --------------------------------------------------------------------------------
## global parameters from config file

shell.prefix("set -euo pipefail; ulimit -s 131072; ") ##possible to plot many samples
FOLDER=config["results_path"]


def get_latest_results_folder():
    """Determine the folder containing the latest results_alleles.tsv file."""
    base_path = os.path.join(FOLDER, "chewBBACA")
    
    # Find all timestamped folders and ensure they're directories
    results_dirs = glob.glob(f"{base_path}/results_20*T*")
    results_dirs = [path for path in results_dirs if os.path.isdir(path)]
    
    # Get the newest folder based on creation time (if any)
    newest_folder = max(results_dirs, key=os.path.getctime, default=None)
    newest_folder_time = os.path.getctime(newest_folder) if newest_folder else 0

    # Path to the results_alleles.tsv file in the base directory
    results_file_path = os.path.join(base_path, "results_alleles.tsv")

    # If the file exists and is newer than the newest folder (or no folders exist), return its folder
    if os.path.isfile(results_file_path) and os.path.getmtime(results_file_path) > newest_folder_time:
        return base_path

    # If a newest folder exists, return its path
    if newest_folder:
        return newest_folder

    # If nothing is found, raise an error
    raise ValueError("No valid file or folder or file found in the given path.")

# Compute the folder containing the latest results_alleles.tsv file
latest_folder = get_latest_results_folder()


## --------------------------------------------------------------------------------
## output file sets


plots = f"{latest_folder}/dendrogram.pdf"


## --------------------------------------------------------------------------------
## rules


rule all:
    input:
        plots

rule dist:
    input:
        f"{latest_folder}/results_alleles.tsv"
    output:
        f"{latest_folder}/dist.mat.tsv"
    threads:
        1
    resources:
        mem_mb=1000,
        slurm_partition=config["partition"],
        slurm_mem_mb=1000,
        slurm_threads=1

    shell:
        """
        cgmlst-dists {input} > {output}
        """

rule plot:
    input:
        f"{latest_folder}/dist.mat.tsv"
    output:
        f"{latest_folder}/dendrogram.pdf"
    threads:
        1
    resources:
        mem_mb=1000,
        slurm_partition=config["partition"],
        slurm_mem_mb=1000,
        slurm_threads=1
    shell:
        """
        automatic_chewbbaca_dendrograms.R --out {output} --dist {input}
        """


