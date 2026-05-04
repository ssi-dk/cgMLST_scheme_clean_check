#!/usr/bin/env python


import os
import argparse
import subprocess
import yaml
from pathlib import Path
import sys
import signal 



def get_unprocessed_assemblies(assemblies_dir, output_dir):
    all_fasta_files = [
        f for f in os.listdir(assemblies_dir)
        if f.endswith(('.fa', '.fasta', '.fas', '.fna'))
    ]
    processed = [os.path.splitext(f)[0] for f in os.listdir(output_dir) if f.endswith('.fa')]

    unprocessed = {}
    for fasta_file in all_fasta_files:
        sample_id = os.path.splitext(fasta_file)[0]  # full stem
        if sample_id not in processed:
            unprocessed[sample_id] = fasta_file
    return unprocessed

def inject_samples_into_config(config_path, sample_map):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    config['samples'] = sample_map
    with open(config_path, 'w') as f:
        yaml.dump(config, f)


def run_incremental(config_path, snakefile_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    assemblies_dir = config['fasta']
    blast_out_dir = os.path.join(config['results_path'], 'blast_gene_call')
    os.makedirs(blast_out_dir, exist_ok=True)

    unprocessed = get_unprocessed_assemblies(assemblies_dir, blast_out_dir)
    sample_names = [os.path.splitext(f)[0] for f in unprocessed]

    if not sample_names:
        print("No new samples to process.")
        return

    print(f"Found {len(sample_names)} new sample(s): {sample_names}")

    inject_samples_into_config(config_path, unprocessed)

    dryrun = f"snakemake --snakefile {snakefile_path} --configfile {config_path} --profile slurm -np"
    subprocess.run(dryrun, shell=True, executable="/bin/bash")

    while True:
        user_input = input("Did the dryrun run without any red errors? [yes/y or no/n]: ").strip().lower()
        if user_input in {"yes", "y"}:
            run_cmd = (
                f"snakemake --snakefile {snakefile_path} --configfile {config_path} --profile slurm")

            proc = None
            try:
                proc = subprocess.Popen(run_cmd, shell=True, executable="/bin/bash", preexec_fn=os.setsid)
                proc.communicate()
            except KeyboardInterrupt:
                print("\n[INTERRUPTED] Terminating all subprocesses.")
                if proc is not None:
                    try:
                        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                    except Exception as e:
                        print(f"[WARNING] Could not kill process group: {e}")
                sys.exit(1)

            if proc.returncode != 0:
                print(f"[ERROR] Snakemake failed with exit code {proc.returncode}")
                return False

            return True

        elif user_input in {"no", "n"}:
            print("Check errors in the config file or snakefile and try again.")
            return False
        else:
            print("Please enter either 'yes' or 'y', or 'no' or 'n'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to config.yml")
    parser.add_argument("--snakefile", required=True, help="Path to Snakefile")
    args = parser.parse_args()

    run_incremental(args.config, args.snakefile)
