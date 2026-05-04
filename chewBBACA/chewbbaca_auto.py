#!/usr/bin/env python

import sys
import argparse
import subprocess
import readline
import os
import glob
import signal
from pathlib import Path


def complete(text, state):
    matches = glob.glob(text + '*')
    try:
        return matches[state] + os.sep if os.path.isdir(matches[state]) else matches[state]
    except IndexError:
        return None

readline.set_completer_delims('')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

def create_config(res, cg_scheme, partition, ass):
    with open(f'{res}/config.yml', 'w') as f:
        f.write(f'results_path: {res}\n')
        f.write(f'scheme: {cg_scheme}\n')
        #f.write(f'blastdb: {b_db}\n')
        f.write(f'partition: {partition}\n')
        f.write(f'fasta: {ass}\n')

def prompt_paths():
    while True:
        results_folder = input("Enter the path to the results folder (. for current directory): ")
        if results_folder and os.path.isdir(results_folder):
            res = os.path.abspath(results_folder)
            break
        elif results_folder:
            res = os.path.abspath(results_folder)
            os.makedirs(res)
            break
        else:
            print("Please enter the path to the results folder (. for current directory): ")

    while True:
        assemblies = input("Enter the path to the folder with the assemblies (. for current directory): ")
        if os.path.isdir(assemblies):
            ass = os.path.abspath(assemblies)
            break
        else:
            print("The path is not valid. Please enter again.")

    while True:
        cgmlst_scheme = input("""Enter the path to the chewBBACA prepared scheme folder:

    Campylobacter: /users/data/Projects/FBI_SOFI/proj/chewBBACA/Resources/dbs_chewBBACA_3/Campy_oxford_cgMLST_scheme_chewBBACA_sizefilter40_st_6p
    E. coli: /dpssi/data/Projects/FBI_SOFI/proj/chewBBACA/Resources/dbs_chewBBACA_3/Escherichia_cgMLSTv1_clean_chew_25p
    C. Difficile: /users/data/Projects/FBI_SOFI/proj/chewBBACA/Resources/dbs_chewBBACA_3/cdiff_cgMLST_ridom_clean_chew_st5p_trn
    Salmonella: /dpssi/data/Projects/FBI_SOFI/proj/chewBBACA/Resources/dbs_chewBBACA_3/enterobase_Salmonella_cgMLSTv2_230207_clean_chew_st38p
    Listeria: /users/data/Projects/FBI_SOFI/proj/chewBBACA/Resources/dbs_chewBBACA_3/listeria_bigsdb_cgmlst1748_chewBBACA_sizefilter40_st_5p

Paste the path for the correct species:
""")

        if os.path.isdir(cgmlst_scheme):
            cg_scheme = os.path.abspath(cgmlst_scheme)
            break
        else:
            print("Invalid path. Try again.")


    while True:
        partition = input("Enter partition (standard, daytime, surveillance, outbreak): ")
        if partition in {"standard", "daytime", "surveillance", "outbreak"}:
            break
        else:
            print("Invalid partition. Try again.")

    return res, cg_scheme, partition, ass

def run_incremental(config_file):
    print("Running BLAST gene calls, please wait...")
    incremental_script = "/users/data/Projects/FBI_SOFI/proj/chewBBACA/Tools/scripts/new_per_locus/new/chewbbaca_incremental.py"
    snakefile_path = "/users/data/Projects/FBI_SOFI/proj/chewBBACA/Tools/scripts/new_per_locus/new/snakefile_genecall.smk"

    run_cmd = [
        "python", incremental_script,
        "--config", config_file,
        "--snakefile", snakefile_path
    ]

    result = None
    try:
        result = subprocess.Popen(run_cmd, preexec_fn=os.setsid)
        result.communicate()
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Terminating all subprocesses.")
        if result is not None:
            try:
                os.killpg(os.getpgid(result.pid), signal.SIGKILL)
            except Exception as e:
                print(f"[WARNING] Could not kill process group: {e}")
        sys.exit(1)

    if result.returncode != 0:
        print(f"\nIncremental step failed with code {result.returncode}")
        sys.exit(result.returncode)

    return True

def run_chewbbaca(config_file):
    chew_snakefile = "/users/data/Projects/FBI_SOFI/proj/chewBBACA/Tools/scripts/new_per_locus/new/snakefile_chewbbaca.smk"
    print("Running dryrun for chewBBACA allele calling...")

    dry = None
    try:
        dry = subprocess.Popen(
            f"snakemake --snakefile {chew_snakefile} --configfile {config_file} --profile slurm -np",
            shell=True, executable="/bin/bash", preexec_fn=os.setsid
        )
        dry.communicate()
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Terminating dryrun.")
        if dry is not None:
            try:
                os.killpg(os.getpgid(dry.pid), signal.SIGKILL)
            except Exception as e:
                print(f"[WARNING] Could not kill dryrun process group: {e}")
        sys.exit(1)

    if dry.returncode != 0:
        print(f"[ERROR] Snakemake dryrun failed with exit code {dry.returncode}")
        return False

    while True:
        snakemake_errors = input("Did the allele calling dryrun run without any errors (red)? [yes/y or no/n]: ")
        if snakemake_errors.lower() in {"yes", "y"}:
            result = None
            try:
                result = subprocess.Popen(
                    f"snakemake --snakefile {chew_snakefile} --configfile {config_file} --profile slurm",
                    shell=True, executable="/bin/bash", preexec_fn=os.setsid
                )
                result.communicate()
            except KeyboardInterrupt:
                print("\n[INTERRUPTED] Terminating allele calling.")
                if result is not None:
                    try:
                        os.killpg(os.getpgid(result.pid), signal.SIGKILL)
                    except Exception as e:
                        print(f"[WARNING] Could not kill process group: {e}")
                sys.exit(1)

            if result.returncode != 0:
                print(f"[ERROR] Snakemake failed with exit code {result.returncode}")
                return False

            return True

        elif snakemake_errors.lower() in {"no", "n"}:
            print("Check errors in the config file or snakefile and try again.")
            return False
        else:
            print("Please enter either 'yes' or 'y', or 'no' or 'n'.")

def run_plot(config_file):
    print("Gene and allele calls completed, now running dryrun of snakemake to complete distance matrix and plot, please wait...")

    plot_snakefile = "/users/data/Projects/FBI_SOFI/proj/chewBBACA/Tools/scripts/new_per_locus/new/Snakefile_chewBBACA_dist_plot.smk"

    dry = None
    try:
        dry = subprocess.Popen(
            f"snakemake --snakefile {plot_snakefile} --configfile {config_file} --profile slurm -np",
            shell=True, executable="/bin/bash", preexec_fn=os.setsid
        )
        dry.communicate()
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Terminating dryrun.")
        if dry is not None:
            try:
                os.killpg(os.getpgid(dry.pid), signal.SIGKILL)
            except Exception as e:
                print(f"[WARNING] Could not kill dryrun process group: {e}")
        sys.exit(1)

    if dry.returncode != 0:
        print(f"[ERROR] Snakemake dryrun failed with exit code {dry.returncode}")
        return False

    while True:
        snakemake_errors = input("Did the plot dryrun run without any errors (red)? [yes/y or no/n]: ")
        if snakemake_errors.lower() in {"yes", "y"}:
            snake = None
            try:
                snake = subprocess.Popen(
                    f"snakemake --snakefile {plot_snakefile} --configfile {config_file} --profile slurm",
                    shell=True, executable="/bin/bash", preexec_fn=os.setsid
                )
                snake.communicate()
            except KeyboardInterrupt:
                print("\n[INTERRUPTED] Terminating plotting step.")
                if snake is not None:
                    try:
                        os.killpg(os.getpgid(snake.pid), signal.SIGKILL)
                    except Exception as e:
                        print(f"[WARNING] Could not kill plot process group: {e}")
                sys.exit(1)

            if snake.returncode != 0:
                print(f"[ERROR] Snakemake failed with exit code {snake.returncode}")
                return False

            return True
        elif snakemake_errors.lower() in {"no", "n"}:
            print("Check errors in the config file or snakefile and try again.")
            return False
        else:
            print("Please enter either 'yes' or 'y', or 'no' or 'n'.")

def main():
    while True:
        first_time = input("Is this the first time running the pipeline for this dataset? [yes/y or no/n]: ")

        if first_time.lower() in {"yes", "y"}:
            res, cg_scheme, partition, ass = prompt_paths()
            create_config(res, cg_scheme, partition, ass)
            config_file = os.path.join(res, "config.yml")
            break

        elif first_time.lower() in {"no", "n"}:
            while True:
                config_path = input("Enter the path to the existing config.yml file: ")
                if os.path.isfile(config_path):
                    config_file = os.path.abspath(config_path)
                    break
                else:
                    print("Invalid config file path. Please provide a valid file.")

            while True:
                res_folder = input("Enter the path to the associated results folder: ")
                if os.path.isdir(res_folder):
                    res = os.path.abspath(res_folder)
                    break
                else:
                    print("Invalid results folder path. Please provide a valid directory.")
            break

        else:
            print("Either enter yes/y or no/n")

    if not run_incremental(config_file):
        sys.exit(1)
    if not run_chewbbaca(config_file):
        sys.exit(1)
    if not run_plot(config_file):
        sys.exit(1)

if __name__ == "__main__":
    main()
