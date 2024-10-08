#!/usr/bin/env python3
"""
GoldPolish-Target: Targeted long read polishing with GoldPolish
Written by Emily Zhang (@ezhang)
"""
import argparse
import os
import shlex
import subprocess
from random import choice
from string import ascii_uppercase
import btllib

def parse_args():
    """Parses Arguments passed by users via CLI"""
    parser = argparse.ArgumentParser(
        description="GoldPolish-Target: Targeted long read polishing with GoldPolish"
    )
    parser.add_argument(
        "-f", "--fasta", help="Input genome fasta file", type=str, required=True
    )
    parser.add_argument(
        "-l",
        "--length",
        help="Length of flanking sequences [64]",
        type=int,
        required=False,
        default=64,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Name of output file",
        type=str,
        required=False,
        default="GoldPolish-Target_out.fa",
    )
    parser.add_argument(
        "-s",
        "--subsample-max-reads-per-10kbp",
        help="s parameter for GoldPolish ",
        default=100,
    )
    parser.add_argument(
        "-x",
        "--mx-max-reads-per-10kbp",
        help="x parameter for GoldPolish",
        default=150
    )
    parser.add_argument("-r", "--reads", help="reads file", required=True)

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--ntLink",
        action="store_true",
        help="Run ntLink to generate read mappings (default).",
    )
    group.add_argument(
        "--minimap2",
        action="store_true",
        help="Run minimap2 to generate read mappings.",
    )

    parser.add_argument(
        "--k_ntlink",
        help="k-mer size used for ntLink mappings (if --ntlink specified)",
        type=int,
        default=88,
    )
    parser.add_argument(
        "--w_ntlink",
        help="Window size used for ntLink mappings (if --ntlink specified)",
        type=int,
        default=1000,
    )
    parser.add_argument(
        "--sensitive", 
        help="Sensitive remapping (if --ntlink specified)", 
        default=True,
        required=True
    )

    parser.add_argument(
        "-b",
        "--bed",
        help="BED file specifying coordinates to polish",
        type=str,
        default="",
    )

    parser.add_argument(
        "-n",
        "--dry-run",
        help="Print out the commands that will be executed",
        action="store_true",
    )
    parser.add_argument(
        "--benchmark",
        action='store_true',
        help="Store benchmarks for each step of the GoldPolish-Target pipeline",
    )
    parser.add_argument(
        "--dev",
        action='store_true',
        help="Keeps intermediate files"
    )
    parser.add_argument("-t", help="Number of threads [48]", type=int, default=48)
    parser.add_argument(
        "-v", "--version", action="version", version="GoldPolish-Target v0.0.1"
    )

    args = parser.parse_args()

    if args.t < 2:
        args.t = 2
        btllib.log_warning("Threads number is too low. Setting to 2.")
    return args


def get_mapping_info(minimap2, ntLink):
    """Determines which mapper to use based on user input"""
    if ntLink:
        mapper = "ntLink"
        s_goldpolish = 100
    elif minimap2:
        mapper = "minimap2"
        s_goldpolish = 40
    else:
        raise ValueError("Unknown mapping tool")
    return (mapper, s_goldpolish)

def split_output_file_name(output_file):
    """Removes suffix from end of output file to create a prefix"""
    if output_file.endswith('.fa'):
        return output_file[:-3]
    elif output_file.endswith('.fasta'):
        return output_file[:-6]
    elif output_file.endswith('.fna'):
        return output_file[:-4]
    return output_file

def cleanup(rand_str):
    files = [f for f in os.listdir() if os.path.isfile(f)]
    for file in files:
        if f'INTERMEDIATE_{rand_str}' in file:
            os.remove(file)

def main():
    "Run goldpolish-target snakemake file"
    args = parse_args()
    base_dir = os.path.dirname(os.path.realpath(__file__))
    mapper, s = get_mapping_info(
        args.minimap2,
        args.ntLink,
    )

    prefix = split_output_file_name(args.output)

    command = (
        f"snakemake -s {base_dir}/goldpolish-target-run-pipeline.smk --cores {args.t} "
        f"target --config fasta={args.fasta} length={args.length} threads={args.t} "
        f"mapper={mapper} prefix={prefix} output={args.output} reads={args.reads} s_param={s} "
        f"x_param={args.mx_max_reads_per_10kbp} sensitive={args.sensitive} "
        f"k_ntlink={args.k_ntlink} w_ntlink={args.w_ntlink} "
    )

    rand_str = ''.join(choice(ascii_uppercase) for i in range(5))
    
    command += "benchmark=True " if args.benchmark else "benchmark=False "
    command += f"delete_intermediates=INTERMEDIATE_{rand_str}. " if not args.dev else ""
    command += f"bed={args.bed}" if args.bed else ""

    if args.dry_run:
        command += " -n "

    print(f"Running {command}", flush=True)

    command = shlex.split(command)

    ret = subprocess.call(command)
    if ret != 0:
        raise subprocess.SubprocessError(
            "GoldPolish-Target failed - check the logs for the error."
        )
    if not args.dev:
        cleanup(rand_str)

if __name__ == "__main__":
    main()
