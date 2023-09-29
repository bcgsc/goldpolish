#!/usr/bin/env snakemake -s
import shutil
import sys
import os

onsuccess:
    shutil.rmtree(".snakemake", ignore_errors=True)

# Read in parameters
fasta = config["f"] if "f" in config else "Must specify input assembly 'fasta'"
reads = config["reads"] if "reads" in config else "Must specify reads"
length = config["l"] if "l" in config else 64
mappings = config["mappings"] if "mappings" in config else ""
agp = config["agp"] if "agp" in config else ""
mapper = config["mapper"] if "mapper" in config else ""
bed = config["b"] if "b" in config else ""
k_ntlink = config["k_ntlink"] if "k_ntlink" in config else 32
w_ntlink = config["w_ntlink"] if "w_ntlink" in config else 1000
prefix = config["p"] if "p" in config else "GoldPolish-Target_out"
s = config["s"] if "s" in config else 100
x = config["x"] if "x" in config else 150
max_threads = config["t"] if "t" in config else 48
benchmark = config["benchmark"] if "benchmark" in config else False
sensitive = config["sensitive"] if "sensitive" in config else True
intermediate = config["delete_intermediates"] if "delete_intermediates" in config else ''

# If want to benchmark, use memusg or /usr/bin/time
benchmark_path = "" 
if benchmark:
    if benchmark_memusg := shutil.which("memusg"):
        benchmark_path = f"{benchmark_memusg} -t"
    elif benchmark_time := shutil.which("time"):
        benchmark_path = f"{benchmark_time} -v"
    else:
        print("WARNING: memusg and time not found in PATH. Benchmarks will not be tallied.")

# Get path to the base directory (where snakemake file is located)
script_path = workflow.basedir

def choose_mapping(Wildcards):
    if config["mapper"] == "ntLink":
        return rules.update_mapping_tsv.output
    elif config["mapper"] == "minimap2":
        return rules.update_mapping_paf.output
    else:
        print("invalid mapping option specified")

rule ntLink_target:
    input: expand("{prefix}.polished.fa", prefix=prefix),
           expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           expand("{prefix}.gaps.goldpolished.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.{intermediate}paf", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink, intermediate=intermediate),
           expand("{prefix}.gaps.fa.k{k_ntlink}.w{w_ntlink}.{intermediate}paf", prefix=prefix, k_ntlink=k_ntlink, w_ntlink=w_ntlink, intermediate=intermediate)

rule minimap2_target:
    input: expand("{prefix}.polished.fa", prefix=prefix),
           expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           expand("{prefix}.gaps.goldpolished.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           expand("{prefix}.gaps.fa.{intermediate}paf", prefix=prefix, intermediate=intermediate),
           expand("{prefix}.unpolished.{intermediate}paf", prefix=prefix, intermediate=intermediate)

rule run_ntLink_pair:
    input: fa=expand("{fasta}", fasta=fasta),
           reads=expand("{reads}", reads=reads)
    output: expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.{intermediate}paf", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink, intermediate=intermediate)
    params: options=expand("sensitive={sensitive} dev=True k={k_ntlink} w={w_ntlink} t={max_threads} paf=True", sensitive=sensitive, k_ntlink=k_ntlink, w_ntlink=w_ntlink, max_threads=max_threads),
            benchmarking=expand("{benchmark_path} -o {prefix}.ntLink_pair.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} ntLink pair target={input.fa} reads={input.reads} {params.options}"

rule run_minimap2:
    input: fa=expand("{fasta}", fasta=fasta),
           reads=expand("{reads}", reads=reads)
    output: expand("{prefix}.unpolished.{intermediate}paf", prefix=prefix, intermediate=intermediate)
    params: options=expand("-t {max_threads}", max_threads=max_threads),
            benchmarking=expand("{benchmark_path} -o {prefix}.minimap2.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} minimap2 {params.options} {input.fa} {input.reads} > {output}"

rule extract_seq:
    input: fa=expand("{fasta}", fasta=fasta),
    output: expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate)
    params: options=expand("-l {length}", length=length),
            path_to_script=expand("{script_path}/goldpolish-target-extract-seq.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {prefix}.extract_seq.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} -f {input.fa} {params.options} -o {output}"

# rule for paf file/minimap2 remapping
rule update_mapping_paf:
    input: gaps=expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           mapping=rules.run_minimap2.output
    output: expand("{prefix}.gaps.fa.{intermediate}paf", prefix=prefix, intermediate=intermediate)
    params: path_to_script=expand("{script_path}/goldpolish-target-update-mapping.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {prefix}.update_mapping.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} -g {input.gaps} -m {input.mapping} -o {output}"

# rule for ntlink pair remapping
rule update_mapping_tsv:
    input: expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           mapping=rules.run_ntLink_pair.output
    output: expand("{prefix}.gaps.fa.k{k_ntlink}.w{w_ntlink}.z1000.{intermediate}paf", prefix=prefix, k_ntlink=k_ntlink, w_ntlink=w_ntlink, intermediate=intermediate)
    params: path_to_script=expand("{script_path}/goldpolish-target-update-mapping.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {prefix}.update_mapping.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} -g {input.gaps} -m {input.mapping} -o {output}"

rule run_goldpolish:
    input: mapping=choose_mapping,
           gaps=expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate)
    output: expand("{prefix}.gaps.goldpolished.{intermediate}fa", prefix=prefix, intermediate=intermediate)
    params: options=expand("-s {s} -x {x} -t {max_threads}", s=s, x=x, max_threads=max_threads),
            benchmarking=expand("{benchmark_path} -o {prefix}.goldpolish.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} goldpolish --mappings {input.mapping} {params.options} {input.gaps} {reads} {output}"

rule run_post_processing:
    input: expand("{prefix}.gaps.goldpolished.{intermediate}fa", prefix=prefix, intermediate=intermediate)
    output: "{prefix}.polished.fa"
    params: options=expand("-f {fasta}", fasta=fasta),
            path_to_script=expand("{script_path}/goldpolish-target-post-processing.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {{prefix}}.post_processing.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} {params.options} -g {input} -o {output}"