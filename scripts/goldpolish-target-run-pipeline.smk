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
max_threads = config["t"] if "t" in config else 48
benchmark = config["benchmark"] if "benchmark" in config else False
sensitive = config["sensitive"] if "sensitive" in config else True

# If want to benchmark, use memusg or /usr/bin/time
# The snakemake benchmarking is quite inaccurate for shorter runtimes
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
           expand("{prefix}.gaps.fa", prefix=prefix),
           expand("{prefix}.gaps.goldpolished.fa", prefix=prefix),
           expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.verbose_mapping.tsv", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink),
           expand("{prefix}.gaps.fa.k{k_ntlink}.w{w_ntlink}.z1000.verbose_mapping.tsv", prefix=prefix, k_ntlink=k_ntlink, w_ntlink=w_ntlink)

rule minimap2_target:
    input: expand("{prefix}.polished.fa", prefix=prefix),
           expand("{prefix}.gaps.fa", prefix=prefix),
           expand("{prefix}.gaps.goldpolished.fa", prefix=prefix),
           expand("{prefix}.gaps.fa.paf", prefix=prefix),
           expand("{prefix}.unpolished.paf", prefix=prefix)

rule run_ntLink_pair:
    input: fa=expand("{fasta}", fasta=fasta),
           reads=expand("{reads}", reads=reads)
    output: expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.verbose_mapping.tsv", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink)
    params: options=expand("sensitive={sensitive} dev=True k={k_ntlink} w={w_ntlink} t={max_threads}", sensitive=sensitive, k_ntlink=k_ntlink, w_ntlink=w_ntlink, max_threads=max_threads),
            benchmarking=expand("{benchmark_path} -o {prefix}.ntLink_pair.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} ntLink pair target={input.fa} reads={input.reads} {params.options}"

rule run_minimap2:
    input: fa=expand("{fasta}", fasta=fasta),
           reads=expand("{reads}", reads=reads)
    output: "{prefix}.unpolished.paf"
    params: options=expand("-x map-ont -t {max_threads}", max_threads=max_threads),
            benchmarking=expand("{benchmark_path} -o {{prefix}}.minimap2.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} minimap2 {params.options} {input.fa} {input.reads} > {output}"

rule extract_seq:
    input: fa=expand("{fasta}", fasta=fasta),
    output: "{prefix}.gaps.fa"
    params: options=expand("-l {length} -o {prefix}.gaps.fa", length=length, prefix=prefix),
            path_to_script=expand("{script_path}/extract_seq.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {{prefix}}.extract_seq.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} -f {input.fa} {params.options}"

# rule for paf file/minimap2 remapping
rule update_mapping_paf:
    input: gaps="{prefix}.gaps.fa",
           mapping=rules.run_minimap2.output
    output: "{prefix}.gaps.fa.paf"
    params: path_to_script=expand("{script_path}/update_mapping.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {{prefix}}.update_mapping.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} -g {input.gaps} -m {input.mapping} -o {output}"

# rule for verbose_mapping file / ntlink pair remapping
rule update_mapping_tsv:
    input: gaps=expand("{prefix}.gaps.fa", prefix=prefix),
           mapping=rules.run_ntLink_pair.output
    output: expand("{prefix}.gaps.fa.k{k_ntlink}.w{w_ntlink}.z1000.verbose_mapping.tsv", prefix=prefix, k_ntlink=k_ntlink, w_ntlink=w_ntlink)
    params: path_to_script=expand("{script_path}/update_mapping.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {prefix}.update_mapping.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} -g {input.gaps} -m {input.mapping} -o {output}"

rule run_goldpolish:
    input: mapping=choose_mapping,
           gaps="{prefix}.gaps.fa"
    output: "{prefix}.gaps.goldpolished.fa"
    params: options=expand("-s {s} -t {max_threads}", s=s, max_threads=max_threads),
            benchmarking=expand("{benchmark_path} -o {{prefix}}.goldpolish.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} goldpolish --mappings {input.mapping} {params.options} {input.gaps} {reads} {output}"

rule run_post_processing:
    input: "{prefix}.gaps.goldpolished.fa"
    output: "{prefix}.polished.fa"
    params: options=expand("-f {fasta} -o {prefix}.polished.fa", fasta=fasta, prefix=prefix),
            path_to_script=expand("{script_path}/post_processing.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {{prefix}}.post_processing.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} {params.options} -g {input}"