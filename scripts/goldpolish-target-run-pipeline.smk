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
bed = config["bed"] if "bed" in config else ""
k_ntlink = config["k_ntlink"] if "k_ntlink" in config else 88
w_ntlink = config["w_ntlink"] if "w_ntlink" in config else 1000
prefix = config["p"] if "p" in config else "GoldPolish-Target_out"
s = config["s"] if "s" in config else 100
x = config["x"] if "x" in config else 150
max_threads = config["t"] if "t" in config else 48
benchmark = config["benchmark"] if "benchmark" in config else False

sensitive = config["sensitive"] if "sensitive" in config else True
intermediate = config["delete_intermediates"] if "delete_intermediates" in config else 'i.'

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

rule target:
    input: expand("{prefix}.polished.fa", prefix=prefix),
           expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           expand("{prefix}.gaps.goldpolished.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           expand("{prefix}.unpolished.{mapper}.{intermediate}paf", prefix=prefix, intermediate=intermediate, mapper=mapper),
           expand("{prefix}.gaps.fa.{intermediate}paf", prefix=prefix, intermediate=intermediate)

rule run_ntLink_pair:
    input: fa=expand("{fasta}", fasta=fasta),
           reads=expand("{reads}", reads=reads)
    output: paf=expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.paf", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink),
            verbose_mapping=expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.verbose_mapping.tsv", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink),
            tsv=expand("{fasta}.k{k_ntlink}.w{w_ntlink}.tsv", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink),
            dot=expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.n1.scaffold.dot", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink)
    params: options=expand("sensitive={sensitive} dev=True k={k_ntlink} w={w_ntlink} t={max_threads} paf=True", sensitive=sensitive, k_ntlink=k_ntlink, w_ntlink=w_ntlink, max_threads=max_threads),
            benchmarking=expand("{benchmark_path} -o {prefix}.ntLink_pair.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} ntLink pair target={input.fa} reads={input.reads} {params.options}"

rule update_ntLink_pair_output_name:
    input: paf=rules.run_ntLink_pair.output.paf,
           verbose_mapping=rules.run_ntLink_pair.output.verbose_mapping,
           tsv=rules.run_ntLink_pair.output.tsv,
           dot=rules.run_ntLink_pair.output.dot
    output: paf=expand("{prefix}.unpolished.ntLink.{intermediate}paf", prefix=prefix, intermediate=intermediate),
            verbose_mapping=expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.verbose_mapping.{intermediate}tsv", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink, intermediate=intermediate),
            tsv=expand("{fasta}.k{k_ntlink}.w{w_ntlink}.{intermediate}tsv", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink, intermediate=intermediate),
            dot=expand("{fasta}.k{k_ntlink}.w{w_ntlink}.z1000.n1.scaffold.{intermediate}dot", fasta=fasta, k_ntlink=k_ntlink, w_ntlink=w_ntlink, intermediate=intermediate)
    shell: "mv {input.paf} {output.paf} && mv {input.verbose_mapping} {output.verbose_mapping} && mv {input.tsv} {output.tsv} && mv {input.dot} {output.dot}"

rule run_minimap2:
    input: fa=expand("{fasta}", fasta=fasta),
           reads=expand("{reads}", reads=reads)
    output: expand("{prefix}.unpolished.minimap2.{intermediate}paf", prefix=prefix, intermediate=intermediate)
    params: options=expand("-t {max_threads}", max_threads=max_threads),
            benchmarking=expand("{benchmark_path} -o {prefix}.minimap2.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} minimap2 {params.options} {input.fa} {input.reads} > {output}"

rule extract_seq:
    input: fa=expand("{fasta}", fasta=fasta),
    output: expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate)
    params: options=expand("-l {length} --bed {bed}", length=length, bed=bed) if bed else expand("-l {length}", length=length),
            path_to_script=expand("{script_path}/goldpolish-target-extract-seq.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {prefix}.extract_seq.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} -f {input.fa} {params.options} -o {output}"

# rule for paf file remapping
rule update_mapping_paf:
    input: gaps=expand("{prefix}.gaps.{intermediate}fa", prefix=prefix, intermediate=intermediate),
           mapping=expand("{prefix}.unpolished.{mapper}.{intermediate}paf", prefix=prefix, intermediate=intermediate, mapper=mapper)
    output: expand("{prefix}.gaps.fa.{intermediate}paf", prefix=prefix, intermediate=intermediate)
    params: path_to_script=expand("{script_path}/goldpolish-target-update-mapping.py", script_path=script_path),
            benchmarking=expand("{benchmark_path} -o {prefix}.update_mapping.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} python {params.path_to_script} -g {input.gaps} -m {input.mapping} -o {output}"

rule run_goldpolish:
    input: mapping=rules.update_mapping_paf.output,
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
