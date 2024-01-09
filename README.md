![Logo](https://github.com/bcgsc/goldpolish/blob/master/goldpolish-logo.png)

# GoldPolish (aka GoldRush-Edit)

GoldPolish (aka GoldRush-Edit) is an efficient draft genome assembly polishing tool that uses long reads for polishing. [ntEdit](https://github.com/bcgsc/ntEdit) polishes the draft assembly and flags additional erroneous regions, then [Sealer](https://github.com/bcgsc/abyss/tree/master/Sealer) fills assembly gaps and erroneous sequence regions flagged by ntEdit. The polisher is adapted from [ntedit_sealer_protocol](https://github.com/bcgsc/ntedit_sealer_protocol/) to use long reads instead of short reads.

## Dependencies

- Build
  * GCC 7+ or Clang 8+ (with OpenMP support)
  * [meson](https://mesonbuild.com/)
  * [ninja](https://ninja-build.org/)
  * [btllib](https://github.com/bcgsc/btllib) v1.6.2+
  * [boost](https://www.boost.org/)

- Run
  * GNU Make
  * Python 3
  * [btllib](https://github.com/bcgsc/btllib) v1.6.2+
  * [ntLink](https://github.com/bcgsc/ntlink) v1.3.5+
  * [minimap2](https://github.com/lh3/minimap2)
  * [snakemake](https://github.com/snakemake/snakemake)

The dependencies can be installed through [Conda](https://docs.conda.io/en/latest/) package manager:
```
conda install -c conda-forge -c bioconda compilers meson ninja boost-cpp btllib ntlink minimap2 snakemake

```

## Citation
If you use GoldPolish in your research, please cite:

Wong J, Coombe L, Nikolić V, Zhang E, Nip KM, Sidhu P, Warren RL and Birol I (2023). Linear time complexity de novo long read genome assembly with GoldRush. Nature Communications, 14(1), 2906. https://doi.org/10.1038/s41467-023-38716-x

## Installation

To build GoldPolish and install it at `$GOLDPOLISH_PREFIX`, run the following commands from within the `goldpolish` directory:
```
meson setup build --buildtype release --prefix $GOLDPOLISH_PREFIX
cd build
ninja install
```

## Usage

To polish a draft assembly named `assembly.fa` with long reads named `reads.fa` and store the results at `assembly-polished.fa`, run the following:
```
goldpolish assembly.fa reads.fa assembly-polished.fa
```

You can run `goldpolish --help` to see the available options:
```
usage: goldpolish [-h] [-k K] [-b BSIZE] [-m SHARED_MEM] [-t THREADS] [-v] [-x MX_MAX_READS_PER_10KBP] [-s SUBSAMPLE_MAX_READS_PER_10KBP]
                  [--ntlink | --minimap2 | --mappings MAPPINGS]
                  seqs_to_polish polishing_seqs output_seqs

positional arguments:
  seqs_to_polish        Sequences to polish.
  polishing_seqs        Sequences to polish with.
  output_seqs           Filename to write polished sequences to.

optional arguments:
  -h, --help            show this help message and exit
  -k K                  k-mer sizes to use for polishing. Example: -k32 -k28 (Default: 32, 28, 24, 20)
  -b BSIZE, --bsize BSIZE
                        Batch size. A batch is how many polished sequences are processed per Bloom filter. (Default: 1)
  -m SHARED_MEM, --shared-mem SHARED_MEM
                        Shared memory path to do polishing in. (Default: /dev/shm)
  -t THREADS, --threads THREADS
                        How many threads to use. (Default: 48)
  -v, --verbose
  -x MX_MAX_READS_PER_10KBP, --mx-max-reads-per-10kbp MX_MAX_READS_PER_10KBP
                        When subsampling, increase the common minimizer count threshold for ntLink mappings until there's at most this many reads per 10kbp of polished sequence.
                        (Default: 150)
  -s SUBSAMPLE_MAX_READS_PER_10KBP, --subsample-max-reads-per-10kbp SUBSAMPLE_MAX_READS_PER_10KBP
                        Random subsampling of mapped reads. For ntLink mappings, this is done after common minimizer subsampling. For minimap2 mappings, only this subsampling is done.
                        By default, 40 if using minimap2 mappings and 100 if using ntLink mappings.
  --ntlink              Run ntLink to generate read mappings (default).
  --minimap2            Run minimap2 to generate read mappings.
  --mappings MAPPINGS   Use provided pre-generated mappings. Accepted formats are PAF, SAM, and *.verbose_mapping.tsv from ntLink.
  --target              Run GoldPolish in targeted mode
  -l LENGTH, --length LENGTH
                        GoldPolish-Target flank length (if --target specified)
  --bed BED             BED file specifying target coordinates (if --target specified)
  --softmask            Target coordinates determined from softmasked regions in the input assembly (if --target specified)
```

## GoldPolish-Target
GoldPolish can be run in a targeted mode, polishing only fragments of the assembly (either by looking for and polishing softmasked sequences or by using target coordinates that are specified in a BED file). To run GoldPolish-Target, use the `--target` flag for the GoldPolish command. Example command running GoldPolish-Target with a BED file:
```
goldpolish --target --bed polishing_coordinates.bed assembly.fa reads.fa assembly-polished.fa 
```  
