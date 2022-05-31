# GoldRush-Edit

GoldRush-Edit is an efficient draft genome assembly polishing tool that uses long reads for polishing. [ntEdit](https://github.com/bcgsc/ntEdit) polishes the draft assembly and flags erroneous regions, then [Sealer](https://github.com/bcgsc/abyss/tree/master/Sealer) fills assembly gaps and erroneous sequence regions flagged by ntEdit. The polisher is adapted from [ntedit_sealer_protocol](https://github.com/bcgsc/ntedit_sealer_protocol/) to use long reads instead of short reads.

## Dependencies

- Build
  * GCC 7+ or Clang 8+ (with OpenMP support)
  * [meson](https://mesonbuild.com/)
  * [ninja](https://ninja-build.org/)
  * [btllib](https://github.com/bcgsc/btllib) v1.4.3+

- Run
  * GNU Make
  * Python 3
  * [btllib](https://github.com/bcgsc/btllib) v1.4.3+
  * [ntLink](https://github.com/bcgsc/ntlink) v1.3.0+
  * [minimap2](https://github.com/lh3/minimap2)

The dependencies can be installed through [Conda](https://docs.conda.io/en/latest/) package manager:
```
conda install -c conda-forge meson ninja 
conda install -c bioconda btllib ntlink minimap2
```

## Installation

To build GoldRush-Edit and install it at `$GOLDRUSH_EDIT_PREFIX`, run the following commands from within the `goldrush-edit` directory:
```
meson setup build --buildtype release --prefix $GOLDRUSH_EDIT_PREFIX
cd build
ninja install
```

## Usage

To polish a draft assembly named `assembly.fa` with long reads named `reads.fa` and store the results at `assembly-polished.fa`, run the following:
```
goldrush-edit assembly.fa reads.fa assembly-polished.fa
```

You can run `goldrush-edit --help` to see the available options:
```
usage: goldrush-edit [-h] [-k K] [-b BSIZE] [-m SHARED_MEM] [-t THREADS] [-v] [-x MX_MAX_READS_PER_10KBP] [-s SUBSAMPLE_MAX_READS_PER_10KBP]
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
```