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

meson, ninja, btllib, and ntLink can be installed through [Conda](https://docs.conda.io/en/latest/) package manager:
```
conda install -c conda-forge meson ninja 
conda install -c bioconda btllib ntlink
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