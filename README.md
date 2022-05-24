# GoldRush-Edit

GoldRush-Edit is an efficient genome polishing tool that uses long reads. [ntEdit](https://github.com/bcgsc/ntEdit) polishes the draft assembly and flags erroneous regions, then [Sealer](https://github.com/bcgsc/abyss/tree/master/Sealer) fills assembly gaps and erroneous sequence regions flagged by ntEdit. The polisher is adapted from [ntedit_sealer_protocol](https://github.com/bcgsc/ntedit_sealer_protocol/) to use long reads instead of short reads.

## Dependencies

- GNU Make
- Python 3
- [btllib](https://github.com/bcgsc/btllib) v1.4.3+

btllib can be installed through [Conda](https://docs.conda.io/en/latest/) package manager:
```
conda install -c bioconda btllib 
```

## Usage

To polish a draft assembly named `assembly.fa` with long reads named `reads.fa` and store the results at `assembly-polished.fa`, run the following:
```
goldrush-edit assembly.fa reads.fa assembly-polished.fa
```