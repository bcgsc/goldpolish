# ntEdit+Sealer Assembly Finishing Protocol

An automated protocol for finishing long-read genome assemblies using long reads. [ntEdit](https://github.com/bcgsc/ntEdit) polishes the draft assembly and flags erroneous regions, then [Sealer](https://github.com/bcgsc/abyss/tree/master/Sealer) fills assembly gaps and erroneous sequence regions flagged by ntEdit. The protocol is implemented as a Makefile pipeline.

## Dependencies

- GNU Make
- Python 3
- [btllib](https://github.com/bcgsc/btllib) v1.4.3+
- [ntEdit](https://github.com/bcgsc/ntEdit) X
- [ABySS](https://github.com/bcgsc/abyss) X (includes Sealer)