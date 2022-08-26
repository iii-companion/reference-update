# reference-update
Automatically build a reference directory structure for a VEuPathDB group for use in Companion pipeline.
## Quick Start
1. Ensure Dependencies are installed.
1. Create a *params_default.config* from the template *params_default.example* and fill out.
1. Run the following on command line (see [Nextflow docs](https://www.nextflow.io/docs/latest/cli.html#run) for OPTIONS):
```
nextflow run reference_update.nf [OPTIONS]
```
1. If successful, fully populated reference directory will be available at **params.REFERENCE_PATH**.

## Dependencies

* AUGUSTUS
* Companion
* EuPathWS
* PorthoMCL
    * Perl 5.x
    * Python 2.x
    * NCBI Blast
    * MCL
