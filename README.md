# Code repository for "Differential co-localisation analysis of multi-sample and multi-condition experiments with `spatialFDA`"
Repository organising the code to generate the results of the `spatialFDA` manuscript. `spatialFDA` is available on [Bioconductor](https://bioconductor.org/packages/release/bioc/html/spatialFDA.html)

The analysis is written as a `snakemake` workflow. In order to run the scripts [install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Requirements

### Dependencies

To run this project you need to have installed `snakemake --version 9.16.3`, `conda` or a similar option like `mamba` and `apptainer`. You can also run this project with `conda` only, but then you have to specify your `python` environments yourself.

### Data

The dataset is available through `Bioconductors` `ExperimentHub` via the package `imcdatasets` and downloaded directly in the workflow. The original data was deposited at https://data.mendeley.com/datasets/cydmwsfztj/2.

## Quick start

In order to run the self-contained reproducible environment please run:

### 1. Clone analysis repo and move to the correct directory

```
git clone https://github.com/mjemons/spatialfda-manuscript.git
cd spatialfda-manuscript/
```
### 2. Run the workflow

```
snakemake --cores <nCores> --sdm conda apptainer
```

If you want to use only the `conda` environment without `apptainer` you might need to specify a correct `python` installation.

It can happen that you run into an `apptainer` error when setting up the environment. In this case please specify a new `tmp` directory via `--apptainer-args`

### Citation

Include citation after bioRxiv upload

## Contact

For questions, please open an issue in this repository
