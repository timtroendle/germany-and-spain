# Decomposing CO2 emissions in Germany and Spain

A Kaya-like decomposition of CO2 emissions in Germany and Spain.

This repository contains the entire scientific project. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml --no-default-packages

## Run the analysis

    snakemake

This will run all analysis steps to reproduce results.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps `build/dag.pdf` run:

    snakemake -f dag

## Be notified of build successes or fails

  As the execution of this workflow may take a while, you can be notified whenever the execution terminates either successfully or unsuccessfully. Notifications are sent by email. To activate notifications, add the email address of the recipient to the configuration key `email`. You can add the key to your configuration file, or you can run the workflow the following way to receive notifications:

      snakemake test --config email=<your-email>

## Repo structure

* `scripts`: contains the Python source code as scripts
* `rules`: contains Snakemake rule definitions
* `envs`: contains execution environments
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)

## Data sources

This project uses and contains data from eurostat, destatis, and Instituto Nacional de Estadísticas.

Source: [eurostat](https://ec.europa.eu/eurostat/web/main/home), 2022

Source: Bevölkerungsstand, [Statistisches Bundesamt (Destatis)](https://www.destatis.de), 2022

Source: [Instituto Nacional de Estadísticas](www.ine.es)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
