# Sustainable Mobility Analysis

This repository contains the R scripts used for the quantitative analyses
presented in the thesis 'Green Travel Preferences in Vienna, Madrid and Tokyo - an international comparison', written at the University of Technology Vienna.

## Scripts

- `01_cluster_analysis.R`: hierarchical cluster analysis for Vienna, Madrid,
  and Tokyo
- `02_sem_model_a.R`: structural equation Model A for Vienna and Madrid
- `03_sem_model_b.R`: structural equation Model B for Tokyo

## Data

The survey datasets are not included in this repository. To run the scripts,
the processed datasets must be stored locally as:

- `data/processed/vienna.csv`
- `data/processed/madrid.csv`
- `data/processed/tokyo.csv`

Generated results are saved in the `outputs` directory.

## Usage

The scripts must be executed from the root directory of the repository.
Required R packages are listed at the beginning of each script.
