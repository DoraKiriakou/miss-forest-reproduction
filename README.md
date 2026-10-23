# Miss Forest Reproduction

This repository contains a reproduction of the main experiment of the missForest algorithm presented in the paper 
> STEKHOVEN, Daniel J.; BÜHLMANN, Peter. MissForest—non-parametric missing value imputation for mixed-type data. Bioinformatics, 2012, 28.1: 112-118.

The algorithm is run on three UCI datasets, one per data type: Parkinson's (continuous), SPECT (categorical), and Muck version 1 (mixed). These three datasets are complete, so values are deleted at random (10, 20, and 30%). Then their missing values are imputed with a baseline method (mean for continuous and most frequent category for categorical), and with missForest. The error of the imputation is measured only on the deleted entries (NRMSE for numeric, PFC for categorical). For each method, dataset, and missing rate the experiment is repeated 50 times as in the paper.

The authors' R package `missForest` is used as is. The code I wrote covers the loading and reading of the datasets, the random deletion, the two error measures, the experiment loop with seeds for reproducibility, and a final summary table.

## Repository Contents

| Path | Description |
| --- | --- |
| `MissForest.R` | The full R script (loading, random deletion, error measures, experiment, summary table) |
| `datasets/` | UCI data files used by the R script |
| `results/` | Output of my full run (50 repetitions): all repetitions per dataset `.rds` files, and the `summary_table.csv` |

## Requirements To Run Code

- R 4.6.1
- R packages: `missForest` 1.6.1, and `this.path` 2.8.0 (they are installed automatically by the script if missing)

## How To Run

**`Note:`** Before running anything, the full experiment took almost 7.5 hours on MacBook M2 Pro. If you want a quicker run, open `MissForest.R` and change on line 24 *reps_musk <- 50* to *reps_musk <- 2* (or 3).

1. Install [**R**](https://cran.r-project.org) and [**RStudio**](https://posit.co/download/rstudio-desktop/) (if you do not have them).
2. Download or clone the repository, and keep the folder structure as it is.
3. Open `MissForest.R` in RStudio.
4. Run the entire file.

No paths need to be edited, if you keep the folder structure unchanged. The script finds its own folder, reads from `datasets/`, and writes on `results/`. The required packages (`missForest`, `this.path`) are installed automatically on the first run, so an internet connection is needed.
