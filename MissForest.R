# =====================================================================
# MissForest reproduction: Parkinson's, SPECT, Musk (mixed)
# =====================================================================

# ---------- Install and import the authors' package ----------
install.packages("missForest") # download the authors' package
install.packages("this.path")
library(missForest) # import the package
library(this.path)


# folder where everything is located in
project_folder <- this.path::this.dir()
datasets <- file.path(project_folder, "datasets") # folder where datasets are stored
res_dir <- file.path(project_folder, "results") # results will be stored here
dir.create(res_dir, showWarnings = FALSE) # create results folder if missing

# Date that this R file runs, use it to change the names of the saving files so there is no overwrite
run_date <- format(Sys.Date(), "%Y-%m-%d") # get a date like 2026-10-20

# Number of repetitions for each rate and dataset
reps_parkinsons <- 50
reps_spect      <- 50
reps_musk       <- 50 # this one takes some hours

# ---------- Load the datasets ----------
# --- Parkinson's ---
pk <- read.csv(file.path(datasets, "parkinsons/parkinsons.data")) # read the file into a dataframe
parkinsons <- pk[, !(names(pk) %in% c("name", "status"))] # drop name and health status columns, keep only continuous cols
print(dim(parkinsons)) # 195 observations and 22 features as described in paper

# SPECT
sp_names <- cat(readLines(file.path(datasets, "spect+heart", "SPECT.names")), sep = "\n") # info about the SPECT dataset
sp_train <- read.csv(file.path(datasets, "spect+heart/SPECT.train"), header = FALSE) # read file into a dataframe
sp_test  <- read.csv(file.path(datasets, "spect+heart/SPECT.test"),  header = FALSE)  # read file into a dataframe
sp <- rbind(sp_train, sp_test) # merge them row-wise to get all 267 observations together
spect <- as.data.frame(lapply(sp[, 2:23], as.factor)) # keep the 22 categorical features, convert each into a category
print(dim(spect)) # 267 observations and 22 features as described in the paper

# Musk, the file is compressed so we unpack it with gzip
musk_info <- cat(readLines(file.path(datasets, "musk+version+1", "clean1.info")), sep = "\n") # info about the musk dataset
mk <- read.csv(file.path(datasets, "musk+version+1", "clean1.data"), header = FALSE) # read file into a dataframe
musk_mixed <- cbind(mk[, 3:168], # columns 3-168 -> the 166 numeric features
                    class = as.factor(mk[, 169])) # column 169 -> musk / non-musk classification as a category, convert into a category
print(dim(musk_mixed)) # 476 observations with 167 features (166 continuous & 1 categorical) as described in the paper

# ---------- Helper functions ----------
# Delete a share of the values completely at random (MCAR)
make_mcar <- function(x, rate) {
  n <- nrow(x) # number of rows
  p <- ncol(x) # number of columns
  mask <- matrix(FALSE, n, p) # a grid of TRUE/FALSE, same shape as the data
  mask[sample(n * p, round(rate * n * p))] <- TRUE   # pick random cells to delete (TRUE = deleted)
  xmis <- x # copy the data
  for (j in seq_len(p)) { # go through the columns one by one
    xmis[mask[, j], j] <- NA # put NA (missing) in the cells of column j that were marked as TRUE (delete them)
  }
  list(xmis = xmis, mask = mask) # return the damaged data and the list of deleted cells
}

# Error for numeric columns (NRMSE), ONLY on the deleted cells
nrmse <- function(ximp, xtrue, mask) {
  nc <- which(sapply(xtrue, is.numeric)) # list of numeric column numbers/positions
  if (length(nc) == 0) return(NA) # no numeric columns, nothing to measure
  t <- unlist(lapply(nc, function(j) xtrue[mask[, j], j]))  # flattened list (single vector) with every deleted true value across all numeric columns
  i <- unlist(lapply(nc, function(j) ximp[mask[, j], j]))   # flattened list (single vector) with every deleted imputed value across all numeric columns
  sqrt(mean((t - i)^2) / var(t)) # squared error / variance of true values, as described in paper
}

# Error for categorical columns (PFC), ONLY on the deleted cells
pfc <- function(ximp, xtrue, mask) {
  cc <- which(!sapply(xtrue, is.numeric)) # list of categorical column numbers/positions
  if (length(cc) == 0) return(NA) # no categorical columns, nothing to measure
  wrong <- sum(unlist(lapply(cc, function(j) # count the deleted cells that were guessed wrong
    as.character(xtrue[mask[, j], j]) != as.character(ximp[mask[, j], j])))) # convert to plain text for the comparison
  wrong / sum(mask[, cc]) # number of wrong guesses / the number of all deleted categorical cells
}

# Simple baseline imputation method to start with | mean for numbers & most common category for categories
mean_mode <- function(x) {
  for (j in seq_along(x)) { # iterate through the columns
    miss <- is.na(x[[j]]) # which cells are missing in this column
    if (is.numeric(x[[j]])) { # if numeric column
      x[miss, j] <- mean(x[[j]], na.rm = TRUE) # replace the missing values with the column's mean
    } else { # if categorical column
      x[miss, j] <- names(which.max(table(x[[j]])))  # replace the missing values with the most frequent category
    }
  }
  x # return the first imputed data using the baseline methods (1st step of algorithm)
}

# ---------- The experiment ----------
run_experiment <- function(name, x, reps, rates = c(0.1, 0.2, 0.3),
                           ntree = 100, maxiter = 10) {
  out <- list() # empty list to collect the results
  for (rate in rates) { # iterate through the different percentages (10%, 20%, 30% missing)
    for (r in seq_len(reps)) { # repeat for each rate and dataset
      set.seed(1000 * r + round(100 * rate)) # fixed seed based on the repetition and rate (reproducible)
      m <- make_mcar(x, rate) # create the damaged data and list of the deleted cells
      
      t0 <- Sys.time() # start the stopwatch
      fit <- missForest(m$xmis, # run MissForest on the damaged data
                        ntree = ntree, # number of trees (from paper 100)
                        maxiter = maxiter, # max iterations of missForest if criterion is not met (we use the default value 10)
                        mtry = floor(sqrt(ncol(x))), # selected variables per split (square root of p from paper)
                        backend = "randomForest") # randomForest retains legacy behavior for compatibility
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))  # seconds MissForest took
      
      bl <- mean_mode(m$xmis) # run the baseline (mean) on the same damaged data
      
      out[[length(out) + 1]] <- data.frame( # store one row per method
        dataset = name, rate = rate, rep = r,
        method  = c("missForest", "mean/mode"),
        NRMSE   = c(nrmse(fit$ximp, x, m$mask), nrmse(bl, x, m$mask)), # numeric error
        PFC     = c(pfc(fit$ximp, x, m$mask),   pfc(bl, x, m$mask)), # categorical error
        seconds = c(secs, NA)) # runtime (only for MissForest)
      
      cat(name, "| rate", rate, "| rep", r, "of", reps, "| took", round(secs), "s\n")  # progress message
    }
  }
  do.call(rbind, out) # merge all rows into one table
}

# ---------- Every Dataset's Run, save every dataset's results ----------
# Run Parkinson's Dataset
res_pk <- run_experiment("Parkinson's", parkinsons, reps = reps_parkinsons) # 50 repetitions of every 3 missing rates
saveRDS(res_pk, file.path(res_dir, paste0("res_parkinsons_", run_date, ".rds"))) # save the results for Parkinson's dataset

# Run SPECT Dataset
res_sp <- run_experiment("SPECT", spect, reps = reps_spect)
saveRDS(res_sp, file.path(res_dir, paste0("res_spect_", run_date, ".rds")))

# Run Musk (mixed) Dataset
res_mk <- run_experiment("Musk (mixed)", musk_mixed, reps = reps_musk)
saveRDS(res_mk, file.path(res_dir, paste0("res_musk_", run_date, ".rds")))

# ---------- Average/Mean for every method and every rate on all 50 repetitions ----------
res <- rbind(res_pk, res_sp, res_mk) # merge the results of all 3 datasets together
groups <- split(res, list(res$dataset, res$rate, res$method), drop = TRUE) # a list of small tables for every dataset/rate/impuation_method combination that exists
summary_table <- do.call(rbind, lapply(groups, function(d) data.frame( # takes each small table from groups & creates a summary row for it
  dataset = d$dataset[1], rate = d$rate[1], method = d$method[1], # in each small table there is the same dataset, rate, method
  NRMSE_mean = mean(d$NRMSE, na.rm = TRUE), NRMSE_sd = sd(d$NRMSE, na.rm = TRUE), # calculate the mean & std NRMSE of all 50 repetitions in the small table
  PFC_mean   = mean(d$PFC,   na.rm = TRUE), PFC_sd   = sd(d$PFC,   na.rm = TRUE), # calculate the mean & std PFC of all 50 repetitions in the small table
  n_runs = nrow(d)))) # number of rows of each small table, it should be 50 for each combination
rownames(summary_table) <- NULL # reset row labels
print(summary_table) # the table to compare with the paper's results
write.csv(summary_table, file.path(res_dir, paste0("summary_table_", run_date, ".csv")), row.names = FALSE)  # save the final table with results

