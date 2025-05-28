# create_fit_predict_hmsc_models_cv

## Description

This function performs cross-validation on a Hierarchical Model of Species Communities (HMSC) using High Performance Computing (HPC). It fits the model on training data and makes predictions on test data for each fold of the cross-validation.

## Usage

```R
create_fit_predict_hmsc_models_cv(
  hM,
  n_folds = 4,
  seed = 123,
  save_dir = getwd(),
  nSamples = NULL,
  thin = NULL,
  nChains = NULL,
  verbose = NULL
)
```

## Arguments

- `hM`: An HMSC model object.
- `n_folds`: Number of folds for cross-validation (default: 4).
- `seed`: Random seed for reproducibility (default: 123).
- `save_dir`: Directory to save output files (default: current working directory).
- `nSamples`: Number of MCMC samples (optional, default: uses `hM$samples`).
- `thin`: Thinning interval for MCMC (optional, default: uses `hM$thin`).
- `nChains`: Number of MCMC chains (optional, default: uses `length(hM$postList)`).
- `verbose`: Verbosity level (optional, default: uses `hM$verbose` or 10 if not set).

## Return Value

A full prediction array combining predictions from all folds.

## Details

The function performs the following steps:

1. Initializes a full prediction array to store results from all folds.
2. Generates fold indices for cross-validation.
3. For each fold:
   a. Prepares training and testing data.
   b. Creates a new HMSC model for the fold.
   c. Initializes MCMC sampling for HPC.
   d. Executes HPC Gibbs sampler using a Python script.
   e. Imports and processes HPC results.
   f. Computes predictions for test data.
   g. Stores predictions in the full prediction array.
   h. Saves fold-specific results (initialization, fitted model, predictions).
4. Saves the full prediction array containing results from all folds.

## Notes

- This function requires a properly set up Conda/pip environment with TensorFlow and other necessary Python dependencies.
- The function uses both R and Python, leveraging HPC capabilities for MCMC sampling.
- Results are saved in a structured directory within the specified `save_dir`.

## Dependencies

- R packages: Hmsc, jsonify, cli, tidyverse, caret, reticulate
- Python environment with TensorFlow and HMSC dependencies

## Example

```R
# Assuming 'hM' is your initial Hmsc model
result <- create_fit_predict_hmsc_models_cv(
  hM,
  n_folds = 4,
  save_dir = "path/to/save/directory"
)
```

This will perform 4-fold cross-validation on the HMSC model, saving results in the specified directory and returning the full prediction array.

## Function
```R
# Load required packages
require(Hmsc)
require(jsonify)
require(cli)
require(tidyverse)
require(caret)
require(reticulate)

#' Create Train-Test Indices for Cross-Validation
#'
#' This function generates indices for train and test sets for k-fold cross-validation.
#'
#' @param n_samples Number of samples in the dataset
#' @param n_folds Number of folds for cross-validation
#' @param seed Random seed for reproducibility
#' @return A list of lists containing train and test indices for each fold
create_train_test_indices <- function(n_samples, n_folds = 4, seed = NULL) {
  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Create folds using caret's createFolds function
  folds <- createFolds(1:n_samples, k = n_folds, list = TRUE, returnTrain = TRUE)
  
  # Initialize list to store train and test indices for each fold
  train_test_indices <- list()
  
  # For each fold, create a list with train and test indices
  for (i in 1:n_folds) {
    train_ind <- folds[[i]]
    # Test indices are all samples not in the training set
    test_ind <- setdiff(1:n_samples, train_ind)
    
    train_test_indices[[i]] <- list(
      train = train_ind,
      test = test_ind
    )
  }
  
  return(train_test_indices)
}

#' Fit and Predict HMSC Models with Cross-Validation
#'
#' This function performs cross-validation on an HMSC model, fitting the model
#' on training data and making predictions on test data for each fold. It uses
#' High Performance Computing (HPC) for model fitting.
#'
#' @param hM An HMSC model object
#' @param n_folds Number of folds for cross-validation
#' @param seed Random seed for reproducibility
#' @param save_dir Directory to save output files
#' @param nSamples Number of MCMC samples (optional)
#' @param thin Thinning interval for MCMC (optional)
#' @param nChains Number of MCMC chains (optional)
#' @param verbose Verbosity level (optional)
#' @return A full prediction array combining predictions from all folds
#'
#' @importFrom reticulate use_condaenv
#' @importFrom cli cli_h1 cli_alert_info cli_alert_success cli_rule
create_fit_predict_hmsc_models_cv <- function(hM,
                                              n_folds = 4,
                                              seed = 123,
                                              save_dir = getwd(),
                                              nSamples = NULL,
                                              thin = NULL,
                                              nChains = NULL,
                                              verbose = NULL) {
  set.seed(seed)
  
  # Helper function to initialize the prediction array
  initialize_prediction_array <- function(hM) {
    nSamples <- sum(sapply(hM$postList, length))
    array(NA, c(hM$ny, hM$ns, nSamples))
  }
  
  # Initialize the full prediction array to store results from all folds
  full_pred_array <- initialize_prediction_array(hM)
  
  # Set default values for MCMC parameters if not provided
  # Use the %||% operator to assign a default value to nSamples if it is NULL.
  # The %||% operator checks if the left-hand side (nSamples) is NULL; 
  # if so, it assigns the right-hand side (hM$samples). Otherwise, nSamples retains its value.
  # This is a concise way to provide a fallback value without using an explicit if statement.
  nSamples <- nSamples %||% hM$samples
  thin <- thin %||% hM$thin
  nChains <- nChains %||% length(hM$postList)
  verbose <- verbose %||% if(is.null(hM$verbose)) 10 else hM$verbose
  
  # # Best to make sure you've set up your Conda/pip venv before running this funciton like below
  # # Set up conda environment
  # use_condaenv("hmsc-hpc")
  # python <- "/home/dansmi/miniconda3/envs/hmsc-hpc/bin/python"
  
  # # Set environment variable to reduce TensorFlow debug output
  # Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3")
  
  # Python code to set up TensorFlow with additional error suppression
  #   tf_setup_code <- "
  # import os
  # import tensorflow as tf
  # 
  # # Suppress TensorFlow logging - Constant NUMA warnings on Ubuntu
  # # https://stackoverflow.com/questions/68128979/tensorflow-repeated-success-messages-and-numa-node-read-warning
  # os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
  # tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
  # 
  # # Check TensorFlow installation
  # print(tf.constant(1))
  # 
  # # Additional debugging information
  # print('TensorFlow version:', tf.__version__)
  # print('GPU available:', tf.test.is_gpu_available())
  # "
  
  # # Execute Python code
  # system2(python, c("-c", shQuote(tf_setup_code)))
  
  # Create folder structure for cross-validation results
  cv_dir <- file.path(save_dir, "CV")
  sapply(c("Init", "Fitted", "Predictions"), function(subdir) {
    dir.create(file.path(cv_dir, subdir), recursive = TRUE, showWarnings = FALSE)
  })

  # Generate fold indices for cross-validation
  fold_ind <- create_train_test_indices(nrow(hM$Y), n_folds = n_folds)
  
  # Loop through each fold
  for (i in 1:n_folds) {
    cli_h1(sprintf("Processing Fold %d of %d", i, n_folds))
    
    # Get train and test indices for this fold
    train <- fold_ind[[i]]$train
    test <- fold_ind[[i]]$test
    
    # Prepare training data for this fold
    Y <- hM$Y[train, , drop = FALSE]
    XData <- hM$XData[train, , drop = FALSE]
    Tr <- hM$Tr
    Phy <- hM$phyloTree
    studyDesign <- droplevels(hM$dfPi[train, , drop = FALSE])
    rLsite <- hM$rL
    distr <- hM$distr
    
    # Create study design data frame for training
    studyDesign <- data.frame(
      sample = factor(paste0("site_", 1:nrow(Y))),
      site = factor(1:nrow(Y))
    )
    rownames(studyDesign) <- studyDesign$sample
    
    # Create random level for the HMSC model
    rlSite <- HmscRandomLevel(units = studyDesign$site)
    
    # Create new HMSC model for this fold
    hMk <- Hmsc(
      Y = Y,
      XData = XData,
      Tr = Tr,
      phyloTree = Phy,
      studyDesign = studyDesign,
      ranLevels = list(site = rlSite),
      distr = distr
    )
    
    # Initialize MCMC sampling for HPC
    cli_alert_info("Initializing for HPC")
    transient <- ceiling(nSamples * thin / 10)
    init_obj <- sampleMcmc(hMk, samples = nSamples, thin = thin,
                           transient = transient, nChains = nChains,
                           verbose = verbose, engine = "HPC")
    
    # Save initialization object
    init_file <- file.path(cv_dir, "Init", sprintf("init_S%d_T%d_Fold%d.rds", nSamples, thin, i))
    saveRDS(to_json(init_obj), file = init_file)
    cli_alert_success("Initialization complete")
    
    # Execute HPC Gibbs sampler
    cli_alert_info("Executing HPC Gibbs sampler")
    post_file <- file.path(cv_dir, "Fitted", sprintf("post_S%d_T%d_Fold%d.rds", nSamples, thin, i))
    
    # Prepare command line arguments for Python script
    python_cmd_args <- paste("-m hmsc.run_gibbs_sampler",
                             "--input", shQuote(init_file),
                             "--output", shQuote(post_file),
                             "--samples", nSamples,
                             "--transient", transient,
                             "--thin", thin,
                             "--verbose", verbose)
    
    # Run Python script for HPC Gibbs sampling
    system2(python, python_cmd_args)
    cli_alert_success("HPC execution complete")
    
    # Import and process HPC results
    cli_alert_info("Importing HPC results")
    importFromHPC <- from_json(readRDS(file = post_file)[[1]])
    postList <- importFromHPC[1:nChains]
    cli_alert_success(sprintf("Fitting time for S%d_T%d: %.1f sec", nSamples, thin, importFromHPC[[nChains + 1]]))
    
    # Post-process HPC results
    cli_alert_info("Post-processing HPC results")
    fitted_model <- importPosteriorFromHPC(hMk, postList, nSamples, thin, transient)
    cli_alert_success("Post-processing complete")
    
    # Save fitted model
    cli_alert_info("Saving fitted model")
    saveRDS(fitted_model, file.path(cv_dir, "Fitted", sprintf("fitted_S%d_T%d_Fold%d.rds", nSamples, thin, i)))
    cli_alert_success("Fitted model saved")
    
    # Prepare testing data
    testing <- list(
      Y = hM$Y[test, , drop = FALSE],
      XData = hM$XData[test, , drop = FALSE],
      studyDesign = droplevels(hM$dfPi[test, , drop = FALSE])
    )
    testing$studyDesign$site <- factor(1:nrow(testing$Y))
    rownames(testing$studyDesign) <- paste0("site_", 1:nrow(testing$Y))
    
    # Compute predictions for test data
    cli_alert_info("Computing predictions")
    preds <- predict(fitted_model, 
                     XData = testing$XData,
                     studyDesign = testing$studyDesign)
    # Ensure predictions are in array format
    preds <- simplify2array(preds)
    
    # Store predictions in the full prediction array
    full_pred_array[test, , ] <- preds
    
    # Save predictions for this fold
    saveRDS(preds, file.path(cv_dir, "Predictions", sprintf("preds_S%d_T%d_Fold%d.rds", nSamples, thin, i)))
    cli_alert_success("Predictions computed and saved")
    
    cli_alert_success(sprintf("Processing complete for Fold %d", i))
    cli_rule()
  }
  
  # Save the full prediction array containing results from all folds
  # saveRDS(full_pred_array, file.path(cv_dir, "full_prediction_array.rds"))
  saveRDS(full_pred_array, file.path(cv_dir, sprintf("PredArray_S%d_T%d.rds", nSamples, thin)))
  cli_alert_success("Full prediction array saved")
  
  cli_alert_success("HMSC model fitting and prediction completed for all folds.")
  
  return(full_pred_array)
}

# Usage example:
# Assuming 'hM' is your initial Hmsc model
# result <- create_fit_predict_hmsc_models_cv(hM,
#                                             n_folds = 4,
#                                             save_dir = "path/to/save/directory")

hM <- readRDS("../sediment/18S/fitTF.rds")

result <- create_fit_predict_hmsc_models_cv(hM,
                                            n_folds = 2,
                                            save_dir = "sediment/18S",
                                            nSamples = 250,
                                            thin = 1000,
                                            nChains = 4,
                                            verbose = 100)
```
