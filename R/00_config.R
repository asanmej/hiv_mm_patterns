# -----------------------------------------------------------------------------
# 00_config.R: Global Project Configuration
#
# Description:
# This script initializes the working environment for the project.

# It handles automated dependency management (installation and loading of 
# required libraries), environment cleanup, standardized definition of relative 
# paths using the 'here' package, and automatic creation of the output directory 
# structure necessary for the entire pipeline.
# ----------------------------------------------------------------------------- 

# 1. List of required packages 
pkg <- c(
  
  ## Directory management
  "here", 
  ## To save the HTML in a different directory
  "rmarkdown",
  
  ## For data manipulation and visualization
  "tidyverse",
  ## To calculate the Hopkins Statistic
  "factoextra",
  ## For efficient Jaccard distance calculation
  "proxy", 
  ## To customize and view dendrograms
  "dendextend",
  ## To evaluate the average silhouette width 
  "cluster",  
  ## Conditional logistic regression
  "survival",
  
  ## Interactive data visualization
  "plotly",
  ## Saving interactive plots as standalone HTML widgets
  "htmlwidgets"
)

# 2. Automatically load and install missing packages
lapply(pkg, function (x){
  if(!require(x, character.only = T)){
    install.packages(x , character.only = T)
  }
})

# 3. Clean the environment by deleting the temporary variable
rm(pkg)  

# 4. Stablish relative paths required for the execution  

## Directories
PATH_PROYECT <- here()

PATH_DATA <- file.path(PATH_PROYECT, "Data")
PATH_R <- file.path(PATH_PROYECT, "R")

PATH_RESULTS <- file.path(PATH_PROYECT, "Results")
PATH_RESULTS_TABLES <- file.path(PATH_RESULTS, "Tables")
PATH_RESULTS_FIGURES <- file.path(PATH_RESULTS, "Figures")

PATH_RESULTS_LCA <- file.path(PATH_RESULTS, "LCA")
PATH_LCA_STRATIFIED <- file.path(PATH_RESULTS_LCA, "stratified")
PATH_LCA_STRATIFIED_FIGURES <- file.path(PATH_LCA_STRATIFIED, "Figures")
PATH_LCA_STRATIFIED_TABLES <- file.path(PATH_LCA_STRATIFIED, "Tables")
PATH_LCA_STRATIFIED_MODELS <- file.path(PATH_LCA_STRATIFIED, "Models")
PATH_LCA_STRATIFIED_PATIENT_DATA <- file.path(PATH_LCA_STRATIFIED,"Patient_data")

PATH_INTERMEDIATE <- file.path(PATH_PROYECT, "Intermediate_data")
PATH_INTERMEDIATE_LCA <- file.path(PATH_INTERMEDIATE, "LCA")

## Output folder are created automatically
dirs <- c(
  PATH_RESULTS,
  PATH_RESULTS_LCA,
  PATH_RESULTS_TABLES,
  PATH_RESULTS_FIGURES,
  PATH_LCA_STRATIFIED,
  PATH_LCA_STRATIFIED_FIGURES,
  PATH_LCA_STRATIFIED_TABLES,
  PATH_LCA_STRATIFIED_MODELS,
  PATH_LCA_STRATIFIED_PATIENT_DATA,
  
  PATH_INTERMEDIATE,
  PATH_INTERMEDIATE_LCA
)

invisible(
  lapply(
    dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

## Original data (.csv)
PATH_BDU_FULL <- file.path(PATH_DATA, "bdu_full.csv")
PATH_COHORT_MATCH_ID <- file.path(PATH_DATA, "cohort_match_id.csv")
PATH_DIAG_BOOL_PREVALENT_FINAL <- file.path(PATH_DATA, "diagnosis_bool_prevalent_final.csv")
