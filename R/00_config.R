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

PATH_INTERMEDIATE <- file.path(PATH_PROYECT, "Intermediate_data")

## Output folder are created automatically
dirs <- c(
  PATH_RESULTS,
  PATH_INTERMEDIATE
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
