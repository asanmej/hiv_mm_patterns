# 00_config.R: 

# 1. List of required packages 
pkg <- c(
  
  ## Directory management
  "here", 
  
  # 1. Hierarchical clustering analysis
  ## For data manipulation and visualization
  "tidyverse",
  ## For efficient Jaccard distance calculation
  "proxy", 
  ## To customize and view dendrograms
  "dendextend",
  ## To evaluate the average silhouette width 
  "cluster",  
  ## Conditional logistic regression (for 03_validation)
  "survival" 
  
  # 2. Data visualization
 
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

PATH_DATA <- file.path(PATH_PROYECT, "Datos")
PATH_R <- file.path(PATH_PROYECT, "R")
PATH_RESULTS <- file.path(PATH_PROYECT, "Results")

## Intermediatr data
PATH_INTERMEDIATE <- file.path(PATH_PROYECT, "Datos_transformados")

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
