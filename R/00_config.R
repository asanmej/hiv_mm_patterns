# -----------------------------------------------------------------------------
# 00_config.R: Global Project Configuration
#
# Description:
#   Single source of truth for packages and PATHS. Every script of the project
#   must build its input/output paths from the PATH_* variables defined here;
#   no script should write file paths by hand.
#
# Folder structure (created automatically):
#
#   Data/                                  (input, read-only)
#   R/                                     (all scripts)
#   Results/
#     Hierarchical/
#       Full_cohort/          one clustering of diseases for the whole cohort
#         <reports .html>     (rendered Rmd 01, 02, 03)
#         Figures/  (OE_heatmaps/<Males|Females>/k_<k>/)
#         Tables/   Models/   Patient_data/
#       Stratified/           clustering repeated in every sex x age x HIV stratum
#         Figures/  (<sex>/<age>/<status>/ and OE_heatmaps/)
#         Tables/   Models/   Patient_data/    (<sex>/<age>/<status>/ inside)
#     LCA/                    Latent Class Analysis (always stratified)
#       Figures/  Tables/  Models/  Patient_data/
#     EFA/                    Exploratory Factor Analysis (always stratified)
#       Figures/  Tables/  Models/  Patient_data/
#
#   To add a new hierarchical attempt (e.g. clustering by sex only) copy the
#   block "Hierarchical / Stratified" below with a new name; the folders are
#   created automatically.
#
# Naming rule:  PATH_<METHOD>[_<ATTEMPT>]_<KIND>
#   KIND = FIGURES | HEATMAPS | TABLES | MODELS | PATIENT_DATA
#     FIGURES      plots of the model itself (scree, tetrachoric, LCA fit...)
#     HEATMAPS     observed/expected heatmaps (same format for all methods)
#     TABLES       tables and disease-to-pattern assignments (csv / xlsx)
#     MODELS       fitted objects (.rds / .RData)
#     PATIENT_DATA patient-level outputs (patient -> pattern/class)
# -----------------------------------------------------------------------------

# 1. Required packages --------------------------------------------------------
pkg <- c(
  "here", "rmarkdown", "knitr", "kableExtra",
  "tidyverse", "data.table",
  "factoextra", "proxy", "dendextend", "cluster", "psych",
  "survival",
  "plotly", "htmlwidgets",
  "writexl"
)

# 2. Check and load the packages -------------------------------------------------
# Packages are NOT installed here: installing while other packages are loaded
# fails on Windows (locked .dll files, "lazy loading failed"). Install them once
# with R/00_install_packages.R in a fresh R session.
missing_pkg <- setdiff(pkg, rownames(installed.packages()))
if (length(missing_pkg) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkg, collapse = ", "), "\n",
    "Run R/00_install_packages.R once, in a fresh R session ",
    "(Session > Restart R, or Rscript from a terminal).",
    call. = FALSE
  )
}
invisible(lapply(pkg, library, character.only = TRUE))
rm(pkg, missing_pkg)

# Patient assignment rules to run and compare (hierarchical and EFA):
#   "count"      pattern where the patient has the most diseases
#   "proportion" pattern whose diseases the patient has in the highest percentage
ASSIGNMENT_RULES <- c("count", "proportion")

# Thresholds shared by the three techniques to decide whether a pattern is worth
# showing to the clinicians (single definition; do not repeat these numbers in
# other scripts):
#   MIN_PATIENTS_CLUSTER  patients assigned to a pattern / cluster / class; below
#                         this the heatmap column is striped (06_heatmap_common.R)
#                         and EFA solutions containing it are not selected
#   MIN_DISEASES_PATTERN  diseases defining a pattern (one disease is not a pattern)
#   EFA_LOADING_CUTOFF    loading from which a disease belongs to an EFA factor
#   K_VIABILITY_RULES     assignment rules under which EVERY EFA pattern must reach
#                         MIN_PATIENTS_CLUSTER for the solution to be selected
MIN_PATIENTS_CLUSTER <- 30
MIN_DISEASES_PATTERN <- 2
EFA_LOADING_CUTOFF   <- 0.20
K_VIABILITY_RULES    <- ASSIGNMENT_RULES

# 3. Paths --------------------------------------------------------------------

## Project root (requires an .Rproj, .git or .here file in the project root)
PATH_PROYECT <- here::here()

PATH_DATA    <- file.path(PATH_PROYECT, "Data")
PATH_R       <- file.path(PATH_PROYECT, "R")
PATH_RESULTS <- file.path(PATH_PROYECT, "Results")

## 3.1 Hierarchical clustering ------------------------------------------------
PATH_HIERARCHICAL <- file.path(PATH_RESULTS, "Hierarchical")

### Attempt 1: full cohort (one clustering of diseases for everybody)
PATH_HIERARCHICAL_FULL              <- file.path(PATH_HIERARCHICAL, "Full_cohort")
PATH_HIERARCHICAL_FULL_FIGURES      <- file.path(PATH_HIERARCHICAL_FULL, "Figures")
PATH_HIERARCHICAL_FULL_HEATMAPS     <- file.path(PATH_HIERARCHICAL_FULL_FIGURES, "OE_heatmaps")
PATH_HIERARCHICAL_FULL_TABLES       <- file.path(PATH_HIERARCHICAL_FULL, "Tables")
PATH_HIERARCHICAL_FULL_MODELS       <- file.path(PATH_HIERARCHICAL_FULL, "Models")
PATH_HIERARCHICAL_FULL_PATIENT_DATA <- file.path(PATH_HIERARCHICAL_FULL, "Patient_data")

### Attempt 2: stratified (clustering repeated in each sex x age x HIV stratum)
PATH_HIERARCHICAL_STRATIFIED              <- file.path(PATH_HIERARCHICAL, "Stratified")
PATH_HIERARCHICAL_STRATIFIED_FIGURES      <- file.path(PATH_HIERARCHICAL_STRATIFIED, "Figures")
PATH_HIERARCHICAL_STRATIFIED_HEATMAPS     <- file.path(PATH_HIERARCHICAL_STRATIFIED_FIGURES, "OE_heatmaps")
PATH_HIERARCHICAL_STRATIFIED_TABLES       <- file.path(PATH_HIERARCHICAL_STRATIFIED, "Tables")
PATH_HIERARCHICAL_STRATIFIED_MODELS       <- file.path(PATH_HIERARCHICAL_STRATIFIED, "Models")
PATH_HIERARCHICAL_STRATIFIED_PATIENT_DATA <- file.path(PATH_HIERARCHICAL_STRATIFIED, "Patient_data")

## 3.2 Latent Class Analysis (LCA) ---------------------------------------------
PATH_LCA              <- file.path(PATH_RESULTS, "LCA")
PATH_LCA_FIGURES      <- file.path(PATH_LCA, "Figures")
PATH_LCA_HEATMAPS     <- file.path(PATH_LCA_FIGURES, "OE_heatmaps")
PATH_LCA_TABLES       <- file.path(PATH_LCA, "Tables")
PATH_LCA_MODELS       <- file.path(PATH_LCA, "Models")
PATH_LCA_PATIENT_DATA <- file.path(PATH_LCA, "Patient_data")

## 3.3 Exploratory Factor Analysis (EFA) ----------------------------------------
PATH_EFA              <- file.path(PATH_RESULTS, "EFA")
PATH_EFA_FIGURES      <- file.path(PATH_EFA, "Figures")
PATH_EFA_HEATMAPS     <- file.path(PATH_EFA_FIGURES, "OE_heatmaps")
PATH_EFA_TABLES       <- file.path(PATH_EFA, "Tables")
PATH_EFA_MODELS       <- file.path(PATH_EFA, "Models")
PATH_EFA_PATIENT_DATA <- file.path(PATH_EFA, "Patient_data")

## 3.4 Output folders are created automatically ----------------------------------
## (every variable PATH_RESULTS / PATH_HIERARCHICAL* / PATH_LCA* / PATH_EFA*
##  is a folder, so a new folder only needs its PATH_ variable above)
output_dirs <- unlist(mget(
  grep("^PATH_(RESULTS|HIERARCHICAL|LCA|EFA)", ls(), value = TRUE)
))
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
rm(output_dirs)

## 3.5 Original data files (.csv) -------------------------------------------------
PATH_BDU_FULL                  <- file.path(PATH_DATA, "bdu_full.csv")
PATH_COHORT_MATCH_ID           <- file.path(PATH_DATA, "cohort_match_id.csv")
PATH_DIAG_BOOL_PREVALENT_FINAL <- file.path(PATH_DATA, "diagnosis_bool_prevalent_final.csv")

## Fail early (and clearly) if an input file is missing
.input_files <- c(PATH_BDU_FULL, PATH_COHORT_MATCH_ID, PATH_DIAG_BOOL_PREVALENT_FINAL)
.missing <- .input_files[!file.exists(.input_files)]
if (length(.missing) > 0) {
  stop("Missing input file(s) in ", PATH_DATA, ":\n  ",
       paste(basename(.missing), collapse = "\n  "), call. = FALSE)
}
rm(.input_files, .missing)
