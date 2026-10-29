# -----------------------------------------------------------------------------
# MAIN PIPELINE EXECUTION
#
# Objective: run the whole comorbidity-pattern analysis (HIV cohort vs matched
#            controls) with four approaches that share the same heatmap format:
#              1. Hierarchical clustering, full cohort
#              2. Hierarchical clustering, stratified (sex x age x HIV)
#              3. Latent Class Analysis (LCA)
#              4. Exploratory Factor Analysis (EFA)
#
# How to run: open the project (.Rproj) and run this file. Every step can be
#             switched on/off in RUN below (e.g. to repeat only the LCA).
#             All paths come from 00_config.R; results go to Results/<method>/.
# -----------------------------------------------------------------------------

# 1. Project and environment configuration ------------------------------------
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
source(here::here("R", "00_config.R"))

# 2. Steps to run -------------------------------------------------------------
RUN <- list(
  hierarchical_full       = TRUE,
  hierarchical_stratified = TRUE,
  lca                     = TRUE,
  efa                     = TRUE
)

# Run a script / report and stop with a clear message if it fails
run_step <- function(label, expr) {
  message("\n", strrep("=", 70), "\n>>> ", label, "\n", strrep("=", 70))
  t0 <- Sys.time()
  tryCatch(
    force(expr),
    error = function(e) stop("Step failed: ", label, "\n", conditionMessage(e), call. = FALSE)
  )
  message("<<< ", label, " finished in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}
r_script <- function(f) source(file.path(PATH_R, f), local = FALSE)
r_report <- function(f, out_dir) {
  rmarkdown::render(file.path(PATH_R, f), output_dir = out_dir, envir = globalenv(), quiet = TRUE)
}

# 3. HIERARCHICAL CLUSTERING - FULL COHORT -> Results/Hierarchical/Full_cohort --
if (RUN$hierarchical_full) {
  run_step("Hierarchical full: clustering of diseases (Jaccard / Ward.D2)",
           r_report("01_hierarchical_clustering.Rmd", PATH_HIERARCHICAL_FULL))
  run_step("Hierarchical full: patient pattern scores",
           r_report("02_patient_pattern_scoring.Rmd", PATH_HIERARCHICAL_FULL))
  run_step("Hierarchical full: epidemiological validation (clogit)",
           r_report("03_epidemiological_validation.Rmd", PATH_HIERARCHICAL_FULL))
  run_step("Hierarchical full: O/E heatmaps", r_script("04_stratified_heatmaps.R"))
}

# 4. HIERARCHICAL CLUSTERING - STRATIFIED -> Results/Hierarchical/Stratified ----
if (RUN$hierarchical_stratified) {
  run_step("Hierarchical stratified: clustering in each stratum",
           r_script("03_hierarchical_clustering_stratified.R"))
  run_step("Hierarchical stratified: patient scores and O/E heatmaps",
           r_script("05_stratified_heatmaps.R"))
}

# 5. LATENT CLASS ANALYSIS -> Results/LCA ---------------------------------------
if (RUN$lca) {
  run_step("LCA: models by stratum", r_script("LCA_HIV.R"))
  run_step("LCA: O/E heatmaps",      r_script("04b_stratified_heatmaps_LCA.R"))
}

# 6. EXPLORATORY FACTOR ANALYSIS -> Results/EFA ---------------------------------
if (RUN$efa) {
  run_step("EFA: factor models by stratum", r_script("3_Exploratory_Factor_Analysis.R"))
  run_step("EFA: O/E heatmaps",             r_script("07_efa_heatmaps.R"))
}

message("\nPipeline finished. Results in: ", PATH_RESULTS)
