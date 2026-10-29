# -----------------------------------------------------------------------------
# 05_stratified_heatmaps.R
#
# Description:
#   Generates interactive heatmaps of disease enrichment within
#   patient-level comorbidity patterns, using the disease patterns
#   obtained from the sex × age stratified hierarchical clustering
#   (03_hierarchical_clustering_stratified.R).
#
#   This script is fully self-contained. It does NOT depend on the
#   original patient pattern scoring script (02) nor on the original
#   heatmap script (04). Patient-level pattern scores and the dominant
#   pattern are recalculated internally, using exactly the same logic
#   as in the original scoring script.
#
#   For each combination of:
#     - Sex (MALE, FEMALE)
#     - Age group (under_45, 45to65, over_65)
#     - HIV status (ctl, hiv)
#     - Number of disease clusters (k = 3, 4, 5, 6, 7)
#
#   the script calculates an Observed/Expected (O/E) ratio for each
#   disease and patient cluster:
#
#       O/E ratio =
#         prevalence of the disease within the patient cluster /
#         prevalence of the disease within the corresponding
#         sex × age × HIV status stratum
#
#   Interpretation and colour conventions are identical to the
#   original 04_stratified_heatmaps.R, so the two sets of heatmaps
#   are directly comparable.
# -----------------------------------------------------------------------------

# 1. PACKAGE MANAGEMENT ------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(plotly)
library(htmlwidgets)

source(file.path(PATH_R, "06_heatmap_common.R"))  # shared heatmap functions, colours, statistics

# 2. GLOBAL CONFIGURATION ----------------------------------------------------

# Candidate numbers of clusters (must match the clustering script)
k_candidates <- 3:7

# Input/output folders: PATH_HIERARCHICAL_STRATIFIED_* (see 00_config.R)
# 3. DATA LOADING ------------------------------------------------------------

# Load baseline demographic and clinical data
bdu_data <- read_csv(PATH_BDU_FULL)

# Load patient-level binary diagnostic matrix
diag_data <- read_csv(PATH_DIAG_BOOL_PREVALENT_FINAL)

# Reconstruct the same sex × age × HIV status classification and
# attach the binary disease indicators from the diagnostic matrix.
# Diseases are required both to compute the O/E ratio and to draw
# the structural cluster borders.
cohort_stratified <- bdu_data %>%
  select(patient_id, sexo, edad, vih_bool) %>%
  mutate(
    sexo_clean = case_when(
      str_detect(tolower(sexo), "h") ~ "MALE",
      str_detect(tolower(sexo), "m") ~ "FEMALE",
      TRUE ~ as.character(sexo)
    ),
    grupo_edad = case_when(
      edad < 45 ~ "under_45",
      edad >= 45 & edad <= 65 ~ "45to65",
      edad > 65 ~ "over_65",
      TRUE ~ NA_character_
    ),
    status_hiv = if_else(as.logical(vih_bool), "hiv", "ctl")
  ) %>%
  inner_join(diag_data, by = "patient_id")


# 4. HELPER: PATIENT PATTERN SCORES AND DOMINANT PATTERN ---------------------

# Compute patient-level proportional scores for all patterns of a
# given clustering solution (k_val) within a single sex × age stratum.
# This reproduces the logic of the original scoring script:
#
#   proportion     = count / n_diseases
#   pattern_weight = n_diseases / n_total_diseases
#   weighted_score = proportion * pattern_weight
compute_patient_pattern_scores <- function(diag_stratum, pat, k_val) {
  
  col_pattern <- paste0("cluster_k", k_val)
  
  pat_k <- pat %>%
    select(enfermedad, enfermedad_col, pattern = all_of(col_pattern)) %>%
    mutate(pattern = as.integer(pattern))
  
  n_total_diseases <- nrow(pat_k)
  
  cluster_sizes <- pat_k %>%
    count(pattern, name = "n_diseases") %>%
    mutate(pattern_weight = n_diseases / n_total_diseases)
  
  # Long format: keep only present diseases
  patient_diagnoses_long <- diag_stratum %>%
    select(patient_id, all_of(pat_k$enfermedad_col)) %>%
    pivot_longer(
      cols      = -patient_id,
      names_to  = "enfermedad_col",
      values_to = "has_disease"
    ) %>%
    filter(has_disease == 1) %>%
    select(patient_id, enfermedad_col)
  
  # Count diseases present per patient and pattern
  patient_pattern_counts <- patient_diagnoses_long %>%
    inner_join(
      pat_k %>% select(enfermedad_col, pattern),
      by = "enfermedad_col",
      relationship = "many-to-many"
    ) %>%
    count(patient_id, pattern, name = "count")
  
  # Complete patient × pattern grid
  all_combinations <- expand_grid(
    patient_id = diag_stratum$patient_id,
    pattern    = cluster_sizes$pattern
  )
  
  # Proportional scores
  scores <- all_combinations %>%
    left_join(
      patient_pattern_counts,
      by = c("patient_id", "pattern"),
      relationship = "one-to-one"
    ) %>%
    left_join(
      cluster_sizes,
      by = "pattern",
      relationship = "many-to-one"
    ) %>%
    mutate(
      count          = replace_na(count, 0L),
      proportion     = count / n_diseases,
      weighted_score = proportion * pattern_weight
    )
  
  scores
}

# Derive the dominant pattern per patient using the same rule as the
# original scoring script.
derive_dominant_pattern <- function(scores, rule = "count") {
  # "count": weighted_score (= number of diseases of the pattern the patient has)
  # "proportion": proportion (= % of the pattern's diseases the patient has)
  scores$.score <- if (rule == "count") scores$weighted_score else scores$proportion
  scores %>%
    group_by(patient_id) %>%
    mutate(
      max_score = max(.score, na.rm = TRUE),
      n_max = sum(abs(.score - max_score) < 1e-12 & max_score > 0, na.rm = TRUE),
      dominant_pattern = case_when(
        max_score == 0 ~ NA_integer_,
        n_max > 1      ~ NA_integer_,
        TRUE           ~ pattern[which.max(.score)]
      )
    ) %>%
    ungroup() %>%
    select(patient_id, dominant_pattern) %>%
    distinct()
}

# 5. COLOUR SCALE, STATISTICS AND PLOTTING: see 06_heatmap_common.R

# 6. PRE-COMPUTE PATTERNS AND PATIENT SCORES PER STRATUM ---------------------

# Twelve sex × age × HIV status strata
strata_full <- expand.grid(
  sexo_val   = c("MALE", "FEMALE"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

# Containers for patterns and scores, keyed by stratum name
# (stratum_name now includes the HIV status)
# Solutions chosen in 03_hierarchical_clustering_stratified.R (3 per stratum)
selected_k_file <- file.path(PATH_HIERARCHICAL_STRATIFIED_TABLES, "k_selection_hierarchical.csv")
if (!file.exists(selected_k_file)) stop("Run 03_hierarchical_clustering_stratified.R first: ", selected_k_file)
selected_k <- read_csv(selected_k_file, show_col_types = FALSE) %>% filter(selected)

all_patterns <- list()
all_scores   <- list()

for (i in seq_len(nrow(strata_full))) {
  
  sexo_actual   <- strata_full$sexo_val[i]
  edad_actual   <- strata_full$edad_val[i]
  status_actual <- strata_full$status_val[i]
  stratum_name  <- paste(sexo_actual, edad_actual, status_actual, sep = "_")
  
  # Load disease-to-cluster assignments for the current stratum.
  # The file is now per sex × age × HIV status, because the clustering
  # was estimated separately within each subgroup.
  pattern_file <- file.path(
    PATH_HIERARCHICAL_STRATIFIED_TABLES,
    sexo_actual,
    edad_actual,
    status_actual,
    "jaccard_disease_patterns_multi_k.csv"
  )
  
  if (!file.exists(pattern_file)) {
    message("Missing pattern file for stratum ", stratum_name, ". Skipping.")
    next
  }
  
  patrones_multi_df <- read_csv(pattern_file, show_col_types = FALSE)
  
  # Recover the disease column name used in the diagnostic matrix
  patrones_multi_df <- patrones_multi_df %>%
    mutate(enfermedad_col = gsub(" ", "_", enfermedad))
  
  # Patients belonging to the current sex × age × HIV status stratum
  ids_stratum <- cohort_stratified %>%
    filter(
      sexo_clean == sexo_actual,
      grupo_edad == edad_actual,
      status_hiv == status_actual
    ) %>%
    pull(patient_id)
  
  diag_stratum <- diag_data %>%
    filter(patient_id %in% ids_stratum) %>%
    select(patient_id, all_of(patrones_multi_df$enfermedad_col))
  
  if (nrow(diag_stratum) == 0) {
    message("No patients in stratum ", stratum_name, ". Skipping.")
    next
  }
  
  # Compute patient pattern scores and dominant pattern for each k
  scores_stratum <- purrr::map_dfr(
    selected_k$k[selected_k$sexo_val == sexo_actual & selected_k$edad_val == edad_actual &
                 selected_k$status_val == status_actual],
    function(k_val) {
    
    scores_k   <- compute_patient_pattern_scores(diag_stratum, patrones_multi_df, k_val)
    dominant_k <- purrr::map_dfr(ASSIGNMENT_RULES, function(r) {
      derive_dominant_pattern(scores_k, r) %>% mutate(rule = r)
    }) %>%
      pivot_wider(names_from = rule, values_from = dominant_pattern,
                  names_prefix = "dominant_pattern_")
    
    scores_k %>%
      left_join(dominant_k, by = "patient_id") %>%
      mutate(
        sexo_val   = sexo_actual,
        edad_val   = edad_actual,
        status_val = status_actual,
        k          = k_val
      )
  })
  
  all_patterns[[stratum_name]] <- patrones_multi_df
  all_scores[[stratum_name]]   <- scores_stratum
}

# 7. STRATA DEFINITION FOR HEATMAPS ------------------------------------------

# One heatmap per sex × age × HIV status × k combination
estratos <- selected_k %>%
  transmute(sexo_val, edad_val, status_val, k_val = k) %>%
  tidyr::crossing(rule = ASSIGNMENT_RULES)

# 8. MAIN HEATMAP LOOP -------------------------------------------------------

rule_results <- purrr::pmap(estratos, function(sexo_val, edad_val, status_val, k_val, rule) {
  
  sexo_sel   <- sexo_val
  edad_sel   <- edad_val
  status_sel <- status_val
  k_sel      <- k_val
  rule_sel   <- rule
  
  stratum_name <- paste(sexo_sel, edad_sel, status_sel, sep = "_")
  
  patrones_multi_df <- all_patterns[[stratum_name]]
  scores_stratum    <- all_scores[[stratum_name]]
  
  # 8.1 Assign each patient to the dominant pattern for the current k
  patient_clusters <- scores_stratum %>%
    filter(k == k_sel, !is.na(.data[[paste0("dominant_pattern_", rule_sel)]])) %>%
    transmute(patient_id, cluster = as.character(.data[[paste0("dominant_pattern_", rule_sel)]])) %>%
    distinct(patient_id, .keep_all = TRUE)
  
  if (nrow(patient_clusters) == 0) {
    message(
      "No patient assignments for ",
      sexo_sel, " | ", edad_sel, " | k=", k_sel
    )
    return(NULL)
  }
  
  # 8.2 Retrieve the disease-to-cluster assignments corresponding to
  # the current hierarchical clustering solution
  col_cluster_name <- paste0("cluster_k", k_sel)
  
  patrones_subset <- patrones_multi_df %>%
    select(
      enfermedad,
      enfermedad_col,
      cluster_enfermedad = all_of(col_cluster_name)
    ) %>%
    mutate(cluster_enfermedad = as.character(cluster_enfermedad))
  
  # 8.3 Check that all diseases are available in the diagnostic matrix
  disease_columns <- patrones_subset$enfermedad_col[
    patrones_subset$enfermedad_col %in% colnames(cohort_stratified)
  ]
  
  if (length(disease_columns) == 0) {
    message(
      "No disease columns were found for: ",
      sexo_sel, " | ", edad_sel, " | ", status_sel, " | k=", k_sel
    )
    return(NULL)
  }
  
  patrones_subset <- patrones_subset %>%
    filter(enfermedad_col %in% disease_columns)
  
  # 8.4 Restrict the cohort to the current sex × age × HIV status stratum
  # and assign each patient to their dominant patient-level pattern
  cohort_stratum <- cohort_stratified %>%
    filter(
      sexo_clean == sexo_sel,
      grupo_edad == edad_sel,
      status_hiv == status_sel
    ) %>%
    left_join(patient_clusters, by = "patient_id")  # unassigned patients stay in the stratum (cluster = NA)
  
  if (nrow(cohort_stratum) == 0) {
    message(
      "No patients for ",
      sexo_sel, " | ", edad_sel, " | ", status_sel, " | k=", k_sel
    )
    return(NULL)
  }
  
  hm_files <- oe_heatmap_files(file.path(PATH_HIERARCHICAL_STRATIFIED_HEATMAPS, rule_sel),
                               paste0("hierarchical_stratified_", rule_sel),
                               sexo_sel, edad_sel, status_sel, k_sel)
  disease_map_plot <- patrones_subset %>%
    transmute(enfermedad_col, enfermedad, cluster_enfermedad, sign = "pos")

  res <- make_oe_heatmap(
    cohort_stratum = cohort_stratum,
    disease_map    = disease_map_plot,
    clusters       = as.character(seq_len(k_sel)),
    title_core     = paste0("Observed/Expected Ratio (k=", k_sel, ") - ",
                            sexo_sel, " | Age: ", edad_sel,
                            " | Status: ", toupper(status_sel),
                            " | Assign: ", rule_sel),
    out_html = hm_files$html,
    out_csv  = hm_files$csv
  )
  if (is.null(res)) return(NULL)
  list(sexo_val = sexo_sel, edad_val = edad_sel, status_val = status_sel,
       k = k_sel, rule = rule_sel, table = res$table,
       assignment = cohort_stratum[, c("patient_id", "cluster")])
})

# Sensitivity of the results to the assignment rule (count vs proportion)
compare_assignment_rules(
  rule_results,
  file.path(PATH_HIERARCHICAL_STRATIFIED_TABLES, "assignment_rule_comparison.csv")
)

# 9. EXPORT PATIENT-LEVEL PATTERN SCORES -------------------------------------

# Export the per-stratum patient scores as a reproducible intermediate
# file inside PATH_HIERARCHICAL_STRATIFIED_PATIENT_DATA.
patient_pattern_scores_stratified <- bind_rows(all_scores) %>%
  select(
    patient_id,
    sexo_val,
    edad_val,
    status_val,
    k,
    pattern,
    count,
    n_diseases,
    proportion,
    weighted_score,
    starts_with("dominant_pattern_")
  ) %>%
  arrange(sexo_val, edad_val, status_val, k, patient_id, pattern)

write_csv(
  patient_pattern_scores_stratified,
  file.path(
    PATH_HIERARCHICAL_STRATIFIED_PATIENT_DATA,
    "patient_pattern_scores_stratified_selected_k.csv"
  )
)