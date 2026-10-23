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
#     - Sex (MEN, WOMEN)
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

# 2. GLOBAL CONFIGURATION ----------------------------------------------------

# Candidate numbers of clusters (must match the clustering script)
k_candidates <- 3:7

# Root directory for stratified clustering outputs
PATH_HCLUST_STRATIFIED <- file.path(
  PATH_INTERMEDIATE,
  "stratified_hierarchical"
)

# Root output directory for the stratified heatmaps
PATH_FIGURES_STRATIFIED <- file.path(
  PATH_RESULTS,
  "stratified_hierarchical",
  "figures"
)

dir.create(
  PATH_FIGURES_STRATIFIED,
  recursive = TRUE,
  showWarnings = FALSE
)

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
      str_detect(tolower(sexo), "h") ~ "MEN",
      str_detect(tolower(sexo), "m") ~ "WOMEN",
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
derive_dominant_pattern <- function(scores) {
  
  scores %>%
    group_by(patient_id) %>%
    mutate(
      max_score = max(weighted_score, na.rm = TRUE),
      n_max = sum(
        weighted_score == max_score & max_score > 0,
        na.rm = TRUE
      ),
      dominant_pattern = case_when(
        max_score == 0 ~ NA_integer_,
        n_max > 1      ~ NA_integer_,
        TRUE           ~ pattern[which.max(weighted_score)]
      )
    ) %>%
    ungroup() %>%
    select(patient_id, dominant_pattern) %>%
    distinct()
}

# 5. FIXED OBSERVED/EXPECTED COLOUR SCALE ------------------------------------

color_0_05  <- "#FFFFFF"
color_05_09 <- "#E8EEF8"
color_09_11 <- "#C7CCD4"
color_11_2  <- "#E89BB7"
color_2_3   <- "#F2D35E"
color_3_4   <- "#78B87A"
color_4_5   <- "#E99A4A"
color_gt5   <- "#B83A4A"

Observed_Expected_colorscale <- list(
  c(0.00 / 5.5, color_0_05),
  c(0.50 / 5.5, color_0_05),
  c(0.50 / 5.5, color_05_09),
  c(0.90 / 5.5, color_05_09),
  c(0.90 / 5.5, color_09_11),
  c(1.10 / 5.5, color_09_11),
  c(1.10 / 5.5, color_11_2),
  c(2.00 / 5.5, color_11_2),
  c(2.00 / 5.5, color_2_3),
  c(3.00 / 5.5, color_2_3),
  c(3.00 / 5.5, color_3_4),
  c(4.00 / 5.5, color_3_4),
  c(4.00 / 5.5, color_4_5),
  c(5.00 / 5.5, color_4_5),
  c(5.00 / 5.5, color_gt5),
  c(5.50 / 5.5, color_gt5)
)

# 6. PRE-COMPUTE PATTERNS AND PATIENT SCORES PER STRATUM ---------------------

# Twelve sex × age × HIV status strata
strata_full <- expand.grid(
  sexo_val   = c("MEN", "WOMEN"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

# Containers for patterns and scores, keyed by stratum name
# (stratum_name now includes the HIV status)
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
    PATH_HCLUST_STRATIFIED,
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
  scores_stratum <- purrr::map_dfr(k_candidates, function(k_val) {
    
    scores_k   <- compute_patient_pattern_scores(diag_stratum, patrones_multi_df, k_val)
    dominant_k <- derive_dominant_pattern(scores_k)
    
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
estratos <- expand.grid(
  sexo_val   = c("MEN", "WOMEN"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  k_val      = k_candidates,
  stringsAsFactors = FALSE
)

# 8. MAIN HEATMAP LOOP -------------------------------------------------------

purrr::pwalk(estratos, function(sexo_val, edad_val, status_val, k_val) {
  
  sexo_sel   <- sexo_val
  edad_sel   <- edad_val
  status_sel <- status_val
  k_sel      <- k_val
  
  stratum_name <- paste(sexo_sel, edad_sel, status_sel, sep = "_")
  
  patrones_multi_df <- all_patterns[[stratum_name]]
  scores_stratum    <- all_scores[[stratum_name]]
  
  # 8.1 Assign each patient to the dominant pattern for the current k
  patient_clusters <- scores_stratum %>%
    filter(k == k_sel, !is.na(dominant_pattern)) %>%
    select(patient_id, cluster = dominant_pattern) %>%
    mutate(cluster = as.character(cluster)) %>%
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
    inner_join(patient_clusters, by = "patient_id")
  
  if (nrow(cohort_stratum) == 0) {
    message(
      "No patients for ",
      sexo_sel, " | ", edad_sel, " | ", status_sel, " | k=", k_sel
    )
    return(NULL)
  }
  
  # 8.5 Long format: one row per patient-disease combination
  diag_long <- cohort_stratum %>%
    select(
      patient_id,
      cluster,
      all_of(disease_columns)
    ) %>%
    pivot_longer(
      cols      = all_of(disease_columns),
      names_to  = "enfermedad_col",
      values_to = "has_disease"
    ) %>%
    mutate(
      has_disease = if_else(
        is.na(has_disease) | has_disease == "0",
        0,
        1
      ),
      has_disease = as.numeric(has_disease)
    ) %>%
    inner_join(
      patrones_subset %>% select(enfermedad_col, enfermedad),
      by = "enfermedad_col"
    )
  
  # 8.6 Expected prevalence within the current stratum
  global_prev <- diag_long %>%
    group_by(enfermedad) %>%
    summarise(
      expected_prev = mean(has_disease, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(expected_prev > 0)
  
  # 8.7 Observed prevalence within each patient cluster
  cluster_prev <- diag_long %>%
    group_by(enfermedad, cluster) %>%
    summarise(
      observed_prev = mean(has_disease, na.rm = TRUE),
      n_patients    = n(),
      n_cases       = sum(has_disease, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 8.8 Observed/Expected enrichment ratio
  heatmap_matrix <- cluster_prev %>%
    left_join(global_prev, by = "enfermedad") %>%
    mutate(
      Observed_Expected_ratio = case_when(
        expected_prev > 0 ~ observed_prev / expected_prev,
        TRUE ~ NA_real_
      )
    ) %>%
    select(
      enfermedad,
      cluster,
      value = Observed_Expected_ratio,
      observed_prev,
      expected_prev,
      n_patients,
      n_cases
    )
  
  # 8.9 Complete the disease × patient-cluster grid (keep NA)
  clusters_expected <- as.character(seq_len(k_sel))
  
  heatmap_matrix <- heatmap_matrix %>%
    complete(
      enfermedad = unique(global_prev$enfermedad),
      cluster    = clusters_expected
    )
  
  # 8.10 Add the structural disease cluster for the black border
  disease_cluster_map <- patrones_subset %>%
    select(enfermedad, cluster_enfermedad) %>%
    distinct()
  
  heatmap_matrix <- heatmap_matrix %>%
    left_join(disease_cluster_map, by = "enfermedad")
  
  # 8.11 Order diseases alphabetically
  sorted_diseases <- sort(unique(heatmap_matrix$enfermedad))
  
  heatmap_matrix <- heatmap_matrix %>%
    mutate(
      cluster_num = match(cluster, clusters_expected),
      disease_num = match(enfermedad, sorted_diseases)
    )
  
  # 8.12 Highlight cells structurally assigned to the same cluster
  highlighted_cells <- heatmap_matrix %>%
    filter(
      !is.na(cluster_enfermedad),
      cluster == cluster_enfermedad
    ) %>%
    distinct(enfermedad, disease_num, cluster, cluster_num)
  
  highlight_shapes <- list()
  if (nrow(highlighted_cells) > 0) {
    highlight_shapes <- lapply(seq_len(nrow(highlighted_cells)), function(j) {
      list(
        type      = "rect",
        xref      = "x",
        yref      = "y",
        x0        = highlighted_cells$cluster_num[j] - 0.5,
        x1        = highlighted_cells$cluster_num[j] + 0.5,
        y0        = highlighted_cells$disease_num[j] - 0.46,
        y1        = highlighted_cells$disease_num[j] + 0.46,
        line      = list(color = "#000000", width = 3),
        fillcolor = "rgba(0,0,0,0)",
        layer     = "above"
      )
    })
  }
  
  # 8.13 Skip if no valid ratio is available
  max_ratio <- max(heatmap_matrix$value, na.rm = TRUE)
  if (!is.finite(max_ratio)) {
    message(
      "No valid ratios for: ",
      sexo_sel, " | ", edad_sel, " | ", status_sel, " | k=", k_sel
    )
    return(NULL)
  }
  
  # 8.14 Build a unique identifier and title for the current heatmap
  tag_name <- paste0(
    sexo_sel, "_", edad_sel, "_", status_sel, "_k", k_sel
  )
  
  file_title <- paste(
    "Observed/Expected Ratio (k=", k_sel, ") -",
    sexo_sel, "| Age:", edad_sel,
    "| Status:", toupper(status_sel)
  )
  
  # 8.15 Format hover values
  heatmap_matrix <- heatmap_matrix %>%
    mutate(
      ratio_text = if_else(
        is.na(value),
        "NA",
        paste0(round(value, 2))
      )
    )
  
  # 8.16 Build the interactive heatmap
  p <- plot_ly(
    data        = heatmap_matrix,
    x           = ~cluster_num,
    y           = ~disease_num,
    z           = ~value,
    zmin        = 0,
    zmax        = 5.5,
    type        = "heatmap",
    colorscale  = Observed_Expected_colorscale,
    hoverongaps = FALSE,
    text = ~paste0(
      "<b>", enfermedad, "</b>",
      "<br><br>",
      "<b>Observed/Expected ratio:</b> ", ratio_text,
      "<br><b>Observed prevalence:</b> ", round(observed_prev, 4),
      "<br><b>Expected prevalence:</b> ", round(expected_prev, 4),
      "<br><b>N patients:</b> ", n_patients,
      "<br><b>Cases:</b> ", n_cases,
      "<br><b>Disease cluster:</b> ",
      ifelse(is.na(cluster_enfermedad), "NA", cluster_enfermedad)
    ),
    hovertemplate = paste0(
      "%{text}",
      "<br><br>",
      "<b>Patient cluster:</b> Cl. %{x}",
      "<extra></extra>"
    ),
    colorbar = list(
      title = list(text = "Observed_Expected ratio", side = "top"),
      tickmode          = "array",
      tickvals          = c(0, 0.5, 0.9, 1.1, 2, 3, 4, 5, 5.25),
      ticktext          = c("0", "0.5", "0.9", "1.1", "2", "3", "4", "5", ">5"),
      len               = 0.8,
      thickness         = 20,
      ticklabelposition = "outside",
      nticks            = 9
    )
  ) %>%
    layout(
      title = list(text = file_title, font = list(size = 11)),
      xaxis = list(
        title      = "",
        tickmode   = "array",
        tickvals   = seq_len(k_sel),
        ticktext   = paste0("Cl. ", clusters_expected),
        range      = c(0.5, k_sel + 0.5),
        fixedrange = TRUE
      ),
      yaxis = list(
        title      = "",
        tickmode   = "array",
        tickvals   = seq_along(sorted_diseases),
        ticktext   = sorted_diseases,
        autorange  = "reversed",
        range      = c(length(sorted_diseases) + 0.5, 0.5),
        tickfont   = list(size = 8),
        fixedrange = TRUE
      ),
      shapes = highlight_shapes,
      margin = list(l = 240, b = 50, r = 100)
    )
  
  # 8.17 Output path: {Males|Females}/k_{k}/{age}/heatmap_....html
  target_dir <- file.path(
    PATH_FIGURES_STRATIFIED,
    if_else(sexo_sel == "MEN", "Males", "Females"),
    paste0("k_", k_sel),
    edad_sel
  )
  
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  
  html_filename <- file.path(
    target_dir,
    paste0("heatmap_interactive_", tag_name, ".html")
  )
  
  saveWidget(p, html_filename, selfcontained = TRUE)
  
  message(
    "Generated: ", tag_name,
    " | N patients = ", nrow(cohort_stratum),
    " | Max O/E ratio = ", round(max_ratio, 2)
  )
})

# 9. EXPORT PATIENT-LEVEL PATTERN SCORES -------------------------------------

# Export the per-stratum patient scores as a reproducible intermediate
# file inside the stratified_hierarchical folder.
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
    dominant_pattern
  ) %>%
  arrange(sexo_val, edad_val, status_val, k, patient_id, pattern)

write_csv(
  patient_pattern_scores_stratified,
  file.path(
    PATH_HCLUST_STRATIFIED,
    "patient_pattern_scores_stratified_k3_k7.csv"
  )
)