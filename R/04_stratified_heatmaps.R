# -----------------------------------------------------------------------------
# 04_stratified_heatmaps.R
#
# Description:
#   Generates interactive heatmaps to evaluate disease enrichment within
#   patient-level comorbidity patterns across demographic and clinical strata.
#
#   For each combination of:
#     - Sex
#     - Age group
#     - HIV status
#     - Number of disease clusters (k = 3, 4, 5, 6, 7)
#
#   the script calculates an Observed/Expected (O/E) ratio for each disease
#   and patient cluster:
#
#       O/E ratio =
#       prevalence of the disease within the patient cluster /
#       prevalence of the disease within the complete stratum
#
#   Patient clusters are defined using the dominant_pattern derived from
#   the patient-level pattern scores.
#
#   Heatmap interpretation:
#     - Color represents the Observed/Expected enrichment ratio
#     - A ratio of 1 indicates that disease prevalence in the patient cluster
#       is equal to the prevalence expected from the corresponding stratum.
#     - Ratios >1 indicate enrichment, whereas ratios <1 indicate depletion.
#     - Black borders identify diseases that are structurally assigned to the
#       corresponding disease cluster by the hierarchical clustering analysis
#
#   The color scale is discrete and fixed across all heatmaps:
#     0-0.5   = white
#     0.5-0.9 = very light blue
#     0.9-1.1 = neutral grey
#     1.1-2   = pink
#     2-3     = yellow
#     3-4     = green
#     4-5     = orange
#     >5      = dark red
#
#   Values above 5 are visually saturated at the upper end of the scale,
#   while their exact Observed/Expected ratio remains available in the
#   interactive tooltip.
#
#   Missing values (NA) indicate that an Observed/Expected ratio could not
#   be estimated, whereas a ratio of 0 indicates that the disease was not
#   observed within the corresponding patient cluster.
# -----------------------------------------------------------------------------

# Load the baseline demographic and clinical data
bdu_data <- read_csv(PATH_BDU_FULL)
# Load the binary disease matrix used in the clustering analysis
diag_data <- read_csv(PATH_DIAG_BOOL_PREVALENT_FINAL)

source(file.path(PATH_R, "06_heatmap_common.R"))  # shared heatmap functions, colours, statistics
# Load disease-to-cluster assignments for all candidate k solutions
patrones_multi_df <- read_csv(file.path(PATH_HIERARCHICAL_FULL_TABLES, "jaccard_disease_patterns_multi_k.csv"))
# Load patient-level pattern scores and dominant pattern assignments
pattern_data <- read_csv(file.path(PATH_HIERARCHICAL_FULL_PATIENT_DATA, "patient_pattern_scores_hierarchical_k3_k7.csv"))


# Define demographic and clinical strata
# Sex, age group and HIV status are used to calculate stratum-specific
# expected disease prevalence
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
      edad > 65 ~ "over_65"
    ),
    status_hiv = if_else(as.logical(vih_bool), "hiv", "ctl")
  ) %>%
  inner_join(diag_data, by = "patient_id")

# Heatmap output folder: PATH_HIERARCHICAL_FULL_HEATMAPS (see 00_config.R)
# Define all combinations of demographic/clinical strata and candidate
# clustering solutions to be evaluated
k_candidatos <- c(3, 4, 5, 6, 7)
estratos <- expand.grid(
  sexo_val   = c("MALE", "FEMALE"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  k_val      = k_candidatos, 
  stringsAsFactors = FALSE
)

# (Colour scale, statistics and plotting now live in 06_heatmap_common.R)

# Generate one heatmap for each demographic/clinical stratum and clustering solution:
pwalk(estratos, function(sexo_val, edad_val, status_val, k_val) {
  
  # Assign each patient to the pattern with the highest representation score.
  # This dominant pattern defines the patient cluster used in the enrichment analysis
  patient_clusters <- pattern_data %>%
    
    filter(k == k_val, !is.na(dominant_pattern)) %>%
    select(
      patient_id,
      cluster = dominant_pattern
    ) %>%
    mutate(cluster = as.character(cluster)) %>%
    distinct(patient_id, .keep_all = TRUE)
  
  if (nrow(patient_clusters) == 0) {
    message(
      "No patient assignments for k = ",
      k_val
    )
    return(NULL)
  }
  
  # Retrieve the disease-to-cluster assignments corresponding to the
  # current hierarchical clustering solution
  col_cluster_name <- paste0(
    "cluster_k",
    k_val
  )
  
  patrones_subset <- patrones_multi_df %>%
    select(enfermedad,cluster_enfermedad = all_of(col_cluster_name)) %>%
    mutate(cluster_enfermedad = as.character(cluster_enfermedad),
           enfermedad_col = gsub(" ","_",enfermedad))
  
  # Check that all diseases included in the clustering analysis are
  # available in the patient-level diagnostic matrix
  disease_columns <- patrones_subset$enfermedad_col[
    patrones_subset$enfermedad_col %in%
      colnames(cohort_stratified)
  ]
  
  if (length(disease_columns) == 0) {
    message(
      "No disease columns were found for: ",
      sexo_val,
      " | ",
      edad_val,
      " | ",
      status_val,
      " | k=",
      k_val
    )
    return(NULL)
  }
  
  # Retain only diseases that are available in the patient-level dataset
  patrones_subset <- patrones_subset %>%
    filter(enfermedad_col %in% disease_columns)
  
  # Restrict the cohort to the current demographic and clinical stratum
  # and assign each patient to their dominant patient-level pattern
  cohort_stratum <- cohort_stratified %>%
    filter(
      sexo_clean == sexo_val,
      grupo_edad == edad_val,
      status_hiv == status_val
    ) %>%
    # Add the patients cluster
    left_join(patient_clusters, by = "patient_id")  # unassigned patients stay in the stratum (cluster = NA)
  
  # Skip the current stratum if no patients are available.
  if (nrow(cohort_stratum) == 0) {
    message(
      "No patients for ",
      sexo_val,
      " | ",
      edad_val,
      " | ",
      status_val,
      " | k=",
      k_val
    )
    return(NULL)
  }
  
  # Convert the disease matrix to long format so that each row represents
  # one patient-disease combination
  # (long format is now built inside make_oe_heatmap)

  hm_files <- oe_heatmap_files(PATH_HIERARCHICAL_FULL_HEATMAPS, "hierarchical_full",
                               sexo_val, edad_val, status_val, k_val)
  # Disease -> pattern map of this k (solid border around the diseases of each pattern)
  disease_map_plot <- patrones_subset %>%
    transmute(enfermedad_col, enfermedad, cluster_enfermedad, sign = "pos")

  make_oe_heatmap(
    cohort_stratum = cohort_stratum,
    disease_map    = disease_map_plot,
    clusters       = as.character(seq_len(k_val)),
    title_core     = paste0("Observed/Expected Ratio (k=", k_val, ") - ",
                            sexo_val, " | Age: ", edad_val,
                            " | Status: ", toupper(status_val)),
    out_html = hm_files$html,
    out_csv  = hm_files$csv
  )
  invisible(NULL)
}
)
