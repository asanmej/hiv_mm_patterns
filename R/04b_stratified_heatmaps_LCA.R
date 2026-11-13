# -----------------------------------------------------------------------------
# 04b_stratified_heatmaps_LCA.R
#
# Description:
#   Generates interactive heatmaps to evaluate disease enrichment within the
#   latent classes obtained by the stratified Latent Class Analysis (LCA)
#   performed in LCA_HIV.R. It is the LCA counterpart of
#   04_stratified_heatmaps.R (hierarchical clustering).
#
#   For each stratum (sex x age group x HIV status) and for each latent-class
#   solution (k) evaluated in that stratum, the script calculates an
#   Observed/Expected (O/E) ratio for each disease and latent class:
#
#       O/E ratio =
#       prevalence of the disease within the latent class /
#       prevalence of the disease within the complete stratum
#
#   Heatmap interpretation:
#     - Color represents the Observed/Expected enrichment ratio
#     - A ratio of 1 indicates that disease prevalence in the latent class
#       is equal to the prevalence expected from the corresponding stratum.
#     - Ratios >1 indicate enrichment, whereas ratios <1 indicate depletion.
#
#   The color scale is discrete and fixed across all heatmaps (identical to
#   the one used in 04_stratified_heatmaps.R):
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
#   while their exact ratio remains available in the interactive tooltip.
#
#   Missing values (NA) indicate that an Observed/Expected ratio could not
#   be estimated (e.g. a latent class with no patients), whereas a ratio of 0
#   indicates that the disease was not observed within that latent class.
#
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(plotly)
library(htmlwidgets)

source(file.path(PATH_R, "06_heatmap_common.R"))  # shared heatmap functions, colours, statistics


#### 1. Load the LCA patient assignments and the binary disease matrix

# Patient-level latent class assignments for every stratum and evaluated k
pattern_data <- read_csv(
  file.path(PATH_LCA_PATIENT_DATA, "patient_pattern_scores_LCA.csv"),
  show_col_types = FALSE
)

# Binary disease matrix (one column per disease)
diag_data <- read_csv(
  PATH_DIAG_BOOL_PREVALENT_FINAL,
  show_col_types = FALSE
)


#### 2. Disease name mapping

# Inside LCA_HIV.R every disease column is renamed "dis_<original name>" and
# then passed through make.names() (which, for example, turns "-" into "."),
# so disease_names saved with each model use that internal naming. This table
# reproduces the same transformation on the original diag_data column names,
# which lets us translate disease_names back to the original columns.
#
# vih_dt, HIV_infection and HIV-related_disease are not disease indicators and
# were excluded before fitting the LCA models (LCA_HIV.R, section 1.4).
non_disease_cols <- c("patient_id", "vih_dt", "HIV_infection", "HIV-related_disease")

disease_name_map <- tibble(
  enfermedad_col = setdiff(colnames(diag_data), non_disease_cols)
) %>%
  mutate(
    # Name as it appears inside the fitted LCA models
    dis_col_clean = make.names(paste0("dis_", enfermedad_col)),
    # Readable label for the heatmap axis (cosmetic only: the original
    # column name above is the one used to read the data)
    enfermedad = trimws(gsub("_+", " ", enfermedad_col))
  )

# The mapping is only valid if every internal name is unique
if (anyDuplicated(disease_name_map$dis_col_clean) > 0) {
  stop(
    "make.names() produced duplicated disease names: ",
    paste(
      unique(disease_name_map$dis_col_clean[duplicated(disease_name_map$dis_col_clean)]),
      collapse = ", "
    )
  )
}


#### 3. Recover the disease list used by the LCA model of each stratum

# The RData files also contain the full fitted models (res), which are large,
# so each file is loaded into a temporary environment and only disease_names
# is kept.
strata_unique <- pattern_data %>%
  distinct(sexo_val, edad_val, status_val)

disease_names_by_stratum <- list()

for (i in seq_len(nrow(strata_unique))) {
  
  sexo_i   <- strata_unique$sexo_val[i]
  edad_i   <- strata_unique$edad_val[i]
  status_i <- strata_unique$status_val[i]
  
  stratum_i <- paste(sexo_i, edad_i, status_i, sep = "_")
  
  rdata_path <- file.path(
    PATH_LCA_MODELS,
    sexo_i,
    edad_i,
    status_i,
    paste0("ResultsLCA_", stratum_i, ".RData")
  )
  
  if (!file.exists(rdata_path)) {
    message("No saved LCA model found for stratum: ", stratum_i)
    next
  }
  
  tmp_env <- new.env()
  load(rdata_path, envir = tmp_env)
  
  disease_names_by_stratum[[stratum_i]] <- as.character(tmp_env$disease_names)
  
  rm(tmp_env)
}


#### 4. Define the strata / k combinations to plot

# Taken from the data instead of assuming a full grid: each stratum retained
# its own k values. Columns are named so that they cannot be confused with the
# columns of pattern_data inside filter().
strata_grid <- pattern_data %>%
  distinct(sexo_val, edad_val, status_val, k) %>%
  arrange(sexo_val, edad_val, status_val, k) %>%
  rename(
    stratum_sex    = sexo_val,
    stratum_age    = edad_val,
    stratum_status = status_val,
    k_val          = k
  )


#### 5. Colour scale, statistics and plotting: see 06_heatmap_common.R

# 6. Generate one heatmap for each stratum and evaluated k

lca_status <- pmap_chr(strata_grid, function(stratum_sex, stratum_age, stratum_status, k_val) {
  
  stratum_name <- paste(stratum_sex, stratum_age, stratum_status, sep = "_")
  
  # 6.1 Diseases used by the LCA model of this stratum 
  model_disease_names <- disease_names_by_stratum[[stratum_name]]
  
  if (is.null(model_disease_names)) {
    message("Skipped (no disease list available): ", stratum_name, " | k = ", k_val)
    return("skipped: no disease list (model .RData not found)")
  }
  
  disease_map <- tibble(dis_col_clean = model_disease_names) %>%
    left_join(disease_name_map, by = "dis_col_clean")
  
  # Every disease of the model must be translated back to a diag_data column.
  # If not, the naming transformation differs from the one assumed above and
  # the problem has to be fixed before trusting any heatmap.
  if (anyNA(disease_map$enfermedad_col)) {
    stop(
      "Diseases of the LCA model could not be mapped to diag_data columns in ",
      stratum_name, ": ",
      paste(disease_map$dis_col_clean[is.na(disease_map$enfermedad_col)], collapse = ", ")
    )
  }
  
  # 6.2 Latent class of each patient in this stratum / k 
  patient_clusters <- pattern_data %>%
    filter(
      sexo_val   == stratum_sex,
      edad_val   == stratum_age,
      status_val == stratum_status,
      k          == k_val,
      !is.na(cluster)
    ) %>%
    transmute(
      patient_id,
      cluster = as.character(cluster)
    ) %>%
    distinct(patient_id, .keep_all = TRUE)
  
  if (nrow(patient_clusters) == 0) {
    message("No patient assignments for ", stratum_name, " | k = ", k_val)
    return("skipped: no patient assignments")
  }
  
  # 6.3 Disease status of the patients of this stratum 
  cohort_stratum <- diag_data %>%
    select(patient_id, all_of(disease_map$enfermedad_col)) %>%
    inner_join(patient_clusters, by = "patient_id")
  
  if (nrow(cohort_stratum) < nrow(patient_clusters)) {
    message(
      "Warning: ",
      nrow(patient_clusters) - nrow(cohort_stratum),
      " patients of ", stratum_name, " | k = ", k_val,
      " were not found in the disease matrix"
    )
  }
  
  if (nrow(cohort_stratum) == 0) {
    message("No patients for ", stratum_name, " | k = ", k_val)
    return("skipped: no patients in disease matrix")
  }
  
  # 6.4-6.10 Statistics (O/E, Fisher vs rest + FDR), marks and HTML output are
  # produced by the shared function, identical for hierarchical / LCA / EFA.
  hm_files <- oe_heatmap_files(PATH_LCA_HEATMAPS, "LCA",
                               stratum_sex, stratum_age, stratum_status, k_val)
  res <- make_oe_heatmap(
    cohort_stratum = cohort_stratum,
    disease_map    = disease_map %>%
      transmute(enfermedad_col, enfermedad, cluster_enfermedad = NA_character_),
    clusters       = as.character(seq_len(k_val)),
    title_core     = paste0("Observed/Expected Ratio (LCA, k=", k_val, ") - ",
                            stratum_sex, " | Age: ", stratum_age,
                            " | Status: ", toupper(stratum_status)),
    out_html = hm_files$html,
    out_csv  = hm_files$csv
  )
  if (is.null(res)) return("skipped: no valid O/E ratios")
  "ok"
}
)

# Summary: tell clearly what was (not) generated
message("\nLCA heatmaps: ", sum(lca_status == "ok"), " generated out of ", length(lca_status), " solutions")
if (any(lca_status != "ok")) print(table(lca_status[lca_status != "ok"]))
if (!any(lca_status == "ok")) {
  stop("No LCA heatmap was generated. Check the messages above (model .RData files in ",
       PATH_LCA_MODELS, " and patient_pattern_scores_LCA.csv in ", PATH_LCA_PATIENT_DATA, ").",
       call. = FALSE)
}
