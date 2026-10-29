# -----------------------------------------------------------------------------
# 07_efa_heatmaps.R
#
# Description:
#   Observed/Expected heatmaps for the Exploratory Factor Analysis (EFA), in
#   exactly the same format as the hierarchical (04_, 05_) and LCA (04b_)
#   heatmaps, so the three methods can be compared side by side.
#
#   Input : loadings exported by 3_Exploratory_Factor_Analysis.R
#           (efa_loadings_<stratum>_<Ctl|Hiv>_2%_<k>factors.csv)
#   Rules (decided by the team):
#     - A disease belongs to a factor when |loading| >= 0.20.
#       A disease can belong to several factors (cross-loadings are kept).
#     - Positive loadings -> solid border; negative loadings -> dashed border
#       (negative loading = disease tends to be absent in that factor).
#     - Patient assignment (same rule as hierarchical): for every factor,
#       count the patient's diseases with a POSITIVE loading >= 0.20;
#       membership = count / sum of counts; the patient goes to the factor
#       with the highest membership; ties or no disease -> not assigned.
#     - O/E, Fisher (cluster vs rest) + FDR and all the marks come from
#       06_heatmap_common.R.
#
#   The EFA was fitted separately in HIV and controls within each sex x age
#   stratum, as in the hierarchical-stratified and LCA analyses.
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

source(file.path(PATH_R, "06_heatmap_common.R"))

# 1. CONFIGURATION ------------------------------------------------------------

EFA_PREV       <- 2          # prevalence threshold (%) used in the EFA file names
K_RANGE        <- 3:7        # same range as the other methods
LOADING_CUTOFF <- 0.20       # |loading| to assign a disease to a factor

# Input loadings: PATH_EFA_TABLES/<stratum>/ ; heatmaps: PATH_EFA_HEATMAPS (00_config.R)
# Labels used by the EFA script for sex and age bands
efa_sex <- c(MALE = "HOMBRE", FEMALE = "MUJER")
efa_age <- c(under_45 = "under 45", `45to65` = "45 to 65", over_65 = "over 65")

# 2. DATA ---------------------------------------------------------------------

bdu_data  <- read_csv(PATH_BDU_FULL)
diag_data <- read_csv(PATH_DIAG_BOOL_PREVALENT_FINAL)

# Same stratum definition as 04_stratified_heatmaps.R
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

estratos <- expand.grid(
  sexo_val   = c("MALE", "FEMALE"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  k_val      = K_RANGE,
  stringsAsFactors = FALSE
)

# 3. ONE HEATMAP PER STRATUM x NUMBER OF FACTORS ------------------------------

efa_assignments <- list()

pwalk(estratos, function(sexo_val, edad_val, status_val, k_val) {
  
  strata_label <- paste0(efa_sex[[sexo_val]], "_", efa_age[[edad_val]])
  hiv_label    <- if (status_val == "hiv") "Hiv" else "Ctl"
  
  loadings_file <- file.path(
    PATH_EFA_TABLES, strata_label,
    sprintf("efa_loadings_%s_%s_%s%%_%dfactors.csv",
            strata_label, hiv_label, EFA_PREV, as.integer(k_val))
  )
  if (!file.exists(loadings_file)) {
    message("Missing EFA loadings: ", loadings_file)
    return(NULL)
  }
  
  # 3.1 Loadings -> disease/factor map ---------------------------------------
  L <- read_csv(loadings_file, show_col_types = FALSE)
  names(L)[1] <- "enfermedad_col"
  factor_cols <- setdiff(names(L), "enfermedad_col")
  
  # Label = column name with spaces (same labelling as the hierarchical files)
  all_dis <- L %>%
    transmute(enfermedad_col, enfermedad = gsub("_", " ", enfermedad_col))
  
  patterns <- L %>%
    pivot_longer(all_of(factor_cols), names_to = "factor", values_to = "loading") %>%
    filter(abs(loading) >= LOADING_CUTOFF) %>%
    mutate(
      cluster_enfermedad = as.character(match(factor, factor_cols)),
      enfermedad = gsub("_", " ", enfermedad_col),
      sign = if_else(loading > 0, "pos", "neg")
    ) %>%
    select(enfermedad_col, enfermedad, cluster_enfermedad, sign, loading)
  
  disease_map <- bind_rows(
    patterns,
    all_dis %>% mutate(cluster_enfermedad = NA_character_, sign = "pos", loading = NA_real_)
  )
  
  # 3.2 Patients of the stratum and their factor (membership rule) -----------
  cohort_stratum <- cohort_stratified %>%
    filter(sexo_clean == sexo_val, grupo_edad == edad_val, status_hiv == status_val)
  
  pos_patterns <- patterns %>%
    filter(sign == "pos") %>%
    transmute(enfermedad_col, pattern = cluster_enfermedad)
  
  if (nrow(cohort_stratum) == 0 || nrow(pos_patterns) == 0) {
    message("Nothing to assign for: ", strata_label, " | ", status_val, " | k=", k_val)
    return(NULL)
  }
  
  assignment <- assign_patients_by_membership(cohort_stratum, pos_patterns)
  cohort_stratum <- cohort_stratum %>%
    left_join(assignment %>% select(patient_id, cluster), by = "patient_id")
  
  # 3.3 Heatmap ---------------------------------------------------------------
  hm_files <- oe_heatmap_files(PATH_EFA_HEATMAPS, "EFA", sexo_val, edad_val, status_val, k_val)

  # Keep the patient -> factor assignment (for concordance between methods)
  assignment$sexo_val   <- sexo_val
  assignment$edad_val   <- edad_val
  assignment$status_val <- status_val
  assignment$k          <- as.integer(k_val)
  efa_assignments[[paste(sexo_val, edad_val, status_val, k_val, sep = "_")]] <<-
    assignment[, c("patient_id", "sexo_val", "edad_val", "status_val", "k", "cluster", "max_membership")]
  make_oe_heatmap(
    cohort_stratum = cohort_stratum,
    disease_map    = disease_map,
    clusters       = as.character(seq_len(k_val)),
    title_core     = paste0("Observed/Expected Ratio (EFA, k=", k_val, ") - ",
                            sexo_val, " | Age: ", edad_val,
                            " | Status: ", toupper(status_val)),
    out_html = hm_files$html,
    out_csv  = hm_files$csv
  )
  invisible(NULL)
})

# 4. EXPORT PATIENT-LEVEL ASSIGNMENTS -------------------------------------------
if (length(efa_assignments) > 0) {
  write_csv(
    bind_rows(efa_assignments),
    file.path(PATH_EFA_PATIENT_DATA, "patient_pattern_assignments_EFA_k3_k7.csv")
  )
}
