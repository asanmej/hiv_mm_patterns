# -----------------------------------------------------------------------------
# 07_efa_heatmaps.R
#
# Description:
#   Observed/Expected heatmaps for the Exploratory Factor Analysis (EFA), in
#   exactly the same format as the hierarchical (05_) and LCA (04b_) heatmaps,
#   so the three methods can be compared side by side.
#
#   Input : loadings exported by 3_Exploratory_Factor_Analysis.R
#           (PATH_EFA_TABLES/<stratum>/efa_loadings_<stratum>_<Ctl|Hiv>_2%_<k>factors.csv)
#           The 3 solutions of each stratum and group are read from the files
#           that exist (they are chosen in the EFA script).
#   Rules (decided by the team):
#     - A disease belongs to a factor when |loading| >= 0.20.
#       A disease can belong to several factors (cross-loadings are kept).
#     - Positive loadings -> solid border; negative loadings -> dashed border
#       (negative loading = disease tends to be absent in that factor).
#     - Patient assignment: only POSITIVE loadings >= 0.20 count. It is done with
#       every rule of ASSIGNMENT_RULES ("count" and "proportion", see
#       06_heatmap_common.R) and the results are compared.
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
LOADING_CUTOFF <- EFA_LOADING_CUTOFF   # defined in 00_config.R: |loading| to assign a disease to a factor

# Labels used by the EFA script for sex and age bands
efa_sex <- c(MALE = "HOMBRE", FEMALE = "MUJER")
efa_age <- c(under_45 = "under 45", `45to65` = "45 to 65", over_65 = "over 65")

# 2. DATA ---------------------------------------------------------------------

bdu_data  <- read_csv(PATH_BDU_FULL)
diag_data <- read_csv(PATH_DIAG_BOOL_PREVALENT_FINAL)

# Same stratum definition as 05_stratified_heatmaps.R
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

# Solutions (number of factors) available for every stratum and group, crossed
# with the assignment rules
estratos_base <- expand.grid(
  sexo_val   = c("MALE", "FEMALE"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

estratos <- pmap_dfr(estratos_base, function(sexo_val, edad_val, status_val) {
  strata_label <- paste0(efa_sex[[sexo_val]], "_", efa_age[[edad_val]])
  hiv_label    <- if (status_val == "hiv") "Hiv" else "Ctl"
  files <- list.files(
    file.path(PATH_EFA_TABLES, strata_label),
    pattern = sprintf("^efa_loadings_.*_%s_%s%%_[0-9]+factors\\.csv$", hiv_label, EFA_PREV)
  )
  if (length(files) == 0) return(NULL)
  tibble(
    sexo_val, edad_val, status_val,
    k_val = sort(as.integer(sub(".*_([0-9]+)factors\\.csv$", "\\1", files)))
  )
}) %>%
  tidyr::crossing(rule = ASSIGNMENT_RULES)

if (nrow(estratos) == 0) stop("No EFA loadings found in ", PATH_EFA_TABLES, ". Run 3_Exploratory_Factor_Analysis.R first.")

# 3. ONE HEATMAP PER STRATUM x SOLUTION x ASSIGNMENT RULE ------------------------

rule_results <- pmap(estratos, function(sexo_val, edad_val, status_val, k_val, rule) {
  
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
  
  # 3.2 Patients of the stratum and their factor ------------------------------
  cohort_stratum <- cohort_stratified %>%
    filter(sexo_clean == sexo_val, grupo_edad == edad_val, status_hiv == status_val)
  
  pos_patterns <- patterns %>%
    filter(sign == "pos") %>%
    transmute(enfermedad_col, pattern = cluster_enfermedad)
  
  if (nrow(cohort_stratum) == 0 || nrow(pos_patterns) == 0) {
    message("Nothing to assign for: ", strata_label, " | ", status_val, " | k=", k_val)
    return(NULL)
  }
  
  assignment <- assign_patients_by_membership(cohort_stratum, pos_patterns, rule = rule)
  cohort_stratum <- cohort_stratum %>%
    left_join(assignment %>% select(patient_id, cluster), by = "patient_id")
  
  # 3.3 Heatmap ---------------------------------------------------------------
  hm_files <- oe_heatmap_files(file.path(PATH_EFA_HEATMAPS, rule), paste0("EFA_", rule),
                               sexo_val, edad_val, status_val, k_val)
  
  res <- make_oe_heatmap(
    cohort_stratum = cohort_stratum,
    disease_map    = disease_map,
    clusters       = as.character(seq_len(k_val)),
    title_core     = paste0("Observed/Expected Ratio (EFA, k=", k_val, ") - ",
                            sexo_val, " | Age: ", edad_val,
                            " | Status: ", toupper(status_val),
                            " | Assign: ", rule),
    out_html = hm_files$html,
    out_csv  = hm_files$csv
  )
  if (is.null(res)) return(NULL)
  
  list(
    sexo_val = sexo_val, edad_val = edad_val, status_val = status_val,
    k = as.integer(k_val), rule = rule, table = res$table,
    assignment = cohort_stratum[, c("patient_id", "cluster")]
  )
})

# 4. EXPORTS ------------------------------------------------------------------

# 4.1 Patient -> factor assignments (for concordance between methods)
efa_assignments <- bind_rows(lapply(Filter(Negate(is.null), rule_results), function(r) {
  r$assignment %>%
    mutate(sexo_val = r$sexo_val, edad_val = r$edad_val, status_val = r$status_val,
           k = r$k, rule = r$rule)
}))
if (nrow(efa_assignments) > 0) {
  write_csv(efa_assignments,
            file.path(PATH_EFA_PATIENT_DATA, "patient_pattern_assignments_EFA.csv"))
}

# 4.2 Sensitivity of the results to the assignment rule (count vs proportion)
compare_assignment_rules(
  rule_results,
  file.path(PATH_EFA_TABLES, "assignment_rule_comparison.csv")
)
