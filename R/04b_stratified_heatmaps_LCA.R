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


#### 1. Load the LCA patient assignments and the binary disease matrix

# Patient-level latent class assignments for every stratum and evaluated k
pattern_data <- read_csv(
  file.path(PATH_INTERMEDIATE_LCA, "patient_pattern_scores_LCA.csv"),
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
    PATH_LCA_STRATIFIED_MODELS,
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


#### 5. Fixed discrete color scale for the Observed/Expected ratio

# Identical to 04_stratified_heatmaps.R so that both sets of heatmaps are
# directly comparable.
color_0_05  <- "#FFFFFF"
color_05_09 <- "#E8EEF8"
color_09_11 <- "#C7CCD4"
color_11_2  <- "#E89BB7"
color_2_3   <- "#F2D35E"
color_3_4   <- "#78B87A"
color_4_5   <- "#E99A4A"
color_gt5   <- "#B83A4A"

Observed_Expected_colorscale <- list(
  
  # 0-0.5
  c(0.00 / 5.5, color_0_05),
  c(0.50 / 5.5, color_0_05),
  
  # 0.5-0.9
  c(0.50 / 5.5, color_05_09),
  c(0.90 / 5.5, color_05_09),
  
  # 0.9-1.1
  c(0.90 / 5.5, color_09_11),
  c(1.10 / 5.5, color_09_11),
  
  # 1.1-2
  c(1.10 / 5.5, color_11_2),
  c(2.00 / 5.5, color_11_2),
  
  # 2-3
  c(2.00 / 5.5, color_2_3),
  c(3.00 / 5.5, color_2_3),
  
  # 3-4
  c(3.00 / 5.5, color_3_4),
  c(4.00 / 5.5, color_3_4),
  
  # 4-5
  c(4.00 / 5.5, color_4_5),
  c(5.00 / 5.5, color_4_5),
  
  # >5
  c(5.00 / 5.5, color_gt5),
  c(5.50 / 5.5, color_gt5)
)


# 6. Generate one heatmap for each stratum and evaluated k

pwalk(strata_grid, function(stratum_sex, stratum_age, stratum_status, k_val) {
  
  stratum_name <- paste(stratum_sex, stratum_age, stratum_status, sep = "_")
  
  # 6.1 Diseases used by the LCA model of this stratum 
  model_disease_names <- disease_names_by_stratum[[stratum_name]]
  
  if (is.null(model_disease_names)) {
    message("Skipped (no disease list available): ", stratum_name, " | k = ", k_val)
    return(NULL)
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
    return(NULL)
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
    return(NULL)
  }
  
  # 6.4 Long format: one row per patient-disease combination 
  diag_long <- cohort_stratum %>%
    pivot_longer(
      cols = all_of(disease_map$enfermedad_col),
      names_to = "enfermedad_col",
      values_to = "has_disease"
    ) %>%
    mutate(
      # 1 = disease present, 0 = disease absent
      has_disease = if_else(is.na(has_disease) | has_disease == 0, 0, 1)
    ) %>%
    # Readable disease label
    left_join(
      disease_map %>% select(enfermedad_col, enfermedad),
      by = "enfermedad_col"
    )
  
  # 6.5 Expected prevalence: whole stratum, regardless of class
  global_prev <- diag_long %>%
    group_by(enfermedad) %>%
    summarise(
      expected_prev = mean(has_disease, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Remove diseases that are completely absent from the stratum
    filter(expected_prev > 0)
  
  # 6.6 Observed prevalence: within each latent class 
  cluster_prev <- diag_long %>%
    group_by(enfermedad, cluster) %>%
    summarise(
      observed_prev = mean(has_disease, na.rm = TRUE),
      n_patients = n(),
      n_cases = sum(has_disease, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 6.7 Observed/Expected ratio 
  # O/E = P(Disease | Latent class) / P(Disease | Stratum)
  heatmap_matrix <- cluster_prev %>%
    inner_join(global_prev, by = "enfermedad") %>%
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
  
  # Complete the disease x latent-class grid. A latent class of the model with
  # no assigned patients would otherwise disappear from the heatmap.
  #   NA = no valid Observed/Expected estimate is available
  #   0  = the disease has an observed prevalence of zero in the class
  clusters_expected <- as.character(seq_len(k_val))
  
  heatmap_matrix <- heatmap_matrix %>%
    complete(
      enfermedad = unique(global_prev$enfermedad),
      cluster = clusters_expected
    )
  
  # Order diseases alphabetically and classes numerically
  sorted_diseases <- sort(unique(heatmap_matrix$enfermedad))
  
  heatmap_matrix <- heatmap_matrix %>%
    mutate(
      cluster_num = match(cluster, clusters_expected),
      disease_num = match(enfermedad, sorted_diseases)
    )
  
  # Size of each latent class (% of the patients of the stratum), shown in the
  # x-axis labels so they can be compared with the class sizes printed in the
  # LCA figures
  cluster_labels <- tibble(cluster = clusters_expected) %>%
    left_join(
      cohort_stratum %>%
        count(cluster) %>%
        mutate(pct = round(100 * n / sum(n))),
      by = "cluster"
    ) %>%
    mutate(
      pct = coalesce(pct, 0),
      label = paste0("Cl. ", cluster, " (", pct, "%)")
    ) %>%
    pull(label)
  
  # Check that valid ratios are available
  max_ratio <- suppressWarnings(max(heatmap_matrix$value, na.rm = TRUE))
  
  if (!is.finite(max_ratio)) {
    message("No valid ratios for: ", stratum_name, " | k = ", k_val)
    return(NULL)
  }
  
  # 6.8 Titles and hover text 
  file_title <- paste0(
    "Observed/Expected Ratio (LCA, k=", k_val, ") - ",
    stratum_sex,
    " | Age: ", stratum_age,
    " | Status: ", toupper(stratum_status),
    " | N = ", nrow(cohort_stratum)
  )
  
  heatmap_matrix <- heatmap_matrix %>%
    mutate(
      ratio_text = if_else(is.na(value), "NA", paste0(round(value, 2)))
    )
  
  # 6.9 Interactive heatmap 
  p <- plot_ly(
    data = heatmap_matrix,
    
    x = ~cluster_num,
    y = ~disease_num,
    
    # The actual Observed/Expected ratio determines the cell colour
    z = ~value,
    
    # Fixed visual scale from 0 to 5.5: values above 5 are displayed in the
    # >5 colour category while their exact values remain in the tooltip
    zmin = 0,
    zmax = 5.5,
    
    type = "heatmap",
    
    colorscale = Observed_Expected_colorscale,
    
    hoverongaps = FALSE,
    
    # Hover text
    text = ~paste0(
      "<b>", enfermedad, "</b>",
      "<br><br>",
      
      "<b>Observed/Expected ratio:</b> ",
      ratio_text,
      
      "<br><b>Observed prevalence:</b> ",
      round(observed_prev, 4),
      
      "<br><b>Expected prevalence:</b> ",
      round(expected_prev, 4),
      
      "<br><b>N patients:</b> ",
      n_patients,
      
      "<br><b>Cases:</b> ",
      n_cases
    ),
    
    hovertemplate =
      paste0(
        "%{text}",
        "<br><br>",
        "<b>Latent class:</b> Cl. %{x}",
        "<extra></extra>"
      ),
    
    # Colorbar with fixed tick positions at the interval boundaries
    colorbar = list(
      title = list(
        text = "Observed_Expected ratio",
        side = "top"
      ),
      
      tickmode = "array",
      
      tickvals = c(0, 0.5, 0.9, 1.1, 2, 3, 4, 5, 5.25),
      
      ticktext = c("0", "0.5", "0.9", "1.1", "2", "3", "4", "5", ">5"),
      
      len = 0.8,
      thickness = 20,
      ticklabelposition = "outside",
      nticks = 9
    )
  ) %>%
    
    layout(
      title = list(
        text = file_title,
        font = list(size = 11)
      ),
      
      xaxis = list(
        title = "",
        tickmode = "array",
        tickvals = seq_along(clusters_expected),
        ticktext = cluster_labels,
        range = c(0.5, length(clusters_expected) + 0.5),
        fixedrange = TRUE
      ),
      
      yaxis = list(
        title = "",
        tickmode = "array",
        tickvals = seq_along(sorted_diseases),
        ticktext = sorted_diseases,
        autorange = "reversed",
        range = c(length(sorted_diseases) + 0.5, 0.5),
        tickfont = list(size = 8),
        fixedrange = TRUE
      ),
      
      margin = list(
        l = 240,
        b = 50,
        r = 100
      )
    )
  
  # 6.10 Save as a self-contained HTML file 
  # The heatmap is saved next to the other LCA figures of the same stratum
  # (metodo1/metodo2/prev_...), following their naming convention.
  target_dir <- file.path(
    PATH_LCA_STRATIFIED_FIGURES,
    stratum_sex,
    stratum_age,
    stratum_status
  )
  
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  
  html_filename <- file.path(
    target_dir,
    paste0("heatmap", k_val, "LC_", stratum_name, ".html")
  )
  
  saveWidget(
    p,
    html_filename,
    selfcontained = TRUE
  )
  
  # Report generation status
  message(
    "Generated: ",
    stratum_name,
    " | k = ", k_val,
    " | N patients = ", nrow(cohort_stratum),
    " | N diseases = ", length(unique(global_prev$enfermedad)),
    " | Max O/E ratio = ", round(max_ratio, 2)
  )
}
)