# -----------------------------------------------------------------------------
# 03_hierarchical_clustering_stratified.R
#
# Description:
#   Stratified hierarchical clustering of disease co-occurrence patterns
#   using Jaccard distance and Ward.D2 linkage.
#
#   The analysis is repeated independently for each of the twelve
#   sex × age × HIV status strata:
#       MEN   | under_45 | ctl
#       MEN   | under_45 | hiv
#       MEN   | 45to65   | ctl
#       MEN   | 45to65   | hiv
#       MEN   | over_65  | ctl
#       MEN   | over_65  | hiv
#       WOMEN | under_45 | ctl
#       WOMEN | under_45 | hiv
#       WOMEN | 45to65   | ctl
#       WOMEN | 45to65   | hiv
#       WOMEN | over_65  | ctl
#       WOMEN | over_65  | hiv
#
#   Unlike the LCA pipeline, where HIV status is also a stratification
#   variable, here the clustering is estimated separately within each
#   HIV status subgroup. As a consequence, disease pattern IDs are
#   NOT directly comparable across ctl and hiv: Pattern 1 in ctl and
#   Pattern 1 in hiv do not necessarily represent the same set of
#   diseases. This is intentional: the objective is to characterise
#   the internal structure of each subgroup on its own terms, and to
#   document how the comorbidity pattern differs between subgroups.
#
#   All thresholds, filters, clustering parameters and candidate
#   solutions are identical to the original unstratified script.
#   This script does not depend on the original pipeline: all outputs
#   are written under PATH_INTERMEDIATE/stratified_hierarchical.
# -----------------------------------------------------------------------------

# 1. PACKAGE MANAGEMENT ------------------------------------------------------

library(readr)
library(dplyr)
library(stringr)
library(psych)
library(plotly)
library(htmlwidgets)

# 2. GLOBAL CONFIGURATION ---------------------------------------

# Twelve sex × age × HIV status strata
strata <- expand.grid(
  sexo_val   = c("MEN", "WOMEN"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

# Candidate numbers of clusters
k_candidates <- 3:7

# Minimum disease prevalence filter (identical to the original script)
min_prev <- 0.02

# Root output directory for the stratified hierarchical clustering
PATH_HCLUST_STRATIFIED <- file.path(
  PATH_INTERMEDIATE,
  "stratified_hierarchical"
)

dir.create(
  PATH_HCLUST_STRATIFIED,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. DATA LOADING ----------------------------------------------------------

# Load prevalent diagnostic data
diag_data <- read_csv(PATH_DIAG_BOOL_PREVALENT_FINAL)

# Load baseline demographic data (needed to define the strata)
bdu_data <- read_csv(PATH_BDU_FULL)

# Reconstruct the sex × age × HIV status classification.
# status_hiv is derived only for downstream use; it is NOT used to
# define the clustering strata.
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
  )

# 4. HELPER: TETRACHORIC CORRELATION HEATMAP ---------------------

# Same discrete tetrachoric correlation heatmap as in the original
# unstratified script. It is returned as a plotly object so that it
# can be saved as a self-contained HTML file per stratum.
build_tetrachoric_heatmap <- function(corr_ordered, title_text) {
  
  zmin <- -0.15
  zmax <-  0.15
  
  breaks <- c(
    -0.15, -0.125, -0.10, -0.075, -0.05, -0.025,
    0,
    0.025, 0.05, 0.075, 0.10, 0.125, 0.15
  )
  
  colors <- c(
    "#67000D", "#A50F15", "#CB181D", "#EF3B2C", "#FC8D59",
    "#FFFFFF", "#FFFFFF",
    "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C"
  )
  
  # Build a truly discrete Plotly colour scale
  colorscale <- list()
  
  for (i in seq_along(colors)) {
    lower_norm <- (breaks[i]     - zmin) / (zmax - zmin)
    upper_norm <- (breaks[i + 1] - zmin) / (zmax - zmin)
    colorscale[[length(colorscale) + 1]] <- list(lower_norm, colors[i])
    colorscale[[length(colorscale) + 1]] <- list(upper_norm, colors[i])
  }
  
  plot_ly(
    x = colnames(corr_ordered),
    y = rownames(corr_ordered),
    z = corr_ordered,
    type        = "heatmap",
    colorscale  = colorscale,
    zmin        = zmin,
    zmax        = zmax,
    hoverongaps = FALSE,
    text        = round(corr_ordered, 3),
    hovertemplate = paste0(
      "Disease 1: %{y}<br>",
      "Disease 2: %{x}<br>",
      "Tetrachoric correlation: %{text}<extra></extra>"
    ),
    colorbar = list(
      title    = "Tetrachoric correlation",
      tickmode = "array",
      tickvals = breaks,
      ticktext = sprintf("%.2f", breaks)
    )
  ) %>%
    layout(
      title = title_text,
      xaxis = list(
        title         = "",
        tickangle     = -45,
        tickfont      = list(size = 8),
        categoryarray = colnames(corr_ordered),
        categoryorder = "array"
      ),
      yaxis = list(
        title         = "",
        tickfont      = list(size = 8),
        autorange     = "reversed",
        categoryarray = rownames(corr_ordered),
        categoryorder = "array"
      ),
      margin = list(l = 150, b = 150)
    )
}

# 5. MAIN STRATIFIED LOOP ----------------------------------------------------

for (i in seq_len(nrow(strata))) {
  
  sexo_actual   <- strata$sexo_val[i]
  edad_actual   <- strata$edad_val[i]
  status_actual <- strata$status_val[i]
  stratum_name  <- paste(sexo_actual, edad_actual, status_actual, sep = "_")
  
  message("\n")
  message("=============================================================")
  message("Processing stratum: ", stratum_name)
  message("=============================================================")
  
  # 5.1 Stratum-specific output folder
  stratum_dir <- file.path(
    PATH_HCLUST_STRATIFIED,
    sexo_actual,
    edad_actual,
    status_actual
  )
  
  dir.create(stratum_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 5.2 Select patients belonging to the current sex × age stratum.
  # HIV status is deliberately NOT used here: both PLWH and controls
  # contribute to the same clustering solution.
  ids_stratum <- cohort_stratified %>%
    filter(
      sexo_clean == sexo_actual,
      grupo_edad == edad_actual,
      status_hiv == status_actual
    ) %>%
    pull(patient_id)
  
  diag_stratum <- diag_data %>%
    filter(patient_id %in% ids_stratum)
  
  message("Number of patients: ", nrow(diag_stratum))
  
  # 5.3 Remove patient identifiers and HIV-related variables that are
  # not used as disease indicators in the clustering analysis
  disease_df <- diag_stratum %>%
    select(-patient_id, -vih_dt, -HIV_infection, -`HIV-related_disease`)
  
  # 5.4 Apply the same minimum prevalence filter of 2%
  prevalences <- colMeans(disease_df)
  frequent_diseases <- names(prevalences[prevalences >= min_prev])
  
  disease_df_filtered <- disease_df %>%
    select(all_of(frequent_diseases))
  
  matrix_disease <- as.matrix(disease_df_filtered)
  
  message("Number of diseases retained: ", ncol(matrix_disease))
  
  # Technical safeguard:
  # After the 2 % prevalence filter, at least two diseases and more
  # than one patient are required for the distance matrix and the
  # hierarchical clustering to be computable. This is a technical
  # minimum, not an analytical threshold.
  if (ncol(matrix_disease) < 2 || nrow(matrix_disease) < 2) {
    message(
      "Stratum skipped: fewer than 2 diseases or patients after ",
      "the 2% prevalence filter."
    )
    next
  }
  
  # 5.5 Tetrachoric correlation heatmap (exploratory visualisation)
  corr_matrix <- psych::tetrachoric(matrix_disease, na.rm = TRUE)$rho
  
  clean_names <- gsub("_", " ", colnames(matrix_disease))
  rownames(corr_matrix) <- clean_names
  colnames(corr_matrix) <- clean_names
  
  # Hierarchical clustering is used only to reorder diseases in the
  # correlation heatmap for visualisation. It is not used to define
  # the Jaccard-based disease patterns.
  row_dist      <- as.dist(1 - corr_matrix)
  row_clust     <- hclust(row_dist, method = "ward.D2")
  ordered_names <- row_clust$labels[row_clust$order]
  
  corr_ordered <- corr_matrix[ordered_names, ordered_names]
  diag(corr_ordered) <- NA
  
  p_corr <- build_tetrachoric_heatmap(
    corr_ordered,
    title_text = paste0(
      "Binary Tetrachoric Correlation Matrix - ",
      stratum_name
    )
  )
  
  saveWidget(
    p_corr,
    file.path(
      stratum_dir,
      paste0("tetrachoric_correlation_", stratum_name, ".html")
    ),
    selfcontained = TRUE
  )
  
  # 5.6 Jaccard dissimilarity between diseases
  # matrix_disease has patients as rows and diseases as columns, so it
  # is transposed before computing pairwise dissimilarities between
  # diseases. method = "binary" returns the Jaccard distance for
  # binary vectors: 1 - |A ∩ B| / |A ∪ B|.
  jaccard_dist <- dist(t(matrix_disease), method = "binary")
  
  # 5.7 Ward.D2 hierarchical clustering
  hc <- hclust(jaccard_dist, method = "ward.D2")
  
  # 5.8 Evaluate candidate solutions k = 3 to 7 and build
  # disease-to-cluster assignments
  cluster_assignments <- data.frame(
    enfermedad = colnames(matrix_disease),
    stringsAsFactors = FALSE
  )
  
  for (k_val in k_candidates) {
    clusters_k <- cutree(hc, k = k_val)
    cluster_assignments[[paste0("cluster_k", k_val)]] <-
      as.integer(clusters_k[colnames(matrix_disease)])
  }
  
  # 5.9 Export disease-to-cluster assignments for the current stratum
  write_csv(
    cluster_assignments,
    file.path(stratum_dir, "jaccard_disease_patterns_multi_k.csv")
  )
  
  # 5.10 Save the hclust object for reproducibility
  saveRDS(
    hc,
    file.path(stratum_dir, "hclust_ward_D2.rds")
  )
  
  message(
    "Exported: ",
    file.path(stratum_dir, "jaccard_disease_patterns_multi_k.csv")
  )
}