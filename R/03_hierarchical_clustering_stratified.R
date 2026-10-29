# -----------------------------------------------------------------------------
# 03_hierarchical_clustering_stratified.R
#
# Description:
#   Stratified hierarchical clustering of disease co-occurrence patterns
#   using Jaccard distance and Ward.D2 linkage.
#
#   The analysis is repeated independently for each of the twelve
#   sex × age × HIV status strata:
#       MALE   | under_45 | ctl
#       MALE   | under_45 | hiv
#       MALE   | 45to65   | ctl
#       MALE   | 45to65   | hiv
#       MALE   | over_65  | ctl
#       MALE   | over_65  | hiv
#       FEMALE | under_45 | ctl
#       FEMALE | under_45 | hiv
#       FEMALE | 45to65   | ctl
#       FEMALE | 45to65   | hiv
#       FEMALE | over_65  | ctl
#       FEMALE | over_65  | hiv
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
#   are written under PATH_HIERARCHICAL_STRATIFIED_* (see 00_config.R).
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
  sexo_val   = c("MALE", "FEMALE"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

# Candidate numbers of clusters
k_candidates <- 3:7

# Minimum disease prevalence filter (identical to the original script)
min_prev <- 0.02

# Output folders are defined in 00_config.R (PATH_HIERARCHICAL_STRATIFIED_*)
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

# Choice of the 3 solutions shown to the clinicians (per stratum):
#   - a solution is VALID if every cluster has at least `min_diseases` diseases
#     (a cluster of a single disease is not a pattern);
#   - among the valid ones, the 3 with the highest mean silhouette are selected
#     (ties: smaller k); if fewer than 3 are valid, the best invalid ones fill in.
select_k_hierarchical <- function(k_stats, n_select = 3, min_diseases = 2) {
  k_stats$valid <- k_stats$min_size >= min_diseases
  ord <- order(!k_stats$valid, -k_stats$mean_silhouette, k_stats$k)
  k_stats$selected <- FALSE
  k_stats$selected[ord[seq_len(min(n_select, nrow(k_stats)))]] <- TRUE
  k_stats
}
k_selection_all <- list()

for (i in seq_len(nrow(strata))) {
  
  sexo_actual   <- strata$sexo_val[i]
  edad_actual   <- strata$edad_val[i]
  status_actual <- strata$status_val[i]
  stratum_name  <- paste(sexo_actual, edad_actual, status_actual, sep = "_")
  
  message("\n")
  message("=============================================================")
  message("Processing stratum: ", stratum_name)
  message("=============================================================")
  
  # 5.1 Stratum-specific output folders (figures / tables / models)
  stratum_figures <- file.path(PATH_HIERARCHICAL_STRATIFIED_FIGURES, sexo_actual, edad_actual, status_actual)
  stratum_tables  <- file.path(PATH_HIERARCHICAL_STRATIFIED_TABLES,  sexo_actual, edad_actual, status_actual)
  stratum_models  <- file.path(PATH_HIERARCHICAL_STRATIFIED_MODELS,  sexo_actual, edad_actual, status_actual)
  invisible(lapply(c(stratum_figures, stratum_tables, stratum_models),
                   dir.create, recursive = TRUE, showWarnings = FALSE))
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
      stratum_figures,
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
  # Diagnostics of every candidate solution and selection of 3
  n_leaves <- length(hc$order)
  k_stats <- do.call(rbind, lapply(k_candidates, function(k_val) {
    cl <- cutree(hc, k = k_val)
    data.frame(
      sexo_val = sexo_actual, edad_val = edad_actual, status_val = status_actual,
      k = k_val,
      n_patients = nrow(diag_stratum),
      n_diseases = n_leaves,
      min_size = min(table(cl)),
      max_size = max(table(cl)),
      mean_silhouette = round(mean(cluster::silhouette(cl, jaccard_dist)[, "sil_width"]), 4),
      # jump in merge height between k+1 -> k and k -> k-1 clusters (information)
      height_gap = round(hc$height[n_leaves - k_val + 1] - hc$height[n_leaves - k_val], 4)
    )
  }))
  k_selection_all[[stratum_name]] <- select_k_hierarchical(k_stats)

  write_csv(
    cluster_assignments,
    file.path(stratum_tables, "jaccard_disease_patterns_multi_k.csv")
  )
  
  # 5.10 Save the hclust object for reproducibility
  saveRDS(
    hc,
    file.path(stratum_models, "hclust_ward_D2.rds")
  )
  
  message(
    "Exported: ",
    file.path(stratum_tables, "jaccard_disease_patterns_multi_k.csv")
  )
}

# Selection of solutions of all strata (read by 05_stratified_heatmaps.R)
k_selection_table <- do.call(rbind, k_selection_all)
write_csv(
  k_selection_table,
  file.path(PATH_HIERARCHICAL_STRATIFIED_TABLES, "k_selection_hierarchical.csv")
)
message("Selected solutions saved in k_selection_hierarchical.csv")
