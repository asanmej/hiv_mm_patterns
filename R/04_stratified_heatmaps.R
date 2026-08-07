# -----------------------------------------------------------------------------
# 04_stratified_heatmaps.R

# Description: 
#   This script performs stratified batch visualization of comorbidity patterns. 
#   By segmenting the cohort into 12 distinct clinical and demographic strata 
#   (combining sex categories, age groups, and HIV status), it computes 
#   disease-to-cluster prevalence matrices and generates individual heatmaps 
#   to evaluate structural variations in multimorbidity across subgroups.
# -----------------------------------------------------------------------------

# Load base datasets and hierarchical pattern definitions
bdu_data <- read_csv(PATH_BDU_FULL)
diag_data <- read_csv(PATH_DIAG_BOOL_PREVALENT_FINAL)
patrones_df <- read_csv(file.path(PATH_INTERMEDIATE, "jaccard_disease_patterns.csv"))

# Define demographic and clinical subgroups
cohort_stratified <- bdu_data %>%
  select(patient_id, sexo, edad, vih_bool) %>%
  mutate(
    # Normalize sex variable
    sexo_clean = case_when(
      str_detect(tolower(sexo), "h") ~ "HOMBRE",
      str_detect(tolower(sexo), "m") ~ "MUJER",
      TRUE ~ as.character(sexo)
    ),
    # Define exact age groups (<45, 45-65, >65)
    grupo_edad = case_when(
      edad < 45 ~ "under_45",
      edad >= 45 & edad <= 65 ~ "45to65",
      edad > 65 ~ "over_65"
    ),
    # Define HIV status label
    status_hiv = if_else(as.logical(vih_bool), "hiv", "ctl")
  ) %>%
  inner_join(diag_data, by = "patient_id")

# Create output directory for figures if it does not exist
PATH_FIGURES <- file.path(PATH_RESULTS, "figures")
if (!dir.exists(PATH_FIGURES)) dir.create(PATH_FIGURES, recursive = TRUE)

# Reshape diagnosis data to long format and merge with pattern assignments
diag_long <- cohort_stratified %>%
  pivot_longer(
    cols = -c(patient_id, sexo, sexo_clean, edad, grupo_edad, vih_bool, status_hiv, vih_dt, HIV_infection, `HIV-related_disease`),
    names_to = "enfermedad",
    values_to = "has_disease"
  ) %>%
  mutate(has_disease = if_else(is.na(has_disease) | has_disease == "0", 0, 1)) %>%
  inner_join(patrones_df, by = "enfermedad")

# Generate the 12 possible combinations (2 Sexes x 3 Age Groups x 2 HIV Statuses)
estratos <- expand.grid(
  sexo_val   = c("HOMBRE", "MUJER"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

# Color palette configuration matching reference standard:
# Soft Red/Pink (-0.3) -> White (0) -> Light Blue (0.3) -> Intense Blue (0.9+)
color_low  <- "#F4BFC4" 
color_mid  <- "#FFFFFF" 
color_high <- "#3C58D0" 

# Iterative rendering and export of stratified heatmaps (PNG format only)
pwalk(estratos, function(sexo_val, edad_val, status_val) {
  
  # Filter data for the specific subgroup stratum
  sub_df <- diag_long %>%
    filter(sexo_clean == sexo_val, grupo_edad == edad_val, status_hiv == status_val)
  
  if (nrow(sub_df) > 0) {
    
    # Calculate percentage prevalence or normalized mean by disease and cluster
    heatmap_matrix <- sub_df %>%
      group_by(enfermedad, cluster) %>%
      summarise(value = mean(has_disease), .groups = "drop") %>%
      complete(enfermedad, cluster = unique(patrones_df$cluster), fill = list(value = 0))
    
    # Define file identifiers and titles
    tag_name <- paste0(sexo_val, "_", edad_val, "_", status_val)
    file_title <- paste("Heatmap -", sexo_val, "| Age:", edad_val, "| Status:", toupper(status_val))
    
    # Render ggplot2 heatmap
    p <- ggplot(heatmap_matrix, aes(
      x = factor(cluster, labels = paste0("Cluster ", sort(unique(cluster)))), 
      y = enfermedad, 
      fill = value)) +
      geom_tile(color = NA) + 
      scale_fill_gradient2(
        low = color_low,
        mid = color_mid,
        high = color_high,
        midpoint = 0,
        limits = c(min(0, min(heatmap_matrix$value)), max(heatmap_matrix$value)),
        name = "value"
      ) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        axis.text.x = element_text(face = "bold", size = 10, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
      ) +
      labs(
        title = file_title,
        x = NULL,
        y = NULL
      )
    
    # Determine target subfolder based on sex (Hombres / Mujeres)
    subfolder_name <- if_else(sexo_val == "HOMBRE", "Hombres", "Mujeres")
    target_dir <- file.path(PATH_FIGURES, subfolder_name)
    
    if (!dir.exists(target_dir)) {
      dir.create(target_dir, recursive = TRUE)
    }
    
    # Save high-resolution PNG version (300 DPI) inside the corresponding subfolder
    png_filename <- file.path(target_dir, paste0("heatmap_", tag_name, ".png"))
    ggsave(png_filename, plot = p, width = 10, height = 12, dpi = 300)
    
    message("Successfully generated PNG: ", tag_name, " -> Subfolder: ", subfolder_name)
  }
})