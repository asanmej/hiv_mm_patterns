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
# Load disease-to-cluster assignments for all candidate k solutions
patrones_multi_df <- read_csv(file.path(PATH_INTERMEDIATE, "jaccard_disease_patterns_multi_k.csv"))
# Load patient-level pattern scores and dominant pattern assignments
pattern_data <- read_csv(file.path(PATH_INTERMEDIATE, "patient_pattern_scores_hierarchical_k3_k7.csv"))


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

PATH_FIGURES <- file.path(PATH_RESULTS, "figures")
if (!dir.exists(PATH_FIGURES)) dir.create(PATH_FIGURES, recursive = TRUE)

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

# Define a fixed discrete color scale for the Observed/Expected ratio.
# The same intervals are used across all heatmaps to ensure comparability
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
    inner_join(
      patient_clusters,
      by = "patient_id"
    )
  
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
  diag_long <- cohort_stratum %>%
    select(
      patient_id,
      cluster,
      all_of(disease_columns)
    ) %>%
    pivot_longer(
      cols = all_of(disease_columns),
      names_to = "enfermedad_col",
      values_to = "has_disease"
    ) %>%
    mutate(
      # Convert disease status to a binary indicator:
      # 1 = disease present, 0 = disease absent
      has_disease = if_else(
        is.na(has_disease) |
          has_disease == "0",
        0,
        1
      ),
      has_disease = as.numeric(
        has_disease
      )
    ) %>%
    # Recover the original disease name after reshaping the diagnostic matrix
    inner_join(
      patrones_subset %>%
        select(
          enfermedad_col,
          enfermedad
        ),
      by = "enfermedad_col"
    )
  
  # Calculate the expected disease prevalence within the current stratum.
  # This represents the proportion of all patients in the stratum with
  # the disease, regardless of their assigned patient cluster.
  global_prev <- diag_long %>%
    group_by(
      enfermedad
    ) %>%
    summarise(
      expected_prev = mean(
        has_disease,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    # Remove diseases that are completely absent from the current stratum
    filter(
      expected_prev > 0
    )
  
  # Calculate observed disease prevalence within each patient cluster:
  # For each disease and patient cluster:
  #   observed prevalence =
  #   number of patients in the cluster with the disease /
  #   total number of patients assigned to the cluster.
  cluster_prev <- diag_long %>%
    group_by(
      enfermedad,
      cluster
    ) %>%
    summarise(
      observed_prev = mean(
        has_disease,
        na.rm = TRUE
      ),
      n_patients = n(),
      n_cases = sum(
        has_disease,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  # Calculate the Observed/Expected enrichment ratio:
  # O/E = P(Disease | Patient cluster) / P(Disease | Stratum)
  # Interpretation:
  #   O/E = 1  -> disease prevalence matches the stratum-level expectation
  #   O/E > 1  -> disease is enriched within the patient cluster
  #   O/E < 1  -> disease is less prevalent than expected
  # The ratio is set to NA when the expected prevalence is zero because
  # enrichment cannot be defined when the disease is absent from the stratum.
  heatmap_matrix <- cluster_prev %>%
    left_join(
      global_prev,
      by = "enfermedad"
    ) %>%
    
    mutate(
      Observed_Expected_ratio = case_when(
        expected_prev > 0 ~
          observed_prev / expected_prev,
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
  
  # Complete the disease x patient-cluster grid:
  # Missing combinations are retained as NA rather than being converted to zero.
  # This distinction is important:
  #   NA = no valid Observed/Expected estimate is available
  #   0  = the disease has an observed prevalence of zero in the cluster
  clusters_expected <- as.character(
    seq_len(k_val)
  )
  
  heatmap_matrix <- heatmap_matrix %>%
    complete(
      enfermedad = unique(
        global_prev$enfermedad
      ),
      cluster = clusters_expected
    )
  
  # Identify the structural disease cluster assigned by hierarchical clustering:
  # The black border is drawn when the patient cluster corresponds to the
  # disease's original hierarchical cluster assignment.
  # The border therefore represents structural cluster membership and is
  # independent of the Observed/Expected enrichment value.
  disease_cluster_map <- patrones_subset %>%
    select(
      enfermedad,
      cluster_enfermedad
    ) %>%
    distinct()
  
  heatmap_matrix <- heatmap_matrix %>%
    left_join(
      disease_cluster_map,
      by = "enfermedad"
    )
  
  # Order diseases and patient clusters 
  sorted_diseases <- sort(
    unique(
      heatmap_matrix$enfermedad
    )
  )
  
  heatmap_matrix <- heatmap_matrix %>%
    mutate(
      cluster_num = match(cluster,clusters_expected),
      disease_num = match(enfermedad, sorted_diseases)
    )
  
  # Identify cells corresponding to the structural disease cluster
  highlighted_cells <- heatmap_matrix %>%
    filter(
      !is.na(cluster_enfermedad),
      cluster == cluster_enfermedad
    ) %>%
    distinct(
      enfermedad,
      disease_num,
      cluster,
      cluster_num
    )
  
  # Create Plotly rectangle shapes used to highlight structurally assigned
  # disease-cluster cells
  highlight_shapes <- list()
  if (
    nrow(highlighted_cells) > 0
  ) {
    highlight_shapes <- lapply(
      seq_len(
        nrow(highlighted_cells)
      ),
      function(i) {
        list(
          type = "rect",
          xref = "x",
          yref = "y",
          x0 =
            highlighted_cells$cluster_num[i] -
            0.5,
          x1 =
            highlighted_cells$cluster_num[i] +
            0.5,
          y0 =
            highlighted_cells$disease_num[i] -
            0.46,
          y1 =
            highlighted_cells$disease_num[i] +
            0.46,
          line = list(
            color = "#000000",
            width = 3
          ),
          fillcolor = "rgba(0,0,0,0)",
          layer = "above"
        )
      }
    )
  }
  
  # Check whether valid Observed/Expected ratios are available for the
  # current stratum and clustering solution
  max_ratio <- max(
    heatmap_matrix$value,
    na.rm = TRUE
  )
  
  if (!is.finite(max_ratio)) {
    message(
      "No valid ratios for: ",
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
  
  # Generate a unique identifier and title for the current heatmap
  tag_name <- paste0(
    sexo_val,
    "_",
    edad_val,
    "_",
    status_val,
    "_k",
    k_val
  )
  
  file_title <- paste(
    "Observed/Expected Ratio (k=",
    k_val,
    ") -",
    sexo_val,
    "| Age:",
    edad_val,
    "| Status:",
    toupper(status_val)
  )
  
  # Prepare formatted values for the interactive hover information
  heatmap_matrix <- heatmap_matrix %>%
    mutate(
      ratio_text = if_else(
        is.na(value),
        "NA",
        
        paste0(round(value, 2))
      )
    )
  
  # Create the interactive heatmap
  p <- plot_ly(
    data = heatmap_matrix,
    
    x = ~cluster_num,
    y = ~disease_num,
    
    # Use the actual Observed/Expected ratio to determine the cell colour
    z = ~value,
    
    # Use a fixed visual scale from 0 to 5.5:
    # Values above 5 are displayed in the >5 colour category while their
    # exact values remain available in the hover information
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
      n_cases,
      
      "<br><b>Disease cluster:</b> ",
      ifelse(
        is.na(cluster_enfermedad),
        "NA",
        cluster_enfermedad
      )
    ),
    
    hovertemplate =
      paste0(
        "%{text}",
        "<br><br>",
        "<b>Patient cluster:</b> Cl. %{x}",
        "<extra></extra>"
      ),
    
    # Configure the colorbar with fixed tick positions corresponding to
    # the predefined Observed/Expected ratio intervals
    colorbar = list(
      title = list(
        text = "Observed_Expected ratio",
        side = "top"
      ),
      
      tickmode = "array",
      
      # Exact interval boundaries displayed on the colorbar
      tickvals = c(
        0,
        0.5,
        0.9,
        1.1,
        2,
        3,
        4,
        5,
        5.25
      ),
      
      ticktext = c(
        "0",
        "0.5",
        "0.9",
        "1.1",
        "2",
        "3",
        "4",
        "5",
        ">5"
      ),
      
      len = 0.8,
      
      # Increase colorbar thickness to improve visual separation between
      # the discrete enrichment intervals
      thickness = 20,
      
      # Position tick labels outside the colorbar
      ticklabelposition = "outside",
      
      # Prevent Plotly from generating automatic tick marks
      nticks = 9
    )
  ) %>%
    
    layout(
      title = list(
        text = file_title,
        font = list(
          size = 11
        )
      ),
      
      xaxis = list(
        title = "",
        tickmode = "array",
        tickvals = seq_len(k_val),
        ticktext = paste0(
          "Cl. ",
          clusters_expected
        ),
        
        range = c(
          0.5,
          k_val + 0.5
        ),
        
        fixedrange = TRUE
      ),
      
      yaxis = list(
        title = "",
        tickmode = "array",
        tickvals = seq_along(
          sorted_diseases
        ),
        ticktext = sorted_diseases,
        
        autorange = "reversed",
        
        range = c(
          length(sorted_diseases) + 0.5,
          0.5
        ),
        
        tickfont = list(
          size = 8
        ),
        
        fixedrange = TRUE
      ),
      
      shapes = highlight_shapes,
      
      margin = list(
        l = 240,
        b = 50,
        r = 100
      )
    )
  
  # Save the interactive heatmap as a self-contained HTML file
  subfolder_name <- file.path(
    if_else(
      sexo_val == "MALE",
      "Males",
      "Females"
    ),
    paste0(
      "k_",
      k_val
    )
  )
  
  target_dir <- file.path(
    PATH_FIGURES,
    subfolder_name
  )
  
  if (!dir.exists(target_dir)) {
    
    dir.create(
      target_dir,
      recursive = TRUE
    )
  }
  
  html_filename <- file.path(
    target_dir,
    paste0(
      "heatmap_interactive_",
      tag_name,
      ".html"
    )
  )
  
  saveWidget(
    p,
    html_filename,
    selfcontained = TRUE
  )
  
  # Report generation status
  message(
    "Generated: ",
    tag_name,
    " | N patients = ",
    nrow(cohort_stratum),
    " | Max O/E ratio = ",
    round(max_ratio, 2)
  )
}
)
