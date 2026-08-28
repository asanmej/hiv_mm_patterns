# 0. PACKAGE MANAGEMENT

## 0.1 CRAN packages
pkg <- c(
  ## Efficient data manipulation
  "data.table",
  ## Writing Excel output files
  "writexl",
  ## Package development and GitHub installation
  "devtools"
)

# 0.2 Install and load CRAN packages
lapply(pkg, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, character.only = TRUE)
    library(x, character.only = TRUE)
  }
})

# 0.3 GitHub packages
# MMLCA is required for the latent class analysis performed in this study.
# It is installed from GitHub because the version used in this analysis is
# distributed through the corresponding GitHub repository.

if (!requireNamespace("MMLCA", quietly = TRUE)) {
  devtools::install_github("ARCbiostat/MMLCA")
}

if (!requireNamespace("ggsankey", quietly = TRUE)) {
  devtools::install_github("davidsjoberg/ggsankey")
}

library(ggsankey)
library(MMLCA)


# 1. DATA IMPORT AND INITIAL PREPARATION

# 1.1 Load patient-level and diagnosis data
patient_raw <- read_csv(PATH_BDU_FULL)
diagnoses_wide <- read_csv(PATH_DIAG_BOOL_PREVALENT_FINAL)

# 1.2 Define demographic and clinical strata
cohort_stratified <- patient_raw %>%
  select(
    patient_id,
    sexo,
    edad,
    vih_bool
  ) %>%
  mutate(
    
    # Standardize sex
    sexo_clean = case_when(
      str_detect(tolower(sexo), "h") ~ "MEN",
      str_detect(tolower(sexo), "m") ~ "WOMEN",
      TRUE ~ as.character(sexo)
    ),
    
    # Define age groups
    grupo_edad = case_when(
      edad < 45 ~ "under_45",
      edad >= 45 & edad <= 65 ~ "45to65",
      edad > 65 ~ "over_65",
      TRUE ~ NA_character_
    ),
    
    # Define HIV/control status
    status_hiv = if_else(as.logical(vih_bool), "hiv","ctl")
  )

# 1.3 Define all candidate strata
k_candidates <- 3:7

strata <- expand.grid(
  sexo_val   = c("MEN", "WOMEN"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

# 1.4 Prepare the disease dataset for latent class analysis

# The diagnosis dataset is already provided in wide format. Therefore, each
# row represents one patient and each disease variable represents whether the
# corresponding condition is present.
# Variables related to HIV identification or other metadata are excluded from
# the LCA input because they are not disease indicators to be used for defining
# multimorbidity patterns.
diseases_bl <- diagnoses_wide %>%
  select(-any_of(c("vih_dt", "HIV_infection", "HIV-related_disease")))

# 1.5 Standardize disease variable names
enfermedades_cols <- setdiff(names(diseases_bl), "patient_id")
if(!all(str_detect(enfermedades_cols, "^dis_"))) {
  names(diseases_bl)[names(diseases_bl) %in% enfermedades_cols] <- paste0("dis_", enfermedades_cols)
}

# 1.6 Save the intermediate diagnosis dataset
write_csv(
  diagnoses_wide,
  file.path(PATH_INTERMEDIATE_LCA, "diagnoses_wide.csv")
)

saveRDS(
  diagnoses_wide,
  file.path(PATH_INTERMEDIATE_LCA, "diagnoses_wide.rds")
)

# 1.7 Restrict the analysis to patients with at least two diseases 
diseases_bl <- diseases_bl %>%
  mutate(n_dis = rowSums(select(., contains("dis")), na.rm = TRUE)) %>%
  filter(n_dis >= 2) %>%
  select(-n_dis)

# 1.8 Save the analysis-ready disease dataset
write_csv(
  diseases_bl,
  file.path(PATH_INTERMEDIATE_LCA, "diseases_bl.csv")
)

saveRDS(
  diseases_bl,
  file.path(PATH_INTERMEDIATE_LCA, "diseases_bl.rds")
)

# 1.9 A reproducible random seed is set
set.seed(202112)

# 2. LATENT CLASS ANALYSIS DATA PREPARATION

# 2.0 Merge disease and stratification information
# Keep only patients present in the multimorbidity analysis dataset
cohort_lca <- cohort_stratified %>%
  inner_join(
    diseases_bl %>%
      select(patient_id),
    by = "patient_id"
  )


# Check the number of patients per stratum
stratum_sizes <- cohort_lca %>%
  count(
    sexo_clean,
    grupo_edad,
    status_hiv,
    name = "n_patients"
  )

print(stratum_sizes)

# 2.1 LCA methodological requirements

# The LCA assumes local independence: conditional on latent class membership,
# the observed disease indicators are assumed to be statistically independent.

# The analysis requires an adequate number of observed indicators (more than 4 
# variables) and a sufficiently large sample size (more than 300 observations)
# to estimate the latent class model reliably.

# Missing values are handled according to the procedures implemented in the
# MMLCA workflow. 

# 2.2 Define analysis parameters:

# Minimum disease prevalence
threshold <- 0.02
# Number of repetitions for each candidate solution
nrep <- 5
# Minimum sample size required to fit an LCA
min_n <- 300

# 2.3 Create containers for results
lca_results <- list()
model_selection_results <- list()
patient_class_results <- list()

# 2.4 Start stratified LCA
for (i in seq_len(nrow(strata))){
  
  # 2.4.0 Define current stratum
  sexo_actual <- strata$sexo_val[i]
  edad_actual <- strata$edad_val[i]
  status_actual <- strata$status_val[i]
  
  stratum_name <- paste(
    sexo_actual,
    edad_actual,
    status_actual,
    sep = "_"
  )
  
  message("\n")
  message("Running LCA for stratum: ", stratum_name)
  message("============================================================")
  
  # 2.4.1 Define output directories for current stratum
  stratum_figures <- file.path(
    PATH_LCA_STRATIFIED_FIGURES,
    sexo_actual,
    edad_actual,
    status_actual
  )
  
  stratum_tables <- file.path(
    PATH_LCA_STRATIFIED_TABLES,
    sexo_actual,
    edad_actual,
    status_actual
  )
  
  stratum_models <- file.path(
    PATH_LCA_STRATIFIED_MODELS,
    sexo_actual,
    edad_actual,
    status_actual
  )
  
  stratum_patient_data <- file.path(
    PATH_LCA_STRATIFIED_PATIENT_DATA,
    sexo_actual,
    edad_actual,
    status_actual
  )
  
  dir.create(
    stratum_figures,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    stratum_tables,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    stratum_models,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    stratum_patient_data,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # 2.4.2 Select patients belonging to current stratum
  ids_stratum <- cohort_lca %>%
    filter(
      sexo_clean == sexo_actual,
      grupo_edad == edad_actual,
      status_hiv == status_actual
    ) %>%
    pull(patient_id)
  
  # Keep only patients belonging to the current stratum
  diseases_stratum <- diseases_bl %>%
    filter(patient_id %in% ids_stratum)
  
  # 2.4.3 Check minimum sample size
  n_stratum <- nrow(diseases_stratum)
  
  message("Number of patients: ",n_stratum)
  
  if (n_stratum < min_n) {
    message("Stratum skipped: fewer than ",min_n," patients.")
    next
  }
  
  # 2.4.4 Preserve the original patient identifiers
  original_ids <- diseases_stratum$patient_id
  
  # 2.4.5 Prepare the dataset for MMLCA
  
  # MMLCA expects the observed disease variables to be supplied in a specific
  # format. The patient identifier is temporarily converted into an "id" variable
  # and removed from the disease matrix used for model fitting.
  
  # The identifier is retained separately in original_ids so that it can be
  # restored after latent class assignment.
  
  df_for_mmlca <- diseases_stratum %>% mutate(id=patient_id) %>% select(-patient_id) 
  
  # 2.4.6 Identify disease indicators
  dis_column_names <- names(df_for_mmlca)[str_detect(names(df_for_mmlca), "^dis_")]
  
  # 2.4.7 Prepare the data using MMLCA
  X <- prepare_data(
    df_for_mmlca,
    dis_cols = dis_column_names
  )
  
  names(X) <- make.names(names(X))
  
  # 2.4.8 Select diseases based on minimum prevalence
  
  # Only diseases with a prevalence of at least 2% are considered for the LCA
  disease_names <- select_conditions(X, threshold = threshold)
  
  message("Number of diseases retained: ",length(disease_names))
  
  # Check whether enough disease indicators remain
  if (length(disease_names) < 5) {
    message("Stratum skipped: fewer than 5 eligible disease indicators.")
    next
  }
  
  # 3.1 Fit LCA models with different numbers of latent classes
  
  # Models containing between 3 and 7 latent classes are fitted using the
  # complete set of eligible patients within the current stratum.
  
  # Five repetitions are used for each number of classes to reduce the risk
  # of selecting a solution resulting from a local optimum of the likelihood
  # function.
  
  # The disease indicators previously selected according to the 2% prevalence
  # threshold are used as the observed variables defining the latent classes.
  res <- select_number_LCA(
    nclasses = k_candidates,
    X = X,
    conditions = disease_names,
    nrep = nrep
  )
  
  # 3.2 Store model-selection results
  lca_results[[stratum_name]] <- res
  model_selection_results[[stratum_name]] <- res$metrics
  
  # 3.3 Save fitted models for current stratum
  save(
    res,
    disease_names,
    original_ids,
    file = file.path(
      stratum_models,
      paste0("ResultsLCA_", stratum_name, ".RData")
    )
  )
  
  # 4. MODEL SELECTION VISUALIZATION
  
  # 4.1 Save model-selection visualization
  ggplot2::ggsave(filename = file.path(
    stratum_figures,
    paste0("comparacion_clases_",stratum_name,".png")
  ),
  plot = res$plot,
  width = 6,
  height = 4,
  dpi = 300
  )
  
  # 4.2 Identify candidate latent class solutions for evaluation
  
  # Three solutions are selected for further evaluation within each stratum:
  # 1) the solution with the minimum AIC,
  # 2) the solution with the minimum entropy,
  # 3) the solution with one additional latent class compared with the
  #    minimum-entropy solution.
  #
  # If two or more criteria identify the same number of classes, the duplicate
  # solution is evaluated only once.
  
  metrics <- as.data.frame(res$metrics)
  
  # Identify AIC column
  aic_column <- names(metrics)[
    str_detect(
      names(metrics),
      regex("^AIC$", ignore_case = TRUE)
    )
  ]
  
  if (length(aic_column) == 0) {
    stop("AIC column could not be identified in res$metrics for ",stratum_name)
  }
  
  # Identify entropy column
  entropy_column <- names(metrics)[
    str_detect(
      names(metrics),
      regex("entropy", ignore_case = TRUE)
    )
  ]
  
  if (length(entropy_column) == 0) {
    stop("Entropy column could not be identified in res$metrics for ",stratum_name)
  }
  
  # Identify number-of-classes column
  nclass_column <- names(metrics)[
    str_detect(
      names(metrics),
      regex("nclass|nclasses|classes|k", ignore_case = TRUE)
    )
  ]
  
  if (length(nclass_column) == 0) {
    stop("Number-of-classes column could not be identified in res$metrics for ",stratum_name)
  }
  
  # Convert number of classes to numeric
  nclass_values <- as.numeric(as.character(metrics[[nclass_column[1]]]))
  
  # Solution with minimum AIC
  k_AIC <- nclass_values[which.min(metrics[[aic_column[1]]])]
  
  # Solution with minimum entropy
  k_entropy <- nclass_values[which.min(metrics[[entropy_column[1]]])]
  
  # Solution with one additional class than the minimum-entropy solution
  k_entropy_plus1 <- k_entropy + 1
  
  # Keep only candidate solutions within the fitted range
  k_to_evaluate <- unique(c(k_AIC, k_entropy, k_entropy_plus1))
  
  k_to_evaluate <- k_to_evaluate[k_to_evaluate %in% k_candidates]
  
  message("Minimum AIC: k = ", k_AIC)
  message("Minimum entropy: k = ", k_entropy)
  message("Minimum entropy + 1: k = ", k_entropy_plus1)
  message("Models selected for evaluation: k = ",paste(k_to_evaluate, collapse = ", "))
  
  # Initialize container for selected solutions in the current stratum
  patient_class_results[[stratum_name]] <- list()
  
  # Evaluate each selected latent class solution
  for (k_eval in k_to_evaluate) {
    
    message("\n")
    message("Evaluating ", k_eval," latent classes for stratum: ", stratum_name)
    message("------------------------------------------------------------")
    
    # Identify the fitted model corresponding to the selected number of classes
    model_index <- which(nclass_values == k_eval)
    
    if (length(model_index) != 1) {
      stop("Could not uniquely identify the model with k = ",k_eval," for ", stratum_name)
    }
    
    modelo_evaluado <- res$obj[[model_index]]
    
    nclass_evaluado <- k_eval
    
    # 5. INTERPRETATION OF THE SELECTED LATENT CLASS SOLUTIONS
    
    # 5.1 Method 1: Observed-to-expected enrichment
    OEx_optimo <- ggOEx(modelo_evaluado , table = FALSE)
    ggplot2::ggsave(filename = file.path(stratum_figures, paste0("metodo1",nclass_evaluado,"LC_",stratum_name,".jpg")), 
                    plot = OEx_optimo, 
                    width = 30, 
                    height = 20, 
                    dpi = 300)
    
    # 5.2 Method 2: Observed-to-expected disease associations
    OE_optimo <- ggOE(modelo_evaluado , table = FALSE)
    ggplot2::ggsave(filename = file.path(stratum_figures, paste0("metodo2",nclass_evaluado,"LC_",stratum_name,".jpg")), 
                    plot = OE_optimo, 
                    width = 28, 
                    height = 22, 
                    dpi = 300)
    
    # 5.3 Method 3: Disease prevalence profiles across latent classes
    prev_optimo <- ggprev_spaghetti(modelo_evaluado )
    
    # 5.4 Customize the prevalence plot
    prev_optimo <- prev_optimo +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.length.x = grid::unit(0.1,"cm"),
        plot.margin = ggplot2::margin(t = 5,r = 5,b = 0,l = 5)
      )
    
    # 5.5 Save the prevalence profile figure
    ggplot2::ggsave(
      filename = file.path(stratum_figures,paste0("prev_",nclass_evaluado,"LC_",stratum_name,".jpg")),
      plot = prev_optimo,
      width = 35,
      height = 15,
      dpi = 300
    )
    
    
    # 6. DISEASE PREVALENCE BY LATENT CLASS
    
    # 6.1 Extract the selected optimal model
    obj <- modelo_evaluado
    
    # 6.2 Determine the number of latent classes
    nclass <- nrow(obj$probs[[1]])
    
    # 6.3 Calculate overall disease prevalence
    E <- apply(obj$y - 1, 2, mean)
    
    # 6.4 Calculate disease prevalence within each latent class
    n <- list()
    for (j in 1:nclass) {
      n[[j]] <- apply(obj$y[obj$predclass == j, ] - 1, 2, mean)
    }
    
    # 6.5 Combine prevalence estimates into a single table
    O <- do.call(cbind, n)
    rownames(O) <- colnames(obj$y)
    
    tabla_prevalencias <- data.frame(
      Disease = colnames(obj$y),
      Overall_prevalence = E * 100,
      O * 100
    ) %>%
      # Format percentages to one decimal place
      mutate(
        across(
          where(is.numeric),
          ~ round(.x, 1)
        )
      )
    
    print(tabla_prevalencias)
    
    # 6.6 Export disease prevalence results
    writexl::write_xlsx(
      tabla_prevalencias,
      file.path(stratum_tables,paste0("tabla_prevalencias_",nclass_evaluado,"LC_",stratum_name,".xlsx"))
    )
    
    # 7. DESCRIPTIVE CHARACTERIZATION OF THE SELECTED LATENT CLASSES
    
    # 7.1 Prepare the dataset for latent class assignment
    X_class <- X
    
    # 7.2 Harmonize disease names between the input data and the fitted model
    # Some disease names are automatically modified during data preparation because
    # certain characters are not handled identically by R and by the MMLCA
    # preprocessing functions.
    
    # 7.3 Assign each patient to a latent class
    mm_pattern <- assign_LCA(modelo_evaluado , X_class)
    
    # 7.4 Inspect the distribution of latent class assignments
    table(mm_pattern)
    prop.table(table(mm_pattern)) * 100
    
    # 7.5 Add latent class membership to the patient-level dataset
    X_class$mm_pattern <- mm_pattern
    X_class$id <- original_ids
    
    # 7.6 Convert latent class membership to a factor
    X_class$mm_pattern <- factor(X_class$mm_pattern, levels = 1:nclass_evaluado)
    
    
    # 8. MERGE LATENT CLASS ASSIGNMENTS WITH PATIENT INFORMATION
    
    # 8.1 Remove duplicate patient records
    patient <- patient_raw %>%
      distinct(patient_id, .keep_all = TRUE)
    
    # 8.2 Restrict patient information to the LCA study population
    ids_keep <- diseases_stratum$patient_id
    patient <- patient %>%
      filter(patient_id %in% ids_keep)
    
    # 8.3 Verify the number of eligible patients
    length(unique(diseases_stratum$patient_id))
    length(unique(patient$patient_id))
    
    
    # 9. DERIVE PATIENT-LEVEL MULTIMORBIDITY VARIABLES
    
    # 9.1 Standardize the patient identifier
    X_clasificado <- X_class %>%
      rename(patient_id = "id")
    
    # 9.2 Identify disease variables
    columnas_enfermedades <- setdiff(
      names(X_clasificado),
      c("patient_id", "mm_pattern")
    )
    
    # 9.3 Convert MMLCA disease coding to binary indicators
    
    # MMLCA uses its own coding for binary disease indicators. The values are
    # converted to a conventional 0/1 representation:
    # 1 -> 0 = absence of disease
    # 2 -> 1 = presence of disease
    X_clasificado <- X_clasificado %>%
      mutate(
        across(
          all_of(columnas_enfermedades),
          ~ ifelse(. == 1, 0, ifelse(. == 2, 1, .))
        )
      )
    
    # 9.4 Calculate the number of diseases per patient
    X_clasificado <- X_clasificado %>%
      mutate(
        num_enfermedades = rowSums(
          across(all_of(columnas_enfermedades)),
          na.rm = TRUE
        )
      )
    
    # 10. MERGE CLINICAL, DEMOGRAPHIC AND LCA INFORMATION
    
    # 10.1 Merge latent class assignments with patient characteristics
    merged_data <- merge(
      X_clasificado,
      patient,
      by = "patient_id",
      all.x = TRUE
    )
    
    # 10.2 Retain variables required for descriptive characterization
    merged_data <- merged_data %>%
      select(
        patient_id,
        sexo,
        nacionalidad,
        zbs_tipo,
        num_enfermedades,
        mm_pattern,
        year
      ) %>%
      distinct()
    
    # 10.3 Inspect the resulting dataset
    summary(merged_data)
    
    
    # 11. DESCRIPTIVE CHARACTERIZATION OF LATENT CLASSES
    
    # 11.1 Calculate age at the end of the study period
    
    # year represents the patient's year of birth. Because the study period ended
    # on 31 December 2022, age at the end of follow-up can subsequently be
    # calculated as 2022 minus the year of birth. Therefore, age is calculated
    # using 2022 as the reference year: age at end of study = 2022 - year of birth
    merged_data <- merged_data %>%
      mutate(
        edad_2022 = 2022 - year
      )
    
    # 11.2 Calculate descriptive statistics by latent class
    tabla_descriptiva <- merged_data %>%
      group_by(mm_pattern) %>%
      summarise(
        n = n(),
        
        Porc_Hombres = mean(
          sexo == "HOMBRE",
          na.rm = TRUE
        ) * 100,
        
        Porc_Mujeres = mean(
          sexo == "MUJER",
          na.rm = TRUE
        ) * 100,
        
        Porc_España = mean(
          nacionalidad == "ESPAÑA",
          na.rm = TRUE
        ) * 100,
        
        Porc_Urbano = mean(
          zbs_tipo == "Urbano",
          na.rm = TRUE
        ) * 100,
        
        Porc_Rural = mean(
          zbs_tipo == "Rural",
          na.rm = TRUE
        ) * 100,
        
        Media_Enfermedades = mean(
          num_enfermedades,
          na.rm = TRUE
        ),
        
        Media_Edad_2022 = mean(
          edad_2022,
          na.rm = TRUE
        ),
        
        .groups = "drop"
      ) %>%
      mutate(
        # Percentages: one decimal
        across(
          c(
            Porc_Hombres,
            Porc_Mujeres,
            Porc_España,
            Porc_Urbano,
            Porc_Rural
          ),
          ~ round(.x, 1)
        ),
        # Continuous variables: three decimals
        across(
          c(
            Media_Enfermedades,
            Media_Edad_2022
          ),
          ~ round(.x, 3)
        ),
        # Sample size: integer
        n = as.integer(n)
      ) %>%
      as.data.frame()
    
    # 11.3 Display the descriptive table
    print(tabla_descriptiva)
    
    # 11.4 Export the descriptive table
    write_xlsx(
      tabla_descriptiva,
      file.path(stratum_tables, paste0("descriptivo_",nclass_evaluado,"LC_",stratum_name,".xlsx"))
    )
    
    # 11.5 Save patient-level class assignments
    patient_class_results[[stratum_name]][[paste0("k_", k_eval)]] <- merged_data
    
    write_csv(
      merged_data,
      file.path(stratum_patient_data,paste0("patient_classes_h",k_eval,"_",stratum_name,".csv"))
    )
    
    # 11.6 Save complete results for current stratum
    saveRDS(
      list(
        stratum = stratum_name,
        k = k_eval,
        sexo = sexo_actual,
        edad = edad_actual,
        status_hiv = status_actual,
        n_patients = n_stratum,
        n_diseases = length(disease_names),
        disease_names = disease_names,
        model = modelo_evaluado ,
        n_classes = nclass_evaluado,
        metrics = res$metrics,
        prevalence = tabla_prevalencias,
        descriptive = tabla_descriptiva,
        patient_data = merged_data
      ),
      file = file.path(stratum_models,
                       paste0("LCA_k",k_eval, "_",stratum_name,".rds")
      )
    )
    
  }
}
# 12. SAVE GLOBAL STRATIFIED RESULTS

saveRDS(
  lca_results,
  file = file.path(PATH_LCA_STRATIFIED_MODELS,"all_LCA_results.rds")
)

saveRDS(
  model_selection_results,
  file = file.path(PATH_LCA_STRATIFIED_MODELS,"all_model_selection_results.rds")
)

saveRDS(
  patient_class_results,
  file = file.path(PATH_LCA_STRATIFIED_PATIENT_DATA, "all_patient_class_results.rds")
)

# 13. COMBINE MODEL SELECTION RESULTS

model_selection_table <- bind_rows(
  lapply(
    names(model_selection_results),
    function(x) {
      df <- as.data.frame(model_selection_results[[x]])
      df$stratum <- x
      df
    }
  )
)

model_selection_table <- model_selection_table %>%
  mutate(
    # Convert numeric columns from character to numeric
    across(
      c(
        nclass,
        log_likelihood,
        df,
        BIC,
        AIC,
        CAIC,
        ABIC,
        likelihood_ratio,
        `Assignment accuracy (%)`,
        Entropy
      ),
      as.numeric
    ),
    
    # Format numerical values
    across(
      c(
        log_likelihood,
        BIC,
        AIC,
        CAIC,
        ABIC,
        likelihood_ratio
      ),
      ~round(.x, 3)
    ),
    
    # Convert proportions to percentages
    `Assignment accuracy (%)` = round(`Assignment accuracy (%)` * 100, 1),
    Entropy = round(Entropy * 100, 1)
  )

write_xlsx(
  model_selection_table,
  file.path(PATH_LCA_STRATIFIED_TABLES,"model_selection_all_strata.xlsx")
)
