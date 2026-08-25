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
patient_raw <- read_csv(file.path(PATH_DATA, "bdu_full.csv"))
diagnoses_wide <- read_csv(file.path(PATH_DATA, "diagnosis_bool_prevalent_final.csv"))

# 1.2 Prepare the disease dataset for latent class analysis

# The diagnosis dataset is already provided in wide format. Therefore, each
# row represents one patient and each disease variable represents whether the
# corresponding condition is present.
# Variables related to HIV identification or other metadata are excluded from
# the LCA input because they are not disease indicators to be used for defining
# multimorbidity patterns.
diseases_bl <- diagnoses_wide %>%
  select(-any_of(c("vih_dt", "HIV_infection", "HIV-related_disease")))

# 1.3 Standardize disease variable names
enfermedades_cols <- setdiff(names(diseases_bl), "patient_id")
if(!all(str_detect(enfermedades_cols, "^dis_"))) {
  names(diseases_bl)[names(diseases_bl) %in% enfermedades_cols] <- paste0("dis_", enfermedades_cols)
}

# 1.4 Save the intermediate diagnosis dataset
write_csv(
  diagnoses_wide,
  file.path(PATH_INTERMEDIATE_LCA, "diagnoses_wide.csv")
)

saveRDS(
  diagnoses_wide,
  file.path(PATH_INTERMEDIATE_LCA, "diagnoses_wide.rds")
)

# 1.5 Restrict the analysis to patients with at least two diseases 
diseases_bl <- diseases_bl %>%
  mutate(n_dis = rowSums(select(., contains("dis")), na.rm = TRUE)) %>%
  filter(n_dis >= 2) %>%
  select(-n_dis)

# 1.6 Save the analysis-ready disease dataset
write_csv(
  diseases_bl,
  file.path(PATH_INTERMEDIATE_LCA, "diseases_bl.csv")
)

saveRDS(
  diseases_bl,
  file.path(PATH_INTERMEDIATE_LCA, "diseases_bl.rds")
)


# 2. LATENT CLASS ANALYSIS DATA PREPARATION

# 2.1 LCA methodological requirements

# The LCA assumes local independence: conditional on latent class membership,
# the observed disease indicators are assumed to be statistically independent.

# The analysis requires an adequate number of observed indicators (more than 4 
# variables) and a sufficiently large sample size (more than 300 observations)
# to estimate the latent class model reliably.

# Missing values are handled according to the procedures implemented in the
# MMLCA workflow. The functions used below determine the eligible disease
# indicators and prepare the data accordingly.

# 2.2 Preserve the original patient identifiers
original_ids <- diseases_bl$patient_id
gc()

# 2.3 Prepare the dataset for MMLCA

# MMLCA expects the observed disease variables to be supplied in a specific
# format. The patient identifier is temporarily converted into an "id" variable
# and removed from the disease matrix used for model fitting.

# The identifier is retained separately in original_ids so that it can be
# restored after latent class assignment.

df_for_mmlca <- diseases_bl %>% mutate(id=patient_id) %>% select(-patient_id) 

# 2.4 Identify disease indicators
dis_column_names <- names(df_for_mmlca)[str_detect(names(df_for_mmlca), "^dis_")]

# 2.5 Prepare the data using MMLCA
X <- prepare_data(
  df_for_mmlca,
  dis_cols = dis_column_names
)

# 2.6 Select diseases based on minimum prevalence

# Only diseases with a prevalence of at least 2% are considered for the LCA.
threshold <- 0.02 # 2% prevalence
disease_names <- select_conditions(X, threshold = threshold)

disease_names

# 3. LATENT CLASS MODEL SELECTION

# 3.1 Define the training and test datasets

# A reproducible random seed is set before splitting the data so that the same
# training/test partition can be reproduced in subsequent analyses.
set.seed(202112)

# Approximately 70% of the observations are used to estimate the LCA models,
# while the remaining 30% are retained as an independent test subset for
# evaluating classification performance.
train <- sample(1:nrow(X), round(nrow(X) * 0.7))
test <- setdiff(1:nrow(X), train)

# 3.2 Fit LCA models with different numbers of latent classes

# Models containing between 3 and 7 latent classes are fitted to evaluate
# alternative multimorbidity structures.

# Five repetitions are used for each number of classes to reduce the risk of
# selecting a solution resulting from a local optimum of the likelihood
# function.

# The disease indicators previously selected according to the 2% prevalence
# threshold are used as the observed variables defining the latent classes.
res <- select_number_LCA(
  nclasses = 3:7,
  X = X[train, ],
  conditions = disease_names,
  nrep = 5
)

# 3.3 Inspect model selection and classification results
res$obj$predclass
res$metrics
res$accuracy_matrix
res$elapsed_time
res$plot
res$obj

# 3.4 Save the fitted LCA models
save(
  res,
  disease_names,
  file = file.path(PATH_INTERMEDIATE_LCA, "ResultsLCA_HIV.RData")
)

# 3.5 Reload the fitted models
load(file.path(PATH_INTERMEDIATE_LCA,"ResultsLCA_HIV.RData"))


# 4. MODEL SELECTION VISUALIZATION

# 4.1 Save the comparison of candidate latent class solutions
ggplot2::ggsave(
  filename = file.path(PATH_RESULTS_FIGURES_LCA,"comparacion_clases.png"),
  plot = res$plot,
  width = 6,
  height = 4,
  dpi = 300
)

# 4.2 Evaluate classification accuracy
plot_traintest <- ggaccuracy_LCA(res)

ggplot2::ggsave(
  filename = file.path(PATH_RESULTS_FIGURES_LCA, "comparacion_traintest.png"),
  plot = plot_traintest,
  width = 6,
  height = 4,
  dpi = 300
)


# 5. INTERPRETATION OF THE THREE-CLASS LATENT SOLUTION

# The three-class solution is selected for the subsequent interpretation based
# primarily on its lower BIC compared with the alternative solutions and its
# very high entropy, indicating well-separated latent profiles.
#
# The selected solution is subsequently characterized using three complementary
# approaches:
#
# 1. Observed/Expected prevalence ratios (O/E)
# 2. Observed/Expected disease prevalence
# 3. Disease prevalence profiles across latent classes

# 5.1 Method 1: Observed-to-expected enrichment
OEx_sol3 <- ggOEx(res$obj[[1]], table = FALSE)
ggplot2::ggsave(filename = file.path(PATH_RESULTS_FIGURES_LCA, "metodo1_3LC.jpg"), 
                plot = OEx_sol3, 
                width = 30, 
                height = 20, 
                dpi = 300)

# 5.2 Method 2: Observed-to-expected disease associations
OE_sol3 <- ggOE(res$obj[[1]], table = FALSE)
ggplot2::ggsave(filename = file.path(PATH_RESULTS_FIGURES_LCA, "metodo2_3LC.jpg"), 
                plot = OE_sol3, 
                width = 28, 
                height = 22, 
                dpi = 300)

# 5.3 Method 3: Disease prevalence profiles across latent classes
prev_sol3 <- ggprev_spaghetti(res$obj[[1]])

# 5.4 Customize the prevalence plot
prev_sol3.2 <- prev_sol3 +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                  angle = 45,
                  hjust = 1,
                  margin = ggplot2::margin(t = -5)
                )
              )

prev_sol3.3 <- prev_sol3 +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.length.x = grid::unit(0.1, "cm"),
    plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5)
  )

# 5.5 Save the prevalence profile figure
ggplot2::ggsave(
  filename = file.path(PATH_RESULTS_FIGURES_LCA,"prev_3LC.jpg"),
  plot = prev_sol3.3,
  width = 35,
  height = 15,
  dpi = 300
)


# 6. DISEASE PREVALENCE BY LATENT CLASS

# 6.1 Extract the selected three-class model
obj <- res$obj[[1]]

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
  Overall_prevalence = E,
  O
)

print(tabla_prevalencias)

# 6.6 Export disease prevalence results
writexl::write_xlsx(
  tabla_prevalencias,
  file.path(PATH_RESULTS_TABLES,"tabla_prevalencias_3LC.xlsx")
)


# 7. DESCRIPTIVE CHARACTERIZATION OF THE THREE LATENT CLASSES

# 7.1 Select the final three-class model
modelo3 <- res$obj[[1]]

# 7.2 Harmonize disease names between the input data and the fitted model
# Some disease names are automatically modified during data preparation because
# certain characters are not handled identically by R and by the MMLCA
# preprocessing functions.

# These two variables are therefore renamed so that their names exactly match
# the disease indicators stored in the fitted LCA model
X <- X %>%
  rename(
    dis_Alcohol.related_disorders = `dis_Alcohol-related_disorders`,
    dis_Other_non.traumatic_joint_disorders = `dis_Other_non-traumatic_joint_disorders`
  )

# 7.3 Assign each patient to a latent class
mm_pattern <- assign_LCA(modelo3, X)

# 7.4 Inspect the distribution of latent class assignments
table(mm_pattern)
prop.table(table(mm_pattern)) * 100

# 7.5 Add latent class membership to the patient-level dataset
X$mm_pattern <- mm_pattern
X$id <- original_ids

# 7.6 Convert latent class membership to a factor
X$mm_pattern <- factor(X$mm_pattern, levels = 1:3)


# 8. MERGE LATENT CLASS ASSIGNMENTS WITH PATIENT INFORMATION

# 8.1 Remove duplicate patient records
patient <- patient_raw %>%
  distinct(patient_id, .keep_all = TRUE)

# 8.2 Restrict patient information to the LCA study population
ids_keep <- diseases_bl$patient_id
patient <- patient %>%
  filter(patient_id %in% ids_keep)

# 8.3 Verify the number of eligible patients
length(unique(diseases_bl$patient_id))
length(unique(patient$patient_id))


# 9. DERIVE PATIENT-LEVEL MULTIMORBIDITY VARIABLES

# 9.1 Standardize the patient identifier
X_clasificado <- X %>%
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
tabla_descriptiva3lc <- merged_data %>%
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
    across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  ) %>%
  as.data.frame()

# 11.3 Display the descriptive table
print(tabla_descriptiva3lc)

# 11.4 Export the descriptive table
write_xlsx(
  tabla_descriptiva3lc,
  file.path(PATH_RESULTS_TABLES, "descriptivo_3LC.xlsx")
)
