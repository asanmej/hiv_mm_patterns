# -----------------------------------------------------------------------------
# MAIN PIPELINE EXECUTION

# Objective: To sequentially execute the entire flow of comorbidity pattern analysis and validation in the HIV cohort compared to controls.

# Execution Flow:
# 1. Load the general project configuration (libraries and paths)
# 2. Render the pathology clustering using Jaccard/Hierarchical
# 3. Render the pattern score calculation per patient
# 4. Render the epidemiological validation using clogit (HIV vs. Controls)
# -----------------------------------------------------------------------------

# 1. Project and Environment Configuration
source("00_config.R")

# 2. Unsupervised clustering of diseases (Hierarchical Clustering)
render(input = "01_hierarchical_clustering.Rmd",output_dir = PATH_RESULTS)

# 3. Mapping and imputation of pattern scores by patient
render(input = "02_patient_pattern_scoring.Rmd", output_dir = PATH_RESULTS)

# 4. Epidemiological validation (Conditional Logistic Regression)
render(input = "03_epidemiological_validation.Rmd", output_dir = PATH_RESULTS)
