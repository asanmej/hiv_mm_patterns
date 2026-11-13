# -----------------------------------------------------------------------------
# 03b_hierarchical_dendrograms.R
#
# Description:
#   Full dendrogram (NO cut) of the stratified hierarchical clustering of every
#   stratum, so that the clinicians can see how the diseases join and decide
#   which number of groups is worth looking at.
#
#   Input : PATH_HIERARCHICAL_STRATIFIED_MODELS/<sex>/<age>/<status>/hclust_ward_D2.rds
#           (created by 03_hierarchical_clustering_stratified.R)
#   Output: PATH_HIERARCHICAL_STRATIFIED_FIGURES/Dendrograms/
#             dendrogram_<stratum>.png         one image per stratum
#             dendrograms_all_strata.pdf       all strata, one page each
#
#   Height of the tree = Ward.D2 distance on Jaccard distances between diseases:
#   the higher two branches join, the more different the groups are.
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(factoextra)

dendro_dir <- file.path(PATH_HIERARCHICAL_STRATIFIED_FIGURES, "Dendrograms")
dir.create(dendro_dir, recursive = TRUE, showWarnings = FALSE)

strata <- expand.grid(
  sexo_val   = c("MALE", "FEMALE"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

# Number of patients of each stratum (shown in the title: with few patients the
# distances between diseases are less reliable)
n_patients_file <- file.path(PATH_HIERARCHICAL_STRATIFIED_TABLES, "k_selection_hierarchical.csv")
n_patients <- if (file.exists(n_patients_file)) {
  readr::read_csv(n_patients_file, show_col_types = FALSE) %>%
    distinct(sexo_val, edad_val, status_val, n_patients)
} else NULL

plots <- list()
for (i in seq_len(nrow(strata))) {
  sexo <- strata$sexo_val[i]; edad <- strata$edad_val[i]; status <- strata$status_val[i]
  stratum <- paste(sexo, edad, status, sep = "_")
  rds <- file.path(PATH_HIERARCHICAL_STRATIFIED_MODELS, sexo, edad, status, "hclust_ward_D2.rds")
  if (!file.exists(rds)) { message("No hclust object for ", stratum); next }
  
  hc <- readRDS(rds)
  hc$labels <- gsub("_", " ", hc$labels)
  
  n_pat <- if (!is.null(n_patients)) {
    n_patients$n_patients[n_patients$sexo_val == sexo & n_patients$edad_val == edad &
                          n_patients$status_val == status][1]
  } else NA
  
  p <- fviz_dend(
    hc, k = NULL, horiz = TRUE, cex = 0.7, rect = FALSE, lwd = 0.5,
    main = paste0("Hierarchical clustering of diseases (Jaccard, Ward.D2) - ", stratum,
                  if (!is.na(n_pat)) paste0(" | N = ", format(n_pat, big.mark = ",")) else ""),
    ylab = "Height (Ward.D2 on Jaccard distance)"
  )
  
  height <- max(6, 0.22 * length(hc$labels))
  ggsave(file.path(dendro_dir, paste0("dendrogram_", stratum, ".png")),
         p, width = 11, height = height, dpi = 150, limitsize = FALSE)
  plots[[stratum]] <- list(plot = p, height = height)
}

# One PDF with all the strata (one page each)
if (length(plots) > 0) {
  pdf(file.path(dendro_dir, "dendrograms_all_strata.pdf"), width = 11, height = 14)
  for (x in plots) print(x$plot)
  dev.off()
  message("Dendrograms saved in: ", dendro_dir)
}
