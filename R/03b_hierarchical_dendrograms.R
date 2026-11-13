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
#
#   Layout: the tree is drawn by hand (no factoextra) so that
#     - the height axis does NOT start at 0: all merges happen far from 0 and
#       the leaves used to hang along a long empty stretch. The axis starts just
#       below the first merge (the axis is marked as truncated in the caption);
#     - there is room on the right for the longest disease name (labels were
#       cut off) and every stratum has the same width.
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

dendro_dir <- file.path(PATH_HIERARCHICAL_STRATIFIED_FIGURES, "Dendrograms")
dir.create(dendro_dir, recursive = TRUE, showWarnings = FALSE)

AXIS_STEP     <- 0.05   # the start of the height axis is rounded down to this step
PLOT_WIDTH_IN <- 7      # width (inches) of the tree itself, without the labels
LABEL_SIZE    <- 3      # size of the disease names (mm)
IN_PER_CHAR   <- 0.065  # approximate width of one character at LABEL_SIZE

#' Segments of a dendrogram drawn horizontally (root on the left, leaves on the right).
#' x = -height (so that the root is on the left), y = position of the leaf.
#' Leaves are drawn from `h_leaf` (the start of the truncated axis), not from 0.
dendrogram_segments <- function(hc, h_leaf) {
  n   <- length(hc$order)
  pos <- integer(n); pos[hc$order] <- seq_len(n)     # position of each observation
  node_y <- numeric(n - 1); node_h <- hc$height
  seg <- vector("list", n - 1)
  for (i in seq_len(n - 1)) {
    ch <- hc$merge[i, ]
    y  <- ifelse(ch < 0, pos[abs(ch)], node_y[pmax(ch, 1)])
    h  <- ifelse(ch < 0, h_leaf,       node_h[pmax(ch, 1)])
    node_y[i] <- mean(y)
    seg[[i]] <- data.frame(
      x = -c(h[1], hc$height[i], hc$height[i]),
      xend = -c(hc$height[i], hc$height[i], h[2]),
      y = c(y[1], y[1], y[2]),
      yend = c(y[1], y[2], y[2]))
  }
  do.call(rbind, seg)
}

# Number of patients of each stratum (shown in the title: with few patients the
# distances between diseases are less reliable)
n_patients_file <- file.path(PATH_HIERARCHICAL_STRATIFIED_TABLES, "k_selection_hierarchical.csv")
n_patients <- if (file.exists(n_patients_file)) {
  readr::read_csv(n_patients_file, show_col_types = FALSE) %>%
    distinct(sexo_val, edad_val, status_val, n_patients)
} else NULL

strata <- expand.grid(
  sexo_val   = c("MALE", "FEMALE"),
  edad_val   = c("under_45", "45to65", "over_65"),
  status_val = c("ctl", "hiv"),
  stringsAsFactors = FALSE
)

# 1. Read all trees first: the label margin is the same for every stratum
trees <- list()
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
  trees[[stratum]] <- list(hc = hc, n_pat = n_pat)
}
if (length(trees) == 0) stop("No hclust objects found in ", PATH_HIERARCHICAL_STRATIFIED_MODELS)

max_chars <- max(vapply(trees, function(t) max(nchar(t$hc$labels)), numeric(1)))
fig_width <- PLOT_WIDTH_IN + 0.3 + max_chars * IN_PER_CHAR

# 2. Draw
make_dendrogram <- function(hc, stratum, n_pat) {
  n      <- length(hc$labels)
  h_lo   <- max(0, floor((min(hc$height) - AXIS_STEP) / AXIS_STEP) * AXIS_STEP)
  h_hi   <- max(hc$height)
  seg    <- dendrogram_segments(hc, h_lo)
  leaves <- data.frame(y = seq_len(n), label = hc$labels[hc$order])
  
  ggplot() +
    geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = 0.4, lineend = "round") +
    geom_text(data = leaves, aes(x = -h_lo, y = y, label = label),
              hjust = 0, nudge_x = 0.004 * (h_hi - h_lo), size = LABEL_SIZE) +
    scale_x_continuous(breaks = -pretty(c(h_lo, h_hi), n = 8),
                       labels = function(b) format(abs(b), nsmall = 1),
                       expand = expansion(mult = c(0.02, 0))) +
    scale_y_continuous(expand = expansion(add = 0.6)) +
    coord_cartesian(xlim = c(-h_hi, -h_lo), clip = "off") +
    labs(title = paste0("Hierarchical clustering of diseases (Jaccard, Ward.D2) - ", stratum,
                        if (!is.na(n_pat)) paste0(" | N = ", format(n_pat, big.mark = ",")) else ""),
         x = "Height (Ward.D2 on Jaccard distance)",
         caption = sprintf("Height axis truncated: it starts at %.2f (first merge at %.2f)",
                           h_lo, min(hc$height))) +
    theme_minimal(base_size = 10) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          panel.grid = element_blank(), panel.grid.major.x = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.margin = margin(8, max_chars * IN_PER_CHAR * 72, 8, 8, unit = "pt"))
}

plots <- list()
for (stratum in names(trees)) {
  t <- trees[[stratum]]
  p <- make_dendrogram(t$hc, stratum, t$n_pat)
  height <- max(5, 0.2 * length(t$hc$labels) + 1.2)
  ggsave(file.path(dendro_dir, paste0("dendrogram_", stratum, ".png")),
         p, width = fig_width, height = height, dpi = 150, limitsize = FALSE)
  plots[[stratum]] <- p
}

# One PDF with all the strata (one page each, same size)
pdf(file.path(dendro_dir, "dendrograms_all_strata.pdf"), width = fig_width,
    height = max(5, 0.2 * max(vapply(trees, function(t) length(t$hc$labels), numeric(1))) + 1.2))
for (p in plots) print(p)
dev.off()
message("Dendrograms saved in: ", dendro_dir)
