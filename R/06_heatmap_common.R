# -----------------------------------------------------------------------------
# 06_heatmap_common.R
#
# Description:
#   Shared functions to draw the Observed/Expected (O/E) heatmaps of the three
#   grouping methods (hierarchical, LCA, EFA), so that all of them use exactly
#   the same statistics, colours and marks and can be compared side by side.
#
#   Usage (from 04_, 04b_, 05_ and 07_ scripts):
#       source(file.path(PATH_R, "06_heatmap_common.R"))   # or the right path
#       make_oe_heatmap(...)
#
# O/E definition (identical for the three methods):
#       O/E = prevalence of the disease in the patient cluster /
#             prevalence of the disease in the WHOLE stratum
#   The whole stratum includes patients that could not be assigned to any
#   cluster (ties / no disease of any pattern), which are reported in the title.
#
# Colour scale (discrete, fixed). Bins of O/E:
#   0-0.5 strong blue | 0.5-0.9 light blue | 0.9-1.1 white (neutral)
#   1.1-2 | 2-3 | 3-4 | 4-5 | 5-6 | >6   (light red -> very strong red)
#   Grey background = no O/E could be estimated (cluster without patients).
#
# Marks:
#   - Diagonal grey stripes over a whole column: cluster with < MIN_PATIENTS
#     patients (default 30). O/E of that cluster is unstable.
#   - "†" in a cell: fewer than MIN_CASES (default 10) patients with the
#     disease inside the cluster.
#   - "ns" in a coloured cell (O/E outside 0.9-1.1): difference NOT significant
#     after FDR (Benjamini-Hochberg q >= 0.05), cluster vs rest of the stratum.
#   - Solid black border: disease structurally assigned to that cluster/factor
#     (hierarchical / EFA positive loading). Dashed border: EFA negative loading.
#
# NOTE ON CIRCULARITY: patients are assigned to a cluster according to the
#   diseases they have, so the bordered cells are enriched by construction.
#   The informative cells are the NON-bordered ones.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(htmlwidgets)
})

# 1. GLOBAL PARAMETERS --------------------------------------------------------

# Lower edge of every O/E bin; OE_ZMAX closes the last bin (">6").
OE_BREAKS <- c(0, 0.5, 0.9, 1.1, 2, 3, 4, 5, 6)
OE_ZMAX   <- 7
# One colour per bin (9 bins): blues, white, reds
OE_COLORS <- c(
  "#2166AC",  # 0   - 0.5
  "#C6DBEF",  # 0.5 - 0.9
  "#FFFFFF",  # 0.9 - 1.1
  "#FEE0D2",  # 1.1 - 2
  "#FCBBA1",  # 2   - 3
  "#FC9272",  # 3   - 4
  "#FB6A4A",  # 4   - 5
  "#DE2D26",  # 5   - 6
  "#99000D"   # > 6
)

MIN_PATIENTS <- MIN_PATIENTS_CLUSTER   # defined in 00_config.R: clusters below this size are striped
MIN_CASES    <- 10    # cells with fewer cases are marked with a dagger
FDR_ALPHA    <- 0.05
NEUTRAL_LOW  <- 0.9
NEUTRAL_HIGH <- 1.1

build_oe_colorscale <- function(breaks = c(OE_BREAKS, OE_ZMAX),
                                colors = OE_COLORS,
                                zmax   = OE_ZMAX) {
  stopifnot(length(colors) == length(breaks) - 1)
  scale <- list()
  for (i in seq_along(colors)) {
    scale[[length(scale) + 1]] <- c(breaks[i] / zmax,     colors[i])
    scale[[length(scale) + 1]] <- c(breaks[i + 1] / zmax, colors[i])
  }
  scale
}

Observed_Expected_colorscale <- build_oe_colorscale()

# 2. HELPERS ------------------------------------------------------------------

# Convert any 0/1, "0"/"1", NA or logical column to an integer 0/1 vector
to_binary <- function(v) {
  if (is.logical(v)) return(as.integer(v %in% TRUE))
  as.integer(!is.na(v) & as.character(v) != "0")
}

fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")

# Diagonal stripes (list of x / y with NA separators) covering one column
hatch_segments <- function(xc, ymin, ymax, s = 4, pitch = 2) {
  xs <- numeric(0)
  ys <- numeric(0)
  for (y0 in seq(ymin - s, ymax, by = pitch)) {
    t_lo <- max(0, (ymin - y0) / s)
    t_hi <- min(1, (ymax - y0) / s)
    if (t_hi <= t_lo) next
    xs <- c(xs, xc - 0.5 + t_lo, xc - 0.5 + t_hi, NA)
    ys <- c(ys, y0 + t_lo * s,   y0 + t_hi * s,   NA)
  }
  list(x = xs, y = ys)
}

# 2b. OUTPUT FILE NAMES (identical for all methods) ----------------------------
# heatmaps_dir : PATH_<METHOD>_HEATMAPS from 00_config.R
# method       : "hierarchical_full", "hierarchical_stratified", "LCA" or "EFA"
# Returns the html and csv paths:
#   <heatmaps_dir>/<Males|Females>/k_<k>/heatmap_<method>_<sex>_<age>_<status>_k<k>.html
oe_heatmap_files <- function(heatmaps_dir, method, sexo, edad, status, k) {
  sex_dir <- if (toupper(sexo) %in% c("MALE", "MEN")) "Males" else "Females"
  tag <- paste(method, sexo, edad, status, paste0("k", k), sep = "_")
  dir <- file.path(heatmaps_dir, sex_dir, paste0("k_", k))
  list(
    html = file.path(dir, paste0("heatmap_", tag, ".html")),
    csv  = file.path(dir, paste0("heatmap_table_", tag, ".csv"))
  )
}

# 3. STATISTICS ---------------------------------------------------------------

# cohort_stratum : data frame with ALL patients of the stratum. Needs the
#                  column `cluster` (character; NA = not assigned) and one 0/1
#                  column per disease.
# dis            : data frame (enfermedad_col, enfermedad), one row per disease.
# clusters       : character vector with the expected cluster labels.
compute_oe_table <- function(cohort_stratum, dis, clusters) {
  
  n_stratum <- nrow(cohort_stratum)
  cl <- cohort_stratum$cluster
  
  X <- sapply(dis$enfermedad_col, function(cn) to_binary(cohort_stratum[[cn]]))
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  colnames(X) <- dis$enfermedad
  
  expected   <- colMeans(X)                 # prevalence in the whole stratum
  total_case <- colSums(X)
  n_cl       <- vapply(clusters, function(k) sum(cl == k, na.rm = TRUE), numeric(1))
  
  grid <- expand.grid(
    enfermedad = dis$enfermedad,
    cluster    = clusters,
    stringsAsFactors = FALSE
  )
  
  res <- lapply(seq_len(nrow(grid)), function(i) {
    d <- grid$enfermedad[i]
    k <- grid$cluster[i]
    n_k <- n_cl[[k]]
    if (n_k == 0) {
      return(data.frame(n_patients = 0, n_cases = NA_real_,
                        observed_prev = NA_real_, p_value = NA_real_))
    }
    cases_k <- sum(X[which(cl == k), d])
    rest_n  <- n_stratum - n_k
    rest_c  <- total_case[[d]] - cases_k
    p <- if (rest_n > 0) {
      fisher.test(matrix(c(cases_k, n_k - cases_k, rest_c, rest_n - rest_c),
                         nrow = 2))$p.value
    } else NA_real_
    data.frame(n_patients = n_k, n_cases = cases_k,
               observed_prev = cases_k / n_k, p_value = p)
  })
  
  tab <- bind_cols(grid, bind_rows(res)) %>%
    mutate(
      expected_prev = unname(expected[enfermedad]),
      value = if_else(!is.na(expected_prev) & expected_prev > 0,
                      observed_prev / expected_prev, NA_real_)
    ) %>%
    filter(expected_prev > 0)          # diseases absent from the stratum
  
  # FDR (Benjamini-Hochberg) over all the cells of this heatmap
  tab$q_value <- NA_real_
  ok <- !is.na(tab$p_value)
  tab$q_value[ok] <- p.adjust(tab$p_value[ok], method = "BH")
  
  tab %>%
    mutate(
      small_cluster = n_patients < MIN_PATIENTS,
      low_cases     = !is.na(n_cases) & n_cases < MIN_CASES,
      colored       = !is.na(value) & (value < NEUTRAL_LOW | value > NEUTRAL_HIGH),
      not_signif    = colored & !is.na(q_value) & q_value >= FDR_ALPHA
    )
}

# 4. HEATMAP ------------------------------------------------------------------

# cohort_stratum : see compute_oe_table()
# disease_map    : data frame with columns
#                    enfermedad_col, enfermedad,
#                    cluster_enfermedad (chr; NA / absent = no pattern, e.g. LCA),
#                    sign (optional: "pos" solid border / "neg" dashed border),
#                    loading (optional, EFA only; shown in the tooltip).
#                  A disease may appear in several rows (EFA cross-loadings).
# clusters       : character vector of cluster labels shown on the x axis
# title_core     : e.g. "Observed/Expected Ratio (k=5) - FEMALE | Age: 45to65 | Status: CTL"
# out_html       : path of the self-contained HTML file
# out_csv        : optional path; the full table (O/E, p, q, N, cases) is saved
make_oe_heatmap <- function(cohort_stratum, disease_map, clusters,
                            title_core, out_html, out_csv = NULL) {
  
  if (!"cluster_enfermedad" %in% names(disease_map)) disease_map$cluster_enfermedad <- NA_character_
  if (!"sign"    %in% names(disease_map)) disease_map$sign    <- "pos"
  if (!"loading" %in% names(disease_map)) disease_map$loading <- NA_real_
  disease_map$cluster_enfermedad <- as.character(disease_map$cluster_enfermedad)
  
  dis <- disease_map %>%
    distinct(enfermedad_col, enfermedad) %>%
    filter(enfermedad_col %in% colnames(cohort_stratum))
  if (nrow(dis) == 0 || nrow(cohort_stratum) == 0) {
    message("Nothing to plot for: ", title_core)
    return(invisible(NULL))
  }
  
  clusters   <- as.character(clusters)
  n_stratum  <- nrow(cohort_stratum)
  n_assigned <- sum(!is.na(cohort_stratum$cluster))
  
  tab <- compute_oe_table(cohort_stratum, dis, clusters)
  if (!any(is.finite(tab$value))) {
    message("No valid ratios for: ", title_core)
    return(invisible(NULL))
  }
  
  sorted_diseases <- sort(unique(tab$enfermedad))
  n_dis <- length(sorted_diseases)
  
  cluster_n <- tab %>% distinct(cluster, n_patients)
  cluster_n <- cluster_n[match(clusters, cluster_n$cluster), ]
  
  tab <- tab %>%
    mutate(
      cluster_num = match(cluster, clusters),
      disease_num = match(enfermedad, sorted_diseases)
    )
  
  # Disease -> pattern text for the tooltip
  pattern_txt <- disease_map %>%
    filter(!is.na(cluster_enfermedad)) %>%
    group_by(enfermedad) %>%
    summarise(
      pattern_text = paste0(
        cluster_enfermedad,
        ifelse(is.na(loading), "", paste0(" (", sprintf("%.2f", loading), ")")),
        collapse = ", "
      ),
      .groups = "drop"
    )
  tab <- tab %>%
    left_join(pattern_txt, by = "enfermedad") %>%
    mutate(
      pattern_text = if_else(is.na(pattern_text), "none", pattern_text),
      ratio_text   = if_else(is.na(value), "NA", sprintf("%.2f", value)),
      q_text       = case_when(
        is.na(q_value) ~ "NA",
        q_value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", q_value)
      ),
      hover = paste0(
        "<b>", enfermedad, "</b><br><br>",
        "<b>Observed/Expected ratio:</b> ", ratio_text,
        "<br><b>Observed prevalence:</b> ", round(observed_prev, 4),
        "<br><b>Expected prevalence (stratum):</b> ", round(expected_prev, 4),
        "<br><b>Cases in cluster:</b> ", n_cases,
        "<br><b>N patients in cluster:</b> ", n_patients,
        "<br><b>FDR q (cluster vs rest):</b> ", q_text,
        "<br><b>Disease pattern(s):</b> ", pattern_text,
        if_else(small_cluster, paste0("<br><i>Cluster < ", MIN_PATIENTS, " patients</i>"), ""),
        if_else(low_cases,     paste0("<br><i>< ", MIN_CASES, " cases in this cell</i>"), "")
      )
    )
  
  # Structural borders (solid = positive, dashed = negative loading)
  borders <- disease_map %>%
    filter(!is.na(cluster_enfermedad)) %>%
    mutate(
      disease_num = match(enfermedad, sorted_diseases),
      cluster_num = match(cluster_enfermedad, clusters)
    ) %>%
    filter(!is.na(disease_num), !is.na(cluster_num)) %>%
    distinct(disease_num, cluster_num, sign)
  
  shapes <- lapply(seq_len(nrow(borders)), function(i) {
    list(
      type = "rect", xref = "x", yref = "y",
      x0 = borders$cluster_num[i] - 0.5, x1 = borders$cluster_num[i] + 0.5,
      y0 = borders$disease_num[i] - 0.46, y1 = borders$disease_num[i] + 0.46,
      line = list(color = "#000000", width = 3,
                  dash = if (identical(borders$sign[i], "neg")) "dash" else "solid"),
      fillcolor = "rgba(0,0,0,0)", layer = "above"
    )
  })
  
  # Cell marks
  marks <- tab %>%
    mutate(mark = case_when(
      low_cases  ~ "\u2020",
      not_signif ~ "ns",
      TRUE       ~ NA_character_
    )) %>%
    filter(!is.na(mark), !small_cluster)
  annotations <- lapply(seq_len(nrow(marks)), function(i) {
    list(x = marks$cluster_num[i], y = marks$disease_num[i],
         xref = "x", yref = "y", text = marks$mark[i], showarrow = FALSE,
         font = list(size = 9, color = "#444444"))
  })
  
  # Caption below the plot
  caption <- paste0(
    "Striped column: cluster < ", MIN_PATIENTS, " patients | \u2020 < ", MIN_CASES,
    " cases in the cell | ns: not significant (FDR q \u2265 ", FDR_ALPHA,
    ", cluster vs rest) | solid border: disease of the pattern",
    " | dashed border: negative EFA loading | grey: not estimable"
  )
  annotations[[length(annotations) + 1]] <- list(
    x = 0, y = 0, xref = "paper", yref = "paper", xanchor = "left",
    yanchor = "top", yshift = -62, showarrow = FALSE, align = "left",
    text = caption, font = list(size = 9, color = "#555555")
  )
  
  # x-axis labels with the size of each cluster
  x_labels <- paste0(
    "Cl. ", clusters, "<br>n=", fmt_n(cluster_n$n_patients),
    " (", sprintf("%.1f", 100 * cluster_n$n_patients / n_stratum), "%)"
  )
  
  title_full <- paste0(
    title_core, " | N = ", fmt_n(n_stratum),
    " (assigned: ", fmt_n(n_assigned), ")"
  )
  
  p <- plot_ly() %>%
    add_trace(
      data = tab,
      x = ~cluster_num, y = ~disease_num, z = ~value,
      type = "heatmap",
      zmin = 0, zmax = OE_ZMAX,
      colorscale = Observed_Expected_colorscale,
      xgap = 1, ygap = 1,
      hoverongaps = FALSE,
      text = ~hover,
      hovertemplate = paste0("%{text}<br><br><b>Patient cluster:</b> Cl. %{x}<extra></extra>"),
      colorbar = list(
        title = list(text = "Observed/Expected ratio", side = "top"),
        tickmode = "array",
        tickvals = c(OE_BREAKS, 6.5),
        ticktext = c("0", "0.5", "0.9", "1.1", "2", "3", "4", "5", "6", ">6"),
        len = 0.8, thickness = 20, ticklabelposition = "outside"
      )
    )
  
  # Diagonal stripes for small clusters (one scatter trace, NA-separated)
  small_cl <- which(cluster_n$n_patients < MIN_PATIENTS)
  if (length(small_cl) > 0) {
    hx <- numeric(0); hy <- numeric(0)
    for (xc in small_cl) {
      h <- hatch_segments(xc, 0.5, n_dis + 0.5)
      hx <- c(hx, h$x); hy <- c(hy, h$y)
    }
    p <- p %>% add_trace(
      x = hx, y = hy, type = "scatter", mode = "lines",
      line = list(color = "rgba(70,70,70,0.65)", width = 1.2),
      hoverinfo = "skip", showlegend = FALSE, inherit = FALSE
    )
  }
  
  p <- p %>% layout(
    title = list(text = title_full, font = list(size = 11)),
    xaxis = list(title = "", tickmode = "array", tickvals = seq_along(clusters),
                 ticktext = x_labels, range = c(0.5, length(clusters) + 0.5),
                 fixedrange = TRUE),
    yaxis = list(title = "", tickmode = "array", tickvals = seq_along(sorted_diseases),
                 ticktext = sorted_diseases, autorange = "reversed",
                 range = c(n_dis + 0.5, 0.5), tickfont = list(size = 8),
                 fixedrange = TRUE),
    plot_bgcolor = "#E5E5E5",
    shapes = shapes,
    annotations = annotations,
    margin = list(l = 240, b = 110, r = 100, t = 60),
    showlegend = FALSE
  )
  
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  saveWidget(p, out_html, selfcontained = TRUE)
  
  if (!is.null(out_csv)) {
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    tab %>%
      select(enfermedad, cluster, n_patients, n_cases, observed_prev,
             expected_prev, oe_ratio = value, p_value, q_value,
             small_cluster, low_cases, not_signif, pattern_text) %>%
      write.csv(out_csv, row.names = FALSE)
  }
  
  message(
    "Generated: ", basename(out_html),
    " | N stratum = ", n_stratum,
    " | assigned = ", n_assigned,
    " | clusters < ", MIN_PATIENTS, " pts: ", length(small_cl),
    " | max O/E = ", round(max(tab$value, na.rm = TRUE), 2)
  )
  invisible(list(table = tab, widget = p))
}

# 5. UNIFIED PATIENT ASSIGNMENT FOR PATTERN-BASED METHODS ---------------------
# Used for hierarchical clustering and EFA (disease-based methods). LCA uses the
# modal posterior class (already a membership probability).
#
# For every patient and pattern j:
#   count_j = number of the patient's diseases that belong to pattern j
#   size_j  = number of diseases of pattern j
# Two assignment rules (see ASSIGNMENT_RULES in 00_config.R):
#   "count"      score_j = count_j            -> pattern where the patient has
#                                                MORE diseases (favours big patterns)
#   "proportion" score_j = count_j / size_j   -> pattern whose diseases the patient
#                                                has in the HIGHEST PERCENTAGE
# The patient goes to the pattern with the highest score; a tie or no disease of
# any pattern leaves the patient unassigned (cluster = NA). Unassigned patients
# stay in the stratum for the O/E calculation.
#
# patterns: data frame (enfermedad_col, pattern) ; a disease can be in several
#           patterns (EFA cross-loadings), it then counts for each of them.
assign_patients_by_membership <- function(cohort, patterns, rule = c("count", "proportion")) {
  rule <- match.arg(rule)
  pats <- sort(unique(as.character(patterns$pattern)))
  size <- vapply(pats, function(p) {
    length(intersect(unique(patterns$enfermedad_col[as.character(patterns$pattern) == p]),
                     colnames(cohort)))
  }, numeric(1))
  counts <- sapply(pats, function(p) {
    cols <- unique(patterns$enfermedad_col[as.character(patterns$pattern) == p])
    cols <- intersect(cols, colnames(cohort))
    if (length(cols) == 0) return(rep(0L, nrow(cohort)))
    m <- sapply(cols, function(cn) to_binary(cohort[[cn]]))
    if (is.null(dim(m))) m <- matrix(m, ncol = 1)
    rowSums(m)
  })
  if (is.null(dim(counts))) counts <- matrix(counts, nrow = nrow(cohort))
  colnames(counts) <- pats
  score <- if (rule == "count") counts else sweep(counts, 2, pmax(size, 1), "/")
  mx    <- apply(score, 1, max)
  n_max <- apply(score, 1, function(r) sum(abs(r - max(r)) < 1e-12))
  cluster <- pats[apply(score, 1, which.max)]
  cluster[mx == 0 | n_max > 1] <- NA_character_
  tibble(patient_id = cohort$patient_id, cluster = cluster, max_score = mx)
}

# 6. COMPARISON OF THE TWO ASSIGNMENT RULES -------------------------------------
# results: list returned by pmap over (stratum, k, rule); every element is NULL or
#   list(sexo_val, edad_val, status_val, k, rule, table, assignment), where
#   table = make_oe_heatmap()$table and assignment = data frame(patient_id, cluster)
#   of ALL the patients of the stratum (cluster = NA when not assigned).
# Output: one row per stratum x k with
#   - how many patients are assigned with each rule,
#   - % of patients assigned with both rules that change cluster,
#   - agreement of the O/E values (Spearman, median |log2 change|) and % of cells
#     that change colour band.
compare_assignment_rules <- function(results, out_csv = NULL) {
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(invisible(NULL))
  key <- vapply(results, function(r) paste(r$sexo_val, r$edad_val, r$status_val, r$k, sep = "|"), "")
  rl  <- vapply(results, function(r) r$rule, "")
  out <- list()
  for (kk in unique(key)) {
    ia <- which(key == kk & rl == "count")
    ib <- which(key == kk & rl == "proportion")
    if (length(ia) != 1 || length(ib) != 1) next
    A <- results[[ia]]; B <- results[[ib]]
    asg   <- merge(A$assignment, B$assignment, by = "patient_id", suffixes = c("_a", "_b"))
    both  <- !is.na(asg$cluster_a) & !is.na(asg$cluster_b)
    cells <- merge(A$table[, c("enfermedad", "cluster", "value")],
                   B$table[, c("enfermedad", "cluster", "value")],
                   by = c("enfermedad", "cluster"), suffixes = c("_a", "_b"))
    ok <- is.finite(cells$value_a) & is.finite(cells$value_b)
    out[[kk]] <- data.frame(
      sexo_val = A$sexo_val, edad_val = A$edad_val, status_val = A$status_val, k = A$k,
      n_stratum = nrow(asg),
      n_assigned_count      = sum(!is.na(asg$cluster_a)),
      n_assigned_proportion = sum(!is.na(asg$cluster_b)),
      pct_assigned_count      = round(100 * mean(!is.na(asg$cluster_a)), 1),
      pct_assigned_proportion = round(100 * mean(!is.na(asg$cluster_b)), 1),
      pct_patients_changing_cluster =
        if (any(both)) round(100 * mean(asg$cluster_a[both] != asg$cluster_b[both]), 1) else NA_real_,
      n_cells_compared = sum(ok),
      oe_spearman = if (sum(ok) > 2) suppressWarnings(round(cor(cells$value_a[ok], cells$value_b[ok], method = "spearman"), 3)) else NA_real_,
      oe_median_abs_log2_change =
        if (any(ok)) round(median(abs(log2((cells$value_b[ok] + 0.01) / (cells$value_a[ok] + 0.01)))), 3) else NA_real_,
      pct_cells_changing_colour_band =
        if (any(ok)) round(100 * mean(findInterval(cells$value_a[ok], OE_BREAKS) != findInterval(cells$value_b[ok], OE_BREAKS)), 1) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, out)
  if (!is.null(out_csv) && !is.null(res)) {
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    write.csv(res, out_csv, row.names = FALSE)
    message("Assignment-rule comparison saved: ", out_csv)
  }
  invisible(res)
}
