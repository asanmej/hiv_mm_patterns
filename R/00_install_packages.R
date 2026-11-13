# -----------------------------------------------------------------------------
# 00_install_packages.R: ONE-TIME package installation
#
# Run it ONCE per computer, in a FRESH R session, before Main.R:
#   - RStudio: Session > Restart R, then run this file (Ctrl+Shift+Enter).
#   - If RStudio still has 'rlang' loaded (the script stops and tells you),
#     run it from a terminal instead:   Rscript R/00_install_packages.R
#     or from the plain "R" (Rgui) program.
#
# Why a separate script: Windows cannot replace a package that is loaded
# (rlang, cli...), and newer packages need newer versions of those. Installing
# inside the analysis session therefore fails with "lazy loading failed".
#
# Keep the CRAN list in sync with the `pkg` vector of 00_config.R.
# -----------------------------------------------------------------------------

cran_packages <- c(
  "here", "rmarkdown", "knitr", "kableExtra",
  "tidyverse", "data.table",
  "factoextra", "proxy", "dendextend", "cluster", "psych",
  "survival",
  "plotly", "htmlwidgets",
  "writexl",
  "remotes"                         # lightweight installer for GitHub packages
)
github_packages <- c(
  MMLCA    = "ARCbiostat/MMLCA",
  ggsankey = "davidsjoberg/ggsankey"
)

# 1. The session must be clean ------------------------------------------------
locked <- intersect(c("rlang", "cli", "vctrs", "ggplot2", "dplyr"), loadedNamespaces())
if (length(locked) > 0) {
  stop(
    "These packages are already loaded in this session: ", paste(locked, collapse = ", "),
    ".\nRestart R (Session > Restart R) and run this script again, or run it with ",
    "Rscript from a terminal.", call. = FALSE
  )
}

# 2. CRAN packages: always Windows BINARIES -----------------------------------
# CRAN publishes source versions before the Windows binaries. Without this
# option R falls back to compiling the newer source (needs Rtools, and breaks
# when a dependency is outdated). The binary version is perfectly fine.
options(install.packages.compile.from.source = "never")
missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing) > 0) {
  message("Installing from CRAN: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
}

# 3. GitHub packages ----------------------------------------------------------
# upgrade = "never": do not touch packages that are already installed.
for (pk in names(github_packages)) {
  if (!requireNamespace(pk, quietly = TRUE)) {
    message("Installing from GitHub: ", github_packages[[pk]])
    tryCatch(
      remotes::install_github(github_packages[[pk]], upgrade = "never"),
      error = function(e) message("  FAILED: ", conditionMessage(e))
    )
  }
}

# 4. Final report ---------------------------------------------------------------
all_pkgs <- c(cran_packages, names(github_packages))
still_missing <- setdiff(all_pkgs, rownames(installed.packages()))
if (length(still_missing) == 0) {
  message("\nAll packages are installed. You can now run Main.R")
} else {
  message("\nStill missing: ", paste(still_missing, collapse = ", "))
  message("If a GitHub package failed while compiling, install Rtools 4.4 ",
          "(https://cran.r-project.org/bin/windows/Rtools/) and run this script again.")
}
