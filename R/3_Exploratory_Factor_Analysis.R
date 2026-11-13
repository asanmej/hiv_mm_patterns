# Author: Alejandro Santos Mejías
# Date last update: 2025-06-23
# Input: Boolean diagnosis datasets
# Output: Exploratory Factor Analysis output 
# Motivational comment: Perceiving reality is embracing all possibilities.
# 
##############
# Acknoledge #
##############
# For concepts: https://m-clark.github.io/posts/2020-04-10-psych-explained/#introduction
# For workflow: https://www.promptcloud.com/blog/exploratory-factor-analysis-in-r/
#               https://rpubs.com/pjmurphy/758265
# Original free source of the psych package author: https://www.personality-project.org/r/book/


# Load the project configuration if the script is run on its own
# (Main.R has already done it). The environment is NOT cleared on purpose.
if (!exists("PATH_R")) source(here::here("R", "00_config.R"))

# assign_patients_by_membership(): same patient -> pattern rule as the heatmaps
source(file.path(PATH_R, "06_heatmap_common.R"))

######################################
##### Exploratory Factor Analysis ####
######################################

################################################################################

# Load the data
diag_bool_prev <- fread(PATH_DIAG_BOOL_PREVALENT_FINAL, encoding = "UTF-8")
bdu <- fread(PATH_BDU_FULL, encoding = "UTF-8")
# cohort <- fread("../Data/cohort_final.csv", encoding = "UTF-8", colClasses = "character")

bdu[, ageband := cut(edad, c(-Inf,44,65,Inf), labels = c("under 45", "45 to 65" , "over 65"))]

# Create parameters for filtering EFA:
stratum <- setDT(expand.grid(sexo = c("HOMBRE", "MUJER"), edad = c("under 45", "45 to 65" , "over 65")))
stratum <- rbind(stratum, list("Global", NA))

prevalenceProportionVec <- 2

# Viability of one EFA solution (see the rule above).
#   faResult   : model returned by psych::fa()
#   groupData  : data.frame of the patients of the group (patient_id + diseases)
# Returns a one-row data.frame with the minimum number of diseases per factor, the
# minimum number of patients per factor under each rule and the verdict.
efaCheckSolution <- function(faResult, groupData) {
  L <- round(faResult$loadings[seq_len(nrow(faResult$loadings)), , drop = FALSE], 3)
  k <- ncol(L)
  idx <- which(L >= EFA_LOADING_CUTOFF, arr.ind = TRUE)
  
  nDis <- tabulate(idx[, "col"], nbins = k)
  minPat <- setNames(rep(0, length(K_VIABILITY_RULES)), K_VIABILITY_RULES)
  
  if (nrow(idx) > 0) {
    patterns <- data.frame(
      enfermedad_col = rownames(L)[idx[, "row"]],
      pattern        = as.character(idx[, "col"]),
      stringsAsFactors = FALSE
    )
    for (r in K_VIABILITY_RULES) {
      a <- assign_patients_by_membership(groupData, patterns, rule = r)
      minPat[[r]] <- min(tabulate(match(a$cluster, as.character(seq_len(k))), nbins = k))
    }
  }
  
  out <- data.frame(
    k = k,
    min_diseases_per_factor = min(nDis),
    stringsAsFactors = FALSE
  )
  for (r in K_VIABILITY_RULES) out[[paste0("min_patients_", r)]] <- minPat[[r]]
  out$viable <- out$min_diseases_per_factor >= MIN_DISEASES_PATTERN &&
    all(minPat >= MIN_PATIENTS_CLUSTER)
  out
}

# Container for the selection table (saved at the end of the script)
efaSelection <- list()
# Solutions exported per stratum and group: 3 at most (like LCA and the
# hierarchical clustering).
#
# Starting point: half of the factors suggested by the parallel analysis
# (rounded down), one below and one above, i.e. the solutions
# floor(n/2) - 1, floor(n/2), floor(n/2) + 1 (never below kFloor, 3 solutions
# kFloor, kFloor + 1, kFloor + 2 in that case).
#
# Cap: a solution is VIABLE only if all its patterns are real patterns, i.e.
#   - every factor has at least MIN_DISEASES_PATTERN diseases with a positive
#     loading >= EFA_LOADING_CUTOFF, and
#   - every factor reaches MIN_PATIENTS_CLUSTER patients under EVERY assignment
#     rule of K_VIABILITY_RULES (so no heatmap column of the solution is striped).
# The number of factors is scanned upwards from kFloor and stops at the first
# solution that is not viable (more factors only make patterns smaller). The
# solutions exported are the three highest viable ones that do not exceed the
# starting point above. The same rule applies to every stratum; strata with few
# patients or few diseases simply end up with fewer factors.
# (All thresholds are defined in 00_config.R.)
kFloor <- 3
# Minimum number of variables per factor (below this the model is not fitted)
minVarPerFactor <- 3
correlationKmo <- 0.5

# Output paths: PATH_EFA_TABLES / PATH_EFA_MODELS / PATH_EFA_FIGURES (00_config.R)
for (strata in 1:nrow(stratum)) {
  
  
  if(stratum[strata, sexo] == "Global"){
    diagSub <- copy(diag_bool_prev)
    strataLabel <- "Global"
    
  }else{
    diagSub <- copy(diag_bool_prev[patient_id %in% bdu[sexo == stratum[strata, sexo] & ageband == stratum[strata, edad], patient_id]])
    strataLabel <- stratum[strata, paste0(sexo, "_", edad)]
  }
  
  
  for (prevalenceProportion in prevalenceProportionVec) {
    
    # Create output folders:
    pathEfaTables <- file.path(PATH_EFA_TABLES, strataLabel)
    pathEfaModels <- file.path(PATH_EFA_MODELS, strataLabel)
    pathEfaPic    <- file.path(PATH_EFA_FIGURES, strataLabel, paste0("Prev_", prevalenceProportion))
    for (d in c(pathEfaTables, pathEfaModels, pathEfaPic)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
    # Filter out none-EFA-appropiate variables, variables with low number of
    # positive cases or with really low correlations with the others.
    # Remove variables with less than 1 positive case:
    cols <- colnames(diagSub)[-c(1:3)]
    
    for (i in 0:1) {
      t <- cols[diagSub[HIV_infection == i, colSums(.SD), .SDcols = cols]>2]
      assign(paste0("cols",ifelse(i == 0, "Ctl", "Hiv")), t)
      
    }
    
    
    # Test idoneity of data for EFA Kaiser-Meyer-Olkin factor adequacy is a measure
    # of the proportion of variance among variables that might be common variance. A
    # higher proportion indicates a higher KMO-value, which means the data is more
    # suited for factor analysis. Rule of Thumbs:
    # * 0.00 to 0.49 unacceptable
    # * 0.50 to 0.59 miserable
    # * 0.60 to 0.69 mediocre
    # * 0.70 to 0.79 middling
    # * 0.80 to 0.89 meritorious
    # * 0.90 to 1.00 marvelous
    for (i in 0:1) {
      
      kmo <- KMO(r = diagSub[HIV_infection == i, get(ifelse(i == 0, "colsCtl", "colsHiv")), with = FALSE])  
      if(i == 0){
        t <- colsCtl[kmo$MSAi > correlationKmo]
      }else{
        t <- colsHiv[kmo$MSAi > correlationKmo]
      }
      
      assign(paste0("cols",ifelse(i == 0, "Ctl", "Hiv")), t)
      
    }
    
    # Remove diseases with low incidence/prevalence. In theory for trustful result
    # at least 10 obs per variable is desired, variables with less than 5 it is
    # recommendable to remove. For prevalence, I will follow previous works,
    # selecting prevalences > 1%, if one diseases index increase prevalence to up to
    # 5%. For incidences, I will select diseases with at least 20 observations.
    for (i in 0:1) {
      
      cols <- ifelse(i == 0, "colsCtl", "colsHiv")
      nObsDis <- diagSub[HIV_infection == i, .(disease = get(cols), N = colSums(.SD)/.N*100), .SDcols = get(cols)]
      t <- nObsDis[N >= prevalenceProportion, c(disease)] 
      assign(paste0("cols",ifelse(i == 0, "Ctl", "Hiv")), t)
    }
    
    # Let's start the analysis:
    for (i in 0:1) {
      t <- ifelse(i == 0, "colsCtl", "colsHiv")
      cols <- get(t)
      
      # There are two possible workflows to elucidate MM patterns for a target disease
      # population:
      #  - 2 Factor analysis, one for the population with and without the target disease.
      #  - 1 Factor analysis for the whole cohort, checking in which factor the loading is greater for the target disease.
      
      # The first approach is followed in this analysis.
      
      # @param fm - factoring method, maximum likelihood chosen for consistency with previous work. minres the default one.
      # @param fa - fa = show eigen values for principal axis factor analysis.
      # @param cor - which correlation should be calculated, in our case tetrachoric (binary data).
      
      png(filename = paste0(pathEfaPic, "/Scree_plot_", gsub("cols", "", t), "_", prevalenceProportion, "_", strataLabel ,".png"), width = 700, height = 500)
      
      nFactors <- suppressMessages(suppressWarnings(fa.parallel(
          x = diagSub[HIV_infection == i, ..cols],
          fm = "ml",
          fa = "fa",
          cor = "tet", 
          main = paste0(gsub("cols", "", t), " Prev ", prevalenceProportion, "% ", strataLabel,":\nParallel Analysis Scree Plot"),
          plot = T
        )))
        
      dev.off()
      
      grpLabel <- gsub("cols", "", t)
      grpData  <- as.data.frame(diagSub[HIV_infection == i])
      
      # Starting point of the solutions: half of the factors suggested by the
      # parallel analysis (rounded down) + 1 is the highest solution considered
      nfactDetected <- nFactors$nfact
      if (is.na(nfactDetected) || nfactDetected <= 2) {
        warning("Optimal estimated factors less or equal 2. Check data!")
        kUpper <- kFloor + 2
      } else {
        kUpper <- max(floor(nfactDetected / 2) + 1, kFloor + 2)
      }
      
      # Scan upwards from kFloor; stop at the first solution that is not viable.
      # Oblique rotation (oblimin): factors are allowed to be correlated.
      efaFits <- list()
      efaRows <- list()
      for (l in kFloor:kUpper) {
        
        if (length(cols) < minVarPerFactor * l) {
          message("Stop scanning at ", l, " factors (", strataLabel, ", ", grpLabel, "): only ", length(cols), " variables")
          break
        }
        
        faResult <- tryCatch(
          suppressMessages(suppressWarnings(
            fa(
              diagSub[HIV_infection == i, ..cols],
              fm = "ml",
              nfactors = l,
              cor = "tet",
              rotate = "oblimin"
            )
          )),
          error = function(e) {
            message("fa() failed with ", l, " factors (", strataLabel, ", ", grpLabel, "): ", conditionMessage(e))
            NULL
          }
        )
        if (is.null(faResult)) break
        
        chk <- efaCheckSolution(faResult, grpData)
        efaFits[[as.character(l)]] <- faResult
        efaRows[[as.character(l)]] <- chk
        
        if (!chk$viable) break
      }
      
      efaScan <- do.call(rbind, efaRows)
      
      # Solutions to export: the three highest viable ones (never below kFloor)
      if (is.null(efaScan)) {
        kEfa <- integer(0)
      } else if (!any(efaScan$viable)) {
        warning("No viable EFA solution in ", strataLabel, " (", grpLabel,
                "): the smallest one is exported anyway, check it.")
        kEfa <- efaScan$k[1]
      } else {
        kTop <- max(efaScan$k[efaScan$viable])
        kEfa <- seq(max(kFloor, kTop - 2), kTop)
      }
      
      message("Parallel analysis suggests ", nfactDetected, " factors (", strataLabel, ", ",
              grpLabel, "); highest solution considered: ", kUpper,
              "; exporting solutions: ", paste(kEfa, collapse = ", "))
      
      # Record the scan of this stratum and group
      if (!is.null(efaScan)) {
        efaScan$stratum        <- strataLabel
        efaScan$group          <- grpLabel
        efaScan$n_patients     <- nrow(grpData)
        efaScan$n_diseases     <- length(cols)
        efaScan$nfact_parallel <- nfactDetected
        efaScan$k_highest_considered <- kUpper
        efaScan$selected       <- efaScan$k %in% kEfa
        efaSelection[[paste(strataLabel, grpLabel, sep = "_")]] <- efaScan
      }
      
      # Remove the solutions of previous runs so that only the current ones remain
      unlink(list.files(pathEfaTables, full.names = TRUE,
                        pattern = paste0("_", grpLabel, "_", prevalenceProportion, "%_[0-9]+factors\\.csv$")))
      unlink(list.files(pathEfaModels, full.names = TRUE,
                        pattern = paste0("^efa_models_", grpLabel, "_", prevalenceProportion, "%_[0-9]+factors\\.rds$")))
      
      for (l in kEfa) {
        
        faResult <- efaFits[[as.character(l)]]
        
        faLoadings <- as.data.frame(round(faResult$loadings[1:dim(faResult$loadings)[1], ], 3))
        setorderv(faLoadings, colnames(faLoadings)[-1], rep(-1, length(colnames(faLoadings)[-1])), na.last = T)
        
        fwrite(faLoadings, paste0(pathEfaTables, "/efa_loadings_", strataLabel, "_", grpLabel, "_", prevalenceProportion,"%_", l, "factors.csv"), encoding = "UTF-8", row.names = T)
        saveRDS(faResult, paste0(pathEfaModels, "/efa_models_", grpLabel, "_", prevalenceProportion,"%_", l, "factors.rds"))
        
      }
      
    }
    
    # Create tile plots of the relevant factor results:
    files <- list.files(path = pathEfaTables, pattern = paste0(prevalenceProportion,"%.*\\.csv$"))  
    
    for(resultFile in files){
      t <- fread(paste0(pathEfaTables, "/", resultFile))
      setorder(t, V1)
      t[, V1 := factor(x = V1, levels = rev(V1))] # Invert factor level for alphabetical 
      temp <- melt(data = t, id.vars = "V1")
      # temp[abs(value) < 0.3 , value := NA] # Comment if all disease loading needed to be seen in the images
      temp[, text := paste0("y: ", V1, "\n", "x: ", variable, "\n", "Value: ", value)]
      
      p <- ggplot(temp, aes(x = variable, y = V1, fill = value, text = text)) + 
        geom_tile() + theme_bw() + scale_fill_gradient2(low = "#e6194B" , mid = "white", high = "#4363d8", na.value = "white", breaks = c(1, 0.9, 0.6, 0.3, 0 , -0.3, -0.6, -0.9,-1)) +
        theme(
          axis.title = element_blank(),
          panel.grid.major = element_blank(), # Remove major grid lines
          panel.grid.minor = element_blank(), # Remove minor grid lines
          panel.background = element_blank(), # Remove background color
          axis.line = element_line(color = "black") # Add axis lines
        )
      
      pp <- ggplotly(p, tooltip = "text")
      fileName <- gsub("csv$", "html", resultFile)
      fileName <- gsub("efa_loadings_", "", fileName)
      fileName <- gsub("%", "prev", fileName)
      
      saveWidget(pp, file=paste0(pathEfaPic, "/", fileName)) 
      
    }
  }
}
# Solutions scanned and selected in every stratum and group
# (viable = every factor has >= MIN_DISEASES_PATTERN diseases and
#  >= MIN_PATIENTS_CLUSTER patients under all the rules of K_VIABILITY_RULES)
if (length(efaSelection) > 0) {
  fwrite(rbindlist(efaSelection, fill = TRUE),
         file.path(PATH_EFA_TABLES, "k_selection_EFA.csv"), encoding = "UTF-8")
}

# TODO reliability statistic such as alpha and omega to check results

