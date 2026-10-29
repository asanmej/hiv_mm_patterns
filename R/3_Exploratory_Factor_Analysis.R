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
# Common rule for the three methods (EFA, hierarchical, LCA): same range of
# solutions. The number suggested by the parallel analysis is still printed
# but it no longer limits the range of factors that are exported.
kRange <- 3:7
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
      
      # If we wanted simple solutions modified nMaxFactors to contain all solutions
      nMaxFactors <- nFactors$nfact
      
      if(nMaxFactors <= 2){
        warning("Optimal estimated factors less or equal 2. Check data!")
        nMaxFactors <- 3
      }else{
        nMaxFactors <- ceiling(nMaxFactors/2)+2
      }
      
      # @param rotate - method to rotate the data for factoring. Take into account:
      # Orthogonal rotations if we assume our factors are NOT correlated. Such as none, varimax, quartimax, bentlerT, equamax, varimin, geominT and bifactor
      # Oblique rotations if we assume our factors are correlated. Such as oblimin, geominQ, Promax or target.rot
      message("Parallel analysis suggests ", nFactors$nfact, " factors; exporting kRange = ", paste(kRange, collapse = ", "))
      for(l in kRange){
        
        if(length(cols) < minVarPerFactor * l){
          message("Skipping ", l, " factors (", strataLabel, ", ", gsub("cols", "", t), "): only ", length(cols), " variables")
          next
        }
        
        faResult <-
          suppressMessages(suppressWarnings(
            fa(
              diagSub[HIV_infection == i, ..cols],
              fm = "ml",
              nfactors = l,
              cor = "tet",
              rotate = "oblimin"
            )
          ))
        
        faLoadings <- as.data.frame(round(faResult$loadings[1:dim(faResult$loadings)[1], ], 3))
        setorderv(faLoadings, colnames(faLoadings)[-1], rep(-1, length(colnames(faLoadings)[-1])), na.last = T)
        
        fwrite(faLoadings, paste0(pathEfaTables, "/efa_loadings_", strataLabel, "_", gsub("cols", "", t), "_", prevalenceProportion,"%_", l, "factors.csv"), encoding = "UTF-8", row.names = T)
        saveRDS(faResult, paste0(pathEfaModels, "/efa_models_", gsub("cols", "", t), "_", prevalenceProportion,"%_", l, "factors.rds"))
        
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
# TODO reliability statistic such as alpha and omega to check results

