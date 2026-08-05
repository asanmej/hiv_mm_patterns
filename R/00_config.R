# 00_config.R: 

# 1. List of required packages 
pkg <- c(
  
  # 1. Hierarchical clustering analysis
  ## For data manipulation and visualization
  "tidyverse",
  ## For efficient Jaccard distance calculation
  "proxy", 
  ## To customize and view dendrograms
  "dendextend",
  ## To evaluate the average silhouette width 
  "cluster"
  
  # 2. Data visualization
 
)

# 2. Automatically load and install missing packages
lapply(pkg, function (x){
  if(!require(x, character.only = T)){
    install.packages(x , character.only = T)
    }
  })

# 3. Clean the environment by deleting the temporary variable
rm(pkg)