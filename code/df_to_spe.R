suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(SpatialExperiment)
})

# df to spe to be able to interface bioc methods 

dfToSpe <- function(df){
  # create a SpatialExperiment object
  spe <- SpatialExperiment(
    colData = df,
    spatialCoordsNames = c("x", "y"),
    sample_id = df$sample_id
  )
  return(spe)
}

df <- readRDS(snakemake@input[["rds"]])

spe <- dfToSpe(df)

saveRDS(spe, snakemake@output[["rds"]])