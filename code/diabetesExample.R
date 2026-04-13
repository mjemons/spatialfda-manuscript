library("SpatialExperiment")
library("dplyr")
library("ggplot2"); theme_set(theme_light())
library("patchwork")

#re-install spatialFDA -> needs to be fixed later
#BiocManager::install("spatialFDA")
remotes::install_github("mjemons/spatialFDA@e70d8210bc46c36616c0ec596ac272631a9c0ec4")

library("spatialFDA")
# load the IMC dataset described in Damond et al. 2019 as SpatialExperiment object
spe <- .loadExample(full = TRUE)

### Intensity Boxplot figure

#code adapted for intensity plot from Elizabeth Purdom
#rename image ID
df <- .speToDf(spe)
selection <- colData(spe)[["cell_type"]] %>% unique() 

ConditionIntensities <- function(df, marks, conditionType, imageId, selection){
  df <- df %>% filter(patient_stage == conditionType)

  dfLs <- base::split(df, df[[imageId]])
  
  intensityDfCellType <- lapply(selection, function(x){
    intensitiesDf <- lapply(dfLs, function(dfSub){
      pp <- .dfToppp(dfSub, marks = marks, continuous = FALSE, window = NULL)
      ppSub <- pp[pp$marks == x, drop = TRUE]
      spatstat.geom::marks(ppSub) <- factor(spatstat.geom::marks(ppSub),
                                                levels = unique(x))
      cellTypeIntensity <- data.frame(intensity = spatstat.geom::intensity(ppSub),
                                      patient_id = unique(dfSub$patient_id),
                                      patient_stage = unique(dfSub$patient_stage),
                                      row.names = NULL)
      return(cellTypeIntensity)
    }) %>% bind_rows()
    intensitiesDf$cellType <- x
    return(intensitiesDf)
  }) %>% bind_rows()
  return(intensityDfCellType)
}

intensityDf <- lapply(list("Non-diabetic", "Onset", "Long-duration"), function(x){
  ConditionIntensities(df = df,
                       marks = "cell_type",
                       conditionType = x,
                       imageId = "image_number",
                       selection = selection)
})%>% bind_rows()

intensityDf$patient_stage <- factor(intensityDf$patient_stage,
                                    levels = c("Non-diabetic", "Onset", "Long-duration"))
p <- ggplot(intensityDf, aes(x=patient_stage,group=patient_id,y=intensity,fill=patient_stage))+ geom_boxplot()+facet_wrap(~cellType,scales="free", ncol = 3)

p

ggsave(snakemake@output[["intensityBoxplot"]], plot = p, width = 10, height = 12)

### run the actual differential co-localisation analysis

colData(spe)[["patient_stage"]] <- factor(colData(spe)[["patient_stage"]])
#relevel to have non-diabetic as the reference category
colData(spe)[["patient_stage"]] <- relevel(colData(spe)[["patient_stage"]],
"Non-diabetic")

#run the spatial statistics inference
res <- crossSpatialInference(
    spe, 
    selection = NULL,
    fun = "Gcross", 
    marks = "cell_type",
    rSeq = seq(0, 100, length.out = 50), 
    correction = "rs",
    sample_id = "patient_id",
    family = gaussian(link = "log"),
    algorithm = "bam",
    image_id = "image_number", 
    condition = "patient_stage",
    ncores = 5
)

names(res)

saveRDS(res, snakemake@output[["rds"]])
