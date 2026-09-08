library("spatialFDA")
library("SpatialExperiment")
library("dplyr")
library("mgcv")

spe <- readRDS(snakemake@input[["rds"]])
comp <- snakemake@wildcards[["comp"]]

spe <- subset(spe, ,condition %in% c("ctrl", comp))

colData(spe)[["condition"]] <- factor(colData(spe)[["condition"]])
#relevel to have ctrl as the reference category
colData(spe)[["condition"]] <- relevel(colData(spe)[["condition"]],
"ctrl")

unique(spe$condition)

#rename image ID
colData(spe)[["image_id"]] <- colData(spe)[["ID"]]

#run the spatial statistics inference
resG <- spatialInference(
    spe, 
    selection = 0, 
    fun = "Gest", 
    marks = "Labels",
    rSeq = seq(0, 100, by = 1), 
    correction = "rs",
    sample_id = "sample_id",
    family = gaussian(link = "log"),
    image_id = "image_id", 
    condition = "condition",
    ncores = 1,
    intensityAdjustment = FALSE
)

mdlG <- resG$mdl

outG <- summary(mdlG, re.test = FALSE)
outG$prob <- unique(colData(spe)$prob)

### run with intensity adjustement ###
resGAdj <- spatialInference(
    spe, 
    selection = 0, 
    fun = "Gest", 
    marks = "Labels",
    rSeq = seq(0, 100, by = 1), 
    correction = "rs",
    sample_id = "sample_id",
    family = gaussian(link = "log"),
    image_id = "image_id", 
    condition = "condition",
    ncores = 1,
    intensityAdjustment = TRUE
)

mdlGAdj <- resGAdj$mdl

outGAdj <- summary(mdlGAdj, re.test = FALSE)
outGAdj$prob <- unique(colData(spe)$prob)

### run wo sandwich correction ###
resGNSW <- spatialInference(
    spe, 
    selection = 0, 
    fun = "Gest", 
    marks = "Labels",
    rSeq = seq(0, 100, by = 1), 
    correction = "rs",
    sample_id = "sample_id",
    family = gaussian(link = "log"),
    image_id = "image_id", 
    condition = "condition",
    ncores = 1,
    sandwich = "none",
    intensityAdjustment = FALSE

)

mdlGNSW <- resGNSW$mdl

outGNSW <- summary(mdlGNSW, re.test = FALSE)
outGNSW$prob <- unique(colData(spe)$prob)


saveRDS(outG, snakemake@output[["rdsG"]])
saveRDS(outGAdj, snakemake@output[["rdsGAdj"]])
saveRDS(outGNSW, snakemake@output[["rdsGNSW"]])
