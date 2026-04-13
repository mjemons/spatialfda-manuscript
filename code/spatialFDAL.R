library("spatialFDA")
library("SpatialExperiment")
library("dplyr")
library("mgcv")

spe <- readRDS(snakemake@input[["rds"]])
comp <- snakemake@wildcards[["comp"]]

spe <- subset(spe, ,condition %in% c("ctrl", comp))

colData(spe)[["condition"]] <- factor(colData(spe)[["condition"]])
#relevel to have non-diabetic as the reference category
colData(spe)[["condition"]] <- relevel(colData(spe)[["condition"]],
"ctrl")

#rename image ID
colData(spe)[["image_id"]] <- colData(spe)[["ID"]]

#run the spatial statistics inference
res <- spatialInference(
    spe, 
    selection = 0, 
    fun = "Lest", 
    marks = "Labels",
    rSeq = seq(0, 100, by = 1), 
    correction = "iso",
    sample_id = "sample_id",
    family = gaussian(link = "log"),
    image_id = "image_id", 
    condition = "condition",
    ncores = 1
)

mdl <- res$mdl

out <- summary(mdl, re.test = FALSE)
out$prob <- unique(colData(spe)$prob)

saveRDS(out, snakemake@output[["rds"]])
