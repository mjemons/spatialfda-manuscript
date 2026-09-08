### coded with the help of chatGPT ###
library("spatialFDA")
library("SpatialExperiment")
library("dplyr")

sce <- readRDS(snakemake@input[["rds"]])
spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    spatialCoords = as.matrix( 
      colData(sce)[, c("Center_X", "Center_Y")])
)
rm(sce)
gc()

colData(spe)$Stage_coarse <- dplyr::case_when(
  colData(spe)$Stage %in% c(1, 2) ~ "I",
  colData(spe)$Stage %in% c(3, 4) ~ "II",
  colData(spe)$Stage %in% c(5, 6) ~ "III",
  colData(spe)$Stage == 7 ~ "IV",
  TRUE ~ NA_character_
)

colData(spe)[["Stage_coarse"]] <- factor(colData(spe)[["Stage_coarse"]])
#relevel to have non-diabetic as the reference category
colData(spe)[["Stage_coarse"]] <- relevel(colData(spe)[["Stage_coarse"]],
"I")

colData(spe)$imageId <- paste0(colData(spe)$ImageNumber,"_", colData(spe)$Patient_ID)


#identify the rMax
p <- rMaxHeuristic(spe,
subsetby = "imageId", marks = "cell_category"
)

#run the spatial statistics inference
resG <- crossSpatialInference(
    spe, 
    selection = NULL,
    fun = "Gcross", 
    marks = "cell_category",
    rSeq = seq(0, 40, length.out = 20), 
    correction = "rs",
    sample_id = "Patient_ID",
    family = gaussian(link = "log"),
    algorithm = "bam",
    image_id = "imageId", 
    condition = "Stage_coarse",
    bs.int = list(bs = "ps", k = 15, m = c(2, 1)),
    bs.yindex = list(bs = "ps", k = 4, m = c(2, 1)),
    gc.level = 1,
    ncores = 10
)

names(resG)

set.seed(123)  # reproducible random selection

# Randomly select one image per patient
selected_images <- as.data.frame(colData(spe)) %>%
  distinct(Patient_ID, imageId) %>%
  group_by(Patient_ID) %>%
  slice_sample(n = 1) %>%
  ungroup()

selected_images

keep <- paste(colData(spe)$Patient_ID, colData(spe)$imageId) %in%
        paste(selected_images$Patient_ID, selected_images$imageId)

spe_one_image <- spe[, keep]

#run the spatial statistics inference
resGFM <- crossSpatialInference(
    spe_one_image, 
    selection = NULL,
    fun = "Gcross", 
    marks = "cell_category",
    rSeq = seq(0, 40, length.out = 20), 
    correction = "rs",
    sample_id = NULL,
    family = gaussian(link = "log"),
    algorithm = "bam",
    image_id = "imageId", 
    condition = "Stage_coarse",
    bs.int = list(bs = "ps", k = 15, m = c(2, 1)),
    bs.yindex = list(bs = "ps", k = 4, m = c(2, 1)),
    ncores = 10
)

out <- list("G" = resG, "GFM" = resGFM)

saveRDS(out, snakemake@output[["rds"]])
