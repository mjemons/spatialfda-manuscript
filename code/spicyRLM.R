library("spicyR")
library("SpatialExperiment")

spe <- readRDS(snakemake@input[["rds"]])
comp <- snakemake@wildcards[["comp"]]

speSub <- subset(spe, ,condition %in% c("ctrl", comp))

#colData(speSub)$cellID <- as.factor(paste0("cell", "_", c(1:ncol(speSub))))
colData(speSub)$condition <- as.factor(colData(speSub)$condition)
colData(speSub)$cellType <- as.factor(colData(speSub)$Labels)
colData(speSub)$imageID <- as.factor(colData(speSub)$ID)
colData(speSub)$sample_id <- as.factor(colData(speSub)$sample_id)

spicyTestPair <- spicy(
  speSub,
  condition = "condition",
  imageID = "imageID",
  cellType = "cellType",
  window = "square",
  from = "0",
  to = "0"
)

out <- try(topPairs(spicyTestPair, n = 1), TRUE)
out$comparison <- comp

saveRDS(out, snakemake@output[["rds"]])




