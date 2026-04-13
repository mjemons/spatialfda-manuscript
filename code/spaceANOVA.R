library("spatstat.geom")
library("spatstat")
library("tidyr")
library("dplyr")
library("fda.usc")
library("ggplot2")
library("cowplot")
library("gridExtra")

set.seed(123)

spe <- readRDS(snakemake@input[["rds"]])
comp <- snakemake@wildcards[["comp"]]

speSub <- subset(spe, ,condition %in% c("ctrl", comp))

colData(speSub)$Group <- as.factor(colData(speSub)$condition)
colData(speSub)$cellType <- as.factor(colData(speSub)$Labels)
colData(speSub)$imageID <- as.factor(colData(speSub)$ID)
colData(speSub)$ID <- as.factor(colData(speSub)$sample_id)
colData(speSub)$x <- spatialCoords(speSub)[,1]
colData(speSub)$y <- spatialCoords(speSub)[,2]
data <- colData(speSub) %>% as.data.frame() %>% dplyr::select(Group, cellType, imageID, ID, x, y)

#coded according to the github repository of spaceANOVA

Final_result = SpaceANOVA::All_in_one(data = data, fixed_r = seq(0, 100, by = 1), Summary_function = "g",  Hard_ths = 10, homogeneous = TRUE, interaction_adjustment = TRUE, perm = TRUE, nPerm = 20, cores = 2)

out <- SpaceANOVA::p_extract(Final_result)
print(out)
saveRDS(out[[1]], snakemake@output[["rdsUni"]])
saveRDS(out[[2]], snakemake@output[["rdsMulti"]])
