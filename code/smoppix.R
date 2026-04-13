library("SpatialExperiment")
library('smoppix')

spe <- readRDS(snakemake@input[["rds"]])
comp <- snakemake@wildcards[["comp"]]

speSub <- subset(spe, ,condition %in% c("ctrl", comp))
colData(speSub)$cellTypes <- as.character(colData(speSub)$Labels)
colData(speSub)$imageID <- as.character(colData(speSub)$ID)
colData(speSub)$sample_id <- as.character(colData(speSub)$sample_id)
colData(speSub)$condition <- as.character(colData(speSub)$condition)

df <- colData(speSub) |> as.data.frame() |> cbind(spatialCoords(speSub))

hypDf <- buildHyperFrame(df,
                          coordVars = c("x", "y"),
                          imageVars = c("condition", "sample_id", "ID"),
                         featureName = "cellTypes"
)

nnObj <- estPis(hypDf,
                pis = c("nn"), null = "background", verbose = FALSE,
                features = c("0", "1", "2", "3")
)

nnObj <- addWeightFunction(nnObj, lowestLevelVar = "ID",
                              pi = "nn")

dfUniNN <- buildDataFrame(nnObj, gene = "0", pi = "nn")

lmeMod <- lmerTest::lmer(pi - 0.5 ~ condition + (1 | sample_id),
                         data = dfUniNN, na.action = na.omit,
                         weights = weight, contrasts = list("condition" = "contr.sum")
)

out <- lmerTest:::get_coefmat(lmeMod) |>
    as.data.frame()
out$comparison <- comp

saveRDS(out, snakemake@output[["rds"]])
