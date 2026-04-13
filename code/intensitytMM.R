library("SpatialExperiment")
library("dplyr")
library("lmerTest")

spe <- readRDS(snakemake@input[["rds"]])
comp <- snakemake@wildcards[["comp"]]

#rename image ID
colData(spe)[["image_id"]] <- colData(spe)[["ID"]]
colData(spe)$Labels <- as.factor(colData(spe)$Labels)

speSub <- subset(spe, ,condition %in% c("ctrl", comp))

#code from https://github.com/mjemons/spatialFDA licensed under GPL-3 to mjemons
.speToDf <- function(spe) {
    df <- data.frame(
        x = SpatialExperiment::spatialCoords(spe)[, 1],
        y = SpatialExperiment::spatialCoords(spe)[, 2]
    )
    df <- cbind(df, colData(spe))
}

#code from https://github.com/mjemons/spatialFDA licensed under GPL-3 to mjemons
.dfToppp <- function(df, marks = NULL, continuous = FALSE, window = NULL) {
    #type checking
    stopifnot(is(df, "data.frame"))
    # this definition of the window is quite conservative
    # - can be set explicitly
    pp <- spatstat.geom::as.ppp(data.frame(x = df$x, y = df$y),
        W = spatstat.geom::owin(
            c(
                base::min(df$x) - 1,
                base::max(df$x) + 1
            ),
            c(
                base::min(df$y) - 1,
                base::max(df$y) + 1
            )
        )
    )
    # set the marks
    if (!continuous) {
        spatstat.geom::marks(pp) <- factor(df[[marks]])
    } else {
        spatstat.geom::marks(pp) <- base::subset(df, select =
                                                   names(df) %in% marks)
    }
    # if window exist, set is as new window and potentially exclude some points
    if (!is.null(window)) {
        pp <- spatstat.geom::as.ppp(spatstat.geom::superimpose(pp, W = window))
    }

    return(pp)
}

df <- .speToDf(speSub)

ConditionIntensities <- function(df, marks, conditionType, cellType, imageId){
  df <- df %>% filter(condition == conditionType)

  dfLs <- base::split(df, df[[imageId]])

  intensitiesDf <- lapply(dfLs, function(dfSub){
    pp <- .dfToppp(dfSub, marks = marks, continuous = FALSE, window = NULL)
    ppSub <- pp[pp$marks == cellType, drop = TRUE]
    spatstat.geom::marks(ppSub) <- factor(spatstat.geom::marks(ppSub),
                                              levels = unique(cellType))
    cellTypeIntensity <- data.frame(intensity=spatstat.geom::intensity(ppSub))
    cellTypeIntensity$condition <- unique(dfSub$condition)
    cellTypeIntensity$sample_id <- unique(dfSub$sample_id)
    return(cellTypeIntensity)
  }) %>% bind_rows()

  return(intensitiesDf)
}

intensitiesControl <- ConditionIntensities(df = df,
					   marks = "Labels",
                                           conditionType = "ctrl",
                                           cellType = "0",
                                           imageId = "image_id")

intensitiesPert <- ConditionIntensities(df = df,
					marks = "Labels",
                                        conditionType = comp,
                                        cellType = "0",
                                        imageId = "image_id")

dfIntensity <- rbind(intensitiesControl, intensitiesPert)

dfIntensity$condition <- factor(dfIntensity$condition, levels = c("ctrl", comp))

lmeMod <- lmerTest::lmer(intensity ~ condition + (1|sample_id), data = dfIntensity)
out <- lmerTest:::get_coefmat(lmeMod) |>
    as.data.frame()

out$comparison <- comp
saveRDS(out, snakemake@output[["rds"]])


