library("spatialFDA")
library("SpatialExperiment")
library("dplyr")
library("ggplot2")
library("patchwork")

### coded with the help of chatGPT ###

#chose one simulation instance for visualisation
spe <- readRDS(snakemake@input[["rds"]])

### adapted from the spatialFDA vignettes ### 
metricResG <- calcMetricPerFov(spe = spe, selection = "0",
                              subsetby = "ID", fun = "Gest", 
                              marks = "Labels",
                              rSeq = seq(0, 100, by = 1), 
                              by = c("condition", "sample_id",
                                     "ID"),
                              ncores = 1)
metricResG$image_id <- metricResG$ID
metricResG$ID <- metricResG$sample_id
# plot metrics
pG <- plotMetricPerFov(metricResG, correction = "rs", x = "r",
                 imageId = "image_id", ID = "sample_id", ncol = 5) +
  ggtitle("G function for the simulated celltype 0")

metricResL <- calcMetricPerFov(spe = spe, selection = "0",
                              subsetby = "ID", fun = "Lest", 
                              marks = "Labels",
                              rSeq = seq(0, 100, by = 1), 
                              by = c("condition", "sample_id",
                                     "ID"),
                              ncores = 1)

metricResL$image_id <- metricResL$ID
metricResL$ID <- metricResL$sample_id

# plot metrics
pL <- plotMetricPerFov(metricResL, correction = "iso", x = "r",
                 imageId = "image_id", ID = "ID", ncol = 5) +
  ggtitle("L function for the simulated celltype 0")


speSubctrl <- subset(spe, ,condition %in% c("ctrl"))

df <- data.frame(spatialCoords(speSubctrl), colData(speSubctrl))

IDLs <- df$ID |> unique()
#plot only the first 10 images for visualisation
df <- df |> subset(ID %in% IDLs[seq(10)])

pCtrl <- ggplot(df, aes(x = x, y = y, color = as.factor(Labels))) +
    geom_point(size = 1) +
    facet_wrap(~ID, scales = "free", ncol = 5) +
    theme(legend.title.size = 20, legend.text.size = 20) +
    xlab("x") +
    ylab("y") +
    labs(color = "cell category") +
    scale_color_manual(
        values = c("0" = "red"),
        na.value = "grey80"
    ) +
    coord_cartesian() +
    theme_light() + 
    ggtitle("Control simulation")

speSubcase <- subset(spe, ,condition %in% c("pert3"))

df <- data.frame(spatialCoords(speSubcase), colData(speSubcase))

IDLs <- df$ID |> unique()
#plot only the first 10 images for visualisation
df <- df |> subset(ID %in% IDLs[seq(10)])

pCase <- ggplot(df, aes(x = x, y = y, color = as.factor(Labels))) +
    geom_point(size = 1) +
    facet_wrap(~ID, scales = "free", ncol = 5) +
    theme(legend.title.size = 20, legend.text.size = 20) +
    xlab("x") +
    ylab("y") +
    labs(color = "cell category") +
    scale_color_manual(
        values = c("0" = "red"),
        na.value = "grey80"
    ) +
    coord_cartesian() +
    theme_light() + 
    ggtitle("Perturbation 3 simulation")

pTotal <- pCtrl/pCase

pCurves <- pG/pL

ggsave(snakemake@output[["pltSims"]], plot = pTotal, width = 15, height = 10)
ggsave(snakemake@output[["pltCurves"]], plot = pCurves, width = 10, height = 10)
