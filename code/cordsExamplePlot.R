### coded with help of ChatGPT ###
library("spatialFDA")
library("dplyr")
library("ggplot2")
library("refund")
library("patchwork")
library("SpatialExperiment")

res <- readRDS(snakemake@input[["rds"]])

resG <- res$G
resGFM <- res$GFM

## adapted from the `spatialFDA` vignette
pG <- plotCrossHeatmap(resG, QCThreshold = 0, QCMetric = "medianMinIntensity") +
  theme(legend.position = "bottom") + guides(shape = "none") +
    theme(legend.position = "bottom") + theme(
      axis.title   = element_text(size = 18),
      axis.text    = element_text(size = 14),
      strip.text   = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      plot.title   = element_text(size = 20)
    )

pGFM <- plotCrossHeatmap(resGFM, QCThreshold = 0, QCMetric = "medianMinIntensity") +
  theme(legend.position = "bottom") + guides(shape = "none") +
    theme(legend.position = "bottom") + theme(
      axis.title   = element_text(size = 18),
      axis.text    = element_text(size = 14),
      strip.text   = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      plot.title   = element_text(size = 20)
    )


##### let's look more closely at tumour-fibroblasts (CAF?) ####
mdl <- resG$Tumour_Fibroblast$mdl
designmat <- resG$Tumour_Fibroblast$designmat

#how does the model fit?
qqnorm(mdl$residuals, pch = 16)
qqline(mdl$residuals)

#look at some functional boxplots 

metricRes <- resG$Tumour_Fibroblast$metricRes
#let's recover the patient_id
metricRes$Patient_ID <- sub("^[^_]+_", "", metricRes$imageId)

#remove NA stage for the plots
metricRes <- metricRes |> subset(!is.na(Stage_coarse))
metricRes$Stage_coarse <- factor(metricRes$Stage_coarse, levels <- c("I", "II", "III", "IV"))
metricRes <- metricRes[order(metricRes$Stage_coarse),]

pdf(snakemake@output[["fbplotTF"]])
par(mfrow = c(2, 2))
  plotFbPlot(metricRes, "r", "rs", aggregateBy = "Stage_coarse", sampleId = "Patient_ID", imageId = "imageId")
  graphics::title(main = "G-function boxplot Tumour -> Fibroblast")
dev.off()

plotLs <- lapply(colnames(designmat), plotMdl,
    mdl = mdl,
    shift = mdl$coefficients[["(Intercept)"]]
)

pTF <- wrap_plots(plotLs) + 
  plot_annotation(title = "fGAMM Tumour -> Fibroblast", 
  theme = theme(plot.title = element_text(size = 15))) & 
  theme(plot.tag = element_text(size = 30)) 

##### let's look more closely at vessels-fibroblasts ####
mdl <- resG$vessel_Fibroblast$mdl
designmat <- resG$vessel_Fibroblast$designmat

#how does the model fit?
qqnorm(mdl$residuals, pch = 16)
qqline(mdl$residuals)

#look at some functional boxplots 

metricRes <- resG$vessel_Fibroblast$metricRes
#let's recover the patient_id
metricRes$Patient_ID <- sub("^[^_]+_", "", metricRes$imageId)

#remove NA stage for the plots
metricRes <- metricRes |> subset(!is.na(Stage_coarse))
metricRes$Stage_coarse <- factor(metricRes$Stage_coarse, levels <- c("I", "II", "III", "IV"))
metricRes <- metricRes[order(metricRes$Stage_coarse),]

pdf(snakemake@output[["fbplotVF"]])
par(mfrow = c(2, 2))
  plotFbPlot(metricRes, "r", "rs", aggregateBy = "Stage_coarse", sampleId = "Patient_ID", imageId = "imageId")
dev.off()

plotLs <- lapply(colnames(designmat), plotMdl,
    mdl = mdl,
    shift = mdl$coefficients[["(Intercept)"]]
)

pVF <- wrap_plots(plotLs) +
  plot_annotation(title = "fGAMM Vessel -> Fibroblast", 
  theme = theme(plot.title = element_text(size = 15))) & 
  theme(plot.tag = element_text(size = 35))


### plot some stage spatial plots ###

sce <- readRDS(snakemake@input[["sce"]])
#sce <- readRDS("data/Cords/SingleCellExperiment\ Objects/sce_all_annotated.rds")
df <- data.frame(colData(sce))
df$Stage_coarse <- dplyr::case_when(
  df$Stage %in% c(1, 2) ~ "I",
  df$Stage %in% c(3, 4) ~ "II",
  df$Stage %in% c(5, 6) ~ "III",
  df$Stage == 7 ~ "IV",
  TRUE ~ NA_character_
)
df$imageId <- paste0(df$ImageNumber,"_", df$Patient_ID, "_", df$Stage_coarse)
rm(sce)
gc()
dfSub <- df %>%
    subset(imageId %in% c("41_175_141_I", "41_87_218_II", "41_86_106_III", "42_87_286_IV"))

p <- ggplot(dfSub, aes(x = Center_X, y = Center_Y, color = cell_category)) +
    geom_point(size= 0.5) +
    facet_wrap(~factor(imageId, labels = c("41_175_141_I" = "I", 
"41_87_218_II" = "II",
"41_86_106_III" = "III",
"42_87_286_IV" = "IV")), ncol = 4) +
    xlab("x") +
    ylab("y") +
    labs(color = "cell category")+
    coord_equal() +
    theme_light() +
    theme(legend.position = "bottom") + theme(
      axis.title   = element_text(size = 18),
      axis.text    = element_text(size = 14),
      strip.text   = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      plot.title   = element_text(size = 20)
    ) + guides(colour = guide_legend(
        title = "Cell categories",
        override.aes = list(size = 5)
    ))

pTotal <- p/pG/(wrap_plots(list(wrap_elements(pTF), wrap_elements(pVF)), 
        widths = c(1,1), ncol = 2)) + 
        plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(size = 25))) & 
        theme(plot.tag = element_text(size = 30))

ggsave(snakemake@output[["pCords"]], plot = pTotal, width = 14, height = 14)
ggsave(snakemake@output[["pCordsSupp"]], plot=pGFM, width = 12, height = 9)
