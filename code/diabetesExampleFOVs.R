library("SpatialExperiment")
library("dplyr")
library("ggplot2"); theme_set(theme_light())
library("patchwork")
library("spatialFDA")
# load the IMC dataset described in Damond et al. 2019 as SpatialExperiment object
spe <- .loadExample(full = TRUE)


### adapted from the spatialFDA vignettes ###
### healthy islets

df <- data.frame(spatialCoords(spe), colData(spe))

df <- df %>%
  subset(cell_type %in% c("alpha", "beta", "delta", "Th", "Tc"))

dfSub <- df %>%
    subset(image_name %in% c("E02", "E03", "E04"))

pHealthy <- ggplot(dfSub, aes(x = cell_x, y = cell_y, color = cell_type)) +
    geom_point(size= 1.5) +
    facet_wrap(~image_name) +
    theme(legend.title.size = 20, legend.text.size = 20) +
    xlab("x") +
    ylab("y") +
    labs(color = "cell category")+
    coord_equal() +
    theme_light() + 
    ggtitle("Non diabetic")

dfSub <- df %>%
    subset(image_name %in% c("A01", "A02", "A03"))

pOnset <- ggplot(dfSub, aes(x = cell_x, y = cell_y, color = cell_type)) +
    geom_point(size= 1.5) +
    facet_wrap(~image_name) +
    theme(legend.title.size = 20, legend.text.size = 20) +
    xlab("x") +
    ylab("y") +
    labs(color = "cell category")+
    coord_equal() +
    theme_light() + 
    ggtitle("Onset diabetes")

dfSub <- df %>%
    subset(image_name %in% c("J13", "J33", "Q07"))

pLongduration <- ggplot(dfSub, aes(x = cell_x, y = cell_y, color = cell_type)) +
    geom_point(size= 1.5) +
    facet_wrap(~image_name) +
    theme(legend.title.size = 20, legend.text.size = 20) +
    xlab("x") +
    ylab("y") +
    labs(color = "cell category")+
    coord_equal() +
    theme_light() + 
    ggtitle("Long duration diabetes")

pTotal <- wrap_plots(list(pHealthy, pOnset, pLongduration), nrow = 3) + 
  plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(size = 25))) & 
  theme(plot.tag = element_text(size = 30)) +
  plot_layout(guides = "collect")

ggsave(snakemake@output[["plt"]], plot = pTotal, width = 15, height = 15)


