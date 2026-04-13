library(tidyverse)
library(plotROC)
library(SpatialExperiment)
library(ggplot2)

results <- readRDS(snakemake@input[["rds"]])

results <- tidyr::pivot_longer(results,cols = c("spicyRLM", "spicyRMM", "spatialFDAL", "spatialFDAG", "intensityMM.Pr...t.."), names_to = "method", values_to = "p.value")

results <- results %>%
  mutate(method = case_when(
    method == "spicyRLM" ~ "spicyR.LM",
    method == "spicyRMM" ~ "spicyR.MM",
    method == "spatialFDAL" ~ "spatialFDA.L",
    method == "spatialFDAG" ~ "spatialFDA.G",
    method == "intensityMM.Pr...t.." ~ "intensity.MM"
  )) 

results$p.value <- unlist(results$p.value)

results <- results %>% filter(lambda %in% c(10,30,50,70,90,100))

results$lambda <- as.factor(results$lambda)
results$method <- factor(results$method, levels = c("intensity.MM", "spatialFDA.G", "spatialFDA.L", "spicyR.LM", "spicyR.MM"))

results <- results %>% dplyr::mutate(is_true = as.numeric(simulation == "Signal"))

colors <- c(spatialFDA.L = "#E31A1C", spatialFDA.G = "#FB9A99", spicyR.LM ="#FDBF6F", spicyR.MM = "#FF7F00", intensity.MM = "#A6CEE3")

g1 <- ggplot(results, aes(m = p.value, d = simulation, colour = method)) + geom_roc(n.cuts = 0, increasing = FALSE, key_glyph = "point", size = 2, linealpha = 0.5) + 
  theme_classic() + facet_wrap(~lambda, nrow = 2) + xlab("FPR") + ylab("TPR") +
  scale_colour_manual(values = colors) + 
  xlim(0,0.25) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5,
                                      hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold"),
		    legend.text = element_text(size = 14),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 14)
        ) +
  guides(colour = guide_legend(override.aes = list(
    shape = 16,
    size = 5,
    linetype = 0 
  )))

ggsave(snakemake@output[["plt"]], plot = g1, width = 10, height = 7)