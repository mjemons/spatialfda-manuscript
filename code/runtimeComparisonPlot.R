library(ggplot2)
library(dplyr)

results <- readRDS(snakemake@input[["rds"]])

results <- results %>%
  mutate(method = case_when(
    method == "SpaceAnova" ~ "SpaceANOVA",
    .default = method
  )) 

colors <- c(smoppix = "#1F78B4", SpaceANOVA = "#33A02C", spatialFDA.L = "#E31A1C", spatialFDA.G = "#FB9A99", spicyR.LM ="#FDBF6F", spicyR.MM = "#FF7F00", mxfda = "#CAB2D6")

p <- ggplot(results, aes(x = nIm, y = elapsed_time_seconds, col = method)) +
  geom_point() +
  geom_line() +
  theme_light() + 
  scale_y_log10() +
  scale_colour_manual(values = colors) +
  ylab("Elapsed time [seconds]") +
  xlab("Number of images") +
  ggtitle("Runtime comparison")+
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
        )

ggsave(snakemake@output[["plt"]], plot = p, width = 7, height = 5)