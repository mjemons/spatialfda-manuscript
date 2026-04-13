library("ggplot2")
library("dplyr")
library("patchwork")

lys <- lapply(snakemake@input[["ls"]], readRDS)
lys_sTable <- lapply(lys, function(elem){
			       df <- as.data.frame(elem) %>% select(c("method", "comp", "p.value", "std", "prop")) 
			         return(df)})

df <- dplyr::bind_rows(lys_sTable)
saveRDS(df, "outs/pvalueDF.rds")
df_plt <- df %>% group_by(method, comp, std, prop) %>% count(p.value < 0.05) %>%           
    mutate(power = prop.table(n)) %>%
    filter(`p.value < 0.05` == TRUE)
labels <- c("0", "10", "20", "30", "40", "50", "60", "70", "80")

#add a mutate step to have shapes by method and then colour by method flavour
df_plt <- df_plt %>% mutate(method_class = case_when(
		method == "spatialFDAG" | method == "spatialFDAL" ~ "spatialFDA",
		method == "spicyRLM" | method == "spicyRMM" ~ "spicyR",
		method == "spaceANOVAUni" | method == "spaceANOVAMulti" ~ "spaceANOVA",
		method == "smoppix" ~ "smoppix",
		method == "intensityMM" ~ "intensityMM",
    method == "mxfdaFM" | method == "mxfdaMM" ~ "mxfda")
 
)

df_plt <- df_plt %>%
  mutate(method = case_when(
    method == "spicyRLM" ~ "spicyR.LM",
    method == "spicyRMM" ~ "spicyR.MM",
    method == "spaceANOVAUni" ~ "SpaceANOVA.Uni",
    method == "spaceANOVAMulti" ~ "SpaceANOVA.Multi",
    method == "spatialFDAL" ~ "spatialFDA.L",
    method == "spatialFDAG" ~ "spatialFDA.G",
    method == "smoppix" ~ "smoppix",
		method == "intensityMM" ~ "intensity.MM",
		method == "mxfdaFM" ~ "mxfda",
		method == "mxfdaMM" ~ "mxfda.MM"
  )) 

colors <- c(intensity.MM = "#A6CEE3", smoppix = "#1F78B4", SpaceANOVA.Uni = "#33A02C", SpaceANOVA.Multi = "#B2DF8A",
	 spatialFDA.L = "#E31A1C", spatialFDA.G = "#FB9A99", spicyR.LM ="#FDBF6F", spicyR.MM = "#FF7F00", mxfda = "#CAB2D6",
	mxfda.MM = "#6A3D9A")

p <- ggplot(df_plt, aes(x = comp, y = power, col = as.factor(method), group = as.factor(method), shape = as.factor(method_class)), linewidth = 1) + 
       	geom_point() + 
	geom_line() +
       	geom_hline(yintercept = 0.05, linetype="dashed") + 	
	theme_light() + 
	xlab("difference in mean probability") + 
	scale_x_discrete(label = labels) +
	scale_y_sqrt() +
	facet_wrap(~prop) +
	scale_colour_manual(values = colors) +
	guides(col=guide_legend(title="method"),
	shape = guide_legend(title="method_class"))

ggplot2::ggsave(snakemake@output[["plt"]], width = 10, height = 10)

df_dens <- df %>% group_by(method, comp, std, prop) %>% filter(comp == "pert1")
df_dens$plotting_id <- paste0(df_dens$std, "|", df_dens$prop)

#inspired from Simone Tiberi distinct paper figure 3
q <- ggplot(df_dens, aes(`p.value`, fill = std)) + 
  geom_histogram(breaks = seq(0, 1, 0.05)) +
  theme_light() +  
  facet_wrap(~ method, scales = "free_y") +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
        guides(fill=guide_legend(title="noise"),
	       col =guide_legend(title="noise"))
  
ggplot2::ggsave("outs/pValueNull.pdf", width = 10, height = 7)

pSupplement <- p/q
pSupplement <- pSupplement + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))

ggplot2::ggsave(snakemake@output[["supplement"]], width = 10, height = 10)

df_plt <- df_plt %>% filter(std == 1)
p2 <- ggplot(df_plt, aes(x = comp, y = power, col = as.factor(method), group = as.factor(method), shape = as.factor(method_class)), linewidth = 1) +
        geom_point() +
        geom_line() +
        geom_hline(yintercept = 0.05, linetype="dashed") +
        theme_light() +
        xlab("difference in mean probability") +
        scale_x_discrete(label = labels) +
        scale_y_sqrt() +
	facet_wrap(~prop) +
	scale_color_brewer(palette = "Paired") +
        guides(col=guide_legend(title="method"),
        shape = guide_legend(title="method_class")) #+
        #theme(legend.position = "bottom")

pTotal <- wrap_plots(list(p2, q), heights = c(1,2), nrow = 2)
pTotal <- pTotal + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20))

ggplot2::ggsave(snakemake@output[["combined"]], plot = pTotal, width = 9, height = 9)	
