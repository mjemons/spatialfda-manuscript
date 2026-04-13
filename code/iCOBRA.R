library("ggplot2")
library("dplyr")
library("tidyr")
library("iCOBRA")
library("patchwork")
set.seed(1234)

lys <- lapply(snakemake@input[["ls"]], readRDS)
lys_sTable <- lapply(lys, function(elem){
			       df <- as.data.frame(elem) %>% select(c("method", "comp", "p.value", "std", "prop")) 
			         return(df)})

df <- dplyr::bind_rows(lys_sTable)

dfTotal <- df %>% filter(comp %in% c("pert1", "pert3"))
dfTotal <- dfTotal %>%
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

dfTotal$p_0 <- dfTotal$prop
nSim = 500

iCobraPlot <- function(dfTotal, overall){
	#df <- dfTotal %>% filter(std == stdVal)
	df <- dfTotal
	df <- df %>% group_by(method, comp, std, p_0) %>% mutate(repetition = rep(1:nSim, length.out = n()))

	df$ID <- paste0(df$comp, df$p_0, df$std, df$repetition)

	df_wide <- df %>% pivot_wider(names_from = "method", values_from = "p.value") %>%
  	mutate(truth = comp == "pert3") %>% as.data.frame()

	rownames(df_wide) <- df_wide$ID

	pval <- df_wide %>% ungroup() %>% select(-c(truth, comp, std, p_0, repetition, ID, prop))

	truth <- df_wide %>% ungroup() %>% select(c("ID", "truth", "p_0", "std"))

	cobraData <- COBRAData(pval = pval, truth = truth)

	cobraData <- calculate_adjp(cobraData)
	
	if(overall){
		cobraperf <- calculate_performance(cobraData, binary_truth = "truth")
	}else{
		cobraperf <- calculate_performance(cobraData, binary_truth = "truth", splv = "p_0")
	}

	colors <- c(intensity.MM = "#A6CEE3", smoppix = "#1F78B4", SpaceANOVA.Uni = "#33A02C", SpaceANOVA.Multi = "#B2DF8A",
	 spatialFDA.L = "#E31A1C", spatialFDA.G = "#FB9A99", spicyR.LM ="#FDBF6F", spicyR.MM = "#FF7F00", mxfda = "#CAB2D6",
	mxfda.MM = "#6A3D9A")
  
  update_geom_defaults("line", list(key_glyph = draw_key_point))
	
	cobraplot <- prepare_data_for_plot(cobraperf,
                                  colorscheme = colors,
                                  facetted = TRUE,
																	incloverall = overall)
	
	cobraplot <- reorder_levels(cobraplot, c("intensity.MM", "smoppix", "spicyR.LM", "spicyR.MM", "SpaceANOVA.Multi",
"SpaceANOVA.Uni", "spatialFDA.G", "spatialFDA.L", "mxfda", "mxfda.MM"))

	return(cobraplot)
}

plotFDPROCCurves <- function(dataSplit){
	fdr_obs <- fdrtpr(dataSplit) |> select(method, thr, FDR, TPR, splitval)
	fdr_obs$FDRMCSE <- sqrt(fdr_obs$FDR * (1-fdr_obs$FDR)/nSim)

	fdr_obs$fdrhi <- fdr_obs$FDR + 1.96 * fdr_obs$FDRMCSE 
	fdr_obs$fdrlo <- fdr_obs$FDR - 1.96 * fdr_obs$FDRMCSE 
	fdr_obs$splitval <- gsub("p_(\\w+)", "p[\\1]", fdr_obs$splitval)
		
	pROC <- plot_roc(dataSplit, linewidth = 1) + xlim(0,0.4) + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
					axis.text.y = element_text(size = 14),
					axis.title.x = element_text(size = 14),
					axis.title.y = element_text(size = 14),
					plot.title = element_text(size = 14, face = "bold"),
					legend.text = element_text(size = 14),
					strip.text.x = element_text(size = 14),
					legend.position = "bottom",
					legend.key = element_blank(),           
					legend.background = element_blank(), 
					legend.box.background = element_blank()
			) +
		guides(colour = guide_legend(override.aes = list(
			shape = 16,
			size = 5,
			stroke = 0,
			alpha = 1   
		)))
		# written by claude.ai
	pROC$data$splitval <- gsub("p_(\\w+)", "p[\\1]", pROC$data$splitval)

	# Now tell ggplot to parse those strings as plotmath
	pROC$facet$params$labeller <- label_parsed

	#change layout of facet_wrap and then add overall on top. 
	pROC$facet$params$ncol <- 3
	pROC$facet$params$nrow <- NULL

	pFDP <- plot_fdrtprcurve(dataSplit, pointsize = 4, linewidth = 1, plottype = "points") +  
			xlab("achieved FDR") + 
			geom_errorbar(
			data        = fdr_obs,
			mapping     = aes(y = TPR,
												xmin = fdrlo,
												xmax = fdrhi,
												colour = method),
			width      = 0.03,
			linewidth   = 0.75,
			inherit.aes = FALSE
			) +
			scale_x_continuous(limits = c(NA, NA)) +
			coord_cartesian(xlim = c(0, 0.25)) +
			facet_wrap(~splitval, labeller = label_parsed) +
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
																				hjust = 1, size = 14),
						axis.text.y = element_text(size = 14),
						axis.title.x = element_text(size = 14),
						axis.title.y = element_text(size = 14),
						plot.title = element_text(size = 14, face = "bold"),
						legend.text = element_text(size = 14),
						strip.text.x = element_text(size = 14),
						legend.position = "bottom",
						legend.key = element_blank(),           
						legend.background = element_blank(), 
						legend.box.background = element_blank()
			) +
			guides(colour = guide_legend(override.aes = list(
			shape = 16,
			size = 5,
			stroke = 0,
			alpha = 1   
			)))
	# written by claude.ai
	pFDP$data$splitval <- gsub("p_(\\w+)", "p[\\1]", pFDP$data$splitval)

	# Now tell ggplot to parse those strings as plotmath
	pFDP$facet$params$labeller <- label_parsed

	#change layout of facet_wrap and then add overall on top. 
	pFDP$facet$params$ncol <- 3
	pFDP$facet$params$nrow <- NULL
	return(c(pFDP, pROC))
}


dataSplit <- iCobraPlot(dfTotal = dfTotal, overall = FALSE)

pLs <- plotFDPROCCurves(dataSplit)

pFDP <- pLs[1]
pROC <- pLs[2]


dataOverall <- iCobraPlot(dfTotal = dfTotal, overall = TRUE)
pLs <- plotFDPROCCurves(dataOverall)

pFDPOverall <- pLs[1]
pROCOverall <- pLs[2]

pFDP <- (pFDPOverall[[1]]/pFDP[[1]]) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(size = 20))) & 
  theme(plot.tag = element_text(size = 20), legend.position = 'bottom')
pROC <- (pROCOverall[[1]]/pROC[[1]]) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(size = 20))) & 
  theme(plot.tag = element_text(size = 20), legend.position = 'bottom')

ggsave(snakemake@output[["plt"]], plot = pFDP, width = 9, height = 9)
ggsave(snakemake@output[["roc"]], plot = pROC, width = 9, height = 9)
