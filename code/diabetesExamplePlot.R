### coded by hand for G function, generalised to plot both G and L by 
### claude.ai

library("dplyr")
library("ggplot2"); theme_set(theme_light())
library("patchwork")
library("spatialFDA")
library("tidyr")
library("spicyR")

out <- readRDS(snakemake@input[["rds"]])
resG <- out$G
resL <- out$L
resSpicyRMM <- out$spicyRMM
resSmoppix <- out$smoppix

## reshape data
coefDf <- resSpicyRMM$coefficient
coefDf$pair <- rownames(coefDf)
coefDfLong <- coefDf %>%
  as.data.frame() %>%
  pivot_longer(
    cols = -c(pair),
    names_to = "coefficient",
    values_to = "Estimate"
  )

pvalDf <- resSpicyRMM$p.value
pvalDf$pair <- rownames(pvalDf)
pvalDfLong <- pvalDf %>%
  as.data.frame() %>%
  pivot_longer(
    cols = -c(pair),
    names_to = "coefficient",
    values_to = "p.value"
  )

resSpicyR_long <- coefDfLong %>%
  dplyr::left_join(
    pvalDfLong,
    by = c("pair", "coefficient")
  ) %>%
  separate(pair, into = c("cell1", "cell2"), sep = "__", remove = FALSE)

resSmoppix <- resSmoppix %>%
  separate(pair, into = c("cell1", "cell2"), sep = "--", remove = FALSE) %>%
  rename(p.value = "Pr(>|t|)") %>%
  mutate(coefficient = recode(
      coefficient,
      "patient_stageOnset" = "conditionOnset",
      "patient_stageLong-duration" = "conditionLong-duration"
    )
  )

buildHeatmap <- function(res) {
  p <- plotCrossHeatmap(res, coefficientsToPlot = c("conditionLong_duration(x)", "conditionOnset(x)"), QCThreshold = 1e-5, QCMetric = "medianMinIntensity")
  p <- p +
    guides(shape = "none") +
    labs(color = "mean coefficient") +
    theme(text = element_text(size = 19), legend.position = "bottom", legend.text = element_text(angle = 45, vjust = 0.1)) +
    facet_wrap(~factor(coefficient, levels = c("conditionOnset(x)", "conditionLong_duration(x)"), labels = c("conditionOnset(x)" = "Onset", "conditionLong_duration(x)" = "Long-Duration")))
  return(p)
}

heatmapCustom <- function(df){
  df[["adj.p-value"]] <- stats::p.adjust(df[["p.value"]], method = "BH")
  df <-  dplyr::filter(df, coefficient != "(Intercept)") 
  p <- ggplot(df, aes(x = .data[["cell1"]], y = .data[["cell2"]])) +
    geom_point(aes(size = -log10(.data[["adj.p-value"]] + 0.001),
                  color = .data[["Estimate"]]))+
    geom_point(aes(size = -log10(.data[["adj.p-value"]] + 0.001)),
              shape = 1,colour = "black")+
    scale_colour_gradient2(midpoint = 0,
                            high = scales::muted("red"),
                            mid = "white",
                            low = scales::muted("blue"),
                            na.value = "black") +
    scale_x_discrete(guide = guide_axis(angle = 50)) +
    guides(shape = "none") +
    labs(color = "mean coefficient") +
    theme(text = element_text(size = 19), legend.position = "bottom", legend.text = element_text(angle = 45, vjust = 0.1)) +
    facet_wrap(~factor(coefficient, levels = c("conditionOnset", "conditionLong-duration"), labels = c("conditionOnset" = "Onset", "conditionLong-duration" = "Long-Duration")))
  return(p)
}

### Main heatmap figure (resG only, as before)

p1 <- buildHeatmap(resG) + ggtitle("Diabetes discoveries spatialFDA.G")

### Supplement: resL heatmap only

p1L <- buildHeatmap(resL) + ggtitle("Diabetes discoveries spatialFDA.L")
pSmoppix <- heatmapCustom(resSmoppix) + ggtitle("Diabetes discoveries smoppix")
pSpicyR <- heatmapCustom(resSpicyR_long) + ggtitle("Diabetes discoveries spicyRMM")
pHeatmapMerged <- p1 / p1L + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))
ggsave(snakemake@output[["heatmapSuppA"]], plot = pHeatmapMerged, width = 12, height = 14)
pHeatmapMerged <- pSmoppix / pSpicyR + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 20))
ggsave(snakemake@output[["heatmapSuppB"]], plot = pHeatmapMerged, width = 12, height = 14)

### Everything below uses resG only

mdlDeltaTh <- resG$delta_Th$mdl

metricResDeltaTh <- resG$delta_Th$metricRes

metricResDeltaTh$ID <- factor(paste0(
    metricResDeltaTh$patient_stage, "|", metricResDeltaTh$patient_id
), levels = c("Non-diabetic|6134","Non-diabetic|6278", "Non-diabetic|6126", "Non-diabetic|6386", "Onset|6362", "Onset|6228", "Onset|6414", "Onset|6380", "Long-duration|6418", "Long-duration|6089", "Long-duration|6180", "Long-duration|6264"))

p4 <- plotMetricPerFov(metricResDeltaTh,
    correction = "rs", x = "r",
    imageId = "image_number", ID = "ID"
) + facet_wrap(~ ID)
colours <- c("#0D0887FF", "#2D0594FF", "#44039EFF", "#5901A5FF", 
"#9512A1FF", "#A72197FF", "#B6308BFF", "#C5407EFF", 
"#DD5E66FF", "#E76E5BFF", "#EF7F4FFF", "#F79044FF")

p4 <- p4 + scale_colour_manual(values = colours) +
   theme(text = element_text(size = 19))

# create a unique ID per row in the dataframe
metricResDeltaTh$ID <- paste0(
    metricResDeltaTh$patient_stage, "x", metricResDeltaTh$patient_id,
    "x", metricResDeltaTh$image_number
)

plotMdlCustom <- function(mdl, predictor, shift = NULL) {
    # type checking
    stopifnot(is(mdl, "pffr"))
    stopifnot(is(predictor, "character"))
    # extract the coefficients from the model
    coef <- coef(mdl)
    if (predictor == "(Intercept)" && !is.null(shift)) {
        #rename as pffr output is without brackets
        predictor = "Intercept"
        coef$smterms[["Intercept(x)"]]$coef$value <-
          exp(coef$smterms[["Intercept(x)"]]$coef$value + shift)
    }
    # get the actual values into a dataframe
    df <- coef$smterms[[paste0(predictor, "(x)")]]$coef
    # plot
    p <- ggplot(df, aes(.data$yindex.vec, .data$value)) +
        geom_line(size = 1) +
        # here, I implement a Wald CI - could be improved
        geom_ribbon(
            data = df,
            aes(ymin = .data$value - 1.96 * .data$se,
                ymax = .data$value + 1.96 * .data$se),
            alpha = 0.3
        ) +
        geom_hline(
            yintercept = 0,
            linetype = "dashed", color = "red", size = 1
        ) +
        ggtitle(predictor) +
        ylab("coefficient estimate") +
        xlab("r") +
        theme_light() + theme(text = element_text(size = 19))
    return(p)
}

plotLs <- lapply(colnames(resG$delta_Th$designmat), plotMdlCustom,
    mdl = resG$delta_Th$mdl,
    shift = resG$delta_Th$mdl$coefficients[["(Intercept)"]]
)

plotLscp <- plotLs

plotLs[[2]] <- plotLscp[[3]] + ggtitle("Onset")
plotLs[[3]] <- plotLscp[[2]] + ggtitle("Long-Duration")

p5 <- patchwork::wrap_plots(plotLs, nrow = 3) + plot_layout(axis_titles = "collect")

pTotal <- p1/(wrap_plots(list(p4, p5), widths = c(2,1), ncol = 2)) + plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(size = 25))) & 
  theme(plot.tag = element_text(size = 30))
pTotal
ggsave(snakemake@output[["heatmap"]], plot = pTotal, width = 16, height = 16)

### QC Figure (resG only)

df <- lapply(resG, function(x){
  if(!is.null(x$mdl))
  rsq <- ((summary(x$mdl)$r.sq))
  else rsq <- NULL
  return(rsq)
}) %>% dplyr::bind_rows() %>% 
  t() %>%
  as.data.frame() %>%
  dplyr::rename("R-sq" = "V1") %>%
  mutate(combination = rownames(.)) %>%
  separate(combination, c("cell1", "cell2"), sep = "_")

#deselect NA columns

p6 <- ggplot(df, aes(x = cell1, y = cell2, fill = `R-sq`)) +
    geom_tile() +
    scale_fill_distiller(direction = 1) + 
    scale_x_discrete(guide = guide_axis(angle = 50)) + 
    theme_light() + 
    ggtitle("Entire model diagnostic") +
    guides(fill=guide_legend(title="adjusted R-sq")) +
    theme(text = element_text(size = 19))
p6

df <- lapply(names(resG), function(x){
  if(!is.null(resG[[x]]$curveFittingQC$residual_standard_errors)){
    rse <- data.frame(rse = resG[[x]]$curveFittingQC$residual_standard_errors)
    rse$coefficient <- resG[[x]]$curveFittingQC$coefficient
  }
  else{
     rse <- NULL
     rse$coefficient <- NULL
  }
  rse$combination <- x
  return(rse)
}) %>% bind_rows() %>%
  separate(combination, c("cell1", "cell2"), sep = "_") %>%
  filter(coefficient %in% c("conditionLong_duration(x)", "conditionOnset(x)"))
p7 <- ggplot(df, aes(x = cell1, y = cell2, fill = rse)) +
    geom_tile() +
    facet_wrap(~factor(coefficient, levels = c("conditionOnset(x)", "conditionLong_duration(x)"), labels = c("conditionOnset(x)" = "Onset", "conditionLong_duration(x)" = "Long-Duration"))) +
    scale_fill_distiller(palette = "RdPu",
                         direction = 1) + 
    scale_x_discrete(guide = guide_axis(angle = 50)) +
    ggtitle("Per Condition diagnostic") +
    theme_light() + 
    theme(text = element_text(size = 19)) +
    guides(fill=guide_legend(title="Residual standard error"))

pTotal <- p6/p7 + plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(size = 25))) & 
  theme(plot.tag = element_text(size = 30))
ggsave(snakemake@output[["qcPlot"]], plot = pTotal, width = 15, height = 15)

## QC figure for delta-Th (resG only)

mdlDeltaTh <- resG$delta_Th$mdl

pdf(snakemake@output[["qcPlotDeltaTh"]])
par(mfrow = c(2, 1))
qqnorm(resid(mdlDeltaTh), pch = 16)
qqline(resid(mdlDeltaTh))
mtext("A", side = 3, line = 1, las = 1, cex = 1.5, adj = 0)
image(cor(resid(mdlDeltaTh)),
      col = colorRampPalette(c("blue", "white", "red"))(20),
      zlim = c(-1, 1))
mtext("B", side = 3, line = 1, las = 1, cex = 1.5, adj = 0)
dev.off()