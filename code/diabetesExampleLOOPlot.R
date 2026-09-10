### coded with the help of chatGPT ###
library("dplyr")
library("ggplot2")
library("patchwork")

df <- readRDS(snakemake@input[["rds"]])

term_stability      <- df$term_stability
jaccard_per_fold    <- df$jaccard_per_fold

pJaccard <- jaccard_per_fold %>%
  mutate(
    left_out = factor(left_out),
    left_out = reorder(left_out, jaccard)
  ) %>%
  ggplot(aes(x = left_out, y = jaccard)) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = "Patient omitted",
    y = "Jaccard similarity with full-cohort discoveries",
    title = "Leave-one-patient-out stability of significant results"
  ) +
  theme_bw(base_size = 12)

sig_terms <- term_stability %>%
  filter(sig_full, !is.na(estimate_full), !is.na(jackknife_se)) %>%
  mutate(
    pair = paste(cell1, "\u2192", cell2),
    pair = reorder(pair, concordance_rate)
  )

pStabilityPairs <- ggplot(
  sig_terms,
  aes(x = concordance_rate, y = pair)
) +
  geom_point(size = 2) +
  facet_wrap(
    ~ coefficient,
    scales = "free_y"
  ) +
  labs(
    x = "Leave-one-patient-out significance concordance",
    y = "cell-type pair",
    title = "Stability of significant and estimable spatial associations"
  ) +
  theme_bw(base_size = 11)

effect_terms <- term_stability %>%
  filter(sig_full, !is.na(estimate_full), !is.na(jackknife_se)) %>%
  mutate(
    pair = paste(cell1, "\u2192", cell2),
    pair = reorder(pair, concordance_rate),
    hi = estimate_full + jackknife_se,
    lo = estimate_full - jackknife_se
  )

pEffectPairs <- ggplot(
  effect_terms,
  aes(x = estimate_full, y = pair)
) +
  geom_errorbar(
    aes(
      xmin = lo,
      xmax = hi
    ),
    orientation = "y",
    width = 0.2,
    linewidth = 0.75
  ) +
  geom_point(size = 2) +
  facet_wrap(
    ~ coefficient,
    scales = "free_y"
  ) +
  labs(
    x = "Full-cohort coefficient estimate",
    y = "cell-type pair",
    title = "Stability of CCoL effect sizes"
  ) +
  theme_bw(base_size = 11)

pTotal <- pJaccard / pStabilityPairs / pEffectPairs + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(size = 20))) & 
  theme(plot.tag = element_text(size = 20), legend.position = 'right')

ggsave(snakemake@output[["plt"]], plot = pTotal, width = 10, height = 13)