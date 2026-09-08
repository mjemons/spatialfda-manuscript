### written with the help of claude.ai and GPT 5.6 ###

library("furrr")
library("purrr")
library("dplyr")
library("SpatialExperiment")

nCores <- 15
plan(multisession, workers = nCores)

library("spatialFDA")
# load the IMC dataset described in Damond et al. 2019 as SpatialExperiment object
spe <- .loadExample(full = TRUE)

colData(spe)[["patient_stage"]] <- factor(colData(spe)[["patient_stage"]])
#relevel to have non-diabetic as the reference category
colData(spe)[["patient_stage"]] <- relevel(colData(spe)[["patient_stage"]],
"Non-diabetic")

patients <- unique(colData(spe)$patient_id)

run_full_matrix <- function(spe_subset) {
  resLs <- spatialFDA::crossSpatialInference(
    spe_subset, 
    selection = NULL,
    fun = "Gcross", 
    marks = "cell_type",
    rSeq = seq(0, 100, length.out = 50), 
    correction = "rs",
    sample_id = "patient_id",
    family = gaussian(link = "log"),
    algorithm = "bam",
    image_id = "image_number", 
    condition = "patient_stage",
    ncores = 1
)

  res <- spatialFDA::extractCrossInferenceData(resLs) %>%
    mutate(q_value = p.adjust(`p-value`, method = "BH"),
           sig = q_value < 0.05)
  return(res)
}

full_result <- run_full_matrix(spe) %>%
  dplyr::rename(estimate_full = mean_coefficient, sig_full = sig)

loo_result <- future_map_dfr(patients, function(pid) {
  spe_sub <- spe[, colData(spe)$patient_id != pid]
  run_full_matrix(spe_sub) %>%
    mutate(left_out = pid)
}, .options = furrr_options(seed = TRUE))

loo_joined <- loo_result %>%
  left_join(
    full_result %>% select(cell1, cell2, coefficient, sig_full, estimate_full),
    by = c("cell1", "cell2", "coefficient")
  )

# Per (pair, coefficient) stability — this is your reportable table
term_stability <- loo_joined %>%
  group_by(cell1, cell2, coefficient, sig_full) %>%
  summarise(
    #full dataset estimate
    estimate_full = dplyr::first(estimate_full),
    #proportion of significant results between full and LOO
    concordance_rate = mean(sig == sig_full),
    #jackknife standard errors
    n_folds = sum(!is.na(mean_coefficient)),
    jackknife_se = if(n_folds == length(patients)){
      sqrt((length(patients) - 1) / length(patients)  *
      sum(
        (mean_coefficient -
           mean(mean_coefficient, na.rm = TRUE))^2,
        na.rm = TRUE
      ))
    }else{
      NA_real_
    },
    .groups = "drop"
  ) %>%
  arrange(desc(sig_full), concordance_rate)

# Jaccard index per fold across the full test set (not just per-pair)
jaccard_per_fold <- loo_joined %>%
  group_by(left_out) %>%
  summarise(jaccard = sum(sig & sig_full) / sum(sig | sig_full))

res <- list(
  "loo_joined" = loo_joined,
  "term_stability" = term_stability,
  "jaccard_per_fold" = jaccard_per_fold
)

saveRDS(res, snakemake@output[["rds"]])