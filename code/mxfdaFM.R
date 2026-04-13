library("SpatialExperiment")
library("mxfda")
library("dplyr")

spe <- readRDS(snakemake@input[["rds"]])
comp <- snakemake@wildcards[["comp"]]

spe <- subset(spe, ,condition %in% c("ctrl", comp))

colData(spe)[["condition"]] <- factor(colData(spe)[["condition"]])
#relevel to have non-diabetic as the reference category
colData(spe)[["condition"]] <- relevel(colData(spe)[["condition"]],
"ctrl")
#rename image ID
colData(spe)[["image_id"]] <- colData(spe)[["ID"]]

meta <- as_tibble(colData(spe)) |> select(-c("X.1", "X", "Labels"))
meta$sample_id <- as.factor(meta$sample_id)
#response needs to be 0<=y<=1 for logistic regression
meta$condition_binary <- ifelse(meta$condition == "ctrl", 0, 1)
meta <- meta |> unique()

spatial <- as_tibble(spatialCoords(spe))
#assumes same ordering of ids
spatial$cellType <- as.character(as_tibble(colData(spe))$Labels)
spatial$image_id <-factor(as_tibble(colData(spe))$image_id)

mxFDAobject = make_mxfda(metadata = meta,
                         spatial = spatial,
                         subject_key = "sample_id",
                         sample_key = "image_id")

mxFDAobject = extract_summary_functions(mxFDAobject,
                                        extract_func = univariate,
                                        summary_func = Kest,
                                        r_vec = seq(0, 100, by = 1),
                                        edge_correction = "iso",
                                        markvar = "cellType",
                                        mark1 = "0")

mxFDAobject <- run_sofr(mxFDAobject, 
                        model_name = "fit_sofr_condition", 
                        formula = condition_binary ~ 1, 
                        family = "binomial",
                        metric = "uni k", r = "r", value = "fundiff",
                        optimizer = "efs")

mdl = extract_model(mxFDAobject, 'uni k', type = 'sofr', model_name = 'fit_sofr_condition')

out <- summary(mdl, re.test = FALSE)
out$prob <- unique(colData(spe)$prob)

saveRDS(out, snakemake@output[["rds"]])