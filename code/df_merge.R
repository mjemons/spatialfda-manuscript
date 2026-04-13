lys <- lapply(snakemake@input[["ls"]], read.csv)
prob <- snakemake@wildcards[["prob"]]

df <- dplyr::bind_rows(lys)
df$prob <- prob
saveRDS(df, snakemake@output[["rds"]])
