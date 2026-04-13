#BiocManager::install("EBImage")

devtools::install_github('sealx017/SpaceANOVA', quiet = TRUE)

# write touch file
x <- data.frame()
write.table(x, file='outs/spaceANOVA_installed', col.names=FALSE)
