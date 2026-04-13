#devtools::install_github("drighelli/SpatialExperiment")
devtools::install_github("sthawinke/smoppix@aa6894c2938bc1847d05eb6c820700f92a53490b")

# write touch file
x <- data.frame()
write.table(x, file='outs/smoppix_installed', col.names=FALSE)
