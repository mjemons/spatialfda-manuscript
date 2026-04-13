remotes::install_github("mjemons/spatialFDA@e70d8210bc46c36616c0ec596ac272631a9c0ec4")
# write touch file
x <- data.frame()
write.table(x, file='outs/spatialFDA_installed', col.names=FALSE)
