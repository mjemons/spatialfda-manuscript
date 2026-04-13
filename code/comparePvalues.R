library("ggplot2")
comp <- snakemake@wildcards[["comp"]]
prob <- snakemake@wildcards[["prob"]]
method <- snakemake@wildcards[["method"]]
std <- snakemake@wildcards[["std"]]
prop <- snakemake@wildcards[["prop"]]

lys <- lapply(snakemake@input[["ls"]], readRDS)

lys_sTable <- lapply(lys, function(elem){
  	
	if(method == "spatialFDAL"){
		df <- as.data.frame(elem$s.table)[paste0('condition', comp, '(x)'),]
		df <- dplyr::rename(df, p.value = `p-value`)
	}
	if(method == "spatialFDAG"){
		df <- as.data.frame(elem$s.table)[paste0('condition', comp, '(x)'),]
		df <- dplyr::rename(df, p.value = `p-value`)
	}	
	if(method == "spicyRMM" || method == "spicyRLM"){
		df <- as.data.frame(elem)
	}
	if(method == "spaceANOVAUni" || method == "spaceANOVAMulti"){
		df <- data.frame(p.value = elem[["0","0"]])
	}
	if(method == "smoppix"){
		df <- data.frame(p.value = elem[[paste0('condition1'), "Pr(>|t|)"]])
	}	
	if(method == "intensityMM"){
		df <- data.frame(p.value = elem[[paste0('condition', comp), "Pr(>|t|)"]])
	}
	if(method == "mxfdaFM" || method == "mxfdaMM"){
		df <- as.data.frame(elem$s.table)["s(xmat.tmat):L.xmat",]
		df <- dplyr::rename(df, p.value = `p-value`)
	}
		df$comp <- comp
  	df$prob <- prob
  	df$method <- method
		df$std <- std
  	df$prop <- prop 	
  	return(df)
})
df <- dplyr::bind_rows(lys_sTable)

saveRDS(df, snakemake@output[["rds"]])

