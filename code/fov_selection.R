#file <- read.csv("outs/1/ctrl/patient1/dataframe-celltype-0.7-patient1-1.csv")
library('dplyr')

seed <- as.numeric(snakemake@wildcards[['rep']])*40 + as.numeric(snakemake@wildcards[['fov']])*5 + as.numeric(snakemake@wildcards[['subfov']]) + as.numeric(snakemake@wildcards[['std']])*3 + as.numeric(snakemake@wildcards[['prop']])*100 + as.numeric(snakemake@wildcards[['prob']])*10

set.seed(seed) 
#subfov <- as.integer(snakemake@wildcards[['subfov']])
subfov <- round(runif(1, min = 1, max = 4))
file <- read.csv(snakemake@input[['file']])

find_fov_edge <- function(max_window, fov, fov_length){
	divisor <- max_window / fov_length
  	row <- ((fov -(fov %% divisor)) / divisor) + 1
    	col <- fov %% divisor
    	if(fov %% divisor == 0){
		row <- ((fov -(fov %% divisor)) / divisor)
        	col = divisor
	 }
	return(list(row = row, col = col))
}

fov_length <- 500
max_window <- 1000
#find edge of window - assumes the window is quadratic
coords <- find_fov_edge(max_window, subfov, fov_length)
x_coord <- fov_length * coords[['col']]
y_coord <- fov_length * coords[['row']]

# subset dataframe
fov_df <- file %>% filter(x>x_coord - fov_length & x<x_coord) %>%
	       filter(y>y_coord - fov_length & y<y_coord)

fov_df$ID <- paste0(fov_df$sample_id, "_", snakemake@wildcards[['fov']], "_", snakemake@wildcards[['subfov']])

write.csv(fov_df, snakemake@output[['file']])
