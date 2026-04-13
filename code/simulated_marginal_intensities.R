### Intensity Boxplot figure
library("SpatialExperiment")
library("dplyr")
library("ggplot2"); theme_set(theme_light())
library("patchwork")
#code adapted for intensity plot from Elizabeth Purdom
#rename image ID

#' Convert SpatialExperiment object to ppp object
#' CODE FROM `spatialFDA`
#'
#' @param df A dataframe with the x and y coordinates from the corresponding
#' SpatialExperiment and the ColData
#' @param marks A vector of marks to be associated with the points, has to be
#' either named 'cell_type' if you want to compare discrete celltypes or else
#' continous gene expression measurements are assumed as marks.
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param window An observation window of the point pattern of class `owin`.
#' @return A ppp object for use with `spatstat` functions
#' @export
#'
#' @examples
#' # retrieve example data from Damond et al. (2019)
#' spe <- .loadExample()
#' speSub <- subset(spe, , image_number == "138")
#' dfSub <- .speToDf(speSub)
#' pp <- .dfToppp(dfSub, marks = "cell_type")
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
.dfToppp <- function(df, marks = NULL, continuous = FALSE, window = NULL) {
    #type checking
    stopifnot(is(df, "data.frame"))
    # this definition of the window is quite conservative
    # - can be set explicitly
    pp <- spatstat.geom::as.ppp(data.frame(x = df$x, y = df$y),
        W = spatstat.geom::owin(
            c(
                base::min(df$x) - 1,
                base::max(df$x) + 1
            ),
            c(
                base::min(df$y) - 1,
                base::max(df$y) + 1
            )
        )
    )
    # set the marks
    if (!continuous) {
        spatstat.geom::marks(pp) <- factor(df[[marks]])
    } else {
        spatstat.geom::marks(pp) <- base::subset(df, select =
                                                   names(df) %in% marks)
    }
    # if window exist, set is as new window and potentially exclude some points
    if (!is.null(window)) {
        pp <- spatstat.geom::as.ppp(spatstat.geom::superimpose(pp, W = window))
    }

    return(pp)
}

lys <- lapply(snakemake@input[["ls"]], readRDS)

df <- dplyr::bind_rows(lys)

selection <- df[["Labels"]] %>% unique() 

ConditionIntensities <- function(df, marks, conditionType, imageId, selection){
  df <- df %>% filter(condition == conditionType)

  dfLs <- base::split(df, df[[imageId]])
  
  intensityDfCellType <- lapply(selection, function(x){
    intensitiesDf <- lapply(dfLs, function(dfSub){
      pp <- .dfToppp(dfSub, marks = marks, continuous = FALSE, window = NULL)
      ppSub <- pp[pp$marks == x, drop = TRUE]
      spatstat.geom::marks(ppSub) <- factor(spatstat.geom::marks(ppSub),
                                                levels = unique(x))
      cellTypeIntensity <- data.frame(intensity = spatstat.geom::intensity(ppSub),
                                      sample_id = unique(dfSub$sample_id),
                                      condition = unique(dfSub$condition),
                                      row.names = NULL)
      return(cellTypeIntensity)
    }) %>% bind_rows()
    intensitiesDf$cellType <- x
    return(intensitiesDf)
  }) %>% bind_rows()
  return(intensityDfCellType)
}

intensityDf <- lapply(list("ctrl", "pert1", "pert3"), function(x){
  ConditionIntensities(df = df,
                       marks = "Labels",
                       conditionType = x,
                       imageId = "ID",
                       selection = selection)
})%>% bind_rows()

intensityDf <- intensityDf %>%
  mutate(conditionTyp = case_when(
    condition == "ctrl"  ~ "H['0,0']^ctrl == 0.5",
    condition == "pert1" ~ "H['0,0']^ctrl == 0.5",
    condition == "pert3" ~ "H['0,0']^pert == 0.3"
  ))

intensityDf <- intensityDf %>%
  mutate(condition = case_when(
    condition == "ctrl" ~ "ref",
    condition == "pert1" ~ "ctrl",
    condition == "pert3" ~ "pert"
  )) 

intensityDf$condition <- factor(intensityDf$condition,
                                    levels = c( "ref", "ctrl", "pert"))
p <- ggplot(intensityDf, aes(x=condition,group=sample_id,y=intensity,fill=conditionTyp))+ geom_boxplot()+facet_wrap(~cellType,scales="free", ncol = 2, nrow = 2) +
  xlab("Interaction probability") +
  scale_fill_discrete(labels = function(x) parse(text = x))

p

ggsave(snakemake@output[["plt"]], plot = p, width = 10, height = 7)