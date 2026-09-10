# optimised for future instead of parallel::mclapply with ChatGPT

## Load packages
library(spicyR)
library(spatstat)
library(tidyr)
library(dplyr)
library(plotROC)
library(SpatialExperiment)
library(spatialFDA)
library(smoppix)
library(mxfda)
library(SpaceANOVA)
library(HDF5Array)

## INITIALISE


seed = 51773
set.seed(seed)
window <- owin(xrange = c(0, 1000),
               yrange = c(0, 1000))

nPatients <- 40 
nIm <- 3
nSim <- 200
nCores <- 10
nsimBoot <- 100
counts <- seq(from = 20, to = 400, by = 10)
Rs <- c(10, 30, 50, 70, 90, 100) #seq(from = 10, to = 100, by = 10)

s1 <- Sys.time()

#code from https://github.com/mjemons/spatialFDA licensed under GPL-3 to mjemons
.speToDf <- function(spe) {
    df <- data.frame(
        x = SpatialExperiment::spatialCoords(spe)[, 1],
        y = SpatialExperiment::spatialCoords(spe)[, 2]
    )
    df <- cbind(df, colData(spe))
}

#code from https://github.com/mjemons/spatialFDA licensed under GPL-3 to mjemons
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

ConditionIntensities <- function(df, marks, conditionType, cellType, imageId, sample_id){
  df <- df %>% filter(condition == conditionType)
  dfLs <- base::split(df, df[[imageId]])
  intensitiesDf <- lapply(dfLs, function(dfSub){
    if(nrow(dfSub)>0){
    pp <- .dfToppp(dfSub, marks = marks, continuous = FALSE, window = NULL)
    ppSub <- pp[pp$marks == cellType, drop = TRUE]
    spatstat.geom::marks(ppSub) <- factor(spatstat.geom::marks(ppSub),
                                              levels = unique(cellType))
    cellTypeIntensity <- data.frame(intensity=spatstat.geom::intensity(ppSub))
    cellTypeIntensity$condition <- unique(dfSub$condition)
    cellTypeIntensity[[sample_id]] <- unique(dfSub[[sample_id]])
    return(cellTypeIntensity)
   }else{
    return(data.frame(intensity = NULL,
		       condition = NULL,
		       condition = NULL))
    }
  }) %>% dplyr::bind_rows()
  return(intensitiesDf)
}

resultsTW = NULL

future::plan(future::multisession, workers = nCores)

for(lam in Rs){
  
  lambda = lam  
  
  ## SIGNAL
  sim  <- function(i, counts, nPatients, nIm, window, lambda){
    
    set.seed(i)
    print(paste0("signal iteration ",lambda," ",i-seed))
    g1 <- rpois(nPatients/2, lambda)
    g2 <- rpois(nPatients/2, lambda + lambda/2)
    adjustSigma = c(g1,g2)+1
    
    x <- c()
    y <- c()
    cellType <- c()
    imageID <- c()
    
    for (p in 1:nPatients) {
      for (j in 1:nIm) {
        sCount1 <- sample(counts,1)
        sCount2 <- sample(counts,1)
        a <- rpoispp(sCount1/1000^2, win = window)
        aDens <- density(a, sigma = adjustSigma[p], kernel = "disc")
        aDens$v <- pmax(aDens$v,0)*sCount2/sCount1
        b <- rpoispp(aDens)
        
        x <- c(x, a$x, b$x)
        y <- c(y, a$y, b$y)
        
        cellType <- c(cellType, rep("A", a$n), rep("B", b$n))
        imageID <- c(imageID, rep(paste(p,j,sep = "_"), a$n+b$n))
      }
    }
    
    imageID <- factor(imageID)
    
    cellExp <- data.frame(
      x = x,
      y = y,
      cellType = factor(cellType),
      imageID = imageID
    )

    phenoData <- data.frame(imageID = unique(imageID),
                            condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
                            subject = rep(1:nPatients, each = nIm))
    
    colDF <- cellExp %>% left_join(phenoData, by = "imageID")

    spe <- SpatialExperiment(
      colData = colDF,
      spatialCoordsNames = c("x", "y"),
      sample_id = as.character(colDF$subject)
    )

    colData(spe)[["condition"]] <- factor(colData(spe)[["condition"]])
    #relevel to have ctrl as the reference category
    colData(spe)[["condition"]] <- relevel(colData(spe)[["condition"]],
    "Group1")
    
    #### spicyR ####
    test.spicyLM <- spicyR::spicy(spe,
                  condition = "condition",
                  from = "B",
                  to = "A",
                  fast = TRUE, Rs = seq(from = 10, to = 100, by = 10),
                  edgeCorrect = FALSE,
		  cores = 1,
                  weights = FALSE,
                  verbose = FALSE)
    
    test.spicyMM <- spicyR::spicy(spe,
                  condition = "condition",
                  from = "B",
                  to = "A",
		  subject = "subject",
                  fast = TRUE, Rs = seq(from = 10, to = 100, by = 10),
                  edgeCorrect = FALSE,
                  cores = 1,
		  weights = FALSE,
                  verbose = FALSE)

    #### spatialFDA ####
    res <- spatialFDA::spatialInference(
                        spe, 
                        selection = c("B", "A"), 
                        fun = "Gcross", 
                        marks = "cellType",
                        rSeq = seq(from = 0, to = 100, by = 10), 
                        correction = "rs",
                        sample_id = "subject",
                        family = gaussian(link = "log"),
                        image_id = "imageID", 
                        condition = "condition",
                        ncores = 1
                    )
    test.spatialFDAG <- summary(res$mdl) 

    dfG <- as.data.frame(test.spatialFDAG$s.table)[paste0('condition', "Group2", '(x)'),]
    dfG <- dplyr::rename(dfG, p.value = `p-value`)

    res <- spatialFDA::spatialInference(
                        spe, 
                        selection = c("B", "A"), 
                        fun = "Lcross", 
                        marks = "cellType",
                        rSeq = seq(from = 0, to = 100, by = 10), 
                        correction = "iso",
                        sample_id = "subject",
                        family = gaussian(link = "log"),
                        image_id = "imageID", 
                        condition = "condition",
                        ncores = 1
                    )
    test.spatialFDAL <- summary(res$mdl)

    dfL <- as.data.frame(test.spatialFDAL$s.table)[paste0('condition', "Group2", '(x)'),]
    dfL <- dplyr::rename(dfL, p.value = `p-value`)

    #### intensity MM ####
    
    dfSpe <- .speToDf(spe)
    intensitiesControl <- ConditionIntensities(df = dfSpe,
          	                                  marks = "cellType",
                                              conditionType = "Group1",
                                              cellType = "A",
                                              sample_id = "subject",
					                                    imageId = "imageID")
   
    intensitiesPert <- ConditionIntensities(df = dfSpe,
                                            marks = "cellType",
                                            conditionType = "Group2",
                                            cellType = "A",
					                                  sample_id = "subject",
                                            imageId = "imageID")
   
    dfIntensity <- rbind(intensitiesControl, intensitiesPert)

    lmeMod <- lmerTest::lmer(intensity ~ condition + (1|subject), data = dfIntensity)
    out1 <- lmerTest:::get_coefmat(lmeMod) |>
        as.data.frame()  
    out1$cellType <- "A"
    
    intensitiesControl <- ConditionIntensities(df = dfSpe,
          	                                  marks = "cellType",
                                              conditionType = "Group1",
                                              cellType = "B",
                                              sample_id = "subject",
					                                    imageId = "imageID")
   
    intensitiesPert <- ConditionIntensities(df = dfSpe,
                                            marks = "cellType",
                                            conditionType = "Group2",
                                            cellType = "B",
					                                  sample_id = "subject",
                                            imageId = "imageID")
   
    dfIntensity <- rbind(intensitiesControl, intensitiesPert)
    dfIntensity$condition <- factor(dfIntensity$condition, levels = c("Group1", "Group2"))

    lmeMod <- lmerTest::lmer(intensity ~ condition + (1|subject), data = dfIntensity)
    out2 <- lmerTest:::get_coefmat(lmeMod) |>
        as.data.frame() 
    out2$cellType <- "B"
    
    #bind the results from both models and only keep the rows with the lower p-value
    total <- rbind(out1["conditionGroup2", ], out2["conditionGroup2",])
    total <- total |> slice_min(.data[["Pr(>|t|)"]], n = 1, with_ties = FALSE)
    rownames(total) <- NULL

    #### smoppix ####
    df <- colData(spe) |> as.data.frame() |> cbind(spatialCoords(spe))

    hypDf <- buildHyperFrame(df,
                              coordVars = c("x", "y"),
                              imageVars = c("condition", "subject", "imageID"),
                            featureName = "cellType"
    )

    nnObj <- estPis(hypDf,
                    pis = c("nnPair"), null = "background", verbose = FALSE,
                    features = c("B", "A")
    )

    nnObj <- addWeightFunction(nnObj, lowestLevelVar = "imageID",
                                  pi = "nnPair")

    dfUniNN <- buildDataFrame(nnObj, gene = "B--A", pi = "nnPair")

    lmeMod <- lmerTest::lmer(pi - 0.5 ~ condition + (1 | subject),
                            data = dfUniNN, na.action = na.omit,
                            weights = weight, contrasts = list("condition" = "contr.treatment")
    )

    smoppix <- lmerTest:::get_coefmat(lmeMod) |>
        as.data.frame()

    #### spaceANOVVA ####
    library("fda.usc")
    library("gridExtra")
    df <- colData(spe)
    df$Group <- as.factor(df$condition)
    df$cellType <- as.factor(df$cellType)
    df$imageID <- as.factor(df$imageID)
    df$ID <- as.factor(df$subject)
    df$x <- spatialCoords(spe)[,1]
    df$y <- spatialCoords(spe)[,2]

    data <- df %>% as.data.frame() %>% dplyr::select(Group, cellType, imageID, ID, x, y)

    #coded according to the github repository of spaceANOVA

    Final_result = SpaceANOVA::All_in_one(data = data, fixed_r = seq(0, 100, by = 1), Summary_function = "g",  Hard_ths = 10, homogeneous = TRUE, interaction_adjustment = TRUE, perm = TRUE, nPerm = 20, cores = 2)

    out <- SpaceANOVA::p_extract(Final_result)
    spaceANOVAUni <- out[[1]]
    spaceANOVAMulti <- out [[2]]
    
    #### mxfda ####

    meta <- as_tibble(colData(spe)) |>
      dplyr::select(subject, imageID, condition, sample_id) |>
      unique()
    meta$subject <- as.factor(meta$subject)
    # response needs to be 0<=y<=1 for logistic regression
    meta$condition_binary <- ifelse(meta$condition == "Group1", 0, 1)

    spatial <- as_tibble(spatialCoords(spe))
    # assumes same ordering of rows between spatialCoords(spe) and colData(spe)
    spatial$cellType <- as.character(as_tibble(colData(spe))$cellType)
    spatial$imageID  <- factor(as_tibble(colData(spe))$imageID)

    mxFDAobject <- make_mxfda(metadata   = meta,
                              spatial   = spatial,
                              subject_key = "subject",
                              sample_key  = "imageID")

    mxFDAobject <- extract_summary_functions(mxFDAobject,
                                              extract_func   = bivariate,
                                              summary_func   = Kcross,
                                              r_vec          = seq(0, 100, by = 10),
                                              edge_correction = "iso",
                                              markvar        = "cellType",
                                              mark1          = "B",
                                              mark2          = "A") 

    mxFDAobject <- run_sofr(mxFDAobject,
                            model_name = "fit_sofr_condition",
                            formula    = condition_binary ~ 1,
                            family     = "binomial",
                            metric     = "bi k", r = "r", value = "fundiff",
                            optimizer  = "efs")

    mdl <- extract_model(mxFDAobject, 'bi k', type = 'sofr', model_name = 'fit_sofr_condition')

    mxfda_outFM <- summary(mdl, re.test = FALSE)

    mxFDAobject <- run_sofr(mxFDAobject,
                            model_name = "fit_sofr_condition",
                            formula    = condition_binary ~ 1 + s(subject, bs = "re"),
                            family     = "binomial",
                            metric     = "bi k", r = "r", value = "fundiff",
                            optimizer  = "efs")

    mdl <- extract_model(mxFDAobject, 'bi k', type = 'sofr', model_name = 'fit_sofr_condition')

    mxfda_outMM <- summary(mdl, re.test = FALSE)

    res = c(spicyRLM = test.spicyLM$p.value[1,"conditionGroup2"],
	          spicyRMM = test.spicyMM$p.value[1,"conditionGroup2"], 
            spatialFDAL = dfL$p.value,
            spatialFDAG = dfG$p.value,
            intensityMM = total[["Pr(>|t|)"]],
            smoppix = smoppix[[paste0('conditionGroup2'), "Pr(>|t|)"]],
            spaceANOVAUni = spaceANOVAUni[["B","A"]],
            spaceANOVAMulti = spaceANOVAMulti[["B","A"]],
            mxfdaFM = as.data.frame(mxfda_outFM$s.table)["s(xmat.tmat):L.xmat","p-value"],
            mxfdaMM = as.data.frame(mxfda_outMM$s.table)["s(xmat.tmat):L.xmat","p-value"]
            )
    return(res)
  }
  
  
  res <- future.apply::future_lapply(as.list(seq_len(nSim)+seed), sim, counts = counts, nPatients = nPatients, nIm = nIm, window = window, lambda = lambda, future.packages = c(
    "tidyr",
    "dplyr",
    "spicyR",
    "spatstat",
    "SpatialExperiment",
    "spatialFDA",
    "smoppix",
    "mxfda",
    "SpaceANOVA",
    "lmerTest"
  ))
  
  res <- do.call('rbind', res)
  
  res <- cbind(data.frame(res), simulation = "Signal", lambda = lam)
  
  resultsTW <- rbind(resultsTW,res)
  
}


s2 <- Sys.time()
s2-s1


## NO SIGNAL

s1 <- Sys.time()

set.seed(seed)
sim  <- function(i, counts, nPatients, nIm, window){
  set.seed(i)
  print(paste0("no signal iteration ",i-seed))
  x <- c()
  y <- c()
  cellType <- c()
  imageID <- c()
  for (p in 1:nPatients) {
    for (j in 1:nIm) {
      sCount1 <- sample(counts,1)
      sCount2 <- sample(counts,1)
      a <- rpoispp(sCount1/1000^2, win = window)
      # aDens <- density(a, adjust = 1000)
      # aDens <- aDens#/sum(aDens)
      b <- rpoispp(sCount2/1000^2, win = window)
      x <- c(x, a$x, b$x)
      y <- c(y, a$y, b$y)
      cellType <- c(cellType, rep("A", a$n), rep("B", b$n))
      imageID <- c(imageID, rep(paste(p,j,sep = "_"), a$n+b$n))
    }
  }

  imageID <- factor(imageID)
  
  cellExp <- data.frame(
    x = x,
    y = y,
    cellType = factor(cellType),
    imageID = imageID
  )
  
  phenoData <- data.frame(imageID = unique(imageID),
                            condition = rep(c("Group1", "Group2"), each = nIm*nPatients/2),
                            subject = rep(1:nPatients, each = nIm))
    
    colDF <- cellExp %>% left_join(phenoData, by = "imageID")

    spe <- SpatialExperiment(
      colData = colDF,
      spatialCoordsNames = c("x", "y"),
      sample_id = as.character(colDF$subject)
    )

    colData(spe)[["condition"]] <- factor(colData(spe)[["condition"]])
    #relevel to have ctrl as the reference category
    colData(spe)[["condition"]] <- relevel(colData(spe)[["condition"]],
    "Group1")
    
    test.spicyLM <- spicyR::spicy(spe,
                  condition = "condition",
                  from = "B",
                  to = "A",
                  fast = TRUE, Rs = seq(from = 10, to = 100, by = 10),
                  edgeCorrect = FALSE,
                  cores = 1,
		  weights = FALSE,
                  verbose = FALSE)
    
    test.spicyMM <- spicyR::spicy(spe,
                  condition = "condition",
                  from = "B",
                  to = "A",
		  subject = "subject",
                  fast = TRUE, Rs = seq(from = 10, to = 100, by = 10),
                  edgeCorrect = FALSE,
                  cores = 1,
		  weights = FALSE,
                  verbose = FALSE)

    res <- spatialFDA::spatialInference(
                        spe, 
                        selection = c("B", "A"), 
                        fun = "Gcross", 
                        marks = "cellType",
                        rSeq = seq(from = 0, to = 100, by = 10), 
                        correction = "rs",
                        sample_id = "subject",
                        family = gaussian(link = "log"),
                        image_id = "imageID", 
                        condition = "condition",
                        ncores = 1
                    )
    test.spatialFDAG <- summary(res$mdl) 

    dfG <- as.data.frame(test.spatialFDAG$s.table)[paste0('condition', "Group2", '(x)'),]
    dfG <- dplyr::rename(dfG, p.value = `p-value`)

    res <- spatialFDA::spatialInference(
                        spe, 
                        selection = c("B", "A"), 
                        fun = "Lcross", 
                        marks = "cellType",
                        rSeq = seq(from = 0, to = 100, by = 10), 
                        correction = "iso",
                        sample_id = "subject",
                        family = gaussian(link = "log"),
                        image_id = "imageID", 
                        condition = "condition",
                        ncores = 1
                    )
    test.spatialFDAL <- summary(res$mdl)

    dfL <- as.data.frame(test.spatialFDAL$s.table)[paste0('condition', "Group2", '(x)'),]
    dfL <- dplyr::rename(dfL, p.value = `p-value`)
  
    dfSpe <- .speToDf(spe)
    intensitiesControl <- ConditionIntensities(df = dfSpe,
          	                                  marks = "cellType",
                                              conditionType = "Group1",
                                              cellType = "A",
                                              sample_id = "subject",
					                                    imageId = "imageID")
   
    intensitiesPert <- ConditionIntensities(df = dfSpe,
                                            marks = "cellType",
                                            conditionType = "Group2",
                                            cellType = "A",
					                                  sample_id = "subject",
                                            imageId = "imageID")
   
    dfIntensity <- rbind(intensitiesControl, intensitiesPert)

    lmeMod <- lmerTest::lmer(intensity ~ condition + (1|subject), data = dfIntensity)
    out1 <- lmerTest:::get_coefmat(lmeMod) |>
        as.data.frame()  
    out1$cellType <- "A"
    
    intensitiesControl <- ConditionIntensities(df = dfSpe,
          	                                  marks = "cellType",
                                              conditionType = "Group1",
                                              cellType = "B",
                                              sample_id = "subject",
					                                    imageId = "imageID")
   
    intensitiesPert <- ConditionIntensities(df = dfSpe,
                                            marks = "cellType",
                                            conditionType = "Group2",
                                            cellType = "B",
					                                  sample_id = "subject",
                                            imageId = "imageID")
   
    dfIntensity <- rbind(intensitiesControl, intensitiesPert)
    dfIntensity$condition <- factor(dfIntensity$condition, levels = c("Group1", "Group2"))

    lmeMod <- lmerTest::lmer(intensity ~ condition + (1|subject), data = dfIntensity)
    out2 <- lmerTest:::get_coefmat(lmeMod) |>
        as.data.frame() 
    out2$cellType <- "B"
    
    #bind the results from both models and only keep the rows with the lower p-value
    total <- rbind(out1["conditionGroup2", ], out2["conditionGroup2",])
    total <- total |> slice_min(.data[["Pr(>|t|)"]], n = 1, with_ties = FALSE)
    rownames(total) <- NULL

    #### smoppix ####
    df <- colData(spe) |> as.data.frame() |> cbind(spatialCoords(spe))

    hypDf <- buildHyperFrame(df,
                              coordVars = c("x", "y"),
                              imageVars = c("condition", "subject", "imageID"),
                            featureName = "cellType"
    )

    nnObj <- estPis(hypDf,
                    pis = c("nnPair"), null = "background", verbose = FALSE
    )

    nnObj <- addWeightFunction(nnObj, lowestLevelVar = "imageID",
                                  pi = "nnPair")

    dfUniNN <- buildDataFrame(nnObj, gene = "B--A", pi = "nnPair")

    lmeMod <- lmerTest::lmer(pi - 0.5 ~ condition + (1 | subject),
                            data = dfUniNN, na.action = na.omit,
                            weights = weight, contrasts = list("condition" = "contr.sum")
    )

    smoppix <- lmerTest:::get_coefmat(lmeMod) |>
        as.data.frame()

    #### spaceANOVVA ####
    df <- colData(spe)
    df$Group <- as.factor(df$condition)
    df$cellType <- as.factor(df$cellType)
    df$imageID <- as.factor(df$imageID)
    df$ID <- as.factor(df$subject)
    df$x <- spatialCoords(spe)[,1]
    df$y <- spatialCoords(spe)[,2]

    data <- df %>% as.data.frame() %>% dplyr::select(Group, cellType, imageID, ID, x, y)

    #coded according to the github repository of spaceANOVA

    Final_result = SpaceANOVA::All_in_one(data = data, fixed_r = seq(0, 100, by = 1), Summary_function = "g",  Hard_ths = 10, homogeneous = TRUE, interaction_adjustment = TRUE, perm = TRUE, nPerm = 20, cores = 2)

    out <- SpaceANOVA::p_extract(Final_result)
    spaceANOVAUni <- out[[1]]
    spaceANOVAMulti <- out [[2]]
    
    #### mxfda ####

    meta <- as_tibble(colData(spe)) |>
      dplyr::select(subject, imageID, condition, sample_id) |>
      unique()
    meta$subject <- as.factor(meta$subject)
    # response needs to be 0<=y<=1 for logistic regression
    meta$condition_binary <- ifelse(meta$condition == "Group1", 0, 1)

    spatial <- as_tibble(spatialCoords(spe))
    # assumes same ordering of rows between spatialCoords(spe) and colData(spe)
    spatial$cellType <- as.character(as_tibble(colData(spe))$cellType)
    spatial$imageID  <- factor(as_tibble(colData(spe))$imageID)

    mxFDAobject <- make_mxfda(metadata   = meta,
                              spatial   = spatial,
                              subject_key = "subject",
                              sample_key  = "imageID")

    mxFDAobject <- extract_summary_functions(mxFDAobject,
                                              extract_func   = bivariate,
                                              summary_func   = Kcross,
                                              r_vec          = seq(0, 100, by = 10),
                                              edge_correction = "iso",
                                              markvar        = "cellType",
                                              mark1          = "B",
                                              mark2          = "A") 

    mxFDAobject <- run_sofr(mxFDAobject,
                            model_name = "fit_sofr_condition",
                            formula    = condition_binary ~ 1,
                            family     = "binomial",
                            metric     = "bi k", r = "r", value = "fundiff",
                            optimizer  = "efs")

    mdl <- extract_model(mxFDAobject, 'bi k', type = 'sofr', model_name = 'fit_sofr_condition')

    mxfda_outFM <- summary(mdl, re.test = FALSE)

    mxFDAobject <- run_sofr(mxFDAobject,
                            model_name = "fit_sofr_condition",
                            formula    = condition_binary ~ 1 + s(subject, bs = "re"),
                            family     = "binomial",
                            metric     = "bi k", r = "r", value = "fundiff",
                            optimizer  = "efs")

    mdl <- extract_model(mxFDAobject, 'bi k', type = 'sofr', model_name = 'fit_sofr_condition')

    mxfda_outMM <- summary(mdl, re.test = FALSE)

    res = c(spicyRLM = test.spicyLM$p.value[1,"conditionGroup2"],
	          spicyRMM = test.spicyMM$p.value[1,"conditionGroup2"], 
            spatialFDAL = dfL$p.value,
            spatialFDAG = dfG$p.value,
            intensityMM = total[["Pr(>|t|)"]],
            smoppix = smoppix[[paste0('condition1'), "Pr(>|t|)"]],
            spaceANOVAUni = spaceANOVAUni[["B","A"]],
            spaceANOVAMulti = spaceANOVAMulti[["B","A"]],
            mxfdaFM = as.data.frame(mxfda_outFM$s.table)["s(xmat.tmat):L.xmat","p-value"],
            mxfdaMM = as.data.frame(mxfda_outMM$s.table)["s(xmat.tmat):L.xmat","p-value"]
            )
  return(res)
  
}

results <- future.apply::future_lapply(as.list(seq_len(nSim)+seed), sim, counts = counts, nPatients = nPatients, nIm = nIm, window = window, future.packages = c(
    "tidyr",
    "dplyr",
    "spicyR",
    "spatstat",
    "SpatialExperiment",
    "spatialFDA",
    "smoppix",
    "mxfda",
    "SpaceANOVA",
    "lmerTest"
  ))

results <- do.call('rbind', results)

s2 <- Sys.time()
s2-s1

resultsFW <- results

resultsFW <- cbind(data.frame(resultsFW), simulation = "noSignal", lambda = lam)



RESFW = resultsFW

resFW <- NULL
for(i in Rs){
  r <- RESFW
  r$lambda <- i
  resFW <- rbind(r, resFW)
}

results <- rbind(resFW,resultsTW)

saveRDS(results, file = snakemake@output[["rds"]])
