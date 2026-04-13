### try and install packages ###
remotes::install_github("mjemons/spatialFDA@e70d8210bc46c36616c0ec596ac272631a9c0ec4")
remotes::install_github("sthawinke/smoppix@aa6894c2938bc1847d05eb6c820700f92a53490b")
devtools::install_github('julia-wrobel/mxfda@0d0d924229d05df16450c2600bb78fc832621c60')
devtools::install_github('sealx017/SpaceANOVA@3efb134b578519ed9c0492c4a12b72b80ece0f9a')
devtools::install_github("joshuaulrich/microbenchmark@5fb8294c6fd9b99db3fc8ae02d4744b0ea4a12e4")
devtools::install_github("moviedo5/fda.usc@bfd3b4861e0d75612be6dc13296152f247e23be7")

library(spicyR)
library(spatstat)
library(tidyverse)
library(plotROC)
library(SpatialExperiment)
library(spatialFDA)
library(smoppix)
library(SpaceANOVA)
library(mxfda)
library(microbenchmark)
library(fda.usc)
library(ggplot2)

## Runtime comparison is adapted from the Canete et al. simulation 
## already used in this analysis in `spicyRsim.R`
seed = 51773
set.seed(seed)
window <- owin(xrange = c(0, 1000),
               yrange = c(0, 1000))

nPatients <- 2 
nCores <- 1
lambda <- 10
counts <- seq(from = 20, to = 300, by = 10)
Rs <- seq(from = 10, to = 100, by = 10)

s1 <- Sys.time()

Images <- c(10, 20, 40, 80, 160, 320, 640)

out <- lapply(Images, function(nIm){
  resDf <- data.frame()
  set.seed(nIm)
  print(paste0("signal iteration ",nIm))
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
  
  ####################### spicyR #######################
  test.spicyLM <- microbenchmark(spicyR::spicy(spe,
                condition = "condition",
                from = "B",
                to = "A",
                fast = TRUE, Rs = Rs,
                edgeCorrect = FALSE,
    cores = 1,
                weights = FALSE,
                verbose = FALSE), times = 10, unit = "s")

  test.spicyMM <- microbenchmark(spicyR::spicy(spe,
                condition = "condition",
                from = "B",
                to = "A",
    subject = "subject",
                fast = TRUE, Rs = Rs,
                edgeCorrect = FALSE,
                cores = 1,
    weights = FALSE,
                verbose = FALSE), times = 10, unit = "s")
  
  ####################### spatialFDA #######################

  resG <- microbenchmark(spatialFDA::spatialInference(
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
                  ), times = 10, unit = "s")

  resL <- microbenchmark(spatialFDA::spatialInference(
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
                  ), times = 10, unit = "s")
  
  ####################### mxfda #######################
  meta <- as_tibble(colData(spe)) |> dplyr::select(c("imageID", "condition", "subject"))
  meta$sample_id <- as.factor(meta$subject)
  #response needs to be 0<=y<=1 for logistic regression
  meta$condition_binary <- ifelse(meta$condition == "Group1", 0, 1)
  meta <- meta |> unique()

  spatial <- as_tibble(spatialCoords(spe))
  #assumes same ordering of ids
  spatial$cellType <- as.character(as_tibble(colData(spe))$cellType)
  spatial$imageID <-factor(as_tibble(colData(spe))$imageID)
  mxFDAobject <- make_mxfda(metadata = meta,
                          spatial = spatial,
                          subject_key = "sample_id",
                          sample_key = "imageID")
  
  mxfda_wrapper <- function(mxFDAobject){
    mxFDAobject <- extract_summary_functions(mxFDAobject,
                                        extract_func = bivariate,
                                        summary_func = Kcross,
                                        r_vec = seq(0, 100, by = 1),
                                        edge_correction = "iso",
                                        markvar = "cellType",
                                        mark1 = "B",
                                        mark2 = "A")

    mxFDAobject <- run_sofr(mxFDAobject, 
                          model_name = "fit_sofr_condition", 
                          formula = condition_binary ~ 1, 
                          family = "binomial",
                          metric = "bi k", r = "r", value = "fundiff",
                          optimizer = "efs")
    return(mxFDAobject)
  }
  resMxfda <- microbenchmark(
    mxfda_wrapper(mxFDAobject), times = 10, unit = "s"
    )
  ####################### smoppix #######################

  df <- colData(spe) |> as.data.frame() |> cbind(spatialCoords(spe))
  hypDf <- buildHyperFrame(df,
                            coordVars = c("x", "y"),
                            imageVars = c("condition", "subject", "imageID"),
                          featureName = "cellType"
  )

  smoppix_wrapper <- function(hypDf){
    nnObj <- estPis(hypDf,
                    pis = c("nnPair"), null = "background", verbose = FALSE,
                    features = c("A", "B")
    )

    nnObj <- addWeightFunction(nnObj, lowestLevelVar = "imageID",
                                  pi = "nnPair")

    dfUniNN <- buildDataFrame(nnObj, gene = "B--A", pi = "nnPair")

    lmeMod <- lmerTest::lmer(pi - 0.5 ~ condition + (1 | subject),
                            data = dfUniNN, na.action = na.omit,
                            weights = weight, contrasts = list("condition" = "contr.sum")
    )
    return(lmeMod)
  }
  resSmoppix <- microbenchmark(
      smoppix_wrapper(hypDf), times = 10, unit = "s"
  )

  ####################### SpaceANOVA #######################

  colData(spe)$Group <- as.factor(colData(spe)$condition)
  colData(spe)$cellType <- as.factor(colData(spe)$cellType)
  colData(spe)$ID <- as.factor(colData(spe)$subject)
  colData(spe)$x <- spatialCoords(spe)[,1]
  colData(spe)$y <- spatialCoords(spe)[,2]
  data <- colData(spe) %>% as.data.frame() %>% dplyr::select(Group, cellType, imageID, ID, x, y)
  
  spaceAnova = microbenchmark(SpaceANOVA::All_in_one(data = data, fixed_r = seq(0, 100, by = 1), Summary_function = "g",  Hard_ths = 10, homogeneous = TRUE, interaction_adjustment = TRUE, perm = TRUE, nPerm = 20, cores = 1), times = 10, unit = "s")


  ####################### aggregate #######################
  resDf <- data.frame(spicyRLM = sum(summary(test.spicyLM)$mean),
                      spicyRMM = sum(summary(test.spicyMM)$mean),
                      spatialFDAG = sum(summary(resG)$mean),
                      spatialFDAL = sum(summary(resL)$mean),
                      SpaceANOVA = sum(summary(spaceAnova)$mean),
                      smoppix = sum(summary(resSmoppix)$mean),
                      mxfda = sum(summary(resMxfda)$mean),
                      nIm = nIm
                      )
  return(resDf)
}) |> bind_rows()

results <- tidyr::pivot_longer(out,cols = c("spicyRLM", "spicyRMM", "spatialFDAL", "spatialFDAG", "SpaceANOVA", "smoppix", "mxfda"), names_to = "method", values_to = "elapsed_time_seconds")

results <- results %>%
  mutate(method = case_when(
    method == "spicyRLM" ~ "spicyR.LM",
    method == "spicyRMM" ~ "spicyR.MM",
    method == "spatialFDAL" ~ "spatialFDA.L",
    method == "spatialFDAG" ~ "spatialFDA.G",
    method == "mxfda" ~ "mxfda",
    method == "SpaceANOVA" ~ "SpaceAnova",
    method == "smoppix" ~ "smoppix"
  )) 

saveRDS(results, snakemake@output[["rds"]])