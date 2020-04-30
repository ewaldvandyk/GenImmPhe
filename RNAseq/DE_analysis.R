library(edgeR)

DE_countDF2DGEList <- function(countDF, rowNameCol = 1, colIgnore = 1:4, group = NULL){
  rowNames <- countDF[[rowNameCol]]
  geneInfo <- countDF[colIgnore]
  countDF <- countDF[-colIgnore]
  rownames(countDF) <- rowNames
  
  edgeRList <- DGEList(counts = countDF, genes = geneInfo, group = group)
  
  
  return(edgeRList)
}

DE_filterLowExpressedGenes <- function(edgeRList, countThresh = 6, minSampsPerGroup = 2){
  
  cpmCounts <- cpm(edgeRList)
  
  minLibSizeSampI <- which(edgeRList$samples$lib.size == min(edgeRList$samples$lib.size))[[1]]
  minCountProfile <- edgeRList$counts[,minLibSizeSampI]
  countThresh = min(minCountProfile[minCountProfile >= countThresh])
  minCountI = which(minCountProfile == countThresh)[[1]]
  cpmProfile <- cpmCounts[,minLibSizeSampI]
  cpmThresh <- cpmProfile[minCountI]
  
  keep <- rowSums(cpmCounts > cpmThresh) >= minSampsPerGroup
  edgeRList <- edgeRList[keep, , keep.lib.sizes=FALSE]
  
  return(edgeRList)
}

DE_TMM_normalize <- function(edgeRList){
  edgeRList <- calcNormFactors(edgeRList)
  # Possibly also correct for GC content
}

DE_init_covMat <- function(countDF, colIgnore, covIDs){
  sampNames <- names(countDF)[-colIgnore]
  numSamps <- length(sampNames)
  covDF <- c()
  for (covName in names(covIDs)){
    print(covName)
    covDF[[covName]] <- factor(rep(covIDs[[covName]][1], numSamps), levels = covIDs[[covName]])
  }
  covDF <- data.frame(covDF)
  rownames(covDF) <- sampNames
  return(covDF)
}

DE_design_2_bioRep <- function(design){
  require(purrr)
  
  binCodes <- t(matrix(map_dbl(.x = 0:(ncol(design)-1), .f = function(x) 2^x)))
  bioRepGroups <- rowSums(design*(matrix(rep(1, nrow(design)))%*%binCodes))
  
  
  return(bioRepGroups)
  
}

  