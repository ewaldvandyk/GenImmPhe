library(edgeR)
library(limma)



voom_normalize <- function(countDF, colIgnore = 1:4, lib.size = NULL, hgnc_keep = NULL){
  # Compute effective library size
  if (is.null(lib.size)){
    lib.size <- raw_2_tmm_effLibSize(countDF = countDF, colIgnore = colIgnore)
  }
  
  # Remove genes with low expression
  I <- filt_low_expr_genes(countDF = countDF, colIgnore = colIgnore, lib.size = lib.size, hgnc_keep = hgnc_keep)
  countDF <- countDF[I,]
  
  # Voom normalize
  countMat <- data.matrix(countDF[, -colIgnore])
  voomStr <- voom(counts = countMat, lib.size = lib.size, plot = TRUE)
  countDF[,-colIgnore] <- voomStr$E
  
  return(countDF)
  
}

raw_2_tmm_effLibSize <- function(countDF = NULL, colIgnore = 1:4){
  
  dataMat <- data.matrix(frame = countDF[,-colIgnore])
  # str(dataMat)
  normFactors <- calcNormFactors(object = dataMat, method = "TMM")
  libSizes <- colSums(dataMat)
  lib.size <- normFactors * libSizes
  return(lib.size)
  
}

filt_low_expr_genes = function(countDF, colIgnore = 1:4, lib.size = NULL, hgnc_keep = NULL){
  numGenes <- nrow(countDF)
  if (is.null(hgnc_keep)){
    keepMask <- rep(FALSE, numGenes)
  } else {
    keepMask <- !is.na(match(countDF$hgnc_symbol, hgnc_keep))
  }
  
  dataMat <- data.matrix(frame = countDF[,-colIgnore])
  geneI <- filterByExpr(y = dataMat, lib.size = lib.size)
  geneI <- geneI | keepMask
  return(geneI)
}

###### Legacy ###### 
raw_2_cpm <- function(countDF, lib.size = NULL, colIgnore = 1:4){
  countMat <- data.matrix(countDF[,-colIgnore])
  if (is.null(lib.size)){
    lib.size <- colSums(countMat)
  }
  numSamps <- length(lib.size)
  norm_multiplier <- 1e6 / lib.size
  for (i in 1:numSamps){
    countMat[,i] <- countMat[,i] * norm_multiplier[i]
    
  }
  countDF[,-colIgnore] <- countMat
  return(countDF)
}

###### Legacy ###### 
cpm_2_log <- function(countDF, colIgnore = 1:4){
  countDF[,-colIgnore] <- log2(countDF[,-colIgnore]+0.5)
  return(countDF)
}





