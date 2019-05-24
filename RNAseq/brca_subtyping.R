library(genefu)
library(varhandle)
library(caret)
source("./RNAseq/TCGA_preproc.R")

train_brca_pam50 <- function(countDF, colIgnore = 1:4, ensg_col = "ensembl_gene_id", entrez_col = "entrezgene"){
  map_cluster2molSubtype <- function(clusterVec, subtypeMap){
    numSubtypes <- nrow(subtypeMap)
    vecLen <- length(clusterVec)
    molSubVec <- rep("", vecLen)
    for (i in 1:numSubtypes){
      molSubVec[clusterVec == subtypeMap$cluster[i]] = subtypeMap$molSubtype[i]
    }
    return(molSubVec)
  }
  
  genefuList <- create_genefu_params(countDF = countDF, colIgnore = colIgnore, ensg_col = ensg_col, entrez_col = entrez_col)
  
  data(pam50.robust)
  pam50_genes <- pam50.robust$centroids.map[, c("probe","EntrezGene.ID")]
  pam50_res <- intrinsic.cluster(data = genefuList$countMat, annot = genefuList$dannot, do.mapping = TRUE, 
                                 std = "robust", rescale.q = 0.05, number.cluster = 5, 
                                 method.cor = "spearman", method.centroids = "mean", intrinsicg = pam50_genes,
                                 verbose = TRUE)
  # Correlate known molecular subtypes with clusters
  molSubtNames <- dimnames(pam50.robust$centroids)[[2]]
  clusterNames <- dimnames(pam50_res$model$centroids)[[2]]
  numSubtypes <- length(molSubtNames)
  corrMat <- cor(x = pam50.robust$centroids, y = pam50_res$model$centroids)
  corrSort <- sort(corrMat, decreasing = TRUE)[1:numSubtypes]
  
  subtypeMap <- data.frame(molSubtype = rep("", numSubtypes), cluster = rep("", numSubtypes), stringsAsFactors = FALSE)
  for (i in 1:numSubtypes) {
    currCorrI <- which(x = corrMat == corrSort[i], arr.ind = TRUE)
    subtypeMap$molSubtype[currCorrI[1]] <- molSubtNames[currCorrI[1]]
    subtypeMap$cluster[currCorrI[1]] <- clusterNames[currCorrI[2]]
  }
  print(subtypeMap)
  heatmap(corrMat)
  
  # Replace default cluster names with molecular subtype names
  dimnames(pam50_res$model$centroids)[[2]] <- map_cluster2molSubtype(clusterVec = dimnames(pam50_res$model$centroids)[[2]], subtypeMap = subtypeMap)
  dimnames(pam50_res$subtype.proba)[[2]] <- map_cluster2molSubtype(clusterVec = dimnames(pam50_res$subtype.proba)[[2]], subtypeMap = subtypeMap)
  dimnames(pam50_res$cor)[[2]] <- map_cluster2molSubtype(clusterVec = dimnames(pam50_res$cor)[[2]], subtypeMap = subtypeMap)
  pam50_res$subtype[1:length(pam50_res$subtype)] <- map_cluster2molSubtype(clusterVec = pam50_res$subtype, subtypeMap = subtypeMap)
  
  #Esnure that correct entrez integers are in model.centroids.map
  entrezDF <- data.frame(probe = countDF$ensembl_gene_id, entrez = countDF$entrezgene, stringsAsFactors = FALSE)
  emptyDF <- data.frame(probe = pam50_res$model$centroids.map[,1], stringsAsFactors = FALSE)
  centroid.map.DF <- merge(x = emptyDF, y = entrezDF, by = "probe", all.x = FALSE, sort = FALSE)
  pam50_res$model$centroids.map[,3] <- centroid.map.DF$entrez
  
  return(pam50_res$model)
  
  
}

classify_brca_molecular_suptypes <- function(countDF, colIgnore = 1:4, ensg_col = "ensembl_gene_id", entrez_col = "entrezgene", sbt.model = NULL){
  # Make sure that counts are log normalized with voom in limma
  # All the columns minus colIgnore should be the sample columns used
  
  
  genefuList <- create_genefu_params(countDF = countDF, colIgnore = colIgnore, ensg_col = ensg_col, entrez_col = entrez_col)
  if (is.null(sbt.model)){
    pam50_pred <- molecular.subtyping(data = genefuList$countMat, annot = genefuList$dannot, sbt.model = "pam50", do.mapping = TRUE)    
  } else {
    pam50_pred <- intrinsic.cluster.predict(data = genefuList$countMat, annot = genefuList$dannot, sbt.model = sbt.model, do.mapping = TRUE, verbose = TRUE)
  }

  pam50StringDF <- data.frame(sampName = names(pam50_pred$subtype), pam50 = pam50_pred$subtype, stringsAsFactors = FALSE)
  pam50ProbDF <- data.frame(sampName = dimnames(pam50_pred$subtype.proba)[[1]], pam50_pred$subtype.proba, stringsAsFactors = FALSE)
  pam50DF <- merge(x = pam50StringDF, y = pam50ProbDF, by = "sampName")
  return(pam50DF)
}

create_genefu_params <- function(countDF, colIgnore = 1:4, ensg_col = "ensembl_gene_id", entrez_col = "entrezgene"){
  ensebl_str_vec <- countDF[,ensg_col]
  # entrez_str_vec <- sprintf("%d", countDF[,entrez_col])
  
  
  countMat <- t(as.matrix(countDF[,-(colIgnore)]))
  colnames(countMat) <- ensebl_str_vec
  
  # dannot <- data.frame(probe = ensebl_str_vec, EntrezGene.ID = entrez_str_vec)
  dannot <- data.frame(probe = ensebl_str_vec, EntrezGene.ID = countDF[,entrez_col])
  
  genefuList <- list(countMat = countMat, dannot = dannot)
  return(genefuList)
  
}

create_sbt_confMat <- function(df1, df2, barc1_field = "sampName", barc2_field = "sampName", 
                               cls1_field = "pam50", cls2_field = "pam50"){
  
  barc1 <- barcodes2short(bar_codes = df1[,barc1_field], id_col = 3)
  barc2 <- barcodes2short(bar_codes = df2[,barc2_field], id_col = 3)
  cls1 <- df1[,cls1_field]
  cls2 <- df2[,cls2_field]
  
  I1 <- !duplicated(barc1)
  I2 <- !duplicated(barc2)
  barc1 <- barc1[I1]
  barc2 <- barc2[I2]
  cls1 <- cls1[I1]
  cls2 <- cls2[I2]
  
  df1 <- data.frame(id = barc1, cls1 = factor(cls1), stringsAsFactors = FALSE)
  df2 <- data.frame(id = barc2, cls2 = factor(cls2), stringsAsFactors = FALSE)
  
  df <- merge(x = df1, y = df2, by = "id", all.x = FALSE, all.y = FALSE)
  
  confList <- confusionMatrix(data = df$cls1, reference = df$cls2)
  return(confList)
}

filt_sbt <- function(countDF, colIgnore = 1:4, pam50_pred, sbt = "Basal"){
  filtSamps <- pam50_pred$sampName[pam50_pred$pam50 == sbt]
  I <- !is.na(match(colnames(countDF), filtSamps))
  I[colIgnore] <- TRUE
  
  countDF <- countDF[,I]
  return(countDF)
}