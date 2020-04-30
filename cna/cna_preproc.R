genImmPhe_path <- getwd()
source(file.path(genImmPhe_path, "utility/io_evd.R"), local = TRUE)
source(file.path(genImmPhe_path, "utility/tcga_proc.R"), local = TRUE)

tcga_cnaSeg_2_geneCN <- function(segFile, gene_names, build = NULL){
  delThresh <- -0.25
  ampThresh <- +0.25
  
  geneLocDF <- get_gene_locs(gene_names = gene_names, build = build)
  segDF <- get_firebrows_segDF(segFile = segFile, sampType = "Primary Solid Tumor")
  
  # Retain only overlapping segments per gene
  numGenes <- nrow(geneLocDF)
  overlapSegList <- list()
  for (geneI in 1:numGenes){
    chromI <- !is.na(match(segDF$Chromosome, geneLocDF$chromosome_name[geneI]))
    startI <- segDF$Start <= geneLocDF$end_position[geneI]
    endI   <- segDF$End >= geneLocDF$start_position[geneI]
    overlapI <- chromI & startI & endI
    overlapSegList[[geneLocDF$hgnc_symbol[geneI]]] <- segDF[overlapI,]
  }
  
  # Initialize copy number matrix
  sampNames <- unique(segDF$Sample)
  numSamps <- length(sampNames)
  geneCnMat <- matrix(0, nrow = numGenes, ncol = numSamps)
  rownames(geneCnMat) <- geneLocDF$hgnc_symbol
  colnames(geneCnMat) <- sampNames
  
  # Select copy number log ratio from overlapping segments per sample per gene
  sampNames <- unique(segDF$Sample)
  numSamps <- length(sampNames)
  for (geneI in 1:numGenes){
    currGeneName <- geneLocDF$hgnc_symbol[geneI]
    currGeneSegDF <- overlapSegList[[currGeneName]]
    
    for (sampI in 1:numSamps){
      currSampName <- sampNames[sampI]
      currSampI <- !is.na(match(currGeneSegDF$Sample, currSampName))
      mult_resolve_flag <- sum(currSampI) > 1
      currGeneSampSegDF <- currGeneSegDF[currSampI,]
      cn_logr <- max_abs(currGeneSampSegDF$Segment_Mean)
      if (!is.na(cn_logr) && cn_logr < ampThresh && cn_logr > delThresh){
        cn_logr <- 0
      }
      geneCnMat[currGeneName, currSampName] <- cn_logr
      # if (mult_resolve_flag){
      #   str(currGeneSampSegDF)
      #   str(cn_logr)
      # }
    }
  }
 geneCnDF <- data.frame(hgnc_symbol = rownames(geneCnMat), geneCnMat, check.names = FALSE, stringsAsFactors = FALSE)
 barcodes <- colnames(geneCnDF)[-1]
 uniSampI <- c(TRUE, filtDupSamples(bar_codes = barcodes, sampType = "DNA"))
 geneCnDF <- geneCnDF[,uniSampI]
 newbarcodes <- colnames(geneCnDF)[-1]
 colnames(geneCnDF)[-1] <- barcodes2short(newbarcodes, id_col = 3)


  return(geneCnDF)
  
}

get_samp_gene_logr <- function(overlapSegs, sampName, resFun = max_abs, zeroRange = c(-0.25, 0.25)){
  
}

max_abs <- function(vec){
  if (length(vec) == 0){
    warning("No overlapping segments detected. NA returned")
    return(NA)
  }
  vec_abs <- abs(vec)
  maxI <- vec_abs == max(vec_abs)
  largest <- unique(vec[maxI])
  if (length(largest) > 1){
    warning("Multiple overlapping segments have the same absolute logRatio with alternating signs. Negative value returned")
    largest <- min(largest)
  }
  return(largest)
}


get_firebrows_segDF <- function(segFile, sampType = "Primary Solid Tumor"){
  segDF <- load_tsv(inFile = segFile)
  I <- filtSampleType(bar_codes = segDF$Sample, sampType = sampType)
  segDF <- segDF[I,]
  
  segDF$Chromosome <- as.character(segDF$Chromosome)
  XI <- !is.na(match(segDF$Chromosome, "23"))
  YI <- !is.na(match(segDF$Chromosome, "24"))
  segDF$Chromosome[XI] <- "X"
  segDF$Chromosome[YI] <- "Y"
  
  return(segDF)
}

get_gene_locs <- function(gene_names = NULL, build = NULL){

  chromFilt <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                 "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
                 "21", "22", "X", "Y", "MT")
  
  require('biomaRt')
  
  if (is.null(build)){
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  } else if (build == "grch37"){
    grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    mart <- useDataset("hsapiens_gene_ensembl", grch37)
  } else {
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  }
  
  G_list <- getBM(filters= "hgnc_symbol", 
                  attributes= c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                  values = gene_names, 
                  mart= mart)
  
  I <- !is.na(match(G_list$chromosome_name, chromFilt))
  G_list <- G_list[I,]
  
  return(G_list)
}