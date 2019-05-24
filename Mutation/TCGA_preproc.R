source("./utility/io_evd.R")
source("./utility/tcga_proc.R")

download_maf <- function(tumor = "BRCA", downloadDir = NULL, 
                         pipelines = c("mutect2", "muse", "varscan2", "somaticsniper")){
  if (is.null(downloadDir)){
    print("Error: Please specify a directory to download to")
    return(NULL)
  }
  library(TCGAbiolinks)
  numPipelines <- length(pipelines)                         
  mafList <- list()
  for (calleri in 1:numPipelines){
    currMaf <- GDCquery_Maf(tumor = "BRCA", directory = downloadDir, save.csv = FALSE, pipelines = pipelines[calleri])
    if (calleri == 1){
      mafDF <- currMaf   
    } else {
      mafDF <- rbind(mafDF, currMaf)
    }
  }
  
  
  return(mafDF)
}

mafDF_2_mutation_table <- function(mafDF){
  samp_names <- barcodes2short(mafDF$Tumor_Sample_Barcode, id_col = 3)
  samp_names_uni <- unique(samp_names)
  numSamps <- length(samp_names_uni)
  gene_names <- unique(mafDF$Hugo_Symbol)
  numGenes <- length(gene_names)
  
  mutMatrix <- structure(.Data = rep("", numGenes*numSamps), .Dim = c(numGenes, numSamps))
  colnames(mutMatrix) <- samp_names_uni
  mutDF <- data.frame(hgnc_symbol = gene_names, mutMatrix, stringsAsFactors = FALSE)
  
  numMut <- nrow(mafDF)
  for (muti in 1:numMut){
    if (muti %% 1000 == 0){
      print(sprintf("%0.2f%%", muti / numMut * 100))
    }
    currHgncSymb <- mafDF$Hugo_Symbol[muti]
    currMutName <- mafDF$HGVSp_Short[muti]
    if (is.na(currMutName)){
      next()
    }
    currSampName <- samp_names[muti]
    storedMutName <- mutDF[mutDF$hgnc_symbol == currHgncSymb, currSampName]
    
    if (storedMutName == ""){
      mutDF[mutDF$hgnc_symbol == currHgncSymb, currSampName] <- currMutName
    }else{
      mutList <- strsplit(x = storedMutName, split = "\\|")[[1]]
      newMutList <- union(mutList, currMutName)
      newMut <- paste(newMutList, collapse = "|")

      mutDF[mutDF$hgnc_symbol == currHgncSymb, currSampName] <- newMut
    }
  }
  return(mutDF)
}

tcga_firebrowse_2_mutation_table <- function(baseDir){
  filesInBaseDir <- dir(baseDir)
  filesI <- grep(pattern = "TCGA", x = dir(baseDir))
  mut_files <- filesInBaseDir[filesI]
  samp_names <- barcodes2short(bar_codes = mut_files, id_col = 3)
  numSamps <- length(samp_names)
  
  gene_names <- get_all_mut_genes(baseDir, mut_files)
  numGenes <- length(gene_names)
  
  mutMatrix <- structure(.Data = rep("NONE", numGenes*numSamps), .Dim = c(numGenes, numSamps))
  colnames(mutMatrix) <- samp_names
  mutDF <- data.frame(hgnc_symbol = gene_names, mutMatrix, stringsAsFactors = FALSE)
  
  for (sampi in 1:numSamps){
    if (sampi %% 100 == 0) {
      print(sampi)
    }
    curr_mut_file <- file.path(baseDir, mut_files[sampi])
    curr_sampName <- samp_names[sampi]
    sampDF <- load_tsv(curr_mut_file)
    numMut <- nrow(sampDF)
    # print(curr_mut_file)
    for (muti in 1:numMut){
      currHgncSymb <- sampDF$Hugo_Symbol[muti]
      currMutName <- sampDF$Protein_Change[muti]
      # print(muti)
      # print(currMutName)
      if (currMutName == "" | is.na(currMutName)){
        currMutName = sampDF$Variant_Classification[muti]
      }
      storedMutName <- mutDF[mutDF$hgnc_symbol == currHgncSymb, curr_sampName]
      if (storedMutName == "NONE"){
        mutDF[mutDF$hgnc_symbol == currHgncSymb, curr_sampName] <- currMutName
      }else{
        mutDF[mutDF$hgnc_symbol == currHgncSymb, curr_sampName] <- sprintf("%s;%s", storedMutName, currMutName)
      }
    }
  }
  return(mutDF)
}

get_all_mut_genes <- function(baseDir, mut_files){
  numSamps <- length(mut_files)
  genes <- c()
  for (sampi in 1:numSamps){
    if (sampi %% 100 == 0) {
      print(sampi)
    }
    curr_mut_file <- file.path(baseDir, mut_files[sampi])
    sampDF <- load_tsv(curr_mut_file)
    genes <- union(genes, sampDF$Hugo_Symbol)
  }
  return(genes)
}