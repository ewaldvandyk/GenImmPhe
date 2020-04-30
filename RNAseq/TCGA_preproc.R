genImmPhe_path <- getwd()
source(file.path(genImmPhe_path, "utility/tcga_proc.R"), local = TRUE)

gdc_TCGA_2_rnaSeqraw <-function(gdc_inDir, sampFile) {
  
  gdc_inDir <- tools::file_path_as_absolute(gdc_inDir)
  sampFile <- tools::file_path_as_absolute(sampFile)
  
  ####### Use only primary tumors:
  # sampsDF <- sampFile2filtDF(sampFile = sampFile,
  #                            fieldName = "Sample.Type",
  #                            valueKeepList = c("Primary Tumor"))
  sampsDF <- sampFile2filtDF(sampFile = sampFile,
                             fieldName = "Sample.Type",
                             valueKeepList = c())
  
  ####### Remove duplicate sample patients:
  # df_dupI <- check_uni_patient(sampsDF$Sample.ID)
  # do_uniPatientStats(base_dir = gdc_inDir, sampsDF, df_dupI, 3)
  # repSelectIs <- c(1, 2, 2, 2, 2, 2)
  # maskI <- get_uni_patient_mask(df_dupI, repSelectIs)
  # sampsDF <- sampsDF[maskI,]
  
  rawData <- extract_samps(gdc_inDir, sampsDF)
  return(rawData)
  # numSamps <- length(sampsDF$File.Name)
  # 
  # 
  # for (sampI in 1:numSamps){
  #   sampID <- sampsDF$Sample.ID[sampI]
  # }
  
  
  # fileI <- 0
  # for (subFolder in dir(gdc_inDir)){
  #   currFolder <- file.path(gdc_inDir, subFolder)
  #   if (!dir.exists(currFolder)) next
  #   gz_file_name <- dir(currFolder, pattern = "*.gz")
  #   gz_file <- file.path(currFolder, gz_file_name)
  #   if (length(gz_file) == 0) next
  #   if (!file.exists(gz_file)) next
  #   fileI <- fileI + 1
  #   fileParts <- unlist(strsplit(gz_file_name, "[.]"))
  #   print(gz_file)
  #   print(fileParts)
  #   print(fileI)
  # }
}

addGeneNames2countData <- function(countDF, g_list = NULL, ambEnsgDF = NULL, ambHgncDF = NULL, ambEntrezDF = NULL){
  if (is.null(g_list)){
    g_list <- ensg_2_geneSymbols(countDF)
  }
  g_list <- rem_g_list_dups(g_list, ambEnsgDF, ambHgncDF, ambEntrezDF)
  
  colnames(countDF)[which(names(countDF) == "ensg_id")] <- "ensembl_gene_id"
  countDF$ensembl_gene_id <- ensgsV_2_ensgs(countDF$ensembl_gene_id)
  
  countDF <- merge(x = g_list, y = countDF, by = "ensembl_gene_id", all.y = TRUE)
  
  I = !is.na(match(countDF$hgnc_symbol, ""))
  countDF$hgnc_symbol[I] <- NA
  
  
  
  return(countDF)
}

ensg_2_geneSymbols <- function(countDF){
  ensgV <- countDF$ensg_id
  ensg <- ensgsV_2_ensgs(ensgV)
  
  library('biomaRt')
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description"), values=ensg,mart= mart)
  return(G_list)
}

rem_g_list_dups <- function(g_list, ambEnsgDF = NULL, ambHgncDF = NULL, ambEntrezDF = NULL){
  if (is.null(ambEnsgDF)){
    I <- duplicated(g_list$ensembl_gene_id, fromLast = FALSE) | duplicated(g_list$ensembl_gene_id, fromLast = TRUE)
    g_list <- g_list[!I,]
  } else {
    g_list <- g_list[!duplicated(g_list$ensembl_gene_id),]
    numAmb <- nrow(ambEnsgDF)
    for (i in 1:numAmb){
      I <- g_list$ensembl_gene_id == ambEnsgDF$ensembl_gene_id[i]
      
      g_list$hgnc_symbol[I] <- ambEnsgDF$hgnc_symbol[i]
      g_list$entrezgene[I] <- ambEnsgDF$entrezgene[i]
      g_list$description[I] <- ambEnsgDF$description[i]
    }
  }
  
  if (is.null(ambHgncDF)){
    I <- duplicated(g_list$hgnc_symbol, fromLast = FALSE) | duplicated(g_list$hgnc_symbol, fromLast = TRUE)
    g_list <- g_list[!I,]
  } else {
    g_list <- g_list[!duplicated(g_list$hgnc_symbol),]
    numAmb <- nrow(ambHgncDF)
    for (i in 1:numAmb){
      I <- g_list$hgnc_symbol == ambHgncDF$hgnc_symbol[i]
      
      g_list$ensembl_gene_id[I] <- ambHgncDF$ensembl_gene_id[i]
      g_list$entrezgene[I] <- ambHgncDF$entrezgene[i]
      g_list$description[I] <- ambHgncDF$description[i]
    }
  }
  
  I <- (duplicated(g_list$entrezgene, fromLast = FALSE) | duplicated(g_list$entrezgene, fromLast = TRUE)) & !is.na(g_list$entrezgene)
  if (is.null(ambEntrezDF)){
    g_list <- g_list[!I,]
  } else {
    is <- which(I)
    
    for (i in is){
      ensg <- g_list$ensembl_gene_id[i]
      I <- which(ambEntrezDF$ensembl_gene_id == ensg)
      if (length(I) == 0){
        I1 <- g_list$ensembl_gene_id == ensg
        entrezTest <- g_list$entrezgene[I1]
        I2 <- !is.na(match(g_list$entrezgene, entrezTest))
        print(g_list[I2,])
        
      }
      g_list$hgnc_symbol[i] <- ambEntrezDF$hgnc_symbol[I]
      g_list$entrezgene[i] <- ambEntrezDF$entrezgene[I]
      g_list$description[i] <- ambEntrezDF$description[I]
    }
    
  }
  
  
  return(g_list)
  
}

ensgsV_2_ensgs <- function(ensgV){
  ensgList <-strsplit(x = ensgV, split = "[.]")
  numEnsg <- length(ensgList)
  ensg <- rep(x = "NONE", numEnsg)
  for (i in 1:numEnsg){
    ensg[i] <- ensgList[[i]][1]
  }
  return(ensg)
}

sampFile2filtDF <- function(sampFile, fieldName, valueKeepList){
  df_all <- read.csv(file = sampFile, 
                     sep = "\t", 
                     as.is = c("File.ID", "File.Name", "Case.ID", "Sample.ID"))
  
  if (length(valueKeepList) == 0){
    df_filt <- df_all
  }else{
    df_filt <- df_all[df_all[[fieldName]] %in% valueKeepList,]  
  }
  
  
}

check_uni_patient <- function(bar_codes){
  id_col <- 3
  num_col <- 4
  
  num_codes <- length(bar_codes)
  barSplits <- unlist(strsplit(bar_codes, "[-]"))[]
  num_splits <- length(barSplits)
  patientIDs <- barSplits[seq(id_col, num_splits, num_col)]
  dupI <- duplicated(patientIDs)
  dupIDs <- unique(patientIDs[dupI])
  numDupIDs <- length(dupIDs)
  df_dupI <- data.frame(matrix(FALSE, nrow = num_codes, ncol = numDupIDs))
  colnames(df_dupI) <- dupIDs
  if (numDupIDs < 1){
    print("No duplicated patient samples found")
    return(0)
  }
  for (dupI in 1:numDupIDs){
    df_dupI[[dupI]] <- patientIDs %in% dupIDs[dupI]
  }
  # x <- data.matrix(df_
  return(df_dupI)
}

load_count_file <- function(dirname, sampsDF, index){
    countFile <- file.path(dirname, sampsDF$File.ID[index], sampsDF$File.Name[index])
    countDF <- read.csv(file = countFile, sep = '\t', header = FALSE, col.names = c("ensg_id", "count"), as.is = "ensg_id")
    ensg_i <- grep(countDF$ensg_id, pattern = "ENSG")
    countDF <- countDF[ensg_i,]
    return(countDF)
}
  
do_uniPatientStats <- function(base_dir, sampsDF, df_dupI, dupNum){
  repSampsDF <- sampsDF[df_dupI[,dupNum],]
  numRep <- nrow(repSampsDF)
  
  
  if (numRep < 2){
    print("No repititions to analyze")
    return(0)
  }
  
  # countTotals <- rep(0, numRep)
  # countZeros <- rep(0, numRep)
  
  countDF <- load_count_file(dirname = base_dir, sampsDF = repSampsDF, index = 1)
  colnames(countDF)[2] <- "count_1"
  
  for (i in 2:numRep){
    currCountDF <- load_count_file(dirname = base_dir, sampsDF = repSampsDF, index = i)
    colnames(currCountDF)[2] <- sprintf("count_%d", i)
    countDF <- merge(x = countDF, y = currCountDF, by = "ensg_id")
  }
  countTotals <- colSums(countDF[,-(1)])
  countZeros <- colSums(countDF[,-(1)] == 0)
  str(countZeros)
  logCountDF <- countDF
  logCountDF[,-(1)] <- log2(countDF[,-(1)]+1)
  
  corrMat <- cor(logCountDF[,-(1)])
  heatmap(corrMat,col = heat.colors(n = 100))
  # plot(logCountDF$count_2, logCountDF$count_3)
  print(corrMat)
  print(countTotals)
  print(countZeros)
}

get_uni_patient_mask <- function(df_dupI, repSelectIs){
  # repSelectIs <- c(1, 2, 2, 2, 2, 2)
  numCols <- ncol(df_dupI)
  numRows <- nrow(df_dupI)
  
  rem_i <- c()
  for (i in 1:numCols){
    curr_i <- which(x = df_dupI[,i])
    curr_rem_i <- curr_i[-(repSelectIs[i])]
    rem_i <- c(rem_i, curr_rem_i)
  }
  
  maskI <- rep(TRUE, numRows)
  maskI[rem_i] <- FALSE
  return(maskI)
}

extract_samps <- function(baseDir, sampsDF){
  numSamps <- nrow(sampsDF)
  str(numSamps)
  currSampName <- sampsDF$Sample.ID[1]
  countDF <- load_count_file(dirname = baseDir, sampsDF = sampsDF, index = 1)
  colnames(countDF)[2] <- currSampName
  str(countDF)
  
  for (i in 2:numSamps){
    if (i %% 10 == 0){
      print(i)
    }
    currSampName <- sampsDF$Sample.ID[i]
    currCountDF <- load_count_file(dirname = baseDir, sampsDF = sampsDF, index = i)
    colnames(currCountDF)[2] <- currSampName
    countDF <- merge(x = countDF, y = currCountDF, by = "ensg_id")
  }
  return(countDF)
}



