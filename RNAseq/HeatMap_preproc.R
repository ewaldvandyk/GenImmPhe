
genHeatMap <- function(df, dataCols, rowNameField){
  df_data <- df[,dataCols]
  uniRowNames <- genUniSymbols(as.character(df[[rowNameField]]))
  rownames(df_data) <- uniRowNames
  heatmap(x = data.matrix(df_data), ylab = "genes", labRow = FALSE)
  return(df_data)
}


filtGenes <- function(df, dataCols, plotFlag = FALSE){
  minNonZeroFrac <- 0.00001  # Minimum fraction of samples with zero reads allowed for gene to be included
  minVar <- 0.5^2
  
  df_noNa <- removeNARows(df, dataCols)
  df_nonZ <- removeZeroRows(df_noNa, dataCols, minNonZeroFrac)
  df_log <- df_nonZ
  df_log[,dataCols] <- log2(df_nonZ[,dataCols]+1)
  df_filt <- filtDEgenes(df_log, dataCols, minVar, plotFlag)
}
  

removeNARows <- function(df, checkCol){
  defMatI <- !is.na(df[,checkCol])
  defI <- apply(defMatI, MARGIN = 1, all)
  df_clean <- df[defI,]
}

removeZeroRows <- function(df, checkCol, minNonZeroFrac){
  zeroLogI <- sign(df[,checkCol])
  nonZeroCount <- rowSums(zeroLogI)
  nonZeroFrac <- nonZeroCount / length(checkCol)
  useI <- nonZeroFrac >= minNonZeroFrac
  df_clean <- df[useI,]
}

filtDEgenes <- function(df, checkCol, minVar, plotFlag = FALSE){
  DE_var <- apply(df[,checkCol], 1, var)
  if (plotFlag){
    hist(sqrt(DE_var), breaks = 1000)  
    tic_seq <- seq(0, max(sqrt(DE_var))+1, by = 0.1)
    axis(1, at = tic_seq, labels = FALSE)
    abline(v = sqrt(minVar))
  }
  useI <- DE_var >= minVar
  df_clean <- df[useI,]
}

genUniSymbols <- function(symbs){
  numSyms <- length(symbs)
  uniSymbs <- symbs
  for (i in 1:numSyms){
    symI <- uniSymbs == symbs[i]
    numRep <- sum(symI)
    if (numRep <= 1) next
    newSyms <- rep_len("", length.out = numRep)
    for (j in 1:numRep){
      newSyms[j] <- sprintf("%s_%d", symbs[i], j-1)
    }
    uniSymbs[symI] <- newSyms
  }
  return(uniSymbs)
}

ad_hoc_anni_csv2df <- function(inFile){
  df <- read.csv(file = inFile)
  df_grps <- ad_hoc_anni_sampGroups(df)
  grp_names <- colnames((df_grps))
  num_grps <- length(grp_names)
  for (grpI in 1:num_grps){
    df[[grp_names[grpI]]] <- apply(df[df_grps[[grp_names[grpI]]]], MARGIN = 1, FUN = mean)
  }
  dataList <- list(df, df_grps)
  
  return(dataList)
}

ad_hoc_anni_sampGroups <- function(df){
  numRepsPerSamp <- 3
  group_names <- c("X1937_CIP2A_KD",
                "X1937_Control",
                "X38_CIP2A_KD",
                "X38_Control",
                "X231_CIP2A_KD",
                "X231_Control",
                "X436_CIP2A_KD",
                "X436_Control",
                "X468_CIP2A_KD",
                "X468_Control")
  
  group_colnames <- group_names
  num_groups <- length(group_names)
  for (i in 1:num_groups){
    group_colnames[i] <- sprintf("%s_Mean", group_names[i])
  }
  sampNames <- colnames(df)[-(1:6)]
  df_sampGroups <- data.frame(matrix("NA", nrow = numRepsPerSamp, ncol = length(sampNames)/3), stringsAsFactors = FALSE)
  colnames(df_sampGroups) <- group_colnames
  for (i in 1:num_groups){
    for (j in 1:numRepsPerSamp){
      df_sampGroups[[group_colnames[i]]][j] <- sprintf("%s_%d", group_names[i], j)
    }
  }
  return(df_sampGroups)
}

ad_hoc_getSampNames <- function(df){
  return(unlist(df, use.names = FALSE))
}

ad_hoc_getGroupNames <- function(df){
  return(colnames(df))
}

exe_HeatMap_preproc_anni <- function(){
  ##### Global paramters
  inFile <- "~/data/rna_seq/BCBL_Srikar/RNAseq_BCBL_Srikar_NormalizedData.csv"
  outFile <- "~/devel/R/RNAseq/csv/RNAseq_BCBL_Srikar_NormalizedData_clean.csv"
  outFileLog <- "~/devel/R/RNAseq/csv/RNAseq_BCBL_Srikar_NormalizedData_clean_log2.csv"
  minNonZeroFrac <- 0.00001  # Minimum fraction of samples with zero reads allowed for gene to be included
  rowNameField <- "ID"
  
  ##### Load and add replicate means
  df_list <- ad_hoc_anni_csv2df(inFile)
  df <- df_list[[1]]
  df_grp <- df_list[[2]]
  ##### Get column names from orignal vs. mean fields
  dataCols <- ad_hoc_getSampNames(df_grp)
  grpCols <- ad_hoc_getGroupNames(df_grp)
  allCols <- colnames(df)[-(1:6)]
  
  ##### Remove NAs and zero only rows from df based on original data
  df_noNa <- removeNARows(df, dataCols)
  df_nonZ <- removeZeroRows(df_noNa, dataCols, minNonZeroFrac)
  df_log <- df_nonZ
  df_log[,allCols] <- log2(df_nonZ[,allCols]+1)

  ##### Write cleaned dataset to csv
  # write.csv(df_nonZ, file = outFile, row.names = FALSE)
  # write.csv(df_log, file = outFileLog, row.names = FALSE)
  # 
  ##### Generate heatmap for all data
  # minVar <- 0.5^2
  # plotVarFlag = TRUE
  # df_DE <- filtDEgenes(df = df_log, minVar = minVar, checkCol = dataCols, plotFlag = plotVarFlag)
  # genHeatMap(df = df_DE, dataCols = dataCols, rowNameField = rowNameField)
  # 
  # ##### Generate heatmap for mean data
  # minVar <- 0.5^2
  # plotVarFlag = TRUE
  # df_DE <- filtDEgenes(df = df_log, minVar = minVar, checkCol = grpCols, plotFlag = plotVarFlag)
  # genHeatMap(df = df_DE, dataCols = grpCols, rowNameField = rowNameField)

  ##### Generate heatmap for cell-line 1937
  minNonZeroFrac <- 1.0
  minVar <- 0.6^2
  plotVarFlag = TRUE
  clCols <- c("X1937_CIP2A_KD_1", "X1937_CIP2A_KD_2", "X1937_CIP2A_KD_3",
              "X1937_Control_1", "X1937_Control_2", "X1937_Control_3")
  
  df_log <- removeZeroRows(df_log, clCols, minNonZeroFrac)
  df_DE <- filtDEgenes(df = df_log, minVar = minVar, checkCol = clCols, plotFlag = plotVarFlag)

  genHeatMap(df = df_DE, dataCols = clCols, rowNameField = rowNameField)

  return(df_log)
}


