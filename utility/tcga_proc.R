

barcodes2short <- function(bar_codes, id_col = 3){
  # Extract field from  ***-***-***-*** barcode
  num_codes <- length(bar_codes)
  barSplits <- unlist(strsplit(bar_codes, "[-]"))[]
  num_splits <- length(barSplits)
  num_col <- floor(num_splits / num_codes)
  short_codes <- barSplits[seq(id_col, num_splits, num_col)]
  return(short_codes)
}


filtSampleType <- function(bar_codes, sampType = "Primary Solid Tumor"){
  stdefs <- getSampleTypeFromBarCode(bar_codes)
  I <- !is.na(match(stdefs, sampType))
  return(I)
}

filtDupSamples <- function(bar_codes, sampType = "DNA"){
  sampTypes <- c("DNA", "RNA")
  
  filtI <- rep(TRUE, length(bar_codes))
  
  samp_names <- barcodes2short(bar_codes = bar_codes)
  dupI <- duplicated(samp_names)
  dup_samp_names <- unique(samp_names[dupI])
  if (length(dup_samp_names) == 0){
    return(filtI)
  }
  
  if (!(sampType %in% sampTypes)){
    stop("sampType not recognized")
  }
  numDupSamps <- length(dup_samp_names)
  for (sampi in 1:numDupSamps){
    currI <- !is.na(match(samp_names, dup_samp_names[sampi]))
    curr_bar_codes <- bar_codes[currI]
    if (sampType == "DNA"){
      curr_top_code <- get_top_dna_dup(curr_bar_codes)
    } else if (sampType == "RNA"){
      curr_top_code <- get_top_rna_dup(curr_bar_codes)
    }
    curr_samps_ignore <- setdiff(curr_bar_codes, curr_top_code)
    currI <- is.na(match(bar_codes, curr_samps_ignore))
    filtI <- filtI & currI
    str(curr_bar_codes)
    str(curr_samps_ignore)
  }
  return(filtI)
}

filtDupRNASamples_interactive <- function(countDF, nonDataFields = NULL, countCutoff = 20){
  bar_codes <- colnames(countDF)[-nonDataFields]
  sampNames <- barcodes2short(bar_codes, id_col = 3)
  dupI <- duplicated(sampNames)
  if (sum(dupI) == 0){
    return(countDF) 
  }
  dupSampNames <- unique(sampNames[dupI])
  
  sampsIgnore <- c()
  for (dupSampi in 1:length(dupSampNames)){
    currSampI <- !is.na(match(sampNames, dupSampNames[dupSampi]))
    curr_barcodes <- bar_codes[currSampI]
    curr_numSamps <- length(curr_barcodes)
    currDF <- countDF[-nonDataFields,curr_barcodes]
    filtI <- rowSums(currDF >= countCutoff) > 0
    currDF <- currDF[filtI,]
    # str(currDF)
    countTotals <- colSums(currDF)
    countZeros <- colSums(currDF  == 0)
    writeLines("")
    writeLines("")
    writeLines("")
    writeLines("")
    print("library sizes:")
    print(countTotals)
    print("num zero reads:")
    print(countZeros)
    
    currDF <- log2(currDF+1)
    corrMat <- cor(currDF)
    print("correlation matrix:")
    print(corrMat)
    promptString <- sprintf("Select a barcode (index integer between 1 and %d): ", curr_numSamps)
    promptStrErr <- sprintf("ERROR: Integer between 1 and %d required. Try again: ", curr_numSamps)
    sampChoicei <- strtoi(readline(prompt = promptString))
    while (is.na(sampChoicei) || sampChoicei < 1 || sampChoicei > curr_numSamps){
      sampChoicei <- strtoi(readline(prompt = promptStrErr))
    }
    sampsIgnore <- c(sampsIgnore, setdiff(curr_barcodes, curr_barcodes[sampChoicei]))
  }
  countDF <- countDF[,-which(names(countDF) %in% sampsIgnore)]
  
  return(countDF)
}

getDupSampleList <- function(bar_codes){
  filtI <- rep(TRUE, length(bar_codes))
  
  samp_names <- barcodes2short(bar_codes = bar_codes, id_col = 3)
  dupI <- duplicated(samp_names)
  dup_samp_names <- unique(samp_names[dupI])
  
  return(dup_samp_names)
}

getSampleTypeTable <- function(){
  sampTypeDF <- data.frame(code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", 
                                    "11", "12", "13", "14", "15", "16", "20", "40", "50", "60", "61", "99"),
                           definition = c("Primary Solid Tumor", "Recurrent Solid Tumor", "Primary Blood Derived Cancer - Peripheral Blood", 
                                          "Recurrent Blood Derived Cancer - Bone Marrow", "Additional - New Primary", "Metastatic",
                                          "Additional Metastatic", "Human Tumor Original Cells", "Primary Blood Derived Cancer - Bone Marrow",
                                          "Blood Derived Normal", "Solid Tissue Normal", "Buccal Cell Normal", "EBV Immortalized Normal",
                                          "Bone Marrow Normal", "sample type 15", "sample type 16", "Control Analyte", 
                                          "Recurrent Blood Derived Cancer - Peripheral Blood", "Cell Lines", "Primary Xenograft Tissue",
                                          "Cell Line Derived Xenograft Tissue", "sample type 99"),
                           letterCode = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM", "15SH",
                                          "16SH", "CELLC", "TRB", "CELL", "XP", "XCL", "99SH"),
                           stringsAsFactors = FALSE)
  return(sampTypeDF)
}

# Private functions
getSampleTypeFromBarCode <- function(bar_codes){
  st.table <- getSampleTypeTable()
  trunc_st <- function(strEntry){return(substr(strEntry, start = 1, stop = 2))}
  sub_st_name <- function(strEntry){return(st.table$definition[st.table$code == strEntry])}
  
  
  st.barcode.loc <- 4
  sts <- barcodes2short(bar_codes, id_col = st.barcode.loc)
  
  sts <- sapply(sts, trunc_st, USE.NAMES = FALSE)
  stdef <- sapply(sts, sub_st_name, USE.NAMES = FALSE)
  
  return(stdef)
  
}

get_top_dna_dup <- function(dupBarCodes){
  # Get analyte and plate numbers
  analytes <- barcodes2short(bar_codes = dupBarCodes, id_col = 5)
  analytes <- sapply(analytes, function(x){substr(x, nchar(x), 1000)}, USE.NAMES = FALSE)
  
  
  # Filter analyte D
  DI <- analytes == "D"
  if (sum(DI) != 0){
    dupBarCodes <- dupBarCodes[DI]
  }
  if (length(dupBarCodes) == 1){
    return(dupBarCodes)
  }
  
  # Filter largest plate number
  plateNumbers <- barcodes2short(bar_codes = dupBarCodes, id_col = 6)
  plateSort <- sort(plateNumbers, decreasing = TRUE)
  topPlateI <- plateNumbers == plateSort[1]
  dupBarCodes <- dupBarCodes[topPlateI]
  if (length(dupBarCodes) == 1){
    return(dupBarCodes)  
  }
  
  # Filter lexicographical
  barcodeSort <- sort(dupBarCodes, decreasing = TRUE)
  topBarCodeI <- dupBarCodes == barcodeSort[1]
  dupBarCodes <- dupBarCodes[topBarCodeI]
  return(dupBarCodes)
}

get_top_rna_dup <- function(dupBarCodes){
  # Get analyte and plate numbers
  
  
  
  # Filter analyte D
  analytes <- barcodes2short(bar_codes = dupBarCodes, id_col = 5)
  analytes <- sapply(analytes, function(x){substr(x, nchar(x), 1000)}, USE.NAMES = FALSE)
  DI <- analytes == "H"
  if (sum(DI) != 0){
    dupBarCodes <- dupBarCodes[DI]
  }
  if (length(dupBarCodes) == 1){
    return(dupBarCodes)
  }
  
  analytes <- barcodes2short(bar_codes = dupBarCodes, id_col = 5)
  analytes <- sapply(analytes, function(x){substr(x, nchar(x), 1000)}, USE.NAMES = FALSE)
  RI <- analytes == "R"
  if (sum(RI) != 0){
    dupBarCodes <- dupBarCodes[RI]
  }
  if (length(dupBarCodes) == 1){
    return(dupBarCodes)
  }
  
  analytes <- barcodes2short(bar_codes = dupBarCodes, id_col = 5)
  analytes <- sapply(analytes, function(x){substr(x, nchar(x), 1000)}, USE.NAMES = FALSE)
  TI <- analytes == "T"
  if (sum(TI) != 0){
    dupBarCodes <- dupBarCodes[TI]
  }
  if (length(dupBarCodes) == 1){
    return(dupBarCodes)
  }
  
  # Filter largest plate number
  plateNumbers <- barcodes2short(bar_codes = dupBarCodes, id_col = 6)
  plateSort <- sort(plateNumbers, decreasing = TRUE)
  topPlateI <- plateNumbers == plateSort[1]
  dupBarCodes <- dupBarCodes[topPlateI]
  if (length(dupBarCodes) == 1){
    return(dupBarCodes)  
  }
  # Filter lexicographical
  barcodeSort <- sort(dupBarCodes, decreasing = TRUE)
  topBarCodeI <- dupBarCodes == barcodeSort[1]
  dupBarCodes <- dupBarCodes[topBarCodeI]
  return(dupBarCodes)
}