load_tsv <- function(inFile, header = TRUE, sep = "\t", quote = "",
                     check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE, colClasses = NA, ...){
  countDF <- NULL
  try(countDF <- read.table(file = inFile, header = header, sep = sep, quote = quote,
                        check.names = check.names, stringsAsFactors = stringsAsFactors, fill = fill, colClasses = colClasses, ...), silent = FALSE)
  return(countDF)
}

save_tsv <- function(countDF, outFile){
  write.table(x = countDF, file = outFile, quote = FALSE, sep = "\t", row.names = FALSE)
}


filt_rows_DF <- function(df, fieldName, values){
  fieldValues <- df[[fieldName]]
  I <- !is.na(match(fieldValues, values))
  df <- df[I,]
  return(df)
}

tsv2xlxs_recursive <- function(inDir, pattern = "*.tsv"){
  require(xlsx)
  inFiles <- dir(inDir, recursive = TRUE, pattern = pattern)
  baseFiles <- get_file_base(inFiles)
  for (i in seq_along(inFiles)){
    currIn <- file.path(inDir, inFiles[i])
    currOut <- file.path(inDir, paste0(baseFiles[i], ".xlsx"))
    inDF <- load_tsv(currIn)
    if (!is.null(inDF)){
      write.xlsx(x = inDF, file = currOut)  
    }
    print(currIn)
    print(currOut)
  }
}

get_file_base <- function(fileList){
  baseList <- fileList
  for (fileI in seq_along(fileList)){
    base <- tools::file_path_sans_ext(baseList[fileI])
    while(base != baseList[fileI] && base != ""){
      baseList[fileI] <- base
      base <- tools::file_path_sans_ext(baseList[fileI])
    }
  }
  return(baseList)
}
