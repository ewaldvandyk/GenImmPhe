load_tsv <- function(inFile, colClasses = NA){
  countDF <- read.table(file = inFile, header = TRUE, sep = "\t", quote = "",
                        check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE, colClasses = colClasses)
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