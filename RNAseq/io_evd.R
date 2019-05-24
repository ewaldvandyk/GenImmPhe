load_tsv <- function(inFile){
  countDF <- read.table(file = inFile, header = TRUE, sep = "\t", quote = "",
                        check.names = FALSE, stringsAsFactors = FALSE)
  return(countDF)
}

save_tsv <- function(countDF, outFile){
  write.table(x = countDF, file = outFile, quote = FALSE, sep = "\t", row.names = FALSE)
}
