
muts_2_binaryMatrix <- function(mutDF, mut_filt){
  # filtMuts should be a data.frame with $hgnc_symbol and $mut_name fields
  numMuts <- nrow(mut_filt)
  mut_strs <- sprintf("%s_%s", mut_filt$hgnc_symbol, mut_filt$mut_name)
  
  sampNames <- colnames(mutDF)[-1]
  numSamps <- length(sampNames)
  
  mutBinMat <- matrix(data = 0, nrow = numMuts, ncol = numSamps)
  colnames(mutBinMat) <- sampNames
  rownames(mutBinMat) <- mut_strs
  
  for (muti in 1:numMuts){
    geneI <- !is.na(match(mutDF$hgnc_symbol, mut_filt$hgnc_symbol[muti]))
    for (sampi in 1:numSamps){
      currMuts <- mutDF[[sampNames[sampi]]]
      if (mut_filt$mut_name[muti] == "ANY"){
        nameI <- !is.na(currMuts) & currMuts != ""
      } else {
        nameI <- grepl(x = currMuts, pattern = mut_filt$mut_name[muti], fixed = TRUE)
      }
      mutI <- geneI & nameI
      mutBinMat[mut_strs[muti], sampNames[sampi]] <- sum(mutI)
    }
  }
  mutBinDF <- data.frame(mutation_names = mut_strs, mutBinMat, stringsAsFactors = FALSE)
  rownames(mutBinDF) <- NULL
  
  return(mutBinDF)
}