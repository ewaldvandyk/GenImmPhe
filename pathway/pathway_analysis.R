source("./utility/io_evd.R")
get_genes_with_pathway_term <- function(pathFile = "~/data/pathways/CTD_genes_pathways.tsv", 
                                        terms = c("AKT", "PI3K", "PIK3CA", "MTOR")){
  pathDF <- load_tsv(inFile = pathFile)
  numEntries <- nrow(pathDF)
  numTerms <- length(terms)
  I <- rep(FALSE, numEntries)
  for (i in 1:numTerms){
    I <- I | grepl(x = pathDF$PathwayName, pattern = terms[i])  
  }
  filtGenes <- unique(pathDF$GeneSymbol[I])

  return(filtGenes)
}