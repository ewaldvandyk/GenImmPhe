source("./utility/io_evd.R")
source("./utility/tcga_proc.R")

rppa_tcga_firebrowse_2_df <- function(inFile, sampType = "Primary Solid Tumor"){
  prot_field_name <- "protein_names"
  
  rppa_df <- load_tsv(inFile)
  cols <- colnames(rppa_df)
  colnames(rppa_df)[1] <- prot_field_name
  barcodes <- cols[-1]
  primTumI <- filtSampleType(bar_codes = barcodes, sampType = "Primary Solid Tumor")
  barcodes <- barcodes[primTumI]
  cols_keep <- c(prot_field_name, barcodes)
  rppa_df <- rppa_df[,cols_keep]
  sampNames <- barcodes2short(bar_codes = barcodes, id_col = 3)
  colnames(rppa_df)[-1] <- sampNames
  
  return(rppa_df)
}

