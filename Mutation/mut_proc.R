library(Biostrings)
##### None standard amino acid table:
# Unusual translations
#   U - Selenocysteine
#   O - Pyrrolysine
# Ambiguous translations
#   X - Could be any aminc acid
#   B - Could be N or D
#   Z - Could be Q or E
#   J - Coulbe be L or I

muts_2_binaryMatrix <- function(mutDF, metaFields = 1, hgncField = 1, mut_filt){
  # filtMuts should be a data.frame with $hgnc_symbol and $mut_name fields
  # Only keep relevant genes in mutDF
  hgncUsed <- unique(mut_filt$hgnc_symbol)
  hgncI <- !is.na(match(mutDF[[hgncField]], hgncUsed))
  mutDF <- mutDF[hgncI,]
  
  #Remove redundent fields in mut_filt
  numMuts <- nrow(mut_filt)
  mut_strs <- sprintf("%s_%s", mut_filt$hgnc_symbol, mut_filt$mut_name)
  dupI <- duplicated(mut_strs)
  mut_filt <- mut_filt[!dupI,]
  mut_strs <- mut_strs[!dupI]
  numMuts <- nrow(mut_filt)
  
  #Get sample names
  sampNames <- colnames(mutDF)[-metaFields]
  numSamps <- length(sampNames)
  
  #Initialize mutBinMat
  mutBinMat <- matrix(data = 0, nrow = numMuts, ncol = numSamps)
  colnames(mutBinMat) <- sampNames
  rownames(mutBinMat) <- mut_strs
  
  supported_muttype_df <- get_supported_protMutTypes()
  for (hgnc in hgncUsed){
    mutTypesUsed <- mut_filt$mut_name[mut_filt$hgnc_symbol == hgnc]
    exactTypes <- setdiff(mutTypesUsed, supported_muttype_df$mutation_types)
    generTypes <- setdiff(mutTypesUsed, exactTypes)
    for (samp in sampNames){
      currMuts <- mutDF[mutDF$hgnc_symbol == hgnc, samp]
      currMuts <- strsplit(x = currMuts, split = "|", fixed = TRUE)[[1]]
      exactMutI <- protExactMutI("", mutation_name_list = exactTypes)
      generMutI <- protGenerMutI("", supported_muttype_df = supported_muttype_df)
      for (mut in currMuts){
        exactMutI <- exactMutI | protExactMutI(mut, mutation_name_list = exactTypes)
        generMutI <- generMutI | protGenerMutI(mut, supported_muttype_df = supported_muttype_df)
      }
      for (mut in exactTypes){
        hgnc_mut <- sprintf("%s_%s", hgnc, mut)
        mutBinMat[hgnc_mut, samp]<- exactMutI[[mut]]
      }
      for (mut in generTypes){
        hgnc_mut <- sprintf("%s_%s", hgnc, mut)
        mutBinMat[hgnc_mut, samp]<- generMutI[[mut]]
      }
    }
  }
  
  # for (muti in seq_along(mut_filt$hgnc_symbol)){
  #   geneI <- !is.na(match(mutDF$hgnc_symbol, mut_filt$hgnc_symbol[muti]))
  #   for (sampi in seq_along(sampNames)){
  #     currMuts <- mutDF[[sampNames[sampi]]]
  #     if (mut_filt$mut_name[muti] == "ANY"){
  #       nameI <- !is.na(currMuts) & currMuts != ""
  #     } else {
  #       nameI <- grepl(x = currMuts, pattern = mut_filt$mut_name[muti], fixed = TRUE)
  #     }
  #     mutI <- geneI & nameI
  #     mutBinMat[mut_strs[muti], sampNames[sampi]] <- sum(mutI)
  #   }
  # }
  mutBinDF <- data.frame(mutation_names = mut_strs, mutBinMat, stringsAsFactors = FALSE)
  rownames(mutBinDF) <- NULL
  
  return(mutBinDF)
}

# Exact mutation test functions
protExactMutI <- function(mutation_name, mutation_name_list){
  mutI <- vector(mode = "logical", length(mutation_name_list))
  names(mutI) <- mutation_name_list
  
  if (is.na(mutation_name) && length(mutation_name) == 1){
    mutI[] <- NA
    return(mutI)
  }
  
  if (length(mutation_name) == 0 || mutation_name == ""){
    return(mutI)
  }
  
  proteinNomen <- is_proteinNomen(mutation_name)
  if (is.na(proteinNomen) || !proteinNomen){
    warning("mutation_name is not in protein nomenclature form. NA returned")
    mutI[] <- NA
    return(mutTypes)
  }
  
  mutI[] <- !is.na(match(mutation_name_list, mutation_name))
  
  return(mutI)
}

# Generic mutation test functions
get_supported_protMutTypes <- function(){
  mutation_types <- c("missense", "nonsense", "initiating_methionine", 
                      "deletion", "insertion", "shortrepeat", "indel", "duplication", 
                      "frameshift", "splice", "silent", "nonsilent", "any")
  fun_names <- sprintf("is_mutType_%s", mutation_types)
  function_list <- mget(fun_names, ifnotfound = NA, envir = as.environment(1))
  names(function_list) <- NULL
  mutation_type_df <- data.frame(mutation_types = mutation_types, test_function = I(function_list), stringsAsFactors = FALSE)
  I <- !is.na(mutation_type_df$test_function)
  mutation_type_df <- mutation_type_df[I,]
  return(mutation_type_df)
}

protGenerMutI <- function(mutation_name, supported_muttype_df = NULL){
  
  if (is.null(supported_muttype_df)){
    supported_muttype_df <- get_supported_protMutTypes()
  }
  
  numTypes <- nrow(supported_muttype_df)
  mutTypes <- vector(mode = "logical", numTypes)
  names(mutTypes) <- supported_muttype_df$mutation_types
  
  if (is.na(mutation_name) && length(mutation_name) == 1){
    mutTypes[] <- NA
    return(mutTypes)
  }
  
  if (length(mutation_name) == 0 || mutation_name == ""){
    return(mutTypes)
  }
  
  proteinNomen <- is_proteinNomen(mutation_name)
  if (is.na(proteinNomen) || !proteinNomen){
    warning("mutation_name is not in protein nomenclature form. NA returned")
    mutTypes[] <- NA
    return(mutTypes)
  }
  mutation_name <- substring(text = mutation_name, first = 3)
  
  for (funi in seq_along(mutTypes)){
    mutTypeTest_fun <- supported_muttype_df$test_function[[funi]]  
    mutTypes[[funi]] <- mutTypeTest_fun(mutation_name = mutation_name)
  }
  return(mutTypes)
  
  
}

is_proteinNomen <- function(mutation_name){
  is_prot <-substr(x = mutation_name, start = 1, stop = 2) == "p."
  return(is_prot)
}

is_mutType_any <- function(mutation_name){
  is_type <- nchar(mutation_name) > 0
  return(is_type)
}

is_mutType_nonsilent <- function(mutation_name){
  is_type <- 
    is_mutType_missense(mutation_name) || is_mutType_nonsense(mutation_name) ||
    is_mutType_initiating_methionine(mutation_name) || is_mutType_deletion(mutation_name) ||
    is_mutType_insertion(mutation_name) || is_mutType_duplication(mutation_name) ||
    is_mutType_shortrepeat(mutation_name) || is_mutType_indel(mutation_name) ||
    is_mutType_frameshift(mutation_name) || is_mutType_splice(mutation_name)
  return(is_type)  
  
}

is_mutType_silent <- function(mutation_name){
  germAA <- substring(text = mutation_name, first = 1, last = 1)
  somaAA <- substring(text = mutation_name, first = nchar(mutation_name))
  pos <- substring(text = mutation_name, first = 2, last = nchar(mutation_name)-1)
  
  is_type <- !grepl("[^0-9]", pos) && germAA %in% AA_STANDARD && somaAA %in% AA_STANDARD && germAA == somaAA
  return(is_type)
}

is_mutType_missense <- function(mutation_name){
  germAA <- substring(text = mutation_name, first = 1, last = 1)
  somaAA <- substring(text = mutation_name, first = nchar(mutation_name))
  pos <- substring(text = mutation_name, first = 2, last = nchar(mutation_name)-1)
  
  is_type <- !grepl("[^0-9]", pos) && germAA %in% AA_STANDARD && somaAA %in% AA_STANDARD && germAA != somaAA
  return(is_type)
}

is_mutType_nonsense <- function(mutation_name){
  germAA <- substring(text = mutation_name, first = 1, last = 1)
  somaAA <- substring(text = mutation_name, first = nchar(mutation_name))
  pos <- substring(text = mutation_name, first = 2, last = nchar(mutation_name)-1)
  
  is_type <- !grepl("[^0-9]", pos) && germAA %in% AA_STANDARD && somaAA == "*"
  return(is_type)
}

is_mutType_initiating_methionine <- function(mutation_name){
  germAA <- substring(text = mutation_name, first = 1, last = 1)
  somaAA <- substring(text = mutation_name, first = nchar(mutation_name))
  pos <- substring(text = mutation_name, first = 2, last = nchar(mutation_name)-1)
  
  is_type <- pos == "1" && germAA == "M"
  return(is_type)
}

is_mutType_deletion <- function(mutation_name){
  is_type <- grepl(x = mutation_name, pattern = "del", fixed = TRUE)
  return(is_type)
}

is_mutType_insertion <- function(mutation_name){
  is_type <- grepl(x = mutation_name, pattern = "ins", fixed = TRUE) || is_mutType_duplication(mutation_name)
  return(is_type)
}

is_mutType_duplication <- function(mutation_name){
  is_type <- grepl(x = mutation_name, pattern = "dup", fixed = TRUE)
  return(is_type)
}

is_mutType_shortrepeat <- function(mutation_name){
  name_letters <- strsplit(x = mutation_name, split = "")[[1]]
  lp_pos <- which(name_letters %in% "(")
  rp_pos <- which(name_letters %in% ")")
  if (length(lp_pos) != 1 || length(rp_pos) != 1){
    return(FALSE)
  }
  if (rp_pos - lp_pos != 2){
    return(FALSE)
  }
  aa <- substring(text = mutation_name, first = lp_pos+1, last = rp_pos-1)
  if (!aa %in% AA_STANDARD){
    return(FALSE)
  }
  return(TRUE)
  
}

is_mutType_indel <- function(mutation_name){
  is_type <- grepl(x = mutation_name, pattern = "delins", fixed = TRUE)
  return(is_type)
}

is_mutType_frameshift <- function(mutation_name){
  is_type <- grepl(x = mutation_name, pattern = "fs", fixed = TRUE)
  return(is_type)
}

is_mutType_splice <- function(mutation_name){
  is_type <- grepl(x = mutation_name, pattern = "splice", fixed = TRUE)
  return(is_type)
}





