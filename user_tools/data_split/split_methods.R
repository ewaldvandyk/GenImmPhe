source("./data_integration/df_processing.R", local = TRUE)
source("./utility/io_evd.R", local = TRUE)

load_df_list <- function(data_sources){
  numInSources <- nrow(data_sources$input_data_sources)
  df_list <- list()
  non_data_fields <- c()
  for (dfi in 1:numInSources){
    currDFname <- load(file.path(data_sources$input_dir, data_sources$input_data_sources$file_names[dfi]))
    df_list[[data_sources$input_data_sources$source_names[dfi]]] <- get(currDFname)
    rm(currDFname)
    curr_non_data_fields <- strsplit(gsub(pattern = " ", 
                                          replacement = "", 
                                          x = data_sources$input_data_sources$field_ids[dfi], 
                                          fixed = TRUE), 
                                     split = "|", 
                                     fixed = TRUE)[[1]]
    non_data_fields <- c(non_data_fields, curr_non_data_fields)
  }
  non_data_fields <- unique(non_data_fields)
  nais <- which(data_sources$input_data_sources$NAs_allowed)
  df_list <- allignDFs(df_list = df_list, 
                       nonDataFieldNames = non_data_fields, 
                       NAsAllowed = nais)
  
  
  return(df_list)
}

save_df_list <- function(df_list, data_sources, parse_struct){
  dir.create(path = data_sources$output_dir)
  num_output_sources <- nrow(data_sources$output_data_sources)
  for (dfi in 1:num_output_sources){
    curr_file <- file.path(data_sources$output_dir, data_sources$output_data_sources$file_names[dfi])
    save_tsv(df_list[[data_sources$output_data_sources$source_names[dfi]]], outFile = curr_file)
  }
  info_file <- file.path(data_sources$output_dir, "data_info.txt")
  info_strs <- get_df_list_stats_string(df_list = df_list, data_sources = data_sources, parse_struct = parse_struct)
  fileConn<-file(info_file)
  writeLines(info_strs, fileConn)
  close(fileConn)
}

get_df_list_stats_string <- function(df_list, data_sources, parse_struct){
 
  sourceNames <- data_sources$output_data_sources$source_names
  numSources <- length(sourceNames)
  
  
  strings <- "Parse_structure:"
  strings <- c(strings, capture.output(print(parse_struct)))
  strings <- c(strings, rep("", 2))
  strings <- c(strings, "Data sources: ")
  strings <- c(strings, capture.output(print(data_sources$output_data_sources)))
  for (dfi in 1:numSources){
    non_data_fields <- strsplit(gsub(pattern = " ", 
                                     replacement = "", 
                                     x = data_sources$output_data_sources$field_ids[dfi], 
                                     fixed = TRUE), 
                                split = "|", 
                                fixed = TRUE)[[1]]
    sampNames <- setdiff(colnames(df_list[[sourceNames[dfi]]]), non_data_fields)
    strings <- c(strings, sprintf("%s:", sourceNames[dfi]))
    strings <- c(strings, sprintf("\tNumber of samples = %d", length(sampNames)))
    strings <- c(strings, "")
  }
  return(strings)
}

get_split_samps <- function(df_list, parseStruct, data_sources){
  numSources <- length(df_list)
  if (numSources == 0){
    return(df_list)
  }
  
  curr_non_data_fields <- strsplit(gsub(pattern = " ", 
                                        replacement = "", 
                                        x = data_sources$input_data_sources$field_ids[1], 
                                        fixed = TRUE), 
                                   split = "|", 
                                   fixed = TRUE)[[1]]
  samp_names <- setdiff(colnames(df_list[[1]]), curr_non_data_fields)
  source_names <- names(df_list)
  parse_phrase <- parse_struct_resolve(parse_struct = parseStruct)
  
  samps_ignore <- c()
  for (samp in samp_names){
    for (dfi in 1:numSources){
      assign(x = source_names[dfi], df_list[[dfi]][[samp]])
    }
    sampValid <- eval(parse(text = parse_phrase))
    if (is.na(sampValid) || !sampValid){
      samps_ignore <- c(samps_ignore, samp)
    }
  }
  
  for (dfi in 1:numSources){
    colNames <- colnames(df_list[[dfi]])
    keepI <- is.na(match(colNames, samps_ignore))
    df_list[[dfi]] <- df_list[[dfi]][,keepI]
  }
  
  return(df_list)
}

parse_struct_resolve <- function(parse_struct){
  paramNames <- names(parse_struct$params)
  parse_phrase <- parse_struct$parsePhrase
  for (paramName in paramNames){
    parse_phrase = gsub(x = parse_phrase, pattern = paramName, replacement = parse_struct$params[[paramName]])
  }
  return(parse_phrase)
}

set_param <- function(df_list, source_name, id_colname, id, data_type){
  data_types <- c("numeric", "character")
  if (!(source_name %in% names(df_list))){
    stop("source_name needs to be one of the following strings: ", sprintf("%s; ", names(df_list)))
  }
  
  df <- df_list[[source_name]]
  if (!(id_colname %in% names(df))){
    stop("id_colname needs to be set to the column name describing the rowname entries in data.frame 'df_list[[source_name]]'")
  }
  
  ids <- df[[id_colname]]
  if (!(id %in% ids)){
    stop("id not found within id_column")
  }
  
  if (!(data_type %in% data_types)){
    stop("data_types need to be either 'numeric' or 'character'")
  }
  
  i <- which(ids == id)
  
  if (length(i) > 1){
    stop("Multiple entries found for id. Make sure that 'id_colname' is set correctly")
  }
  if (data_type == data_types[2]){
    return(sprintf("%s[%d]", source_name, i))
  } else {
    return(sprintf("as.numeric(%s[%d])", source_name, i))
  }
}

select_output_sources <- function(input_data_sources, source_name_subset, output_file_type = ".tsv"){
  if (length(source_name_subset) == 0 || !all(source_name_subset %in% input_data_sources$source_names)){
    stop("source_name_subset has to be a proper subset of input_data_sources$source_names")
  }
  I <- !is.na(match(input_data_sources$source_names, source_name_subset))
  output_data_sources <- input_data_sources[I,]
  
  num_sources <- length(source_name_subset)
  for (si in 1:num_sources){
    dotsplits <- strsplit(x = output_data_sources$file_names[si], split = ".", fixed = TRUE)[[1]] 
    if (length(dotsplits) == 1){
      output_file <- paste(dotsplits, output_file_type, sep = "")
    } else {
      output_file <- paste(dotsplits[-length(dotsplits)], collapse = '.')
      output_file <-  paste(output_file, output_file_type, sep = "")
    }
    output_data_sources$file_names[si] <- output_file
  }
  
  return(output_data_sources)
}


