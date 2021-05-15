library(sticky)

genImmPhe_path <- getwd()
source(file.path(genImmPhe_path, "utility/io_evd.R"), local = TRUE)
source(file.path(genImmPhe_path, "utility/plot_scripts.R"), local = TRUE)
source(file.path(genImmPhe_path, "data_integration/df_processing.R"), local = TRUE)

load_df_list <- function(data_sources){
  require(purrr)
  numInSources <- nrow(data_sources$input_data_sources)
  df_list <- list()
  non_data_fields <- c()
  for (dfi in seq_along(data_sources$input_data_sources$source_names)){
    currDfFile <- file.path(data_sources$input_dir, data_sources$input_data_sources$file_names[dfi])
    currDfExt <- tools::file_ext(currDfFile)
    if (currDfExt == "Rda"){
      currDF <- readRDS(currDfFile)
    } else if ((currDfExt == "tsv") || (currDfExt == "txt")){
      currDF <- load_tsv(currDfFile)  
    }
    keepI <- map_lgl(currDF, function(x) sum(is.na(x)) != length(x))
    currDF <- currDF[,keepI]
    df_list[[data_sources$input_data_sources$source_names[dfi]]] <- currDF
    rm(currDF)
    curr_non_data_fields <- strsplit(gsub(pattern = " ", 
                                          replacement = "", 
                                          x = data_sources$input_data_sources$field_ids[dfi], 
                                          fixed = TRUE), 
                                     split = "|", 
                                     fixed = TRUE)[[1]]
    non_data_fields <- c(non_data_fields, curr_non_data_fields)
  }
  non_data_fields <- unique(non_data_fields)
  
  NAsAllowed <- setdiff(seq_along(df_list), which(names(df_list) == data_sources$heatmap_source$source_name))
  df_list <- allignDFs(df_list = df_list, 
                       nonDataFieldNames = non_data_fields, 
                       NAsAllowed = NAsAllowed)
  
  
  return(df_list)
}

set_sideSamp <- function(df_list , source_name,
                         id_colname, id,
                         data_type = NA, zero_ref = "median", colPallet = NULL, saturation_frac = 0.0, lowColor=NULL, highColor=NULL, naColor = "black", firstLevel=NULL){
  
  zero_ref_choices <- c("median", "mean")
  data_type_choices <- c("numeric", "logical", "factor")
  
  if (!data_type %in% data_type_choices){
    stop(sprintf("data_type must be one of the following strings: '%s'", paste(data_type_choices, collapse = "', '")))
  }
  
  if (!(is.numeric(zero_ref) || is.character(zero_ref)) || length(zero_ref) != 1){
    stop("zero_ref parameter has to be either a numeric scalar or a string scalar")
  }
  if (is.character(zero_ref) && !(zero_ref %in% zero_ref_choices)){
    stop(sprintf("zero_ref must be one of the following strings: '%s'", paste(zero_ref_choices, collapse = "', '")))
  }
  
  if (!(source_name %in% names(df_list))){
    stop(sprintf("source_name must be one of the following strings: '%s'", paste(names(df_list), collapse = "', '")))
  }
  
  sourceDF <- df_list[[source_name]]
  
  if (!(id_colname %in% names(sourceDF))){
    stop("id_colname not found")
  }
  
  i <- which(sourceDF[[id_colname]] == id)
  if (length(i) == 0){
    stop("id not found")
  }
  if (length(i) > 1){
    stop("Multiple occurrences of 'id' found")
  }
  
  info <- list()
  info$location$source_name <- source_name
  info$location$id_colname  <- id_colname
  info$location$idi <- i
  info$transform$type <- data_type
  info$transform$zero_ref <- zero_ref
  info$transform$colPallet = colPallet
  info$transform$saturation_frac <- saturation_frac
  info$transform$lowColor <- lowColor
  info$transform$highColor <- highColor
  info$transform$naColor <- naColor
  info$transform$firstLevel <- firstLevel
  return(info)
  
  
}

create_heatmap <- function(df_list, sideSampList, data_sources, colPallet = NULL, ...){
  numSources <- length(df_list)
  if (numSources == 0){
    stop("No data source specified")
  }
  
  curr_non_data_fields <- strsplit(gsub(pattern = " ", 
                                        replacement = "", 
                                        x = data_sources$input_data_sources$field_ids[1], 
                                        fixed = TRUE), 
                                   split = "|", 
                                   fixed = TRUE)[[1]]
  samp_names <- setdiff(colnames(df_list[[1]]), curr_non_data_fields)
  
  heatmapDF <- df_list[[data_sources$heatmap_source$source_name]]
  if (!is.null(data_sources$heatmap_source$ids)){
    I <- !is.na(match(heatmapDF[[data_sources$heatmap_source$id_colname]], data_sources$heatmap_source$ids))
    heatmapDF <- heatmapDF[I,]
  }
  rownameField <- data_sources$heatmap_source$id_colname
  dataCols <- samp_names
  
  colSideAnn <- df_list2sideAnn(df_list, samp_names, sideSampList)

  if (!is.null(data_sources$output_dir)){
    dir.create(data_sources$output_dir, showWarnings = FALSE)
    heatmapFileName <- sprintf("heatmap_%s_%s", data_sources$heatmap_source$source_name, data_sources$heatmap_source$tag)
    jpeg(file.path(data_sources$output_dir, paste0(heatmapFileName, ".jpg")), width = 210, height = 297,units = "mm", res = 72*4)
  }
  # str(heatmapDF)
  heatData <- genHeatMap(df = heatmapDF, rownameField = data_sources$heatmap_source$id_colname, dataFields = samp_names, 
             colSideAnn = colSideAnn, colPallet, ...)
  if (!is.null(data_sources$output_dir)){
    dev.off()
  }
  
  if (!is.null(data_sources$output_dir)){
    
    save(file = file.path(data_sources$output_dir, paste0(heatmapFileName, ".Rda")), heatData)
  }
  
  
  
  
  return(heatData)
}

df_list2sideAnn <- function(df_list, samp_names, sideSampList){
  colSideAnn <- list()
  for (ann in names(sideSampList)){
    colSideAnn[[ann]] <- getSideAnnVec(df_list, samp_names, sideSampList[[ann]])
  }
  return(as.data.frame(colSideAnn))
}

getSideAnnVec <- function(df_list, samp_names, sideSampInfo){
  require(purrr)
  sourceDF <- df_list[[sideSampInfo$location$source_name]]
  dfSlice <- sourceDF[sideSampInfo$location$idi, samp_names]
  if (sideSampInfo$transform$type == "factor"){
    vec <- factor(map_chr(dfSlice, as.character))
    if (!is.null(sideSampInfo$transform$firstLevel)){
      vec <- relevel(x = vec, ref = sideSampInfo$transform$firstLevel)  
    }
  } else if (sideSampInfo$transform$type == "logical"){
    vec <- map_lgl(dfSlice, as.logical)
  } else if (sideSampInfo$transform$type == "numeric"){
    vec <- map_dbl(dfSlice, as.numeric)
    if (sideSampInfo$transform$zero_ref == "mean"){
      zero_ref <- mean(vec, na.rm = T)
    } else if (sideSampInfo$transform$zero_ref == "median"){
      zero_ref <- median(vec, na.rm = T)
    } else{
      zero_ref <- sideSampInfo$transform$zero_ref
    }
    vecSort <- sort(x = vec, decreasing = FALSE, na.last = NA)
    leftMax <- max(vecSort[vecSort < zero_ref])
    rightMin <- min(vecSort[vecSort > zero_ref])
    leftSat <- vecSort[max(round(length(vecSort)*sideSampInfo$transform$saturation_frac/2),1)]
    rightSat <- vecSort[round(length(vecSort)*(1-sideSampInfo$transform$saturation_frac/2))]
    leftSat <- min(leftMax, leftSat) - zero_ref
    rightSat <- max(rightMin, rightSat) - zero_ref
    vec <- vec - zero_ref
    vec[!is.na(vec) & vec < 0] <- -vec[!is.na(vec) & vec < 0]/leftSat
    vec[!is.na(vec) & vec > 0] <- vec[!is.na(vec) & vec > 0]/rightSat
    attr(vec, "colPallet") <- sideSampInfo$transform$colPallet

  }
  attr(vec, "naColor") <- sideSampInfo$transform$naColor
  return(sticky(vec))
}
