transpose_df <- function(df, rowNameCol, colTypeName = "row_name"){
  colNames <- colnames(df)
  colKeep <- colNames != rowNameCol
  
  rowNames <- df[[rowNameCol]]
  df <- df[,colKeep]
  rownames(df) <- rowNames
  tMat <- t(df)

  t_df <- data.frame(tMat, stringsAsFactors = FALSE)
  t_df[[colTypeName]] <- rownames(tMat)
  t_df <- t_df[, c(ncol(t_df), 1:(ncol(t_df)-1))]
  return(t_df)
}

order_df_cols <- function(df, nonSortCols = NULL, 
                          orderColVec = NULL){
  if (is.null(orderColVec)){
    return(df)
  }
  if (is.null(nonSortCols)){
    df_sort <- df
    df_nonSort <- data.frame(matrix(data = NA, nrow=nrow(df), ncol=0))
  } else {
    df_sort <- df[,-nonSortCols, drop = FALSE]
    df_nonSort <- df[,nonSortCols, drop = FALSE]
  }
  
  colNames <- colnames(df_sort)
  numRow <- nrow(df_sort)
  missingCols <- setdiff(orderColVec, colNames)
  for (newCol in missingCols){
    df_sort[[newCol]] <- rep(NA, numRow)
  }
  colNames <- colnames(df_sort)
  I <- match(orderColVec, colNames)

  df_sort <- df_sort[,I, drop = FALSE]

  return(cbind(df_nonSort, df_sort))  
 
}

get_df_numeric_colI <- function(df){
  fieldTypes <- sapply(df, class)
  dataFields <- names(fieldTypes[fieldTypes == "numeric"])
  I <- !is.na(match(colnames(df), dataFields))
  return(I)
}