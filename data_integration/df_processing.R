source("./utility/data_transform.R")

allignDFs <- function(df_list, nonDataFieldNames, NAsAllowed = 1:length(df_list)){
  # NAsAllowed is vector of list indexes showing which data frames are allowed to have NA columns 
  numDFs <- length(df_list)
  notNAsAllowed <- setdiff(1:numDFs, NAsAllowed)
  
  dataField_list <- list()
  colVec <- c()
  for (dfi in 1:numDFs){
    currDF <- df_list[[dfi]]
    currDFColNames <- colnames(currDF)
    dataField_list[[dfi]] <- setdiff(currDFColNames, nonDataFieldNames)
    colVec <- union(colVec, dataField_list[[dfi]])
  }
  for (dfi in notNAsAllowed){
    colVec <- intersect(colVec, dataField_list[[dfi]])
  }
  colVec <- sort(colVec, decreasing = FALSE)
  
  for (dfi in 1:numDFs) {
    currDF <- df_list[[dfi]]
    currDFColNames <- colnames(currDF)
    nonSortCols <- which(!is.na(match(currDFColNames, nonDataFieldNames)))
    
    currDF <- order_df_cols(df = currDF, nonSortCols = nonSortCols, orderColVec = colVec)


    df_list[[dfi]] <- currDF
  }
  

  return(df_list)
}