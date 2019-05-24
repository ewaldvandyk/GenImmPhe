plot_cip2a_expression <- function(countDF, colIgnore = 1:4, low_thresh = 10, high_thresh = 20){
  I <- !is.na(match(countDF$hgnc_symbol, "CIP2A"))
  
  x <- t(data.matrix(countDF[I, -colIgnore]))

  plot(sort(x))
  hist(x, breaks = 1000)
}

split_low_high_cip2a <- function(countDF, colIgnore = 1:4, low_thresh = 12.2, high_thresh = 25.7){
  I <- !is.na(match(countDF$hgnc_symbol, "CIP2A"))
  x <- t(data.matrix(countDF[I, -colIgnore]))
  I_low = x <= low_thresh
  I_high = x >= high_thresh
  
  samps <- dimnames(x)[[1]]
  samps_low <- samps[I_low]
  samps_high <- samps[I_high]
  sampList <- list(low = samps_low, high = samps_high)
  return(sampList)
}

split_cibersort_results <- function(ciberfile, samps){
  ciberDF <- load_tsv(inFile = ciberfile)
  I_low <- !is.na(match(ciberDF$`Input Sample`, samps$low))
  I_high <- !is.na(match(ciberDF$`Input Sample`, samps$high))
  
  ciberLowDF <- ciberDF[I_low,]
  ciberHighDF <- ciberDF[I_high,]
  str(ciberLowDF)
  str(ciberHighDF)
  ciberDF <- list(low = ciberLowDF, high = ciberHighDF)
  return(ciberDF)
}

add_cip2A_expr_2_ciberData <- function(ciberDF, countDF, colIgnore = 1:4){
  names(ciberDF)[1] <- "sampName"
  
  cip2AI <- !is.na(match(countDF$hgnc_symbol, "CIP2A"))
  
  cip2A_matrix <- t(data.matrix(countDF[cip2AI, -colIgnore]))
  colnames(cip2A_matrix)[1] <- "cip2a_exp"
  cip2aDF <- data.frame(sampName = rownames(cip2A_matrix), cip2a_exp = cip2A_matrix, stringsAsFactors = FALSE)
  str(cip2aDF)
  
  ciberCIp2aDF <- merge(x = cip2aDF, y = ciberDF, by = "sampName")
  
  return(ciberCIp2aDF)
}