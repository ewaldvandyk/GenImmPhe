library(heatmap3)
source("./utility/data_transform.R")

genHeatMap <- function(df, rownameField = NULL, dataFields = NULL,
                       colPallet = gen3ColPallet(nCol = 1024, satLeft = 0.45, satRight = 0.55),
                       colSideColors = NULL, 
                       maxRowSizeToLabel = 40, maxColSizeToLabel = 40){
  
  numRows <- nrow(df)
  if (is.null(dataFields)){
    dataFields <- colnames(df)[get_df_numeric_colI(df)]
  }
  if (is.null(rownameField)){
    rowNames <- sprintf("%d", 1:numRows)
  } else {
    rowNames <- df[[rownameField]]
  }
  
  dataMat <- data.matrix(df[,dataFields])
  rownames(dataMat) <- rowNames

  nrows <- dim(dataMat)[1]
  ncols <- dim(dataMat)[2]
  if (nrows <= maxRowSizeToLabel){
    labRow <- rownames(dataMat)
  }else{
    labRow <- NA
  }
  if (ncols <= maxColSizeToLabel){
    labCol <- colnames(dataMat)
  }else{
    labCol <- NA
  }
  
  if (is.null(colSideColors)){
    heatmap3(x = dataMat,
             balanceColor = TRUE, showColDendro = TRUE,
             showRowDendro = TRUE, col=colPallet, labRow = labRow, labCol = labCol)
    
  } else {
    # colSideColors <- create_binary_sideColors(binaryDF = colSideBinDF, orderColList = colnames(dataMat))
    heatmap3(x = dataMat,
             balanceColor = TRUE, showColDendro = TRUE,
             showRowDendro = TRUE, col=colPallet, ColSideColors = colSideColors, 
             labRow = labRow, labCol = labCol)
    
  }
  
  
  return(dataMat)
}



gen3ColPallet <- function(nCol = 1024, lowCol = "green", midColor = "black", highColor = "red", satLeft = 0, satRight = 1){
  interpStart <- round(satLeft*nCol)+1
  interpEnd <- round(satRight*nCol)
  nInterp <- interpEnd - interpStart + 1
  colPallet <- colorRampPalette(c(lowCol, midColor, highColor))(nInterp)
  colPallet <- c(rep(colPallet[1], interpStart-1), colPallet, rep(colPallet[nInterp], nCol  - interpEnd))
  return(colPallet)
}


create_binary_sideColors <- function(binaryDF, nonDataCols = 1,rowNameCol = 1, orderColList = NULL, 
                                     true_color = "red", false_color = "black", na_color = "gray"){
  
  # colNames <- colnames(binaryDF)
  # numRow <- nrow(binaryDF)
  # if (is.null(orderColList)){
  #   numCol <- ncol(binaryDF)
  #   I <- 1:numCol
  # } else {
  #   missingCols <- setdiff(orderColList, colNames)
  #   for (newCol in missingCols){
  #     binaryDF[[newCol]] <- rep(2, numRow) # Use "2" to represent NA
  #   }
  #   colNames <- colnames(binaryDF)
  #   I <- c(1,match(orderColList, colNames))
  # }
  # binaryDF <- binaryDF[,I]
  rowNames <- binaryDF[[rowNameCol]]
  binaryDF <- binaryDF[,-nonDataCols]
  binaryDF <- order_df_cols(df = binaryDF, orderColVec = orderColList)
  colNames <- colnames(binaryDF)
  numRows <- length(rowNames)
  numCols <- length(colNames)
  sideColors <- matrix(data = na_color, nrow = numRows, ncol = numCols)
  rownames(sideColors) <- rowNames
  colnames(sideColors) <- colNames
  for (i in 1:numCols){
    currColBool <- binaryDF[[i]] == 1
    sideColors[currColBool,i] <- true_color
    currColBool <- binaryDF[[i]] == 0
    sideColors[currColBool,i] <- false_color
  }
  sideColors <- t(sideColors)
  
  return(sideColors)
}

create_cont_sideColors <- function(contDF, nonDataCols = 1, rowNameCol = 1, orderColList = NULL, 
                                   low_col = "green", neut_col = "black", high_col = "red", na_color = "gray", 
                                   neut_level = 0, saturate_low = -1, saturate_high = +1){
  # Allign columns
  rowNames <- contDF[[rowNameCol]]
  contDF <- contDF[,-nonDataCols]
  contDF <- order_df_cols(df = contDF, orderColVec = orderColList)
  colNames <- colnames(contDF)
  numRows <- length(rowNames)
  numCols <- length(colNames)
  
  # Map colors
  nColors <- 1025
  colScale <- gen3ColPallet(nCol = nColors, lowCol = low_col, midColor = neut_col, highColor = high_col, 
                            satLeft = 0, satRight = 1)
  na_code <- colorRampPalette(na_color)(1)
  colScale <- c(colScale, na_code)
  
  valueMat <- data.matrix(contDF)
  valiMat <- matrix(data = 0, nrow = numRows, ncol = numCols)
  naI <- is.na(valueMat)
  lowI <- !is.na(valueMat) & valueMat <= neut_level
  satLowI <- !is.na(valueMat) & valueMat < saturate_low
  higI <- !is.na(valueMat) & valueMat > neut_level
  satHighI <- !is.na(valueMat) & valueMat > saturate_high

  valiMat[lowI] <- (valueMat[lowI] - saturate_low) / (neut_level - saturate_low) * (nColors-1)/2 + 1
  valiMat[higI] <- (valueMat[higI] - neut_level) / (saturate_high - neut_level) * (nColors-1)/2 + (nColors-1)/2 + 1
  valiMat[satLowI] <- 1
  valiMat[satHighI] <- nColors
  valiMat[naI] <- nColors+1
  valiMat <- round(valiMat)
  
  sideColors <- colScale[valiMat]
  dim(sideColors) <- dim(valiMat)
  rownames(sideColors) <- rowNames
  colnames(sideColors) <- colNames
  
  sideColors <- t(sideColors)
  return(sideColors)
  
}


# Generic non-core functions
showCols <- function(cl=colors(), bg = "grey",
                     cex = 0.75, rot = 30) {
  # Generic function to show the different types of colors available in R
  
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm))
}


test_heatmap_complex <- function(dataList){
  heatmap3(x = dataList$rnormData,
           ColSideColors=dataList$ColSideColors,
           showRowDendro=FALSE,
           colorCell=dataList$colorCell,
           highlightCell=dataList$highlightCell)
  
  result<-heatmap3(x = dataList$rnormData,
                   ColSideCut=1.2,
                   ColSideAnn=dataList$ColSideAnn,
                   ColSideFun=function(x)showAnn(x),
                   ColSideWidth=0.8,
                   RowSideColors=dataList$RowSideColors,
                   col=colorRampPalette(c("green","black", "red"))(1024),
                   RowAxisColors=1,
                   legendfun=function()showLegend(legend=c("Low","High"),
                                                  col=c("chartreuse4","firebrick")),
                   verbose=TRUE)
  return(result)
}

test_generate_heatmap_data <- function(){
  set.seed(123456789)
  rnormData<-matrix(rnorm(1000), 40, 25)
  rnormData[1:15, seq(6, 25, 2)] = rnormData[1:15, seq(6, 25, 2)] + 2
  rnormData[16:40, seq(7, 25, 2)] = rnormData[16:40, seq(7, 25, 2)] + 4
  colnames(rnormData)<-c(paste("Control", 1:5, sep = ""), paste(c("TrtA", "TrtB"),
                                                                rep(1:10,each=2), sep = ""))
  rownames(rnormData)<-paste("Probe", 1:40, sep = "")
  ColSideColors<-cbind(Group1=c(rep("steelblue2",5), rep(c("brown1", "mediumpurple2"),10)),
                       Group2=sample(c("steelblue2","brown1", "mediumpurple2"),25,replace=TRUE))
  colorCell<-data.frame(row=c(1,3,5),col=c(2,4,6),color=c("green4","black","orange2"),
                        stringsAsFactors=FALSE)
  highlightCell<-data.frame(row=c(2,4,6),col=c(1,3,5),color=c("black","green4","orange2"),
                            lwd=1:3,stringsAsFactors=FALSE)
  ColSideAnn<-data.frame(Information=rnorm(25),Group=c(rep("Control",5), rep(c("TrtA", "TrtB"),10)))
  row.names(ColSideAnn)<-colnames(rnormData)
  RowSideColors<-colorRampPalette(c("chartreuse4", "white", "firebrick"))(40)
  
  heatmapData <- list()
  heatmapData$rnormData <- rnormData
  heatmapData$ColSideColors <- ColSideColors
  heatmapData$colorCell <- colorCell
  heatmapData$highlightCell <- highlightCell
  heatmapData$ColSideAnn <- ColSideAnn
  heatmapData$RowSideColors <- RowSideColors
  
  return(heatmapData)
}
