library(heatmap3)

genImmPhe_path <- getwd()
source(file.path(genImmPhe_path, "utility/data_transform.R"), local = TRUE)

genHeatMap <- function(df, rownameField = NULL, dataFields = NULL,
                       colPallet = NULL,
                       colSideAnn = NULL, 
                       maxRowSizeToLabel = 40, maxColSizeToLabel = 40, ...){
  if (is.null(colPallet)){
    colPallet <- gen3ColPallet(nCol = 1024, satLeft = 0.3, satRight = 0.7)
  }
  numRows <- nrow(df)
  if (is.null(dataFields)){
    dataFields <- colnames(df)[get_df_numeric_colI(df)]
  }
  if (is.null(rownameField)){
    rowNames <- sprintf("%d", 1:numRows)
  } else {
    rowNames <- df[[rownameField]]
  }
  # print(rowNames)
  
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
  
  # str(colSideAnn)
  # str(dataMat)
  if (is.null(colSideAnn)){
   clustRes <- heatmap3(x = dataMat,
             balanceColor = TRUE, showColDendro = TRUE,
             showRowDendro = TRUE, col=colPallet, labRow = labRow, labCol = labCol, 
             na.rm = TRUE, ...)
    
  } else {
    clustRes <- heatmap3(x = dataMat, 
             showColDendro = TRUE, showRowDendro = TRUE,
             labRow = labRow, labCol = labCol,
             balanceColor = TRUE, col=colPallet, 
             ColSideAnn = colSideAnn, ColSideFun = showColSideCols, ColSideWidth = getColSideColStripLength(colSideAnn),
             na.rm = TRUE, ...)
    
  }
  
  return(clustRes)
}

genHeatMapOld <- function(df, rownameField = NULL, dataFields = NULL,
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
    heatmap3(x = dataMat,
             balanceColor = TRUE, showColDendro = TRUE,
             showRowDendro = TRUE, col=colPallet, ColSideColors = colSideColors, 
             labRow = labRow, labCol = labCol)
    
  }
  
  
  return(dataMat)
}

gen3ColPallet <- function(nCol = 1024, lowCol = "navy", midColor = "white", highColor = "firebrick3", satLeft = 0, satRight = 1){
  interpStart <- round(satLeft*nCol)+1
  interpEnd <- round(satRight*nCol)
  nInterp <- interpEnd - interpStart + 1
  colPallet <- colorRampPalette(c(lowCol, midColor, highColor))(nInterp)
  colPallet <- c(rep(colPallet[1], interpStart-1), colPallet, rep(colPallet[nInterp], nCol  - interpEnd))
  return(colPallet)
}

lgl2Color <- function(lglVec, true_color = "red", false_color = "black", na_color = "grey"){
  vecLen <- length(lglVec)
  colVec <- rep(na_color, vecLen)
  colVec[lglVec] <- true_color
  colVec[!lglVec] <- false_color
  return(colVec)
}

num2Color <- function(numVec, low_col = "green", neut_col = "black", high_col = "red", na_color = "grey", 
                    neut_level = 0, saturate_low = -1, saturate_high = +1){
  # Map colors
  nColors <- 1025
  colScale <- gen3ColPallet(nCol = nColors, lowCol = low_col, midColor = neut_col, highColor = high_col, 
                            satLeft = 0, satRight = 1)
  na_code <- colorRampPalette(na_color)(1)
  colScale <- c(colScale, na_code)
  
  naI <- is.na(numVec)
  lowI <- !is.na(numVec) & numVec <= neut_level
  satLowI <- !is.na(numVec) & numVec < saturate_low
  higI <- !is.na(numVec) & numVec > neut_level
  satHighI <- !is.na(numVec) & numVec > saturate_high
  
  valiVec <- rep(0, length(numVec))
  valiVec[lowI] <- (numVec[lowI] - saturate_low) / (neut_level - saturate_low) * (nColors-1)/2 + 1
  valiVec[higI] <- (numVec[higI] - neut_level) / (saturate_high - neut_level) * (nColors-1)/2 + (nColors-1)/2 + 1
  valiVec[satLowI] <- 1
  valiVec[satHighI] <- nColors
  valiVec[naI] <- nColors+1
  valiVec <- round(valiVec)
  
  colVec <- colScale[valiVec]
}

fac2Color <- function(facVec, brewer_pal = "Set1", na_color = "grey"){
  require(RColorBrewer)
  fac_levels <- levels(facVec)
  numLevels <- length(fac_levels)
  
  tryCatch(
    warning = function(cnd){
      stop("Number of factor levels exceeds the chosen pallet size", call. = FALSE)
    },
    pal <- brewer.pal(n = max(3, numLevels), name = brewer_pal)
  )

  colVec <- pal[facVec]
  colVec[is.na(colVec)] <- na_color
  names(pal) <- fac_levels
  return(list(colorVec = colVec, levelColors = pal))
}

showColSideCols <- function(annData){
  require(purrr)
  fieldTypes <- map_chr(.x = annData, .f = class)
  lglI <- fieldTypes == "logical"
  numI <- fieldTypes == "numeric"
  facI <- fieldTypes == "factor"
  
  dataLgl <- annData[,lglI, drop = FALSE]
  dataNum <- annData[,numI, drop = FALSE]
  dataFac <- annData[,facI, drop = FALSE]
  
  binWidth <- 1/nrow(annData)
  halfBinWidth <- binWidth/2
  LeftBound <- -halfBinWidth
  RightBound <- 1 + halfBinWidth 
  
  numLines <- ncol(dataLgl) + ncol(dataNum) + 2*ncol(dataFac)
  plot(c(LeftBound, RightBound), c(0, numLines), 
       type = "n", xaxt = "n", yaxt = "n", xlab = "", 
       ylab = "", bty = "n", axes = FALSE, xaxs = "i")
  lines(x = c(LeftBound, LeftBound, RightBound, RightBound, LeftBound), 
        y = c(0, numLines, numLines, 0, 0))
  
  xleft <- seq(from = LeftBound, to = 1-halfBinWidth, length.out = nrow(annData))
  xright <- xleft + binWidth
  
  offset <- 0
  if (ncol(dataLgl) != 0){
    mtext(side = 2, at = seq_along(dataLgl) + offset - 0.5, 
          text = sprintf("%s ", colnames(dataLgl)), las = 1)
    for (lglVec in dataLgl){
      colorVec <- lgl2Color(lglVec)
      
      rect(xleft = xleft, xright = xright, ybottom = offset, ytop = offset+1, col = colorVec, border = colorVec)
      lines(x = c(LeftBound, RightBound), y = c(offset+1, offset+1))
      offset <- offset+1
    }
  }
  
  if (ncol(dataNum) != 0){
    mtext(side = 2, at = seq_along(dataNum) + offset - 0.5,
          text = sprintf("%s ", colnames(dataNum)), las = 1)
    for (numVec in dataNum){
      colorVec <- num2Color(numVec)
      rect(xleft = xleft, xright = xright, ybottom = offset, ytop = offset+1, col = colorVec, border = colorVec)
      lines(x = c(LeftBound, RightBound), y = c(offset+1, offset+1))
      offset <- offset + 1
    }
  }
  
  legLeft <- 1/20
  legColFrac <- 0.5
  legColHight <- 0.5
  legxBias <- -0.1*legLeft
  legyBias <- -0.1
  
  legColWidth <- legLeft*legColFrac
  legLeftWidth <- legLeft - legColWidth
  legColHightRes <- (1 - legColHight)/2
  
  if (ncol(dataFac) != 0){
    mtext(side = 2, at = 2*seq_along(dataFac) + offset - 1,
          text = sprintf("%s ", colnames(dataFac)), las = 1)
    for (facVec in dataFac){
      pal_info <- fac2Color(facVec)
      colorVec <- pal_info$colorVec
      legColVec <- pal_info$levelColors
      rect(xleft = xleft, xright = xright, ybottom = offset, ytop = offset+1, col = colorVec, border = colorVec)
      levelNames <- levels(facVec)
      numLevels <- length(levelNames)
      FacWidth <- 1.0 / numLevels
      xLegTextPos <- (seq_along(levelNames)-1)*FacWidth + legColWidth +legLeftWidth
      yLegTextPos <- rep(x = offset + 1.5, times = numLevels)
      xLegLeft <- xLegTextPos - legColWidth + legxBias
      xLegRight <- xLegTextPos + legxBias
      yLegBottom <- yLegTextPos - legColHightRes + legyBias
      yLegTop <- yLegTextPos + legColHightRes + legyBias
      text(x = xLegTextPos, y = yLegTextPos, labels = levelNames, adj = 0)
      rect(xleft = xLegLeft, xright = xLegRight, ybottom = yLegBottom, ytop = yLegTop, col = legColVec, border = NA)
      
      lines(x = c(LeftBound, RightBound), y = c(offset+2, offset+2))
      offset <- offset + 2
    }
  }
  
  # offset <- numLines
  return(c(0.0, offset))
}

getColSideColStripLength <- function(annData){
  lenMult <- 0.2
  
  require(purrr)
  fieldTypes <- map_chr(.x = annData, .f = class)
  lglI <- fieldTypes == "logical"
  numI <- fieldTypes == "numeric"
  facI <- fieldTypes == "factor"
  numLines <- sum(lglI) + sum(numI) + 2*sum(facI)
  lenMult*numLines
}

# Generic non-core functions
test_heatmap_complex <- function(dataList){
  heatmap3(x = dataList$rnormData,
           ColSideColors=dataList$ColSideColors,
           showRowDendro=FALSE,
           colorCell=dataList$colorCell,
           highlightCell=dataList$highlightCell)
  
  result<-heatmap3(x = dataList$rnormData,
                   ColSideAnn=dataList$ColSideAnn,
                   ColSideFun=function(x)showColSideCols(x),
                   ColSideWidth=getColSideColStripLength(dataList$ColSideAnn),
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
