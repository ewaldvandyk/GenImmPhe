#################  Input parameters #################  
xlsxIn  <- "~/data/kevin/2020_06_17/RNAseq dataset for ewald.xlsx"
xlsxOut <- "~/data/kevin/2020_06_17/RNAseq dataset for kevin.xlsx"

minCountThresh <- 30
maxSampsBelowThresh <- 0

NonCountColumnPatterns <- c("logFC", "adj.P", "gene_id")
group1_pattern <- "_WT_"
group2_pattern <- "_KEP_"

###################################################  


# Install (if needed) and load required packages
requiredPackages <- c("purrr", "openxlsx")
for (package in requiredPackages){
  if (!package %in% rownames(installed.packages())){
    print(paste0("Installing ", package, "..."))
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = T)
}

# Load all xlsx sheets into a list of data frames
sheetNames <- getSheetNames(xlsxIn)
df_list <- map(seq_along(sheetNames), function(sheetIndex) readWorkbook(xlsxFile = xlsxIn, sheet = sheetIndex))
names(df_list) <- sheetNames

for (sheet in sheetNames){
  numCols <- ncol(df_list[[sheet]])
  
  #Figure out which columns contain counts and store it in countColI
  metaColI <- rep(F, numCols)
  for (pttr in NonCountColumnPatterns){
    metaColI <- metaColI | grepl(pattern = pttr, x = names(df_list[[sheet]]))
  }
  countColI <- !metaColI
  
  #Figure out which columns contain group 1 and 2 counts
  group1I <- countColI & grepl(pattern = group1_pattern, x = names(df_list[[sheet]]))
  group2I <- countColI & grepl(pattern = group2_pattern, x = names(df_list[[sheet]]))
  
  # Split data frames into groups
  group1DF <- df_list[[sheet]][,group1I,drop=F]
  group2DF <- df_list[[sheet]][,group2I,drop=F]
  
  #Apply thresholding
  group1numBelow <- rowSums(group1DF < minCountThresh)
  group2numBelow <- rowSums(group2DF < minCountThresh)
  keepI <- (group1numBelow <= maxSampsBelowThresh) | (group2numBelow <= maxSampsBelowThresh)
  
  df_list[[sheet]] <- df_list[[sheet]][keepI,,drop=F]
}

# Create and save new excel file
wb <- createWorkbook()
for (sheet in sheetNames){
  addWorksheet(wb, sheetName = sheet)
  writeData(wb, sheet = sheet, x = df_list[[sheet]])
}
saveWorkbook(wb = wb, file = xlsxOut)



