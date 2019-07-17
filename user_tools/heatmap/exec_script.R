setwd("~/devel/R/GenImmPhe")
source("./user_tools/heatmap/setup.R")
source("./user_tools/heatmap/heatmap_methods.R")

# Setup input parameters and load data sources into memory
data_sources <- setup_data_sources(cluster_ids = c("Akt_pS473", "Akt_pT308", "mTOR_pS2448", "p70S6K_pT389", "4E-BP1_pS65"))
df_list <- load_df_list(data_sources)
sideSampList <- setup_sideSamp(df_list)

# Create heatmap
heatmapData <- create_heatmap(df_list = df_list, sideSampList = sideSampList, data_sources = data_sources)

