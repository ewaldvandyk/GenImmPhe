setwd("~/devel/R/GenImmPhe")
source("./user_tools/data_split/setup.R")
source("./user_tools/data_split/split_methods.R")

# Setup input parameters and load data sources into memory
data_sources <- setup_data_sources()
df_list <- load_df_list(data_sources)
parse_struct <- setup_split(df_list)

# Split samples based on parse_struct
df_list <- get_split_samps(df_list = df_list, parseStruct = parse_struct, data_sources = data_sources)

# Save filtered data sources
save_df_list(df_list, data_sources, parse_struct)