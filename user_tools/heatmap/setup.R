source("./user_tools/heatmap/heatmap_methods.R", local = TRUE)

setup_data_sources <- function(cluster_ids = NULL){
  input_dir <- "/Volumes/Ewald/Antoinette/TCGA/lumA_pten_nonsilent/"
  output_dir <- file.path(input_dir, "plots") #"~/temp/plots/" 
  input_data_sources  <- 
    data.frame(source_names = 
                 c("MOLECULAR_SUBTYPES",
                   "MUTATIONS",
                   "CNAS",
                   "RPPA"),
               file_names = 
                 c("molecular_subtypes.tsv",
                   "mut_pik_akt_pten.tsv",
                   "cna_pten.tsv",
                   "rppa_pi3k.tsv"), 
               field_ids = 
                 c("molecular_subtype",
                   "mutation_names",
                   "hgnc_symbol",
                   "protein_names"), stringsAsFactors = FALSE)
  
  heatmap_source <- list(source_name = "RPPA", 
                         id_colname = "protein_names",
                         ids = cluster_ids, 
                         tag = "akt_mtor_p70S6K_EBP1")
  
  data_sources <-  list(input_dir=input_dir, output_dir = output_dir, 
                       input_data_sources=input_data_sources, heatmap_source = heatmap_source)
  return(data_sources)
}


setup_sideSamp <-function(df_list){
  sideSampList <- list()
  
  ##############  Setup sample side colours for clustering   ##############
  # To see available data_sources, console: str(df_list, max.level = 1)
  # To see available id_colnames, console: names(df_list$<data_source>)[1:4]  (e.g. names(df_list$pam50DF)[1:4])
  # To see available ids, console: df_list$<data_source>$<id_colname> (e.g. df_list$subtype$molecular_subtype)
  sideSampList$pam50       <- set_sideSamp(df_list = df_list, source_name = "MOLECULAR_SUBTYPES", 
                                       id_colname = "molecular_subtype", id = "pam50", 
                                       data_type = "factor")
  
  sideSampList$PIK3CA_E545K  <- set_sideSamp(df_list = df_list, source_name = "MUTATIONS", 
                                       id_colname = "mutation_names", id = "PIK3CA_p.E545K", 
                                       data_type = "logical")
  sideSampList$PIK3CA_H1047R <- set_sideSamp(df_list = df_list, source_name = "MUTATIONS",
                                       id_colname = "mutation_names", id = "PIK3CA_p.H1047R", 
                                       data_type = "logical")
  sideSampList$PIK3CA_nonsilent  <- set_sideSamp(df_list = df_list, source_name = "MUTATIONS", 
                                           id_colname = "mutation_names", id = "PIK3CA_nonsilent", 
                                           data_type = "logical")
  sideSampList$AKT1_E17K     <- set_sideSamp(df_list = df_list, source_name = "MUTATIONS",
                                       id_colname = "mutation_names", id = "AKT1_p.E17K", 
                                       data_type = "logical")
  sideSampList$AKT1_nonsilent  <- set_sideSamp(df_list = df_list, source_name = "MUTATIONS", 
                                         id_colname = "mutation_names", id = "AKT1_nonsilent", 
                                         data_type = "logical")
  sideSampList$PTEN_nonsilent <- set_sideSamp(df_list = df_list, source_name = "MUTATIONS",
                                        id_colname = "mutation_names", id = "PTEN_nonsilent", 
                                        data_type = "logical")
  sideSampList$PTEN_cna      <- set_sideSamp(df_list = df_list, source_name = "CNAS",
                                       id_colname = "hgnc_symbol", id = "PTEN",
                                       data_type = "numeric", zero_ref = 0,  saturation_frac = 0.05)
  
  return(sideSampList)
}



