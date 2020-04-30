genImmPhe_path <- getwd()
source(file.path(genImmPhe_path, "user_tools/data_split/split_methods.R"), local = TRUE)

setup_data_sources <- function(){
  input_dir <- "~/data/pipeline_interm/2019_05_TCGA_BRCA"
  input_data_sources  <- 
    data.frame(source_names = 
                 c("MOLECULAR_SUBTYPES",
                   "MUTATIONS",
                   "CNAS",
                   "RNASEQ",
                   "RPPA"),
               file_names = 
                 c("molecular_subtypes.Rda",
                   "mut_pik_akt_pten.Rda",
                   "cna_pten.Rda",
                   "voomRNAseq.Rda",
                   "rppa_pi3k.Rda"), 
               field_ids = 
                 c("molecular_subtype",
                   "mutation_names",
                   "hgnc_symbol",
                   "ensembl_gene_id | entrezgene | hgnc_symbol | description",
                   "protein_names"), 
               NAs_allowed =
                 c(TRUE,
                   TRUE,
                   TRUE,
                   TRUE,
                   TRUE), stringsAsFactors = FALSE)
  
  output_dir <- "/Volumes/Ewald/antoinette/TCGA/lumA_pten_loss"
  output_data_sources <- select_output_sources(input_data_sources = input_data_sources, 
                                               source_name_subset = 
                                                 c("MOLECULAR_SUBTYPES", 
                                                   "MUTATIONS", 
                                                   "CNAS", 
                                                   "RNASEQ",
                                                   "RPPA"),
                                               output_file_type = ".tsv")
 
  data_sources <- list(input_dir=input_dir, output_dir=output_dir,
                input_data_sources=input_data_sources, output_data_sources=output_data_sources)
  return(data_sources)
}

setup_TCGA_default_params <- function(df_list){
  paramList <- list()
  
  ##############  Setup parameters to use in logical expression   ##############
  # To see available data_sources, console: str(df_list, max.level = 1)
  # To see available id_colnames, console: names(df_list$<data_source>)[1:4]  (e.g. names(df_list$pam50DF)[1:4])
  # To see available ids, console: df_list$<data_source>$<id_colname> (e.g. df_list$subtype$molecular_subtype)
  paramList$SUBTYPE       <- set_param(df_list = df_list, source_name = "MOLECULAR_SUBTYPES", 
                                       id_colname = "molecular_subtype", id = "pam50", 
                                       data_type = "character")
  paramList$LUMA          <- set_param(df_list = df_list, source_name = "MOLECULAR_SUBTYPES",
                                       id_colname = "molecular_subtype", id = "LumA",
                                       data_type = "numeric")
  paramList$PIK3CA_E545K  <- set_param(df_list = df_list, source_name = "MUTATIONS", 
                                       id_colname = "mutation_names", id = "PIK3CA_p.E545K", 
                                       data_type = "numeric")
  paramList$PIK3CA_H1047R <- set_param(df_list = df_list, source_name = "MUTATIONS",
                                       id_colname = "mutation_names", id = "PIK3CA_p.H1047R", 
                                       data_type = "numeric")
  paramList$PIK3CA_nonsilent  <- set_param(df_list = df_list, source_name = "MUTATIONS", 
                                           id_colname = "mutation_names", id = "PIK3CA_nonsilent", 
                                           data_type = "numeric")
  paramList$AKT1_E17K     <- set_param(df_list = df_list, source_name = "MUTATIONS",
                                       id_colname = "mutation_names", id = "AKT1_p.E17K", 
                                       data_type = "numeric")
  paramList$AKT1_nonsilent  <- set_param(df_list = df_list, source_name = "MUTATIONS", 
                                         id_colname = "mutation_names", id = "AKT1_nonsilent", 
                                         data_type = "numeric")
  paramList$PTEN_nonsilent <- set_param(df_list = df_list, source_name = "MUTATIONS",
                                        id_colname = "mutation_names", id = "PTEN_nonsilent", 
                                        data_type = "numeric")
  paramList$PTEN_CNA      <- set_param(df_list = df_list, source_name = "CNAS",
                                       id_colname = "hgnc_symbol", id = "PTEN",
                                       data_type = "numeric")
  return(paramList)
}

setup_TCGA_default_splits <- function(){
  ##############  Set logical expression ##############
  # Basic unary relations
  # is.na(<x>)  : is NA value
  #
  # Basic binary relations
  # <x> == <y>  : <x> and <y> have the same value
  # <x> <  <y>  : <x> is less than <y>
  # <x> <= <y>  : <x> is less than or equal to <y>
  #
  # Basic unary logical operators:
  # !<X>        :  NOT <X>
  #
  #Basic binary logical operators:
  # <X> && <Y>  : <X> AND <Y> - Both <X> and <Y> have to be true
  # <X> || <Y>  : <X>  OR <Y> - Either <X> or <Y> or both have to be true
  
  ##############  Set split condition ##############
  split <- list()
  # Example 1: choose all samples
  split$all <- "TRUE"
  # Example 2: choose all Luminal A samples
  split$lumA <- "LUMA >= 0.5"
  # Example 3: Choose all lumninal A samples that are wild type in PIK3CA, AKT1 and PTEN
  split$lumA_pik_akt_pten_WT <- "LUMA >= 0.5 && PIK3CA_nonsilent == 0 && AKT1_nonsilent == 0 && PTEN_nonsilent == 0 && PTEN_CNA > -0.25"
  # Example 4: Choose all luminal A samples with different specific mutations
  split$lumA_pik3_E545K <- "LUMA >= 0.5 && PIK3CA_E545K == 1"
  split$lumA_pik3_H1047R <- "LUMA >= 0.5 && PIK3CA_H1047R == 1"
  split$lumA_akt1_E17K <- "LUMA >= 0.5 && AKT1_E17K == 1"
  split$lumA_pten_nonsilent <- "LUMA >= 0.5 && PTEN_nonsilent == 1"
  # Example 5: Choose all luminal A samples with copy number loss in PTEN
  split$lumA_pten_loss <- "LUMA >= 0.5 && PTEN_CNA <= -0.25"
  return(split)
}






