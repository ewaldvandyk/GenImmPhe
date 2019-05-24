source("./user_tools/data_split/split_methods.R", local = TRUE)

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
  
  output_dir <- "~/temp/wild_types"
  output_data_sources <- select_output_sources(input_data_sources = input_data_sources, 
                                               source_name_subset = 
                                                 c("RNASEQ",
                                                   "RPPA"),
                                               output_file_type = ".tsv")
 
  data_sources <- list(input_dir=input_dir, output_dir=output_dir,
                input_data_sources=input_data_sources, output_data_sources=output_data_sources)
  return(data_sources)
}


setup_split <-function(df_list){
  paramList <- list()
  
  ##############  Setup parameters to use in logical expression   ##############
  # To see available data_sources, console: str(df_list, max.level = 1)
  # To see available id_colnames, console: names(df_list$<data_source>)[1:4]  (e.g. names(df_list$pam50DF)[1:4])
  # To see available ids, console: df_list$<data_source>$<id_colname> (e.g. df_list$subtype$molecular_subtype)
  paramList$SUBTYPE       <- set_param(df_list = df_list, source_name = "MOLECULAR_SUBTYPES", 
                                       id_colname = "molecular_subtype", id = "pam50", 
                                       data_type = "character")
  paramList$PIK3CA_E545K  <- set_param(df_list = df_list, source_name = "MUTATIONS", 
                                       id_colname = "mutation_names", id = "PIK3CA_p.E545K", 
                                       data_type = "numeric")
  paramList$PIK3CA_H1047R <- set_param(df_list = df_list, source_name = "MUTATIONS",
                                          id_colname = "mutation_names", id = "PIK3CA_p.H1047R", 
                                          data_type = "numeric")
  paramList$AKT1_E17K     <- set_param(df_list = df_list, source_name = "MUTATIONS",
                                       id_colname = "mutation_names", id = "AKT1_p.E17K", 
                                       data_type = "numeric")
  paramList$PTEN_MUT      <- set_param(df_list = df_list, source_name = "MUTATIONS",
                                       id_colname = "mutation_names", id = "PTEN_ANY", 
                                       data_type = "numeric")
 
  paramList$PTEN_CNA      <- set_param(df_list = df_list, source_name = "CNAS",
                                       id_colname = "hgnc_symbol", id = "PTEN",
                                       data_type = "numeric")
  
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
  
  #split_condition <- "SUBTYPE == 'LumA' && PIK3CA_E545K == 1 && (PTEN_MUT == 1 || PTEN_CNA <= -0.25)"
  split_condition <- "SUBTYPE == 'LumA' && PIK3CA_E545K == 0 && PIK3CA_H1047R == 0 && AKT1_E17K == 0 && PTEN_MUT == 0 && PTEN_CNA > -0.25"
  #split_condition <- "TRUE"
  
  ##############  Set logical expression   ##############
  parseStruct <- list(params = paramList, parsePhrase = split_condition)
  return(parseStruct)
}



