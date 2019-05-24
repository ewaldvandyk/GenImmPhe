source("./utility/io_evd.R")
source("./rppa/rppa_preproc.R")

####################### RPPA processing #######################

#Files required for processing
rppa_file <- "~/data/TCGA/BRCA/Protein/BRCA.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt"


#Intermediate data files for loading if needed
rppa_all_DF_file <- "~/data/pipeline_interm/2019_05_TCGA_BRCA/rppa_all.Rda"
rppa_pi3k_DF_file <- "~/data/pipeline_interm/2019_05_TCGA_BRCA/rppa_pi3k.Rda"

#Load TCGA rppa data
rppaDF <- rppa_tcga_firebrowse_2_df(inFile = rppa_file, sampType = "Primary Solid Tumor") # Saved in rppa_all_DF_file

#Filter proteins
akt_prots <- c("Akt", "Akt_pS473", "Akt_pT308")
apop_prots <- c("FOXO3a", "FOXO3a_pS318_S321", "Bad_pS112")
mtor_prots <- c("mTOR", "mTOR_pS2448")
other_prots <- c("IRS1", "PDK1", "PDK1_pS241", "4E-BP1", "4E-BP1_pS65", "4E-BP1_pT37_T46", "4E-BP1_pT70")
pi3k_prots <- c("PI3K-p110-alpha", "PI3K-p85")
pten_prots <- c("PTEN")
semi_related_prots <- c("p70S6K_pT389")

pi3k_related_prots <- c(akt_prots, apop_prots, mtor_prots, other_prots, pi3k_prots, pten_prots, semi_related_prots)

rppaDF <- filt_rows_DF(df = rppaDF, fieldName = "protein_names", values = pi3k_related_prots) #Saved in rppaDF_file