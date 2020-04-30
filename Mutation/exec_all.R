setwd("~/devel/R/GenImmPhe/")
source("./utility/io_evd.R")
source("./Mutation/TCGA_preproc.R")
source("./Mutation/mut_proc.R")

######################## Mutation calling ####################################

#Files required for processing
mut_filt_file <- "./Mutation/data/mutation_list.txt"

#Intermediate data files for loading if needed
interm_folder <- "~/data/pipeline_interm/2019_05_TCGA_BRCA"
mut_all_file <- file.path(interm_folder, "mut_all.Rda")
mut_pik_akt_pten_file <- file.path(interm_folder, "mut_pik_akt_pten.Rda")

#Load TCGA mutations
mafDF <- download_maf(tumor = "BRCA", downloadDir = "~/temp/", pipelines = c("mutect2", "muse", "varscan2", "somaticsniper"))
mutDF <- mafDF_2_mutation_table(mafDF = mafDF) # Saved in mut_all_file

#Create binary matrix with limited list of mutations
mut_filt <- load_tsv(mut_filt_file)
mutBinDF <- muts_2_binaryMatrix(mutDF = mutDF, mut_filt = mut_filt) # Saved in mut_pik_akt_pten_file