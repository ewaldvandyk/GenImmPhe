setwd("~/devel/R/GenImmPhe")
source("./utility/plot_scripts.R")
source("./data_integration/df_processing.R")

######################### Heatmap of RPPA data #########################

#Files required for processing
data_dir <- "~/data/pipeline_interm/2019_05_TCGA_BRCA"
mut_file <- file.path(data_dir, "mut_pik_akt_pten.Rda")
cna_file <- file.path(data_dir, "cna_pten.Rda")
rppa_file <-file.path(data_dir, "rppa_pi3k.Rda") 

#Load required data
load(cna_file)
load(mut_file)
load(rppa_file)

#Create plotting column side colours
orderColList <- colnames(rppaDF)[-1]
colSideColorsMUT <- create_binary_sideColors(binaryDF = mutBinDF, orderColList = orderColList)
colSideColorsCNA <- create_cont_sideColors(contDF = cnaDF, orderColList = orderColList)
colSideColors <- cbind(colSideColorsMUT, colSideColorsCNA)

#Create heatmap
genHeatMap(df = rppaDF, rownameField = "protein_names", 
           colSideColors = colSideColors)

######################### Integrate data frames #########################

#Files required for processing
data_dir <- "~/data/pipeline_interm/2019_05_TCGA_BRCA"
sub_file <- file.path(data_dir, "molecular_subtypes.Rda")
mut_file <- file.path(data_dir, "mut_pik_akt_pten.Rda")
cna_file <- file.path(data_dir, "cna_pten.Rda")
rna_file <- file.path(data_dir, "voomRNAseq.Rda")
rppa_file <-file.path(data_dir, "rppa_pi3k.Rda") 

#Load required data
load(sub_file)
load(mut_file)
load(cna_file)
load(rna_file)
load(rppa_file)


# Allign data frames
df_list <- list(subtype = pam50DF, mut = mutBinDF, cna = cnaDF, rna = voomRNASeqDF, rppa = rppaDF)
nonDataFields <- c("molecular_subtype", "mutation_names", "hgnc_symbol", "ensembl_gene_id", "entrezgene",
                   "description", "protein_names")
df_alligned_list <- allignDFs(df_list = df_list, nonDataFieldNames = nonDataFields, NAsAllowed = c(2,3,5))

