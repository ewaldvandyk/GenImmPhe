source("./utility/io_evd.R")
source("./utility/tcga_proc.R")
source("./RNAseq/TCGA_preproc.R")
source("./RNAseq/norm_raw_data.R")
source("./RNAseq/brca_subtyping.R")
source("./utility/data_transform.R")

####################### RNA-seq preprocessing #######################

#Files required for processing
gdc_inDir = "~/data/TCGA/BRCA/geneExp/RNAseq/data/" # Directory containing individual sample files with raw read counts
sampFile <- "~/data/TCGA/BRCA/geneExp/RNAseq/gdc_sample_sheet.2018-11-28.tsv" # File relating sample file names to bar-codes
perplexEsngFile <- "~/data/TCGA/BRCA/geneExp/RNAseq/dataOut/perplexed_ensg.tsv" # File mapping ensg ids to unique hgnc and entrez for genes with multiple ensg entries in biomart
perplexHgncFile <- "~/data/TCGA/BRCA/geneExp/RNAseq/dataOut/perplexed_hgnc.tsv" # File mapping ensg ids to unique hgnc and entrez for genes with multiple hgnc entries in biomart
perplexEntrezFile <- "~/data/TCGA/BRCA/geneExp/RNAseq/dataOut/perplexed_entrez.tsv" #File mapping ensg ids to unique hgnc and entrez for genes with multiple entrez entries in biomart

#Intermediate data files for loading if needed
interm_folder <- "~/data/pipeline_interm/2019_05_TCGA_BRCA"
rawRNAseq_basic_file <- file.path(interm_folder, "rawCountsRNAseq_basic.Rda")
bioMart_genes_file <- file.path(interm_folder, "bioMart_genes.Rda")
rawRNAseq_file <- file.path(interm_folder, "rawCountsRNAseq.Rda")
rawRNAseq_final_file <- file.path(interm_folder, "rawCountsRNAseq_final.Rda")

# RNA-seq extraction from TCGA sample files
rawRNAseqDF <- gdc_TCGA_2_rnaSeqraw(gdc_inDir = gdc_inDir, sampFile = sampFile) # Saved in rawRNAseq_basic_file

# Add Hgnc symbols, entrez ids and descriptions from ensgs
g_list <- ensg_2_geneSymbols(rawRNAseqDF) # Saved in rawRNAseq_basic_file
perplexEnsgDF <- load_tsv(inFile = perplexEsngFile)
perplexHgncDF <- load_tsv(inFile = perplexHgncFile)
perplexEntrezDF <- load_tsv(inFile = perplexEntrezFile)
rawRNAseqDF <- addGeneNames2countData(countDF = rawRNAseqDF, g_list = g_list,
                                      ambEnsgDF = perplexEnsgDF, 
                                      ambHgncDF = perplexHgncDF, 
                                      ambEntrezDF = perplexEntrezDF) ## Saved in rawRNAseq_file

#Sample filters 
#Filter samples for primary tumors
bar_codes <- colnames(rawRNAseqDF)[-(1:4)]
primaryI <- c(rep(TRUE, 4), filtSampleType(bar_codes = bar_codes, sampType = "Primary Solid Tumor"))
rawRNAseqDF <- rawRNAseqDF[,primaryI]

#Filter duplicate samples based on interactive read count comparisons and remove non sample id tags from bar_codes
rawRNAseqDF <- filtDupRNASamples_interactive(countDF = rawRNAseqDF, nonDataFields = 1:4)
colnames(rawRNAseqDF)[-(1:4)] <- barcodes2short(colnames(rawRNAseqDF)[-(1:4)], id_col = 3) # Saved in rawRNAseq_final_file

####################### RNA-seq normalization #######################

#Files required for processing
pam50_genes_file <- "~/data/TCGA/BRCA/geneExp/RNAseq/dataOut/pam50_genes.tsv"

#Intermediate data files for loading if needed
interm_folder <- "~/data/pipeline_interm/2019_05_TCGA_BRCA"
voomRNAseq_file <- "~/data/pipeline_interm/2019_05_TCGA_BRCA/voomRNAseq.Rda"
molecular_subtype_file <- "~/data/pipeline_interm/2019_05_TCGA_BRCA/molecular_subtypes.Rda"

#Molecular subtyping
pam50_genes <- load_tsv(pam50_genes_file)$genes
voomRNASeqDF <- voom_normalize(countDF = rawRNAseqDF, colIgnore = 1:4, hgnc_keep = pam50_genes) # Saved in voomRNAseq_file
pam50DF <- classify_brca_molecular_suptypes(countDF = voomRNASeqDF, colIgnore = 1:4, ensg_col = "ensembl_gene_id", entrez_col = "entrezgene")
pam50DF <- transpose_df(df = pam50DF, rowNameCol = "sampName", colTypeName = "molecular_subtype") # Saved in molecular_subtype_file


