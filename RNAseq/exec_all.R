setwd("~/devel/R/GenImmPhe/")
source("./utility/io_evd.R")
source("./utility/tcga_proc.R")
source("./utility/data_transform.R")
source("./RNAseq/TCGA_preproc.R")
source("./RNAseq/norm_raw_data.R")
source("./RNAseq/brca_subtyping.R")
source("./RNAseq/DE_analysis.R")


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

##### DE expression with edgeR
setwd(dir = "~/devel/R/GenImmPhe/")
source("./RNAseq/DE_analysis.R")
library(edgeR)
outFolder <- "~/analysis/DE_pathway/Anni/CIP2A_CL/"
load("~/data/rna_seq/BCBL_Srikar/CountData.Rda")

#Setup design matrix
covIDs <- list(CIP2A = c("KD", "Control"), 
               CL = c("HCC1937", "HCC38", "M231", "M436", "M468"))
covDF <- DE_init_covMat(countDF = countDF, colIgnore = 1:6, covIDs = covIDs)
covDF$CIP2A[-c(1:3, 7:9, 13:15, 19:21, 25:27, 31:33)] <- "Control"
covDF$CL[1:6]   <- "HCC1937"
covDF$CL[7:12]  <- "HCC38"
covDF$CL[13:18] <- "M231"
covDF$CL[19:24] <- "M436"
covDF$CL[25:30] <- "M468"

design <- model.matrix(~ 0 + covDF$CIP2A + covDF$CL)
group <- DE_design_2_bioRep(design)

# Setup contrast
contrast <- c(-1,+1, 0, 0, 0, 0, 0)

# Filter and normalize data
edgeRList <- DE_countDF2DGEList(countDF = countDF, rowNameCol = "ID", colIgnore = 1:6, group = group)
edgeRListFilt <- DE_filterLowExpressedGenes(edgeRList, countThresh = 6, minSampsPerGroup = 2)
edgeRListNorm <- calcNormFactors(edgeRListFilt)

# Estimate disrpersion 
edgeRListNorm <- estimateDisp(edgeRListNorm, design)

#Fit model
fit <- glmQLFit(edgeRListNorm, design)

#DE analysis
qlf <- glmQLFTest(fit, contrast=contrast)
qlf_list <- topTags(qlf, n = 1000, sort.by = "p")
tr  <- glmTreat(fit, contrast = contrast, lfc = 1.0)
outFile <- file.path(outFolder, "DE_genes.tsv")
save_tsv(countDF = qlf_list$table, outFile = outFile)

# Go enrichment
go_DF <- edgeR::goana.DGEExact(de = qlf, geneid = qlf$genes$EntrezIDs, species = "Hs")
go_list <- topGO(results = go_DF, number = Inf)
outFile <- file.path(outFolder, "go_enrichment.tsv")
save_tsv(countDF = go_list, outFile = outFile)


# KEGG pathway analysis
kegg_DF <- edgeR::kegga.DGEExact(de = qlf, geneid = qlf$genes$EntrezIDs, species = "Hs")
kegg_list <- topKEGG(results = kegg_DF, number = Inf)
outFile <- file.path(outFolder, "kegg_pathway.tsv")
save_tsv(countDF = kegg_list, outFile = outFile)

# Gene set testing
romer_DF <- edgeR::romer.DGEList(y = edgeRListNorm, index = c(1,2,3))
