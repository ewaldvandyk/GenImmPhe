source("./cna/cna_preproc.R")

####################### CNA #######################
#Files required for processing
segFile <- "~/data/TCGA/BRCA/cna/gdac.broadinstitute.org_BRCA.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"

#Intermediate data files for loading if needed
cnaDF_file <- "~/data/pipeline_interm/2019_05_TCGA_BRCA/cna_pten.Rda"

#Load segfile and compute log ratios for fixed list of genes
gene_list <- c("PTEN")
cnaDF <- tcga_cnaSeg_2_geneCN(segFile = segFile, gene_names = gene_list, build = "grch37") # Saved in cnaDF_file
