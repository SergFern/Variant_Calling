#!/usr/bin/Rscript

library(tidyverse)
library(vcfR)
library(xlsx)
source('/home/bioinfo/biosoft/R_scripts/stop_quietly.R')

database_dir <- '/home/bioinfo/FARMA/cipic_info_genes'

######################## BASIC HELP MESSAGE ##########################
# Print if --help is given as an argument
# Print if argument checks fail

files <- commandArgs(trailingOnly=TRUE)
#debugging
#files <- c("my-results/annotation/HCOL10.gatk.norm.decomp.snpeff.snpsift.annot.vcf", "my-results/annotation/HCOL.ID.genes2.tsv")

# Help reminder message:
help_message <- 'USAGE:\n\nallele_def_id.R [VCF_FILE] [TSV_FILE]\n\nVCF_FILE is a variant calling file with extension *.vcf\nTSV_FILE is a file with extension *.tsv, direct output of FARMA.nf process "extract_info"\n\nParameters:\n\n--help to display this message.\n\n'
help_variants <- c("--help", "-help")
payload <- paste0("|",paste0(files,collapse = '|'),"|")

if(any(str_detect(payload, pattern = help_variants))){
  cat(help_message)
  stop_quietly()
}else if(isFALSE(str_detect(payload, pattern = '.vcf'))){
  cat(help_message)
  stop('\nVCF_FILE is missing\n')
}else if(isFALSE(str_detect(payload, pattern = '.tsv'))){
  cat(help_message)
  stop('\nTSV_FILE is missing\n')
}

######################## DDBB ##########################
# TODO: Introduce DDBB as a parameter

diplo_gene_list <- read_lines(paste0(database_dir,'/Diplotype-Phenotype/available_genes.list'))
allele_def_gene_list <- read_lines(paste0(database_dir,'/Allele_definition/available_genes.list'))
allele_func_gene_list <- read_lines(paste0(database_dir,'/Allele_functionality/available_genes.list'))

######################## INPUT 1 ######################## 

vcf_file <- grep(files, pattern = ".vcf", value = TRUE)
if(!file.exists(vcf_file)){
  
  stop(paste(vcf_file, '¡Error! file not found. Exiting.'))
  
}

cat("\n######################\n")
cat(paste("## Procesing file:",vcf_file))
cat("\n######################\n\n")

## VCF file load not really used, commented out until I feel comfortable enough to delete it.
# vcf_file <- read.vcfR('/home/bioinfo/Variant_Calling/my-results/annotation/HCOL10.gatk.norm.decomp.snpeff.snpsift.annot.vcf')
# 
# # FOR loop for every sample in a directory? Or would it be handled by nextflow sample by sample?
# vcf_file.df <- vcfR2tidy(vcf_file)
cat("\n######################\n")

######################## INPUT 2 ######################## 
# TSV_FILE is a file with extension *.tsv, direct output of FARMA.nf process "extract_info"
#####.

tsv_file <- grep(files, pattern = ".tsv", value = TRUE)
extracted_Data <- read_tsv(file = tsv_file, col_names = TRUE)

genes <- unique(extracted_Data$`ANN[0].GENE`)
cat(paste("# Genes found in sample and present in database:", paste(genes, collapse = ', ')))

######################## QUERY ##########################
# TODO: Not all genes have entries with rsIds on all columns, we probably have to use REF and ALT or nucleotide change (87G>C)

allele_def_data <- dir(paste0(database_dir,'/Allele_definition'), full.names = TRUE, pattern = 'allele')
cat("\n######################\n")
cat("\n")

for(gene in pull(unique(extracted_Data[4]))){
  cat("\n######################\n")
  cat(paste('## Variants for gene', gene, ':\n'))
  #Prepare reference table
  def_table <- tibble(read.xlsx2(grep(allele_def_data, pattern = gene, value = TRUE), 1))
  rsID_row <- which(def_table[[1]] == 'rsID')
  def_table.data <- tibble(read.xlsx2(grep(allele_def_data, pattern = gene, value = TRUE), sheetIndex = 1, startRow = rsID_row+1))
    # clean data
    colnames(def_table.data) <- str_remove_all(colnames(def_table.data), pattern = '\\.')
  def_table.header <- tibble(read.xlsx2(grep(allele_def_data, pattern = gene, value = TRUE), sheetIndex = 1, endRow = rsID_row+1))
  #reference <- which(def_table[[1]] == 'Reference')
  #def_table <- def_table %>% slice(reference:nrow(def_table)) %>% filter(rsID != "")
  # Find rsIds of our sample
  rsID_list <- extracted_Data %>% filter(`ANN[0].GENE` == gene) %>% select(ID) %>% pull()
  cat(paste('#',rsID_list, collapse = ', '))
  cat("\n")
  # Process the rsIDs found in def_table to simplify query 
  binari_info_sample <- def_table.data %>% select(all_of(rsID_list)) %>% transmute_all(~if_else(. == "",true = 0, false = 1))
    # Include rsID column
    binari_info_sample <- cbind(def_table.data %>% select(rsID), binari_info_sample)
  binari_info_table <- def_table.data %>% select(-rsID) %>% transmute_all(~if_else(. == "",true = 0, false = 1))
    # Include rsID column
    binari_info_table <- cbind(def_table.data %>% select(rsID), binari_info_table)
  # Select from reference table only alleles of interest
  def_table_selected <- binari_info_sample %>% filter(across(c(2,ncol(binari_info_sample)), ~. == 1 ))
  cat(paste('## Possible alleles:', paste(def_table_selected$rsID, collapse = ', ')))
  cat("\n\n")
  # Check all required variants for allele are present
  def_table_alleles <- tibble(binari_info_table[match(def_table_selected$rsID, binari_info_table$rsID),])
    # Alleles
    alleles <- def_table_alleles %>% select(rsID) %>% pull()
    # All rsIDs in definition table
    rsID_list_def <- colnames(binari_info_table %>% select(-1))
    
    # Fail counter
    fail = 0
  for(i in 1:length(alleles)){
  
    allele_variants <- def_table_alleles %>% slice(i) %>% select(-1) %>% unlist(., use.names = FALSE)
    allele_variants <- rsID_list_def[allele_variants == 1]
    similarity_coeff <- sum(rsID_list %in% allele_variants)/length(allele_variants)
    if(similarity_coeff == 1){
      # If true include in report
      cat("#Match:\n")
      cat("#------------\n")
      cat(paste('#Variants\tAllele\tGene\n'))
      cat(paste(paste(rsID_list, collapse = '|'),'\t',alleles[i],'\t',gene,'\n'))
      cat("#-------------\n\n")
      
    }else if(similarity_coeff < 0.5){
      # If not all variants are found, which are missing?
      fail = fail + 1
      missing_variants_for_allele = allele_variants[!allele_variants %in% rsID_list]
      if(fail == length(alleles)){
      cat(paste('# No allele match found for the combination of variants.\n\n'))}
    }else if(similarity_coeff >= 0.5){
      
      cat(paste0('#Partial match with allele ', alleles[i], ' with ', similarity_coeff*100, '% common variants.\n'))
      
    }
     
  }
  
}

######################## OUTPUT #########################

