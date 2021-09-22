#!/usr/bin/Rscript

library(tidyverse)
library(vcfR)
library(xlsx)

# Find home directory:
home <- paste0(str_split(getwd(),'/')[[1]][1:3], collapse = '/')
source(paste0(home,'/biosoft/R_scripts/aid_functions.R'))

database_dir <- paste0(home,'/FARMA/cipic_info_genes')

######################## BASIC HELP MESSAGE ##########################
# Print if --help is given as an argument
# Print if argument checks fails

files <- commandArgs(trailingOnly=TRUE)
#debugging
# files <- c("~/Variant_Calling/output/FARMA/HCOL10.gatk.norm.decomp.snpeff.snpsift.genes.vcf", 
           # "/home/bioinfo/Variant_Calling/output/FARMA/HCOL10.gatk.norm.decomp.snpeff.snpsift.genes.ID.tsv")

# Help reminder message:
help_message <- 'USAGE:\n\nallele_def_id.R [VCF_FILE] [TSV_FILE] > output.log\n\nVCF_FILE is a variant calling file with extension *.vcf\nTSV_FILE is a file with extension *.tsv, direct output of FARMA.nf process "extract_info"\n\nParameters:\n\n--help to display this message.\n\n'
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

######################## INPUT 1 ######################## 

vcf_file <- grep(files, pattern = ".vcf", value = TRUE)
if(!file.exists(vcf_file)){
  
  stop(paste(vcf_file, '¡Error! file not found. Exiting.'))
  
}

cat("\n######################\n")
cat(paste("## Procesing file:",vcf_file))
cat("\n######################\n\n")

######################## INPUT 2 ######################## 
# TSV_FILE is a file with extension *.tsv, direct output of FARMA.nf process "extract_info"
#####.

tsv_file <- grep(files, pattern = ".tsv", value = TRUE)
extracted_Data <- read_tsv(file = tsv_file, col_names = TRUE)

# Clean extracted data (Column `ANN[*].GENE` can have many GENE_NAME per row, select only one.)
genes_lst_multiple <- str_split(extracted_Data$`ANN[*].GENE`,pattern = ',')
# Need a list of available genes
list_of_available_genes <- read_lines(paste0(home,'/FARMA/cipic_info_genes/Allele_definition/available_genes.list'))
genes_lst_single <- c()
for(i in 1:length(genes_lst_multiple)){
  
  gene_available <- genes_lst_multiple[[i]][genes_lst_multiple[[i]] %in% list_of_available_genes][1]
  genes_lst_single <- c(genes_lst_single, gene_available)
  
}
extracted_Data <- extracted_Data %>% mutate(`ANN[*].GENE` = genes_lst_single)
# Substitute tsv file with corrected gene data

write_tsv(path = tsv_file, x = extracted_Data, col_names = TRUE)

# Build HGVS column
extracted_Data <- extracted_Data %>% mutate(HGVS = paste0('g.',POS,REF,'>',ALT))

genes <- unique(extracted_Data$`ANN[*].GENE`)
cat(paste("# Genes found in sample and present in database:", paste(genes, collapse = ', ')))

######################## REFERENCE TABLE ##########################

allele_def_data <- dir(paste0(database_dir,'/Allele_definition'), full.names = TRUE, pattern = 'allele')
cat("\n######################\n")
cat("\n")

for(gene in unique(extracted_Data$`ANN[*].GENE`)){
  cat("\n######################\n")
  cat(paste('## Variants for gene', gene, ':\n'))
  #Prepare reference table DATA  ####
  def_table <- tibble(read.xlsx2(grep(allele_def_data, pattern = gene, value = TRUE), 1))
    # clean data
    def_table <- remove_empty(def_table)
  rsID_row <- which(def_table[[1]] == 'rsID')
  def_table.data <- tibble(read.xlsx2(grep(allele_def_data, pattern = gene, value = TRUE), sheetIndex = 1, startRow = rsID_row+1))
    # clean data
    colnames(def_table.data) <- str_remove_all(colnames(def_table.data), pattern = '\\.')
    # remove empty columns
    def_table.data <- remove_empty(def_table.data)
    # remove empty rows
    def_table.data <- def_table.data %>% filter(rsID != '')
    #Prepare reference table HEADER ####
    def_table.header <- tibble(read.xlsx2(grep(allele_def_data, pattern = gene, value = TRUE), sheetIndex = 1, endRow = rsID_row+1))  %>% transmute_all(.funs = ~as.character(.))
    # clean header
    def_table.header <- remove_empty(def_table.header)
    
    #Substitute def_table.data non-rsID data for HGVS
      # Index non-rsids if present
    index_nonrsIDs <- which(!str_detect(colnames(def_table.data), 'rs'))
    if(length(index_nonrsIDs)!=0){
      HGVS_db <- def_table.header[3,index_nonrsIDs]
      colnames(def_table.data)[index_nonrsIDs] <- HGVS_db
    }
    
    #Remove columns with Structural Variantion information
    def_table.data <- def_table.data %>% select(-which(str_detect(def_table.header, pattern = 'Structural Variation')))
  


######################## QUERY ##########################
    
  # Find queries of our sample for the present gene
  extracted_Data.filtered <- extracted_Data %>% filter(`ANN[0].GENE` == gene)
  # Build query column
  extracted_Data.filtered <- extracted_Data.filtered %>% mutate(list_of_queries = if_else(condition = str_detect(ID,'rs')|ID %in% colnames(def_table.data),true = ID, false = HGVS))
  
  #Building query data structure
  rsID_list <- extracted_Data.filtered %>% select(ID) %>% pull()
  HGVS <- paste0('g.',extracted_Data.filtered$POS, extracted_Data.filtered$REF,'>', extracted_Data.filtered$ALT)
  #I want to use a rsID if it is available AND it is found in the DB, if not i want to use HGVS
  query.tb <- tibble()
  query.tb <- tibble(rsID_list,HGVS) %>% mutate(list_of_queries = if_else(condition = str_detect(rsID_list,'rs')|rsID_list %in% colnames(def_table.data),true = rsID_list, false = HGVS)) %>%
    bind_cols(POS = extracted_Data.filtered$POS,REF = extracted_Data.filtered$REF,ALT = extracted_Data.filtered$ALT)
  
  cat(paste('#',extracted_Data.filtered$list_of_queries, collapse = ', '))
  cat("\n")
  # Process the rsIDs found in def_table to simplify query
  # rsIDs must be present though
  
  #########################
  # column ID can be: rsID present in DB, rsID NOT present in DB, a dot (.)
  #########################
  
  binari_info_sample <- def_table.data %>% select(all_of(extracted_Data.filtered$list_of_queries)) %>% transmute_all(~if_else(. == "",true = 0, false = 1))
    # Include rsID column
    binari_info_sample <- cbind(def_table.data %>% select(rsID), binari_info_sample)
  binari_info_table <- def_table.data %>% select(-rsID) %>% transmute_all(~if_else(. == "",true = 0, false = 1))
    # Include rsID column
    binari_info_table <- cbind(def_table.data %>% select(rsID), binari_info_table)
  # Select from reference table only alleles of interest
  def_table_selected <- binari_info_sample %>% filter_all(any_vars(. == 1))
  
  entries_of_interest <- match(def_table_selected$rsID, binari_info_table$rsID)
  
  cat(paste('## Possible alleles:', paste(def_table_selected$rsID, collapse = ', ')))
  cat("\n\n")
  # Check all required variants for allele are present
    def_table_alleles <- tibble(binari_info_table[entries_of_interest,])
    # required alleles: Alleles
    alleles <- def_table_alleles %>% select(rsID) %>% pull()
    # All rsIDs in definition table
    rsID_list_def <- colnames(binari_info_table %>% select(-1))
    
    # Fail counter
    fail = 0
    match = 0
  for(i in 2:length(alleles)){
  
    allele_variants <- def_table_alleles %>% slice(i) %>% select(-1) %>% unlist(., use.names = FALSE)
    allele_variants <- rsID_list_def[allele_variants == 1]

    similarity_coeff <- sum(extracted_Data.filtered$list_of_queries %in% allele_variants)/length(allele_variants)
    if(similarity_coeff == 1){
      # If true include in log
      # Was there already a match? DO not print #Match: again
      if(match == 0){
        
      cat("#Match:\n")
      cat("#------------\n")
      cat(paste('#Variants|Allele|Gene|AF\n'))
      
      }
      cat(paste0(paste0(extracted_Data.filtered$list_of_queries[extracted_Data.filtered$list_of_queries %in% allele_variants], collapse = ','),'|',alleles[i],'|',gene,'|',extracted_Data.filtered[extracted_Data.filtered$list_of_queries %in% allele_variants,]$AF,'\n'))
      match = 1
      
    }else if(similarity_coeff < 0.5){
      # If not all variants are found, which are missing?
      fail = fail + 1
      missing_variants_for_allele = allele_variants[!allele_variants %in% extracted_Data.filtered$list_of_queries]
      if(fail == length(alleles)){
      cat(paste('# No allele match found for the combination of variants.\n\n'))}
    }else if(similarity_coeff >= 0.5){
      
      cat(paste0('#Partial match with allele ', alleles[i], ' with ', similarity_coeff*100, '% common variants.\n'))
      
    }
  }
    if(match == 1){
      cat("#-------------\n\n")
    }
    match = 0
}

######################## OUTPUT #########################
# Redirect output with ">"