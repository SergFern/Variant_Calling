#!/usr/bin/Rscript

library(tidyverse)
library(xlsx)

# Find home directory:
home <- paste0(str_split(getwd(),'/')[[1]][1:3], collapse = '/')
source(paste0(home,'/biosoft/R_scripts/aid_functions.R'))

database_dir <- paste0(home,'/FARMA/cipic_info_genes')
# TODO: write output to .log
######################## BASIC HELP MESSAGE ##########################
# Print if --help is given as an argument
# Print if argument checks fail

files <- commandArgs(trailingOnly=TRUE)
# Debugging
# files <- c("/home/bioinfo/Variant_Calling/results/FARMA/HCOL10.gatk.norm.decomp.snpeff.snpsift.set.genes.alleleDef.log")

file_prefix <- str_split(str_split(files, '/')[[1]][length(str_split(files, '/')[[1]])], '\\.')[[1]][1]

# Help reminder message:
help_message <- 'USAGE:\n\nmatch_allele.R [log FILE]\n\nlog FILE is a file produced by allele_def_id.R with extension *.log\n\nParameters:\n\n--help to display this message.\n\n'
help_variants <- c("--help", "-help")
payload <- paste0("|",paste0(files,collapse = '|'),"|")

if(any(str_detect(payload, pattern = help_variants))){
  cat(help_message)
  stop_quietly()
}else if(isFALSE(str_detect(payload, pattern = '.log'))){
  cat(help_message)
  stop('\nlog_FILE (*.log) is missing\n')
}

######################## DDBB ##########################

allele_func_gene_list <- read_lines(paste0(database_dir,'/Allele_functionality/available_genes.list'))

######################## INPUT 1 ######################## 

# minor checks
log_file <- grep(files, pattern = ".log", value = TRUE)
if(!file.exists(log_file)){
  stop(paste(log_file, '¡Error! file not found. Exiting.'))
}else if(length(log_file) > 1){
  stop('¡Error! only one file procesable at a time.')
}

log_file.lines <- read_lines(log_file, skip_empty_rows = TRUE)
if(length(log_file.lines) == 0){
  stop('¡Error! Empty log file.')
  }

######################### cleanup #######################
#Matches DO NOT start their line with a #
matches <- grep(log_file.lines, pattern = '#', fixed = TRUE, value = TRUE, invert = TRUE)

data <- str_split(str_trim(matches),n = 4, pattern = "\\|")

######################### OPERATIONS ######################### 

if(length(matches) == 0){#no matches
  stop('No mathes found in log file.')
}else{
  n = 0
  tibble_list <- NULL
  for(i in data){ # skips first line - header
    n = n + 1
    # Identify data of interest in log file OBSOLETE
    variants <- i[1]
    allele <- i[length(i)-2]
    AF <- i[length(i)]
    gene <- i[length(i)-1]
    
    ######################### REFERENCE TABLE ######################### 
    # 1st check we have gene in DDBB
    
    # Identify ddbb file
    func_file <- dir(paste0(database_dir,'/Allele_functionality'),full.names = TRUE)[str_detect(dir(paste0(database_dir,'/Allele_functionality'),full.names = TRUE), pattern = gene)]
    func_file.data <- tibble(read.xlsx2(func_file, sheetIndex = 1,startRow = 2, header = TRUE, colClasses = 'character'))
 
    ######################### QUERY ######################### 
    
    # Match our prerequisites to data
    allele_location <- match(allele, func_file.data[,1] %>% pull())
    func_file.data.query <- func_file.data[allele_location,]
    functional_info.index <- grep(colnames(func_file.data.query), pattern = glob2rx('*Clinical*Function*Required*'))
    clinical_function <- func_file.data.query[,functional_info.index] %>% pull()
    
    strenth_of_evidence.index <- grep(colnames(func_file.data.query), pattern = glob2rx('*Strength*of*Evidence*'))
    strenth_of_evidence <- func_file.data.query[,strenth_of_evidence.index] %>% pull()
    
    # Store and fuse multiple matches
    tmp_store <- tibble(Gene = gene, Allele = allele, Clinical_function = clinical_function, Variants = variants, AF = AF, Strength_of_Evidence = strenth_of_evidence)
    tibble_list <- append(tibble_list, list(tmp_store))
    
  }
  
  results <- bind_rows(tibble_list)
  
}

######################### SAVE RESULTS ######################### 

write_tsv(results, col_names = TRUE, file = paste0(getwd(),'/',file_prefix,'.allele_func_results.tsv'))

#.log is a plain text file, we basically have to look for matching characters.

