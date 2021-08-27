#!/usr/bin/Rscript

#TODO: make a nextflow process and module

source('/home/bioinfo/biosoft/R_scripts/stop_quietly.R')
library(tidyverse)
library(xlsx)

database_dir <- '/home/bioinfo/FARMA/cipic_info_genes'

######################## BASIC HELP MESSAGE ##########################
# Print if --help is given as an argument
# Print if argument checks fail

files <- commandArgs(trailingOnly=TRUE)
# Debugging
# files <- c("../Variant_Calling/results/results/FARMA/HCOL.gatk.norm.decomp.snpeff.annot.set.set.genes.alleleDef.report")

file_prefix <- str_split(str_split(files, '/')[[1]][length(str_split(files, '/')[[1]])], '\\.')[[1]][1]

# Help reminder message:
help_message <- 'USAGE:\n\nmatch_allele.R [REPORT FILE]\n\nREPORT FILE is a file produced by allele_def_id.R with extension *.report\n\nParameters:\n\n--help to display this message.\n\n'
help_variants <- c("--help", "-help")
payload <- paste0("|",paste0(files,collapse = '|'),"|")

if(any(str_detect(payload, pattern = help_variants))){
  cat(help_message)
  stop_quietly()
}else if(isFALSE(str_detect(payload, pattern = '.report'))){
  cat(help_message)
  stop('\nVCF_FILE is missing\n')
}

######################## DDBB ##########################
# TODO: Introduce DDBB as a parameter

diplo_gene_list <- read_lines(paste0(database_dir,'/Diplotype-Phenotype/available_genes.list'))
allele_def_gene_list <- read_lines(paste0(database_dir,'/Allele_definition/available_genes.list'))
allele_func_gene_list <- read_lines(paste0(database_dir,'/Allele_functionality/available_genes.list'))

######################## INPUT 1 ######################## 

# minor checks
report_file <- grep(files, pattern = ".report", value = TRUE)
if(!file.exists(report_file)){
  stop(paste(report_file, '¡Error! file not found. Exiting.'))
}else if(length(report_file) > 1){
  stop('¡Error! only one file procesable at a time.')
}

report_file.lines <- read_lines("../Variant_Calling/results/results/FARMA/HCOL10.gatk.norm.decomp.snpeff.snpsift.annot.set.set.genes.alleleDef.report",
                                skip_empty_rows = TRUE)
if(length(report_file.lines) == 0){
  stop('¡Error! Empty report file.')
  }

######################### cleanup #######################

matches <- grep(report_file.lines, pattern = '#', fixed = TRUE, value = TRUE, invert = TRUE)

data <- str_split(str_trim(str_replace(matches, "\t ", ""))," ")

######################### OPERATIONS ######################### 

if(length(matches) == 0){#no matches
  stop('No mathes found in report file.')
}else{
  n = 0
  tibble_list <- NULL
  for(i in data){ # skips first line - header
    n = n + 1
    # Identify data of interest in report file
    gene <- i[length(i)]
    allele <- i[length(i)-1]
    variants <- i[1]
    
    ######################### QUERY ######################### 
    # 1st check we have gene in DDBB
    if(!gene %in% allele_func_gene_list){break}
    
    # Identify ddbb file
    func_file <- grep(dir(paste0(database_dir,'/Allele_functionality/'),full.names = TRUE), pattern = gene, value = TRUE)
    func_file.data <- read.xlsx2(func_file, sheetIndex = 1,startRow = 2, header = TRUE)
    
    # Match our prerequisites to data
    allele_location <- match(allele, func_file.data[,1])
    func_file.data.query <- func_file.data[2,]
    functional_info.index <- grep(colnames(func_file.data.query), pattern = glob2rx('*Clinical*Functional*'))
    clinical_function <- func_file.data.query[,functional_info.index]
    

    tmp_store <- tibble(Gene = gene, Allele = allele, Clinical_function = clinical_function, Variants = variants)
    
    tibble_list <- append(tibble_list, list(tmp_store))
    
  }
  
  results <- bind_rows(tibble_list)
  
}

######################### SAVE RESULTS ######################### 

write_tsv(results, col_names = TRUE, file = paste0(getwd(),'/',file_prefix,'.allele_func_results.tsv'))

#.report is a plain text file, we basically have to look for matching characters.

