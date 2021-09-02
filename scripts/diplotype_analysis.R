#!/usr/bin/Rscript

source('/home/bioinfo/biosoft/R_scripts/stop_quietly.R')
library(tidyverse)
library(xlsx)

database_dir <- '/home/bioinfo/FARMA/cipic_info_genes/Diplotype-Phenotype'

########################## functions #########################

iscontained <- function(a,b){return(a %in% b)}

######################## BASIC HELP MESSAGE ##########################
# Print if --help is given as an argument
# Print if argument checks fail

files <- commandArgs(trailingOnly=TRUE)
# Debugging
# files <- c("../Variant_Calling/results/FARMA/HCOL10.gatk.norm.decomp.snpeff.snpsift.annot.set.set.genes.alleleDef.report")
files <- c("../Variant_Calling/results/FARMA/HCOL10.allele_func_results.tsv")


file_prefix <- str_split(str_split(files, '/')[[1]][length(str_split(files, '/')[[1]])], '\\.')[[1]][1]

# Help reminder message:
# TODO: edit help message
help_message <- 'USAGE:\n\ndiplotype_test.R [FILE]\n\FILE is a file produced by allele_def_id.R with extension *.report\n\nParameters:\n\n--help to display this message.\n\n'
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

diplo_gene_list <- read_lines(paste0(database_dir,'/Diplotype-Phenotype/available_genes.list'))
allele_def_gene_list <- read_lines(paste0(database_dir,'/Allele_definition/available_genes.list'))
allele_func_gene_list <- read_lines(paste0(database_dir,'/Allele_functionality/available_genes.list'))

######################## INPUT 1 ######################## 
# Input allele func tsv

# minor checks
tsv_file <- grep(files, pattern = ".tsv", value = TRUE)
if(!file.exists(tsv_file)){
  stop(paste(tsv_file, '¡Error! file not found. Exiting.'))
}else if(length(tsv_file) > 1){
  stop('¡Error! only one file procesable at a time.')
}

tsv_file.data <- read_tsv(tsv_file, col_names = TRUE)
if(length(tsv_file.data) == 0){
  stop('¡Error! Empty report file.')
}

######################## QUERY ######################## 

# Identify genes to be queried

genes <- tsv_file.data %>% select(Gene)

# There's no need to check diplotype DB if we don't have multiple alleles present.
# Introduce an empty line
write_lines(file = tsv_file, x = '\n', append = TRUE)

for(gene in unique(genes)){
  
  alleles <- tsv_file.data %>% filter(Gene == gene) %>% select(Allele)
  af <- tsv_file.data %>% filter(Gene == gene) %>% select(AF)

cat(paste0("\n#Checking diplotype effect for gene: ",gene,"\n"))
cat("#################\n")
if(!nrow(alleles) > 1 && AF == 0.5){
  
  cat(paste("#Not enough alleles for dyplotype effect.\n"))                                                                                                                                                                                                                                                                                                                                                
  
}else{
  #Case for a gene we have two alleles
  
  ######################## read DDBB ########################
  
  DiplotypeDB <- read.xlsx2(paste0(database_dir,'/',gene,"_Diplotype_Phenotype_Table.xlsx"), sheetIndex = 1,header = TRUE)
  # clean whitespaces
  DiplotypeDB <- DiplotypeDB %>% mutate(CYP2C9.Diplotype = str_replace_all(CYP2C9.Diplotype, pattern = ' ', replacement = ''))
  # Separate Allele pairs
  # DiplotypeDB <- DiplotypeDB %>% separate(CYP2C9.Diplotype, sep = '/', into = c('CYP2C9.Diplotype.Allele1', 'CYP2C9.Diplotype.Allele2'))
  
  ######################## query DB ########################
  
  if(nrow(alleles) == 1 && AF == 1){
    #Same allele twice
    query_allele <- paste(alleles[[1]][1], alleles[[1]][1], sep = '/')
  }else{
    #Two different alleles
    query_allele <- paste(alleles[[1]][1], alleles[[1]][2], sep = '/')
  }
  
    DiplotypeDB.query <- DiplotypeDB %>% filter(CYP2C9.Diplotype == query_allele)
    cat(paste0('#Allele combination found:\n\n'))
    for(i in c(1:4)){
      cat(paste0(colnames(DiplotypeDB.query)[i],'\t'))
    }
    cat("\n")
    for(i in c(1:4)){
      cat(paste0(DiplotypeDB.query[[i]],'\t'))
    }
    write_tsv(file = tsv_file, x = DiplotypeDB.query, append = TRUE, col_names = TRUE)
}
}
cat(paste0('#Saved in: ',tsv_file,'\n'))
