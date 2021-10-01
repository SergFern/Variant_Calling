#!/usr/bin/Rscript

library(tidyverse)
library(xlsx)

# Find home directory:
home <- paste0(str_split(getwd(),'/')[[1]][1:3], collapse = '/')
source(paste0(home,'/biosoft/R_scripts/aid_functions.R'))

database_dir <- paste0(home,'/FARMA/cipic_info_genes')

######################## BASIC HELP MESSAGE ##########################
# Print if --help is given as an argument
# Print if argument checks fail

files <- commandArgs(trailingOnly=TRUE)
# Debugging
# files <- c("../Variant_Calling/results/FARMA/HCOL10.allele_func_results.tsv")

# Help reminder message:
# TODO: edit help message
help_message <- 'USAGE:\n\ndiplotype_test.R [TSV_FILE]\nTSV_FILE is a file produced by match_allele.R with extension *.tsv\n\nParameters:\n\n--help to display this message.\n\n'
help_variants <- c("--help", "-help")
payload <- paste0("|",paste0(files,collapse = '|'),"|")

if(any(str_detect(payload, pattern = help_variants))){
  cat(help_message)
  stop_quietly()
}else if(isFALSE(str_detect(payload, pattern = '.tsv'))){
  cat(help_message)
  stop('\nTSV FILE is missing\n')
}

######################## INPUT 1 ######################## 
# Input allele func tsv

# minor checks
tsv_file <- grep(files, pattern = ".tsv", value = TRUE)
if(!file.exists(tsv_file)){
  stop(paste(tsv_file, '¡Error! file not found. Exiting.'))
}else if(length(tsv_file) > 1){
  stop('¡Error! only one file procesable at a time.')
}

tsv_file.data <- read_tsv(tsv_file, col_names = TRUE, comment = "#")
file_as_lines <- read_lines(tsv_file)
if(length(tsv_file.data) == 0){
  stop('\n######### ERROR ###########\n#Empty tsv file.\n')
}

######################## QUERY ######################## 

# Identify genes to be queried

genes <- tsv_file.data %>% select(Gene)

# There's no need to check diplotype DB if we don't have multiple alleles present.
# Introduce an empty line
# tmp_file <- tempfile()

# write_lines(file = tmp_file, x = tsv_file.data)
write_lines(file = tsv_file, x = '\n########## Dyplotype analysis results ##########\n', append = TRUE)
cat('\n##################################################')
cat('\n########## Dyplotype analysis results ############\n')
cat('##################################################\n')
for(gene in unique(genes %>% pull())){
  alleles <- tsv_file.data %>% filter(Gene == gene) %>% select(Allele)
  af <- tsv_file.data %>% filter(Gene == gene) %>% select(AF)
  
  cat(paste0("\n#Checking diplotype effect for gene: ",gene,"\n"))
  if(!nrow(alleles) > 1 && af == 0.5){
    
    cat("#Not enough alleles for dyplotype effect.\n")
    
  }else{
    #Case for a gene we have two alleles
    
    ######################## read DDBB ########################
    
    DiplotypeDB <- read.xlsx2(paste0(database_dir,'/Diplotype-Phenotype/',gene,"_Diplotype_Phenotype_Table.xlsx"), sheetIndex = 1, header = TRUE)
    # clean whitespaces
    DiplotypeDB <- DiplotypeDB %>% mutate(across(ends_with("Diplotype"), ~str_replace_all(., pattern = ' ', replacement = '')))
    # Separate Allele pairs
    # DiplotypeDB <- DiplotypeDB %>% separate(CYP2C9.Diplotype, sep = '/', into = c('CYP2C9.Diplotype.Allele1', 'CYP2C9.Diplotype.Allele2'))
    
    ######################## query DB ########################
    
    if(nrow(alleles) == 1 && af == 1){
      #Same allele twice
      query_allele <- paste(alleles[[1]][1], alleles[[1]][1], sep = '/')
    }else{
      #Two different alleles
      #TODO: Despite there being two alleles, one of the variants might have AF = 1
      query_allele <- paste(alleles[[1]][1], alleles[[1]][2], sep = '/')
    }
    
      DiplotypeDB.query <- DiplotypeDB %>% filter(across(ends_with("Diplotype"), ~. == query_allele))
      
      if(nrow(DiplotypeDB.query) != 0){
        
        cat(paste0('#Allele combination found for alleles: ',query_allele,'\n\n'))
        for(i in c(1:4)){
          cat(paste0(colnames(DiplotypeDB.query)[i],'\t'))
        }
        cat("\n")
        for(i in c(1:4)){
          cat(paste0(DiplotypeDB.query[[i]],'\t'))
        }
        write_tsv(file = tsv_file, x = DiplotypeDB.query, append = TRUE, col_names = TRUE)
      }else{
        
        cat(paste0('#No Diplotype effect for allele combination: ',query_allele,'\n\n'))
        
      }
  }
  cat("\n#--------------\n")
}
# file.copy(from = tmp_file,to = tsv_file,overwrite = TRUE)
cat(paste0('#Saved in: ',tsv_file,'\n'))
