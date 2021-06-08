#!/usr/bin/env nextflow

process Indexing_custom_genome {
    tag "Indexes supplied reference FASTA file (Samtools)"
    label 'big_mem'

    publishDir "data/custom_reference", mode: 'copy'


    input:
    file reference_file from ch_reference

    output:
    file("${reference_file[0]}")
    file '*.fai'
    
    """
    samtools faidx ${reference_file[0]}
    """
    
}

process Building_genome_dictionary {
    tag "Creating Dictionary file (GATK)"
    label 'big_mem'

    publishDir "data/custom_reference", mode: 'copy'

    input:
    file reference_file from ch_reference

    output:
    file '*.dict'

    script:
    
    """
    gatk CreateSequenceDictionary -R ${reference_file[0]}
    """
    
  }

  // Extract file name:
  custom_reference = file("$params.genome")
  prefixRef = custom_reference.name.take(custom_reference.name.lastIndexOf('.'))

  process Custom_genome_indexing {
    tag "Indexes reference file using the specified aligner"
    label 'big_mem'

    publishDir "data/custom_reference", mode: 'copy'

    input:
    file reference_file from ch_reference

    output:
    file "*"

    script:

    if(params.aln == 'bwa'){

    """
      bwa index ${reference_file[0]}
    """

    }else if(params.aln == 'bowtie2'){

    """
      bowtie2-build ${reference_file[0]} ${prefixRef}.bowtie2
    """
    }
  }