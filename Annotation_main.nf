#!/usr/bin/env nextflow

//TODO: --help message Annotation main

nextflow.enable.dsl=2

//TODO: send process Normalization and Decomposition to main Variant Calling pipelines.

def helpMessage() {
    log.info"""
================================================================
V A R I A N T  A N N O T A T O R  - I R Y C I S    v 1
================================================================

    Usage:
    The typical command for running the pipeline is as follows:
    ./Annotation_main.nf [OPTIONS]


    Options:

    --VCF_files [FILE]                    Input vcf files. Required.
    --genome <GRCh37 | GRCh38 | [FILE]>   Reference genome to undergo the maping. Options: GRCh37, GRCh38, [/path/to/reference.fasta] (default: $params.genome)
    --genome_annot <GRCh37.87>            Annotation reference for snpEff. default: GRCh37.87
    --snpEffdb [PATH]                     Path to snpEff db. default: "$params.BBDD/snpEffdb"
    --dbSNP_annot [FILE]                  File used to identify dbSNP IDs. default: "$params.dbSNP_dir/$params.genome/All_20180423.vcf.gz"

    --outdir [PATH]                       The output directory where the results will be saved. default: "$baseDir/$params.outdir" 

    """.stripIndent()
}


// Show help message if --help specified
if (params.help){
    helpMessage()
    exit 0
}

log.info """\

================================================================
A N N O T A T O R  - snpEff    v 0.1
================================================================
VCF_files               : $params.VCF_files
genome reference        : $params.seqRef
outdir                  : $params.outdir

#### snpEff ####

genome_annot            : $params.genome_annot
snpEffdb                : $params.snpEffdb

#### snpSift ####

dbSNP annotating db     : $params.dbSNP_annot
rsID set                : $params.rsID_list_FARMA
================================================================
"""

include { VCF_Normalization; VCF_Decomposition; snpEff; snpSift_annotate } from './modules/Annotation.nf'


workflow {

    data = channel.fromPath(params.VCF_files, checkIfExists: true)

    VCF_Normalization(data)
    VCF_Decomposition(VCF_Normalization.out)
    snpEff(VCF_Decomposition.out)
    snpSift_annotate(snpEff.out)

}