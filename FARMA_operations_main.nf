#!/usr/bin/env nextflow

//TODO: make --chain_path optional

nextflow.enable.dsl=2

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

    --rsID set  [FILE]                    A list of rsIDs to filter entries with. All IDs present in file will be kept.
    --farmaDB [PATH]                      Directory where FARMA DB are stored, default: "${params.BBDD}/FARMA/cipic_info_genes"
    --chain_path [FILE]                   File for genome reference liftover if necesary: default: "$params.chain_path"

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

rsID set                : $params.rsID_list_FARMA
farmaDB                 : $params.farmaDB
chain_path              : $params.chain_path
================================================================
"""

include { snpSift_filter_rsID; snpSift_filter_def_genes; extract_info; allele_def; match_alleles; diplotype_analysis } from './modules/FARMA.nf'


workflow {

    data = channel.fromPath(params.VCF_files, checkIfExists: true)

    snpSift_filter_rsID(data)
    snpSift_filter_def_genes(snpSift_filter_rsID.out)
    extract_info(snpSift_filter_def_genes.out)
    allele_def(extract_info.out)
    match_alleles(allele_def.out)
    diplotype_analysis(match_alleles.out)

}