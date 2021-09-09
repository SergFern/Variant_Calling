#!/usr/bin/env nextflow

//TODO: --help message Annotation main

nextflow.enable.dsl=2

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

Annotating db           : $params.dbSNP_annot
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