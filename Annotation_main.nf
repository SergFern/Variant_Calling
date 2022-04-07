#!/usr/bin/env nextflow

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

include { VCF_Normalization; VCF_Decomposition; snpEff; snpSift_annotate; snpSift_filter } from './modules/Annotation.nf'


workflow {

    data = channel.fromPath(params.VCF_files, checkIfExists: true)

    VCF_Normalization(data)
    VCF_Decomposition(VCF_Normalization.out)
    snpEff(VCF_Decomposition.out)
    snpSift_annotate(snpEff.out)
    //snpSift_filter(snpSift_annotate.out)


}

/* COMMENT SECTION
TODO: modular workflows
TODO: add generic label to VariantFIltering and jointVariantCalling script
TODO: investigate - specific config file
TODO: separate profiles configuration from general .config file
TODO: Log files with process versions, container status, etc
*/