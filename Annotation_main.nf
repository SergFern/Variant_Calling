#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\

================================================================
A N N O T A T O R  - snpEff    v 0.1
================================================================
VCF_files               : $params.VCF_files
genome reference        : $params.seqRef
outdir                  : $params.outdir
================================================================
"""

include { VCF_Normalization; VCF_Decomposition; Annotation } from './modules/Annotation.nf'


workflow {

    data = channel.fromPath(params.VCF_files, checkIfExists: true)

    VCF_Normalization(data)
    VCF_Decomposition(VCF_Normalization.out)
    //Annotation(VCF_Decomposition.out)

}

/* COMMENT SECTION
TODO: use container bcftools
TODO: integrated workflow
TODO: define snpEff process as module
TODO: investigate - specific config file
TODO: separate profiles configuration from general .config file
*/