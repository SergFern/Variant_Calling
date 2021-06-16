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

include { VCF_Normalization; VCF_Decomposition } from './modules/Annotation.nf'


workflow {

    data = channel.fromPath(params.VCF_files, checkIfExists: true)

    VCF_Normalization(data)
    VCF_Decomposition(VCF_Normalization.out)

}

/* COMMENT SECTION
TODO: bcftools as container
TODO: integrated workflow
TODO: SampleId to identify different files - check GenomeMapper processes
TODO: define snpEff process as module
TODO: decide - different modules for bcftools and snpeff or integration?
TODO: investigate - specific config file
TODO: separate profiles configuration from general .config file
*/