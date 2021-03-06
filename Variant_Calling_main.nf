#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\

================================================================
A N N O T A T O R  - snpEff    v 0.1
================================================================
BAM_files_location      : $params.BAM_files_location
genome reference        : $params.seqRef
outdir                  : $params.outdir
================================================================
"""

include { Variant_Calling_single } from './modules/Variant_Calling.nf'

workflow {

    data = channel.fromFilePairs("$params.BAM_files_location/*.{bam,bai}", checkIfExists: true).take( params.dev ? params.number_of_inputs : -1 )
    data.view()
    Variant_Calling_single(data)

}