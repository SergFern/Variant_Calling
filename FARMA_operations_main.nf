#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\

================================================================
A N N O T A T O R  - snpEff    v 0.1
================================================================
VCF_files               : $params.VCF_files
genome reference        : $params.seqRef
outdir                  : $params.outdir

rsID set                : $params.rsID_list_FARMA
farmaDB                 : $params.farmaDB
================================================================
"""

include { snpSift_filter_rsID, snpSift_filter_def_genes} from './modules/FARMA.nf'


workflow {

    data = channel.fromPath(params.VCF_files, checkIfExists: true)
    snpSift_filter_rsID(data)
    snpSift_filter_def_genes(snpSift_filter_rsID.out)

}