#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process VCF_Normalization {
    publishDir = "results"

    input:
        path data
    output:
        file('*.vcf')
    script:
        """
        bcftools norm -m-both -o ${data.baseName}.norm.vcf $data
        """
}

process VCF_Decomposition {
    publishDir = "results"

    input:
        path data
    output:
        file('*.vcf')
    script:
        """
        bcftools norm -f $params.seqRef -o ${data.baseName}.decomp.vcf $data
        """
}

/* COMMENT SECTION

*/