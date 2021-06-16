#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process VCF_Normalization {
    publishDir = "results"
    label 'bcftools'

    input:
        path vcf
    output:
        file('*')
    script:
        """
        bcftools norm -m-both -o ${vcf.baseName}.norm.vcf $vcf
        """
}

process VCF_Decomposition {
    publishDir = "results"
    label 'bcftools'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        bcftools norm -f $params.seqRef -o ${vcf.baseName}.decomp.vcf $vcf
        """
}

process Annotation {
    publishDir = "results"

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        snpEff -v -dataDir $params.snpEffdb GRCh37.75 $vcf > ${vcf.baseName}.annot.vcf
        """

}

/* 

COMMENT SECTION

*/