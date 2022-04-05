#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process VCF_Normalization {
    publishDir = "$params.outdir/vcf_cleanup"
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
    publishDir = "$params.outdir/vcf_cleanup"
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

process snpEff {
    //publishDir = "$params.outdir/annotation"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        snpEff -v -dataDir $params.snpEffdb $params.genome_annot $vcf > ${vcf.baseName}.snpeff.vcf
        """

}

process snpSift_annotate {
    publishDir = "$params.outdir/annotation"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        snpSift annotate -id $params.dbSNP_annot $vcf > ${vcf.baseName}.snpsift.vcf
        """

}