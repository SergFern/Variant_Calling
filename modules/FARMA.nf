#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ################## FILTERS ########################

process snpSift_filter_rsID {
    publishDir = "$params.outdir/results/annotation"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        cat $vcf | snpSift filter --set $params.rsID_list_FARMA "ID in SET[0]" > ${vcf.baseName}.set.vcf
        """

}

process snpSift_filter_def_genes {
    publishDir = "$params.outdir/results/annotation"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        cat $vcf | snpSift filter --set $params.farmaDB/Allele_definition/available_genes.list "ANN[*].GENE in SET[0]" > ${vcf.baseName}.set.genes.vcf
        """
}
process snpSift_filter_func_genes {
    publishDir = "$params.outdir/results/annotation"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        cat $vcf | snpSift filter --set $params.farmaDB/Allele_functionality/available_genes.list "ANN[*].GENE in SET[0]" > ${vcf.baseName}.set.genes.vcf
        """
}
process snpSift_filter_dip_genes {
    publishDir = "$params.outdir/results/annotation"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        cat $vcf | snpSift filter --set $params.farmaDB/Diplotype-Phenotype/available_genes.list "ANN[*].GENE in SET[0]" > ${vcf.baseName}.set.genes.vcf
        """
}

// ################## Data Manipulation ########################

process extract_info {
    publishDir = "$params.outdir/results/FARMA/"
    label 'snpEffect'

    input:
        path vcf
    output:
        tuple file('*tsv'), file(vcf)
    script:
        """
        snpSift extractFields $vcf CHROM POS ID ANN[0].GENE > ${vcf.baseName}.ID.genes.tsv
        """
}

// ################## R scripts ########################

process match_alleles {
    publishDir = "$params.outdir/results/FARMA"
    label 'R'

    input:
        tuple file(tsv), file(vcf)
    output:
        file('*.report')
    script:
    """
    ~/Variant_Calling/scripts/allele_def_id.R $vcf $tsv > ${vcf.baseName}.alleleDef.report
    """

}