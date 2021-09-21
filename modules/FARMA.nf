#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*

process liftover {
    publishDir = "$params.outdir/annotation"
    label 'liftover'

    input:
        path vcf
    output:
        file('*.vcf')

    when:
        $params.genome == "GRCh37"
    script:
        """
        java -Xmx8g -jar /usr/picard/picard.jar LiftoverVcf INPUT=$vcf OUTPUT=${vcf.baseName}.liftover.vcf CHAIN=$params.chain_path R=$params.seqRef --WRITE_ORIGINAL_POSITION
        """

}


*/

// ################## FILTERS ########################

process snpSift_filter_rsID {
    publishDir = "$params.outdir/FARMA"
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
    publishDir = "$params.outdir/FARMA"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        cat $vcf | snpSift filter --set $params.farmaDB/Allele_definition/available_genes.list "ANN[*].GENE in SET[0]" > ${vcf.baseName}.genes.vcf
        """
}
/*
process snpSift_filter_func_genes {
    publishDir = "$params.outdir/FARMA"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        cat $vcf | snpSift filter --set $params.farmaDB/Allele_functionality/available_genes.list "ANN[*].GENE in SET[0]" > ${vcf.baseName}.genes.vcf
        """
}
process snpSift_filter_dip_genes {
    publishDir = "$params.outdir/FARMA"
    label 'snpEffect'

    input:
        path vcf
    output:
        file('*.vcf')
    script:
        """
        cat $vcf | snpSift filter --set $params.farmaDB/Diplotype-Phenotype/available_genes.list "ANN[*].GENE in SET[0]" > ${vcf.baseName}.genes.vcf
        """
}
*/
// ################## Data Manipulation ########################

process extract_info {
    publishDir = "$params.outdir/FARMA/"
    label 'snpEffect'

    input:
        path vcf
    output:
        tuple file('*tsv'), file(vcf)
    script:
        """
        snpSift extractFields $vcf CHROM POS REF ALT ID ANN[0].GENE AF > ${vcf.baseName}.ID.tsv
        """
}

// ################## R scripts ########################

process allele_def {
    publishDir = "$params.outdir/FARMA"
    label 'R'

    input:
        tuple file(tsv), file(vcf)
    output:
        file('*.log')
    script:
    """
    ~/Variant_Calling/scripts/allele_def_id.R $vcf $tsv > ${vcf.baseName}.alleleDef.log
    """

}

process match_alleles {
    publishDir = "$params.outdir/FARMA"
    label 'R'

    input:
        file(log_file)
    output:
        tuple file('*.tsv'), file(log_file)
    script:
    """
    ~/Variant_Calling/scripts/match_allele.R $log_file
    """

}

process diplotype_analysis {
    publishDir = "$params.outdir/FARMA"
    label 'R'

    input:
        tuple file(tsv), file(log_file)
    output:
        tuple file(tsv), file(log_file)
    script:
    """
    ~/Variant_Calling/scripts/diplotype_analysis.R $tsv >> $log_file
    """

}