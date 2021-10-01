#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process liftover {
    publishDir = "$params.outdir/annotation"
    label 'picard'

    input:
        path vcf
    output:
        tuple file('*liftover.vcf'), file('.*processes.log')

    when:
        params.genome == "GRCh37"
    script:
        """
        java -Xmx8g -jar /usr/picard/picard.jar LiftoverVcf INPUT=$vcf OUTPUT=${vcf.baseName}.liftover.vcf CHAIN=$params.chain_path R=$params.seqRef_GRCh38 WRITE_ORIGINAL_POSITION=true REJECT=${vcf.baseName}.liftover_rejected.vcf &> .${vcf.baseName}.FARMA_processes.log
        """

}

// ################## FILTERS ########################
/*
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
}*/

process snpSift_filter_def_genes {
    publishDir = "$params.outdir/FARMA"
    label 'snpEffect'

    input:
        tuple file(vcf), file(processes)
    output:
        tuple file('*.vcf'), file(processes)
    script:
        """
        cat $vcf | snpSift filter --set $params.farmaDB/Allele_definition/available_genes.list "ANN[*].GENE in SET[0]" > ${vcf.baseName}.genes.vcf 2>> $processes
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
        tuple file(vcf), file(processes)
    output:
        tuple file('*tsv'), file(vcf), file(processes)

    script:
        """
        snpSift extractFields $vcf -s ',' CHROM POS REF ALT ID ANN[*].GENE AF > ${vcf.baseName}.ID.tsv 2>> $processes
        """
}

// ################## R scripts ########################

process allele_def {
    publishDir = "$params.outdir/FARMA"
    label 'R'

    input:
        tuple file(tsv), file(vcf), file(processes)
    output:
        tuple file('*.log'), file(processes)
    script:
    """
    ~/Variant_Calling/scripts/allele_def_id.R $vcf $tsv > ${vcf.baseName}.alleleDef.log 2>> $processes
    """

}

process match_alleles {
    publishDir = "$params.outdir/FARMA"
    label 'R'

    input:
        tuple file(log_file), file(processes)
    output:
        tuple file('*.tsv'), file(log_file), file(processes)
    script:
    """
    ~/Variant_Calling/scripts/match_allele.R $log_file 2>> $processes
    """

}

process diplotype_analysis {
    publishDir = "$params.outdir/FARMA"
    label 'R'

    input:
        tuple file(tsv), file(log_file), file(processes)
    output:
        tuple file(tsv), file(log_file), file(processes)
    script:
    """
    ~/Variant_Calling/scripts/diplotype_analysis.R $tsv >> $log_file 2>> $processes
    """

}