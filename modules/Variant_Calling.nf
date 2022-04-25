#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def min_alt_fraction_var = params.min_alt_fraction == '' ? 0.2:"${params.min_alt_fraction}"
def seqRef = file(params.seqRef)
def region_interval = params.region_intervals != 'NO_FILE' ? "-L ${params.region_intervals} -ip 100 ":''

process Variant_Calling_single {
    tag "Variant calling using selected Variant Caller (GATK, freebayes, varscan)"
    label 'big_mem'
    label 'generic'
    publishDir "$params.outdir/vcf_raw", mode: 'copy'

    input:
        tuple val(SampleId), file(bam_file)
    output:
        file('*vcf')

    script:

    if(params.vc == 'gatk'){
        """
        gatk HaplotypeCaller --native-pair-hmm-threads ${params.threads} ${params.rmDups_GATK} ${region_interval} -I ${bam_file[1]} -O ${SampleId}.${params.vc}.vcf -R ${seqRef} ${params.vcOpts}
        """
    }else if(params.vc == 'freebayes'){
        """
        freebayes --min-alternate-fraction ${min_alt_fraction_var} -f ${seqRef} ${bam_file[1]} > ${SampleId}.${params.vc}.vcf
        """
    }else if(params.vc == 'varscan'){
        """
        samtools mpileup -B -f ${seqRef} ${bam_file[1]} | varscan mpileup2cns --variants --output-vcf 1 > ${SampleId}.${params.vc}.vcf
        """
    }
}