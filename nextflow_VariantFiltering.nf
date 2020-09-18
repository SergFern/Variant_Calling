#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ==================================================================
    ${workflow.manifest.name}  ~  version ${workflow.manifest.version}
    ==================================================================

    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run [OPTIONS]
    Options:
      
      --indir                          The input directory, all vcf files in this directory will be processed. (default: $params.outdir/raw_variant_calling)
      --outdir                         The output directory where the results will be saved (default: $params.outdir)
      

    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

log.info """\

================================================================
V A R I A N T  C A L L E R  - I R Y C I S    v 0.9-Custom Ref
================================================================
read_directory       : ./$params.VCF_files
results              : ./$params.outdir/curated_variant_calling_files
================================================================
"""

Channel
  .fromFilePairs(params.VCF_files, checkifExists : true, size : 2)
  .take( params.dev ? params.number_of_inputs : -1 ) //TESTING: should only run partial data
  .set{ch_vcf}

//------------------
if(params.genome != "GRCh37" && params.genome != "GRCh38"){

custom_reference = file("${params.genome}") 
prefixRef = custom_reference.name.take(custom_reference.name.lastIndexOf('.'))

}
// Defines regions for the variant calling to focus on.--------
def region_interval = params.region_intervals != 'NO_FILE' ? "-L ${params.region_intervals} -ip 100 ":''

//---------------------------------------------Variant Calling----------------------------------------

process SelectVariants {
  tag "Selects SNV or indels"
  publishDir "$params.outdir/raw_variant_calling_files/intermediate_vcfs"

  input:
    set sampleId, file(vcf_file) from ch_vcf
  output:
    set sampleId, file("*.snp.vcf"), file("*.snp.vcf.idx") into ch_SNV
    set sampleId, file("*.indel.vcf"), file("*.indel.vcf.idx") into ch_indels
  script:
  """
  gatk SelectVariants -V ${vcf_file[0]} -select-type SNP -O ${sampleId}.${params.vc}.snp.vcf
  gatk SelectVariants -V ${vcf_file[0]} -select-type INDEL -O ${sampleId}.${params.vc}.indel.vcf
  """
}

process VariantFiltration {
  tag "Filter both SNP and indel variants"
  publishDir "$params.outdir/raw_variant_calling_files/intermediate_vcfs"

  input:
    set sampleId, file(snp_file) from ch_SNV
    set sampleId, file(indel_file) from ch_indels
  output:
    set sampleId, file("*.snp.filtered.vcf") into ch_SNV_filtered
    set sampleId, file("*.indel.filtered.vcf") into ch_indel_filtered
  script:
  """
  gatk VariantFiltration -V ${snp_file[0]} -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR > 3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' -O ${sampleId}.${params.vc}.snp.filtered.vcf
  gatk VariantFiltration -V ${indel_file[0]} -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' -O ${sampleId}.${params.vc}.indel.filtered.vcf

  """
}

process Sorting_Variants {
  tag "Sort Variants in preparation for the MergeVcfs process"
  publishDir "$params.outdir/raw_variant_calling_files/intermediate_vcfs"

  input:
    set sampleId, file(snp_file) from ch_SNV_filtered
    set sampleId, file(indel_file) from ch_indel_filtered
  output:
    set sampleId, file("*.snp.filtered.sorted.vcf") into ch_SNV_sorted
    set sampleId, file("*.indel.filtered.sorted.vcf") into ch_indel_sorted

  script:
  """
  gatk SortVcf -I ${snp_file[0]} -R ${params.seqRef} -O ${sampleId}.${params.vc}.snp.filtered.sorted.vcf
  gatk SortVcf -I ${indel_file[0]} -R ${params.seqRef} -O ${sampleId}.${params.vc}.indel.filtered.sorted.vcf
  """
}

process MergeVcfs {
  tag "Merges both snp and indel vcfs, filtered and sorted"
  publishDir "$params.outdir/raw_variant_calling_files", mode: 'copy'

  input:
    set sampleId, file(snp_file) from ch_SNV_sorted
    set sampleId, file(indel_file) from ch_indel_sorted
  output:
    set sampleId, file("*.curated.vcf") into ch_curated

  script:
  """
  gatk MergeVcfs -I ${snp_file[0]} -I ${indel_file[0]} -R ${params.seqRef} -O ${sampleId}.${params.vc}.curated.vcf

  """
}

process IndexFeatureFile {
  tag "Indexes final curated vcf"
  publishDir "$params.outdir/curated_variant_calling_files", mode: 'copy'

  input:
    set sampleId, file(vcf_file) from ch_curated
  output:
    set sampleId, file("*.vcf.idx")

  script:
  """
  gatk IndexFeatureFile --input ${vcf_file[0]}
  """
}