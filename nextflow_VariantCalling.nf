#!/usr/bin/env nextflow

log.info """\

================================================================
V A R I A N T  C A L L E R  - I R Y C I S    v 0.6
================================================================
genome               : $params.genome
region               : $params.region_intervals
read_directory       : ./$params.indir
paired               : $params.paired
reads                : $params.reads
adapters             : $params.adapter_file
variant_caller       : $params.vc
remove_duplicates    : $params.remove_duplicates
results              : ./$params.outdir
================================================================
"""

Channel
  .fromFilePairs(params.maped_reads, checkifExists : true)
  .set{ch_recalculated_bam}

def region_interval = params.region_intervals != 'NO_FILE' ? "-L ${params.region_intervals} -ip 100 ":''

//---------------------------------------------Variant Calling----------------------------------------

process VariantCalling {
  tag "Variant calling using selected Variant Caller (GATK, freebayes, varscan)"
  publishDir "$params.outdir/variant_calling_files", mode: 'copy'

  input:
    set sampleId, file(bam_file) from ch_recalculated_bam
  output:
    set sampleId, file("*vcf"), file("*.vcf.idx") optional true into ch_vcf

  def GVCF = params.GVCFmode == 'false' ? "":"-ERC GVCF"

  print


  script:

    if(params.vc == 'gatk'){
    """
    gatk HaplotypeCaller --native-pair-hmm-threads ${params.threads} ${params.rmDups_GATK} ${region_interval} ${params.vcOpts} -I ${bam_file[1]} -O ${sampleId}.${params.vc}.vcf -R ${params.seqRef}
    """
    }else if(params.vc == 'freebayes'){
    """
    freebayes -f ${params.seqRef} ${params.vcOpts} ${bam_file[0]} > ${sampleId}.${params.vc}.vcf
    """
    }else if(params.vc == 'varscan'){
    """
    samtools mpileup -B -f ${params.seqRef} ${bam_file[0]} | varscan mpileup2cns --variants --output-vcf 1 > ${sampleId}.${params.vc}.vcf ${GVCF}
    """
    }
  }

process SelectVariants {
  tag "Selects SNV or indels"
  publishDir "$params.outdir/variant_calling_files"

  input:
    set sampleId, file(vcf_file), file(vcf_idx_file) from ch_vcf
  output:
    set sampleId, file("*.snp.vcf"), file("*.snp.vcf.idx") into ch_SNV
    set sampleId, file("*.indel.vcf"), file("*.indel.vcf.idx") into ch_indels

  script:
  """
  gatk SelectVariants -V ${vcf_file[0]} -select-type SNP -O ${sampleId}.${params.vc}.snp.vcf
  gatk SelectVariants -V ${vcf_file[0]} -select-type SNP -O ${sampleId}.${params.vc}.indel.vcf
  """
}

process VariantFiltration {
  tag "Filter both SNP and indel variants"
  publishDir "$params.outdir/variant_calling_files"

  input:
    set sampleId, file(snp_file), file(snp_idx_file) from ch_SNV
    set sampleId, file(indel_file), file(indel_idx_file) from ch_indels
  output:
    set sampleId, file("*.snp.filtered.vcf") into ch_SNV_filtered
    set sampleId, file("*.indel.filtered.vcf") into ch_indel_filtered
  script:
  """
  gatk VariantFiltration -V ${snp_file[0]} -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'SOR > 3.0' --filter-name 'SOR3' -filter 'FS > 60.0' --filter-name 'FS60' -filter 'MQ < 40.0' --filter-name 'MQ40' -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' -O ${sampleId}.${params.vc}.snp.filtered.vcf
  gatk VariantFiltration -V ${indel_file[0]} -filter 'QD < 2.0' --filter-name 'QD2' -filter 'QUAL < 30.0' --filter-name 'QUAL30' -filter 'FS > 200.0' --filter-name 'FS200' -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' -O ${sampleId}.${params.vc}.indel.filtered.vcf

  """
}

process MergeVcfs {
  tag "Merges both snp and indel vcfs"
  publishDir "$params.outdir/variant_calling_files", mode: 'copy'

  input:
    set sampleId, file(snp_file), file(snp_idx_file) from ch_SNV_filtered
    set sampleId, file(indel_file), file(indel_idx_file) from ch_indel_filtered
  output:
    set sampleId, file("*.curated.vcf") into ch_curated

  script:
  """
  gatk MergeVcfs -I ${snp_file[0]} -I ${indel_file[0]} -O ${sampleId}.${params.vc}.curated.vcf

  """
}

process IndexFeatureFile {
  tag "Indexes final curated vcf"
  publishDir "$params.outdir/variant_calling_files", mode: 'copy'

  input:
    set sampleId, file(vcf_file) from ch_curated
  output:
    set sampleId, file("*.vcf.idx")

  script:
  """
  gatk IndexFeatureFile --input ${vcf_file[0]}
  """
}