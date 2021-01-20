#!/usr/bin/env nextflow

// Remove single_mode Variant_calling

def helpMessage() {
    log.info"""
================================================================
V A R I A N T  C A L L E R  - I R Y C I S    v 1.2
================================================================

    Usage:
    The typical command for running the pipeline is as follows:
    ./nextflow_jointVariantCalling.nf [OPTIONS]


    Options:
      --genome <GRCh37 | GRCh38 | [FILE]>   Reference genome to undergo the maping. Options: GRCh37, GRCh38, [/path/to/reference.fasta] (default: GRCh37)
      --region_intervals [BED FILE]         Specific genomic region in bed format (without chr) to constrict mapping and variant calling. Necessary for Whole Exome Sequencing and Panels. (default: NO_FILE)
      --dbSNP [FILE]                        Automatically provided when selecting GRCh37 or GRCh38, if using a custom reference and not provided base recalibration will not be performed.

      --vc <freebayes | gatk | varscan>     Variant caller to use for variant calling. Options: gatk, freebayes, varscan (default: freebayes)
      --GVCFmode <true | false>             This flag indicates all samples provided come from the same experiment and thus should be called together. As of Oct-2020 only compatible with freebayes. (default: true)
      --common_id [STRING]                  Id by which to identify all samples as coming from the same experiment. Assumed to be leading the file name. (default: first two characters of file name are used as experiment identifier)

      --min_alt_fraction                    Freebayes specific option, minimumn threshold at which allele frequency is considered real. (default: 0.2)

      --remove_duplicates <true | false>    Remove marked as duplicated reads. Options: true, false (default: false)
      --indir [DIR]                         The input directory, all .bam files in this directory will be processed. (default: "$baseDir/$params.outdir/alignment")
      --outdir [DIR]                        The output directory where the results will be saved (default: "$baseDir/$params.outdir")
      

    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}


log.info """\

================================================================
V A R I A N T  C A L L E R  - I R Y C I S    v 1.2
================================================================
genome               : $params.genome
region               : $params.region_intervals
dbSNP                : $params.dbSNP

variant_caller       : $params.vc
remove_duplicates    : $params.remove_duplicates

min_alt_fraction     : $params.min_alt_fraction

read_directory       : ./$params.indir
vcf_directory        : ./$params.outdir/alignment
results              : ./$params.outdir/raw_variant_calling_files
===============================================================
"""

if(params.vc == 'gatk' && params.GVCFmode != 'false'){
  println("gatk right now is incompatible with GVCFmode")
  exit(0)
}else if(params.vc == 'varscan' && params.GVCFmode != 'false'){
  println("varscan right now is incompatible with GVCFmode")
  exit(0)
}

if(params.genome != "GRCh37" && params.genome != "GRCh38"){
  custom_reference = file("${params.genome}") 
  prefixRef = custom_reference.name.take(custom_reference.name.lastIndexOf('.'))
}


def ploidy = params.ploidy != 'no' || params.ploidy == 'yes' && params.ploidy.getClass() == java.lang.Integer ? "--ploidy ${params.ploidy} ":''
def reference = params.genome != "GRCh37" && params.genome != "GRCh38" ? "${params.working_dir}/data/custom_reference/${prefixRef}.fasta": params.indexRef
def seqRef = file(params.seqRef)

if(!params.GVCFmode){

    def samplePrefix = params.common_id != '*' ? "${params.common_id}":""

      Channel
        .fromFilePairs("${params.outdir}/alignment/${samplePrefix}*{${bam_name},${bam_name}.bai}", size : 2, type: "any")
        .set{ch_variant_calling}
  
    process Variant_Calling_single {
      tag "Variant calling using selected Variant Caller (GATK, freebayes, varscan)"
      label 'med_mem'
      //publishDir "$params.outdir/raw_variant_calling_files", mode: 'copy'

      //def reference = params.genome != "GRCh37" && params.genome != "GRCh38" ? "${params.working_dir}/${params.indir}/reference/${prefixRef}.fasta": params.indexRef

      input:
      // Imma need to generate a tsv with all bam files paths listed so I can feed all of them to freebayes and other variant callers together
        set sampleId, file(bam_file),file(bai_file) from ch_variant_calling //.combine(ch_variant_calling2)

      output:

        set sampleId, file("*.vcf") into ch_vcf

      script:

        if(params.vc == 'gatk'){

        """
        gatk HaplotypeCaller --native-pair-hmm-threads ${params.threads} ${params.rmDups_GATK} ${region_interval} -I ${bam_file[0]} -O ${sampleId}.${params.vc}.vcf -R ${seqRef} ${params.vcOpts}
        """
        }else if(params.vc == 'freebayes'){

        """
        freebayes ${ploidy} -f ${seqRef} ${bam_file[0]} > ${sampleId}.${params.vc}.vcf
        """
        }else if(params.vc == 'varscan'){
        """
        samtools mpileup -B -f ${seqRef} ${bam_file[0]} | varscan mpileup2cns --variants --output-vcf 1 > ${sampleId}.${params.vc}.vcf
        
        """
        }
      }

    }else{

      // Gather all alignment files from the same experiment to process together

      def samplePrefix = params.common_id != '*' ? "${params.common_id}":""
      def bam_name = params.dbSNP == 'NO_FILE' ? "sort.bam":"bqsr.bam"

      if(params.common_id == '*'){
          Channel
            .fromFilePairs("${params.outdir}/alignment/${samplePrefix}*{${bam_name},${bam_name}.bai}", size : 2, type: "any", checkIfExists: true)
            .map{it -> [it[0][0,1], it[1]]}
            .groupTuple()
            .set{ch_variant_calling}

      }else{
          Channel
            .fromFilePairs("${params.outdir}/alignment/${samplePrefix}*{${bam_name},${bam_name}.bai}", size : 2, type: "any")
            .map{it -> [params.common_id, it[1]]}
            .groupTuple()
            .set{ch_variant_calling}
      }

      process Variant_Calling_batch {
        tag "Variant calling using selected Variant Caller (GATK, freebayes, varscan)"
        label 'med_mem'
        //publishDir "$params.outdir/raw_variant_calling_files", mode: 'copy'

        //def reference = params.genome != "GRCh37" && params.genome != "GRCh38" ? "${params.working_dir}/${params.indir}/reference/${prefixRef}.fasta": params.indexRef

        input:
          set expId, val(data) from ch_variant_calling //.combine(ch_variant_calling2)

        output:
          set expId, file('*vcf') into ch_vcf_sample

        script:

          String bams = data.flatten().collate(1,2).flatten().join(" ")

          if(params.vc == 'gatk'){

            println("gatk right now is incompatible with GVCFmode")

          }else if(params.vc == 'freebayes'){

            def GVCF = "--gvcf"
            def min_alt_fraction_var = params.min_alt_fraction == '' ? 0.2:"${params.min_alt_fraction}"
          """
          freebayes ${ploidy} --min-alternate-fraction ${min_alt_fraction_var} ${GVCF} -f ${seqRef} ${bams} > ${expId}.${params.vc}.vcf
          cat ${expId}.${params.vc}.vcf | grep -v '<\\*>' > ${expId}.${params.vc}.clean.vcf
          """
          }else if(params.vc == 'varscan'){

          println("varscan right now is incompatible with GVCFmode")

          }
        }


      process Sample_isolation {
        tag "From a gVCF isolate different samples with their particular mutations"
        publishDir "$params.outdir/raw_variant_calling_files/samples_isolates", mode: 'copy'
      
        input:
          set sampleId, file(vcf_file), file(idx_file) from ch_vcf_sample
        output:
          file("*.vcf") into ch_vcf
      
        script:
        """
        $baseDir/scripts/extract_all_samples ${vcf_file[0]}
        """
      }
    }
    
    process Sample_indexing {
    tag "index vcf files of isolated samples"
    publishDir "$params.outdir/raw_variant_calling_files", mode: 'copy'
  
    input:
      file(vcf_file) from ch_vcf.flatten()
    output:
      file("${vcf_file}")
      file('*.vcf.idx')
  
    script:
    """
    gatk IndexFeatureFile --input ${vcf_file}
    """
  }
