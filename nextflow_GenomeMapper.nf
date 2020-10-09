#!/usr/bin/env nextflow

//TODO:
//Provide an error if adapter filtering was usuccessful

def helpMessage() {
    log.info"""

	================================================================
	V A R I A N T  C A L L E R  - I R Y C I S    v 1.2
	================================================================


    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run [OPTIONS]
    Options:
      --genome <GRCh37|GRCh38|[FILE]>     Reference genome to undergo the maping. Options: GRCh37, GRCh38, [/path/to/reference.fasta] (default: GRCh37)
      --adapter_file [FILE]               Adapter file to trimm reads by. (Trimmomatic adapters provided in $baseDir/adapter)
      --region_intervals [BED FILE]       Specific genomic region in bed format (without chr) to constrict mapping and variant calling. Necessary for Whole Exome Sequencing and Panels. (default: NO_FILE)
      --dbSNP [FILE]                      Automatically provided when selecting GRCh37 or GRCh38, if a custom reference is used and a custom_dbSNP is not provided base recalibration will not be performed.

      --paired                            Execute pipleine in single-end or paired-end mode. If "--paired true" then all fastq files in $params.indir will be processed as samples from the same experiment.
                                          If "--paired false" a csv with the single-end files path and their IDs will be used to identify the fasq files. Options: true, false (default: true)

      --reads [GLOB]                      Path to paired-end reads or single-end csv file  (default: "$baseDir/data/*R{1,2}*.fastq.gz (if paired) "$baseDir/data/*.csv" (if single-end))
                                          CSV format: SampleID, [path/to/read].fastq
      --aln <bwa|bowtie2>                 Aligner chosen to map reads to the reference. Options: bwa, bowtie2 (default: bwa)
      --vc  <freebayes|gatk|varscan> Variant caller to use for variant calling. Overrideed by --skip_variant_calling Options: gatk, freebayes, varscan (default:gatk)
                                          By default the variant caller will execute in single sample mode. For joint variant calling use jointVariantCalling.nf pipeline.
      --common_id [STRING]                Id by which to identify all samples as coming from the same experiment. Assumed to be leading the file name. (default: first two characters of file name are used as experiment identifier)

      --skip_variant_calling              Skips variant calling process entirely, only perform the alignment (default:false)
      --remove_duplicates                 Remove marked as duplicated reads. Options: true, false (default: false)
      --min_alt_fraction                  Freebayes specific option, minimumn threshold at which allele frequency is considered real. (default: 0.2)
      --indir [DIR]                       The input directory, all fastq files or csv files in this directory will be processed. (default: "data")
      --outdir [DIR]                      The output directory where the results will be saved (default: "my-results")
      

    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

if(params.paired){

  //Reads must be read twice, both as a trupple and as an array
   Channel
      .fromFilePairs(params.reads, size: 2,  checkIfExists: true, flat:true)
      .take( params.dev ? params.number_of_inputs : -1 ) //TESTING: should only run partial data
      .into{ [ch_pre_samples, ch_FlowCell_lane] }


   //Make ch_pre_samples a tuple to process SampleId in processes more easily and to include it in ReadGroup (RG) info.
   ch_pre_samples.map{it -> new Tuple(it[0].split("_")[0],it[1,2])}.set{ ch_samples_with_id }

   //Identify fastq data for Read group definition "ID". RG:ID = {flowcell}.{lane}.{uniqueId} [must be unique despite documentation stating otherwise]
   ch_FlowCell_lane.splitFastq(record: true, pe: true, limit: 1).map{it -> new Tuple(it[0].split("_")[0], it[1].readHeader.split(":")[2,3,9].join("."))}.set{ch_RG_ID}

}else{

   Channel
      .fromPath(params.reads, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> [row.sampleId, [row.read]] }
      .take( params.dev ? params.number_of_inputs : -1 ) //TESTING: should only run partial data
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into{ [ch_samples, ch_sampleName, ch_FlowCell_lane] }

   ch_FlowCell_lane.splitFastq(record: true, pe: true, limit: 1).map{it -> new Tuple(it[0].split("_")[0], it[1].readHeader.split(":")[2,3,9].join("."))}.set{ch_RG_ID}

 }
//Fuse Ids and samples to manage SampleId as a tuple together with the samples they identify.
ch_RG_ID.concat(ch_samples_with_id).groupTuple().map{ it -> [[it[0],it[1][0]],it[1][1]] }.set{ ch_samples }


ch_dbSNP = file(params.dbSNP)
ch_adapter = file(params.adapter_trimm)

def region_interval = params.region_intervals != 'NO_FILE' ? "-L ${params.region_intervals} -ip 100 ":''
def ploidy = params.ploidy != 'no' || params.ploidy == 'yes' && params.ploidy.getClass() == java.lang.Integer ? "--ploidy ${params.ploidy} ":''

log.info """\

================================================================
V A R I A N T  C A L L E R  - I R Y C I S    v 1.2
================================================================
genome               : $params.genome
reads                : $params.reads
adapters             : $params.adapter_file

region_intervals     : $params.region_intervals
dbSNP                : ${params.dbSNP}

paired               : $params.paired
aligner              : $params.aln
variant_caller       : $params.vc
remove_duplicates    : $params.remove_duplicates
ploidy               : $params.ploidy

read_directory       : $params.indir
results              : $params.outdir
===============================================================
"""

//------------------------------------------------------------Trimming-------------------------------------------------

if(params.genome != "GRCh37" && params.genome != "GRCh38"){


  ch_reference = file(params.genome, checkIfExists: true)
  //(ch_reference_idx, ch_reference_dic) = Channel.from(ch_reference).into(2) apparently there is no need for this?

  process Indexing_custom_genome {
    tag "Indexes supplied reference FASTA file (Samtools)"
    label 'med_mem'

    publishDir "$params.indir/custom_reference", mode: 'copy'


    input:
    file reference_file from ch_reference

    output:
    file("${reference_file[0]}")
    file '*.fai'
    
    """
    samtools faidx ${reference_file[0]}
    """
    
  }

   process Building_genome_dictionary {
    tag "Creating Dictionary file (GATK)"
    label 'med_mem'

    publishDir "$params.indir/custom_reference", mode: 'copy'

    input:
    file reference_file from ch_reference

    output:
    file '*.dict'

    script:
    
    """
    gatk CreateSequenceDictionary -R ${reference_file[0]}
    """
    
  }

  custom_reference = file("$params.genome")    
  custom_reference.copyTo("$params.indir/reference/")

  prefixRef = custom_reference.name.take(custom_reference.name.lastIndexOf('.'))

  process Custom_genome_indexing {
    tag "Indexes reference file using the specified aligner"
    label 'med_mem'

    publishDir "$params.indir/custom_reference", mode: 'copy'

    input:
    file reference_file from ch_reference

    output:
    file "*"

    script:

    if(params.aln == 'bwa'){

    """
      bwa index ${reference_file[0]}
    """

    }else if(params.aln == 'bowtie2'){

    """
      bowtie2-build ${reference_file[0]} ${prefixRef}.bowtie2
    """
    }
  }
}


process FASTQ_Trimming {

  tag "Quality checks and trimms reads using a trimmomatic"
  label 'med_mem'

  publishDir "$params.outdir/trimmed_data"

   input:
   set sampleId, file(samples) from ch_samples
   file adapter from ch_adapter

   output:
   set sampleId, file('*.fastq.gz') into ch_alignment

   script:
   //Makes filtering with an adapter an option without using an annoying amount of different conditionals and scripts.
   def adapter_trimm = adapter.name != 'NO_FILE' ? "ILLUMINACLIP:$adapter:2:30:10" : ''

   //WARNING: Trimmomatic does not mind if the adapter file is not found
   //it will continue without processing it and without a warning.
   if(params.paired){
      // def adapter_trimm does not work with template command.
      //template 'trimmomatic/trimmomatic_PE_adapter_test'

      """
      trimmomatic PE -threads ${params.threads} ${samples[0]} ${samples[1]} ${sampleId[0]}_R1.fastq.gz bad_1 ${sampleId[0]}_R2.fastq.gz bad_2 ${adapter_trimm} SLIDINGWINDOW:15:${params.minqual} MINLEN:${params.minlen}
      """
    //  }
   }
   else{

      """
      trimmomatic SE -threads ${params.threads} ${samples[0]} ${sampleId[0]}_trimmed.fastq.gz ${adapter_trimm} SLIDINGWINDOW:15:${params.minqual} MINLEN:${params.minlen}
      """
   }
  }

//------------------------------------------------------------Alignment----------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

def reference = params.genome != "GRCh37" && params.genome != "GRCh38" ? "${params.working_dir}/${params.indir}/custom_reference/${prefixRef}.fasta": file("${params.indexRef}")
def bowtie2ref = params.genome != "GRCh37" && params.genome != "GRCh38" ? "${params.working_dir}/${params.indir}/custom_reference/${prefixRef}.bowtie2": "${params.indexRef}"
  
process Alignment {

  tag "Aligns reads to a reference genome using bwa or bowtie2"
  label 'med_mem'
   
  publishDir "$params.outdir/alignment"

   input:
   set sampleId, file(fastq_file) from ch_alignment

   output:
   set sampleId, file('*.sam') into ch_sam_to_bam

   script:

   if(params.paired){

      if(params.aln == 'bwa'){

      RG = "@RG\\tID:${sampleId[1]}\\tSM:${sampleId[0]}\\tPL:${params.platform}"

      """
      bwa mem ${params.mappingOptions}-o ${sampleId[0]}.bwa.sam -R "${RG}" -t ${params.threads} ${reference} ${fastq_file[0]} ${fastq_file[1]} 
      """

      }else if(params.aln == 'bowtie2'){

      RG = " --rg-id ${sampleId[1]} --rg SM:${sampleId[0]} --rg PL:${params.platform}"

      """
      bowtie2 -p ${params.threads} -x ${bowtie2ref} -1 ${fastq_file[0]} -2 ${fastq_file[1]} -S ${sampleId[0]}.bowtie2.sam ${RG} ${params.mappingOptions}
      """
      }else if(params.aln == 'novoalign'){

      """
      novoalign --version
      """
      }
    }else{

      if(params.aln == 'bwa'){

      RG = "@RG\\tID:${sampleId[1]}\\tSM:${sampleId[0]}\\tPL:${params.platform}"

      """
      bwa mem ${params.mappingOptions} -o ${sampleId[0]}.bwa.sam -R "${RG}" -t ${params.threads} ${reference} ${fastq_file[0]}
      """

      }else if(params.aln == 'bowtie2'){

      RG = " --rg-id ${sampleId[1]} --rg SM:${sampleId[0]} --rg PL:${params.platform}"

      """
      bowtie2 -p ${params.threads} -x ${bowtie2ref} -1 ${fastq_file[0]} -S ${sampleId[0]}.bowtie2.sam ${RG} ${params.mappingOptions}
      """
      }else if(params.aln == 'novoalign'){

      """
      novoalign --version
      """
      }


    }
}

process SAM_to_BAM{

  tag "Converts SAM file to BAM file using samtools view"
  label 'med_mem'

  publishDir "$params.outdir/alignment"

  input:
  set sampleId, file(sam_file) from ch_sam_to_bam

  output:
  set sampleId, file('*') into ch_bam_sorting

  script:

  """
    samtools view --threads ${params.threads} -b -o ${sampleId[0]}.${params.aln}.bam ${sam_file[0]}
  """
  }

process BAM_sorting{

  tag "Sorts BAM file using Samtools"
  label 'med_mem'

  publishDir "$params.outdir/alignment", mode: 'copy'

  input:
  set sampleId, file(bam_file) from ch_bam_sorting

  output:
  //set sampleId, file('*') into ch_remove_duplicates
  set sampleId, file('*') into ch_bam_final

  script:

  """
  samtools sort --threads ${params.threads} -o ${sampleId[0]}.${params.aln}.sort.bam ${bam_file[0]}
  """
  }

//---------------------------------------------------Mark Duplicates or not-----------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

if(params.remove_duplicates){
  ch_remove_duplicates = Channel.create()
  ch_bam_final.into(ch_remove_duplicates)
}else{
  ch_index_bam = Channel.create()
  ch_bam_final.into(ch_index_bam)
}

if(params.remove_duplicates){

  process Remove_duplicates{
    tag "Remove duplicates if 'remove_duplicates' is true using MarkDuplicates"
    label 'med_mem'

    publishDir "$params.outdir/alignment"

     input:
     set sampleId, file(bam_file) from ch_remove_duplicates

     output:
     set sampleId, file('*.bam') into ch_index_bam

     //when:
     //params.remove_duplicates

     script:

     """
     gatk MarkDuplicates -I ${bam_file[0]} -M ${sampleId[0]}.metrix.dups -O ${sampleId[0]}.${params.aln}.sort.rmdups.bam
     """

      }
  }

process BAM_file_indexing{

  tag "Indexing BAM file"
  label 'med_mem'

  publishDir "$params.outdir/alignment", mode: 'copy'

  input:
  set sampleId, file(bam_file) from ch_index_bam

  output:
  set sampleId,  file("${bam_file[0]}"), file('*.bai') into ch_bamFilesForBaseRecalibration

  script:

   """
   samtools index -@ ${params.threads} ${bam_file[0]}
   """

  }

//---------------------------------------------------Recalibration-----------------------------------------------
//---------------------------------------------------------------------------------------------------------------
/* Is this still relevant?
(ch_gather_bams) = ( params.dbSNP == 'NO_FILE'
                    ? [ ch_bamFilesForBaseRecalibration ]
                    : [ Channel.empty() ] )*/



if(params.dbSNP != 'NO_FILE'){

  process BaseRecalibrator {
    tag "Calculate Base recalibration table"
    label 'med_mem'

    publishDir "$params.outdir/alignment"
    

    input:
      set sampleId, file(bam_files), file(bai_file) from ch_bamFilesForBaseRecalibration
      //set sampleId, file(dbSNP) from ch_dbSNP
    output:
      set sampleId,  file("${bam_files[0]}"), file("${bai_file[0]}"), file("*table") into ch_bamFilesForApplyBQSR

    script:

    """
    gatk BaseRecalibrator ${region_interval} -I ${bam_files[0]} -known-sites ${params.dbSNP} -output ${sampleId[0]}.BQSR.table -reference ${params.seqRef}
    """
    }

  process ApplyBQSR {
    tag "Apply previously recalibrated table"
    label 'med_mem'
    publishDir "$params.outdir/alignment/final", mode: 'copy'

    input:
      set sampleId,file(bam),file(bai),file(bqsr) from ch_bamFilesForApplyBQSR
      //set sampleId, file(varBQSR) from ch_BQSR

    output:
      set sampleId, file('*.bam'), file('*.bai') into ch_variant_calling

    def rmdups = params.remove_duplicates == true ? ".rmdups":""
    
    script:
      """
      gatk ApplyBQSR --bqsr-recal-file ${bqsr[0]} -I ${bam[0]} ${region_interval} -O ${sampleId[0]}.${params.aln}.sort${rmdups}.bqsr.bam

      """

    }  
}else{ch_bamFilesForBaseRecalibration.set{ch_variant_calling}}


if(params.skip_variant_calling){}else{

  process Variant_Calling_single {
    tag "Variant calling using selected Variant Caller (GATK, freebayes, varscan)"
    label 'med_mem'
    //publishDir "$params.outdir/raw_variant_calling_files", mode: 'copy'

    input:
      set sampleId, file(bam_file),file(bai_file) from ch_variant_calling //.combine(ch_variant_calling2)
    output:
      set sampleId, file('*vcf') into ch_vcf
    //reference = file(params.seqRef)

    script:

      if(params.vc == 'gatk'){

      """
      gatk HaplotypeCaller --native-pair-hmm-threads ${params.threads} ${params.rmDups_GATK} ${region_interval} -I ${bam_file[0]} -O ${sampleId[0]}.${params.vc}.vcf -R ${reference} ${params.vcOpts}
      """
      }else if(params.vc == 'freebayes'){

      """
      freebayes --min-alternate-fraction ${params.min_alt_fraction} -f ${reference} ${bam_file[0]} > ${sampleId[0]}.${params.vc}.vcf
      """
      }else if(params.vc == 'varscan'){
      """
      samtools mpileup -B -f ${reference} ${bam_file[0]} | varscan mpileup2cns --variants --output-vcf 1 > ${sampleId[0]}.${params.vc}.vcf
      
      """
      }
    }  

  process VCF_indexing {
    tag "Indexes vcf files generated by Variant_Calling"
    label 'med_mem'

    publishDir "$params.outdir/raw_variant_calling_files", mode: 'copy'

    input:
      set sampleId, file(vcf_file) from ch_vcf
    output:
      set sampleId, file("${vcf_file[0]}"), file('*.vcf.idx')

    script:
    """
    gatk IndexFeatureFile --input ${vcf_file[0]}
    """
  }
}
 //skiping main conditional key

/*
 * result.view { it.trim() }
 */
