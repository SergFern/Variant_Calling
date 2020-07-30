#!/usr/bin/env nextflow

//TODO:
//meter mapping interval
//Provide an error if adapter filtering was usuccessful



if(params.paired){

   Channel
      .fromFilePairs(params.reads)
      .set{ ch_samples }

}else{

   Channel
      .fromPath(params.reads, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> [row.sampleId, [row.read]] }
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .set{ ch_samples }

}

ch_dbSNP = file(params.dbSNP)
ch_adapter = file(params.adapter_trimm)

def region_interval = params.region_intervals != 'NO_FILE' ? "-L ${params.region_intervals} -ip 100 ":''

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
aligner              : $params.aln
remove_duplicates    : $params.remove_duplicates
results              : ./$params.outdir
================================================================
"""

//------------------------------------------------------------Trimming-------------------------------------------------

//if(!params.skipMain){ //-Uncomment to skip main section

if(params.genome != "GRCh37" || params.genome != "GRCh38"){

  ch_reference = file(params.genome, checkIfExists: true)

  process indexing_custom_genome {
    tag "Indexes supplied reference file"

    publishDir "$params.outdir/reference", mode: 'copy'

    input:
    file reference from ch_reference

    output:
    file '*.idx'

  }

}


process trimming {

  tag "Quality checks and trimms reads using a trimmomatic"

  publishDir "$params.outdir/trimmed_data"

   input:
   set sampleId, samples from ch_samples
   file adapter from ch_adapter

   output:
   file '*.fastq.gz'
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
      trimmomatic PE ${samples[0]} ${samples[1]} ${sampleId}_R1.fastq.gz bad_1 ${sampleId}_R2.fastq.gz bad_2 ${adapter_trimm} SLIDINGWINDOW:15:${params.minqual} MINLEN:${params.minlen}
      """
    //  }
   }
   else{

      """
      trimmomatic SE ${samples[0]} ${sampleId}_trimmed.fastq.gz ${adapter_trimm} SLIDINGWINDOW:15:${params.minqual} MINLEN:${params.minlen}
      """
   }
//}
}

//------------------------------------------------------------Alignment----------------------------------------------------

process alignment {

  tag "Aligns reads to a reference genome using bwa or bowtie2"
   
   publishDir "$params.outdir/alignment"

   input:
   set sampleId, file(fastq_file) from ch_alignment
   //file indexRef from alignment_index_ch

   output:
   file '*.sam'
   set sampleId, file('*.sam') into ch_sam_to_bam

   script:

   if(params.paired){

      if(params.aln == 'bwa'){

      """
      bwa mem ${params.mappingOptions} -o ${sampleId}.bwa.sam -R "${params.RG}" -t ${params.threads} ${params.indexRef} ${fastq_file[0]} ${fastq_file[1]} 
      """

      }else if(params.aln == 'bowtie2'){

      """
      bowtie2 ${params.mappingOptions} -p ${params.threads} -x ${params.indexRef} -1 ${fastq_file[0]} -2 ${fastq_file[1]} -S ${sampleId}.bowtie2.sam ${params.RG}
      """
      }else if(params.aln == 'novoalign'){

      """
      novoalign --version
      """
      }
   }else{

      if(params.aln == 'bwa'){

      """
      bwa mem ${params.mappingOptions} -o ${sampleId}.bwa.sam -t ${params.threads} ${params.indexRef} ${fastq_file}
      """
      }else if(params.aln == 'bowtie2'){

      """
      bowtie2 ${params.mappingOptions} -p ${params.threads} -x ${params.indexRef} -q ${fastq_file} -S ${sampleId}.bowtie2.sam
      """
      }else if(params.aln == 'novoalign'){

      """
      novoalign --version
      """
      }
  } 
}

process sam_to_bam{

  tag "Converts SAM file to BAM file using samtools view"

  input:
  set sampleId, file(sam_file) from ch_sam_to_bam

  output:
  set sampleId, file('*') into ch_bam_sorting

  script:

  """
    samtools view --threads ${params.threads} -b -o ${sampleId}.${params.aln}.bam ${sam_file[0]}
  """
}

process bam_sorting{

  tag "Sorts BAM file using Samtools"

  publishDir "$params.outdir/alignment"

  input:
  set sampleId, file(bam_file) from ch_bam_sorting

  output:
  //set sampleId, file('*') into ch_remove_duplicates
  set sampleId, file('*') into ch_bam_final

  script:

  """
  samtools sort --threads ${params.threads} -o ${sampleId}.${params.aln}.sort.bam ${bam_file[0]}
  """
}

//---------------------------------------------------Mark Duplicates or not-----------------------------------------------


if(params.remove_duplicates){
  ch_remove_duplicates = Channel.create()
  ch_bam_final.into(ch_remove_duplicates)
}else{
  ch_index_bam = Channel.create()
  ch_bam_final.into(ch_index_bam)
}

if(params.remove_duplicates){

process remove_duplicates{
  tag "Remove duplicates if 'remove_duplicates' is true using MarkDuplicates"

  publishDir "$params.outdir/alignment"

   input:
   set sampleId, file(bam_file) from ch_remove_duplicates

   output:
   set sampleId, file('*.bam') into ch_index_bam

   //when:
   //params.remove_duplicates

   script:

   """
   gatk MarkDuplicates -I ${bam_file[0]} -M ${sampleId}.metrix.dups -O ${sampleId}.${params.aln}.sort.rmdups.bam
   """

   }
}

process indexing{

  tag "Indexing BAM file"

  publishDir "$params.outdir/alignment"

  input:
  set sampleId, file(bam_file) from ch_index_bam

  output:
  set sampleId,  file("${bam_file[0]}"), file('*.bai') into ch_bamFilesForBaseRecalibration
  set sampleId,  file("${bam_file[0]}"), file('*.bai') into ch_bamFilesForApplyBQSR

  script:

   """
   samtools index -@ ${params.threads} ${bam_file[0]}
   """

}

//(ch_bamFilesForBaseRecalibration,
 //ch_bamFilesForApplyBQSR) = ch_indexed_bam.separate(2) { x -> [ x, x ] }

//---------------------------------------------------Recalibration-----------------------------------------------

process BaseRecalibrator {
  tag "Calculate Base recalibration table"
  publishDir "$params.outdir/alignment"
  

  input:
    set sampleId, file(bam_files), file(dbSNP) from ch_bamFilesForBaseRecalibration
    //set sampleId, file(dbSNP) from ch_dbSNP
  output:
    file("*table") into ch_BQSR

  script:

  """
  gatk BaseRecalibrator ${region_interval} -I ${bam_files[0]} -known-sites ${params.dbSNP} -output BQSR.table -reference ${params.seqRef}
  """
}

process ApplyBQSR {
  tag "Apply previously recalibrated table"
  publishDir "$params.outdir/alignment", mode: 'copy'

  input:
    set sampleId,file(bam),file(bai),file(bqsr) from ch_bamFilesForApplyBQSR.combine(ch_BQSR)
    //set sampleId, file(varBQSR) from ch_BQSR

  output:
    file("*.bam")
    file("*.bai")

  def rmdups = params.remove_duplicates == true ? ".rmdups":"" //We define a variable rmdups to mark files that which duplicates were removed.
  

  script:
    """
    gatk ApplyBQSR --bqsr-recal-file ${bqsr[0]} -I ${bam[0]} ${region_interval} -O ${sampleId}.${params.aln}.sort${rmdups}.bqsr.bam
    """
}

 //skiping main conditional key

/*
 * result.view { it.trim() }
 */
