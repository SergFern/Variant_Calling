#!/usr/bin/env nextflow

//TODO:
//Add adapter filtering - done
//Provide an error if adapter filtering was usuccessful

if(params.paired){

   Channel
      .fromFilePairs(params.reads)
      .set{ samples_ch }

}else{

   Channel
      .fromPath(params.reads, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> [row.sampleId, [row.read]] }
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .set { samples_ch }

}

adapter_ch = file(params.adapter_trimm)

//------------------------------------------------------------Trimming-------------------------------------------------

process trimmomatic{

   //Channel
   //   .fromPath(params.adapter_trimm, checkIfExists: true)
   //   .set{ adapter_ch }

      publishDir "$params.outdir/trimmed_data"

   input:
   set sampleId, samples from samples_ch
   file adapter from adapter_ch

   output:
   file '*.fastq.gz'
   set sampleId, file('*.fastq.gz') into alignment_ch

   script:
   //def adapter_trimm to make filtering with an adapter an option without using an annoying amount of different conditionals and scripts.
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
}

//------------------------------------------------------------Alignment----------------------------------------------------

process alignment {
   
 publishDir "$params.outdir/alignment"

 input:
 set sampleId, file(fastq_file) from alignment_ch
 //file indexRef from alignment_index_ch

 output:
 file '*.sam'
 set sampleId, file('*.sam') into sam_to_bam_ch

 script:

 //def indexRef = file(params.indexRef)

 if(params.paired){

    if(params.aln == 'bwa'){

      """
      bwa mem ${params.mappingOptions} -o ${sampleId}.bwa.sam -t ${params.threads} ${params.indexRef} ${fastq_file[0]} ${fastq_file[1]}
      """

    }else if(params.aln == 'bowtie2'){

      """
      bowtie2 ${params.mappingOptions} -p ${params.threads} -x ${params.indexRef} -1 ${fastq_file[0]} -2 ${fastq_file[1]} -S ${sampleId}.bowtie2.sam
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

  input:
  set sampleId, file(sam_file) from sam_to_bam_ch

  output:
  set sampleId, file("*${rmdups}.sort.bam"), file('*.bam.bai') into bam_out_ch

  def rmdups = params.remove_duplicates == true ? ".rmdups":""

  script:

  """
    samtools view --threads ${params.threads} -b -o ${sampleId}.${params.aln}.bam ${sam_file[0]}
    samtools sort --threads ${params.threads} -o ${sampleId}.${params.aln}.sort.bam ${sampleId}.${params.aln}.bam
  """

    if(params.remove_duplicates){

      """
        gatk MarkDuplicates -I ${sampleId}.${params.aln}.sort.bam -M ${sampleId}.metrix.dups -O ${sampleId}.${params.aln}${rmdups}.sort.bam        
      """
    }
}

process indexing{
  tag "Indexing BAM files"

  publishDir "$params.outdir/alignment"

  input:
    set sampleId, file(bam_file) from bam_out_ch
  output:
    set sampleId, file("*.bam"), file("*.bam.bai") into bam_final_ch

  script:

  """
  samtools index -@ ${params.threads} ${sampleId}.${params.aln}${rmdups}.sort.bam
  """
}


//---------------------------------------------------Recalibration-----------------------------------------------


/*
 * result.view { it.trim() }
 */

