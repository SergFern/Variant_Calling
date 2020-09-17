#!/usr/bin/env nextflow

//TODO:
//meter mapping interval
//Provide an error if adapter filtering was usuccessful

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
      --genome                         Reference genome to undergo the maping. Options: GRCh37, GRCh38, [/path/to/reference.fasta] (default: GRCh37)
      --adapter_file                   Adapter file to trimm reads by. (Trimmomatic adapters provided in $baseDir/adapter)
      --region_intervals               Specific genomic region in bed format (without chr) to constrict mapping and variant calling. Necessary for Whole Exome Sequencing and Panels. (default: NO_FILE)

      --paired                         Execute pipleine in single-end or paired-end mode. If "--paired true" then all fastq files in $params.indir will be processed as samples from the same experiment.
                                       If "--paired false" a csv with the single-end files path and their IDs will be used to identify the fasq files. Options: true, false (default: true)

      --reads                          Path to paired-end reads or single-end csv file  (default: "$baseDir/$params.indir/*R{1,2}*.fastq.gz (if paired) "$baseDir/$params.indir/*.csv" (if single-end))
                                       CSV format: SampleID, [path/to/read].fastq
      --aln                            Aligner chosen to map reads to the reference. Options: bwa, bowtie2 (default: bwa)
      --vc                             Variant caller to use for variant calling. Options: gatk, freebayes, varscan (default:gatk)
      --remove_duplicates              Remove marked as duplicated reads. Options: true, false (default: false)
      --indir                          The input directory, all fastq files or csv files in this directory will be processed.
      --outdir                         The output directory where the results will be saved (default: $params.outdir)
      

    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

if(params.paired){

   Channel
      .fromFilePairs(params.reads, size: 2)
      .take( params.dev ? params.number_of_inputs : -1 ) //TESTING: should only run partial data
      .set{ ch_samples }

}else{

   Channel
      .fromPath(params.reads, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> [row.sampleId, [row.read]] }
      .take( params.dev ? params.number_of_inputs : -1 ) //TESTING: should only run partial data
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .set{ ch_samples }

 }


ch_dbSNP = file(params.dbSNP)
ch_adapter = file(params.adapter_trimm)

def region_interval = params.region_intervals != 'NO_FILE' ? "-L ${params.region_intervals} -ip 100 ":''
def ploidy = params.ploidy != 'no' || params.ploidy == 'yes' && params.ploidy.getClass() == java.lang.Integer ? "--ploidy ${params.ploidy} ":''

log.info """\

================================================================
V A R I A N T  C A L L E R  - I R Y C I S    v 0.9-Custom Ref
================================================================
genome               : $params.genome
reads                : $params.reads
adapters             : $params.adapter_file

region               : $params.region_intervals

paired               : $params.paired
aligner              : $params.aln
variant_caller       : $params.vc
remove_duplicates    : $params.remove_duplicates
ploidy               : $params.ploidy

read_directory       : ./$params.indir
results              : ./$params.outdir
================================================================
"""

//------------------------------------------------------------Trimming-------------------------------------------------

if(params.genome != "GRCh37" && params.genome != "GRCh38"){

  ch_reference = file(params.genome, checkIfExists: true)
  //(ch_reference_idx, ch_reference_dic) = Channel.from(ch_reference).into(2) apparently there is no need for this?

  process Indexing_custom_genome {
    tag "Indexes supplied reference FASTA file (Samtools)"

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
  }

//------------------------------------------------------------Alignment----------------------------------------------------
def reference = params.genome != "GRCh37" && params.genome != "GRCh38" ? "${params.working_dir}/${params.indir}/custom_reference/${prefixRef}": params.indexRef
  
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
        bwa mem ${params.mappingOptions} -o ${sampleId}.bwa.sam -R "${params.RG}" -t ${params.threads} ${reference}.fasta ${fastq_file[0]} ${fastq_file[1]} 
        """

        }else if(params.aln == 'bowtie2'){

        """
        bowtie2 -p ${params.threads} -x ${reference}.bowtie2 -1 ${fastq_file[0]} -2 ${fastq_file[1]} -S ${sampleId}.bowtie2.sam ${params.RG} ${params.mappingOptions}
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

  publishDir "$params.outdir/alignment"

  input:
  set sampleId, file(sam_file) from ch_sam_to_bam

  output:
  set sampleId, file('*') into ch_bam_sorting

  script:

  """
    samtools view --threads ${params.threads} -b -o ${sampleId}.${params.aln}.bam ${sam_file[0]}
  """
  }

process BAM_sorting{

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

process Remove_duplicates{
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

process BAM_file_indexing{

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

(ch_variant_calling) = ( params.dbSNP == 'NO_FILE'
                    ? [ ch_bamFilesForBaseRecalibration ]
                    : [ Channel.empty() ] )

if(params.dbSNP != 'NO_FILE'){

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
    gatk BaseRecalibrator ${region_interval} -I ${bam_files[0]} -known-sites ${params.dbSNP} -output BQSR.table -reference ${reference}.fasta
    """
    }

  process ApplyBQSR {
    tag "Apply previously recalibrated table"
    publishDir "$params.outdir/alignment", mode: 'copy'

    input:
      set sampleId,file(bam),file(bai),file(bqsr) from ch_bamFilesForApplyBQSR.combine(ch_BQSR)
      //set sampleId, file(varBQSR) from ch_BQSR

    output:
      set sampleId, file('*.bam'), file('*.bai') into ch_variant_calling //ch_variant_calling2

    def rmdups = params.remove_duplicates == true ? ".rmdups":"" //We define a variable rmdups to mark files that which duplicates were removed.
    

    script:
      """
      gatk ApplyBQSR --bqsr-recal-file ${bqsr[0]} -I ${bam[0]} ${region_interval} -O ${sampleId}.${params.aln}.sort${rmdups}.bqsr.bam
      """
    }  
}

process Variant_Calling {
  tag "Variant calling using selected Variant Caller (GATK, freebayes, varscan)"
  //publishDir "$params.outdir/raw_variant_calling_files", mode: 'copy'

  //def reference = params.genome != "GRCh37" && params.genome != "GRCh38" ? "${params.working_dir}/${params.indir}/reference/${prefixRef}.fasta": params.indexRef

  input:
  // Imma need to generate a tsv with all bam files paths listed so I can feed all of them to freebayes and other variant callers together
    set sampleId, file(bam_file),file(bai_file) from ch_variant_calling //.combine(ch_variant_calling2)

  output:
    set sampleId, file('*vcf') into ch_vcf

  def GVCF = params.GVCFmode == 'false' ? "":"-ERC GVCF"
  //reference = file(params.seqRef)


  script:
    /*
    println("${bam_file[0]}")
    println("${bam_file[1]}")
    */

    if(params.vc == 'gatk'){
    """
    gatk HaplotypeCaller --native-pair-hmm-threads ${params.threads} ${params.rmDups_GATK} ${region_interval} -I ${bam_file[0]} -O ${sampleId}.${params.vc}.vcf -R ${reference}.fasta ${GVCF} ${params.vcOpts}
    """
    }else if(params.vc == 'freebayes'){

    """
    freebayes ${ploidy}${params.vcOpts} -f ${reference}.fasta ${bam_file[0]} > ${sampleId}.${params.vc}.vcf
    """
    }else if(params.vc == 'varscan'){
    """
    samtools mpileup -B -f ${reference}.fasta ${bam_file[0]} | varscan mpileup2cns --variants --output-vcf 1 > ${sampleId}.${params.vc}.vcf
    
    """
    }
  }

process VCF_indexing {
  tag "Indexes vcf files generated by Variant_Calling"

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
 //skiping main conditional key

/*
 * result.view { it.trim() }
 */