#!/usr/bin/env nextflow

//TODO:
//meter mapping interval
//Provide an error if adapter filtering was usuccessful

params.outdir = "$baseDir/results"
params.reads = "$baseDir/reads/*R{1,2}*.fastq.gz"

outdir = params.outdir

def helpMessage() {
    log.info"""
    ==================================================================
    ${workflow.manifest.name}  ~  version ${workflow.manifest.version}
    ==================================================================

    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run  -profile standard --outdir /output/path
    Options:
      --reads                          Help message for "reads" (default: "$baseDir/reads/*R{1,2}*.fastq.gz")
    Other options:
      --outdir                      The output directory where the results will be saved (default: $outdir)
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -profile                      Configuration profile to use. [standard, other_profiles] (default 'standard')
    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}


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
aligner              : $params.aln
remove_duplicates    : $params.remove_duplicates
results              : ./$params.outdir
================================================================
"""

//------------------------------------------------------------Trimming-------------------------------------------------

if(params.genome != "GRCh37" && params.genome != "GRCh38"){

  ch_reference = file(params.genome, checkIfExists: true)
  (ch_reference_idx, ch_reference_dic) = Channel.from(ch_reference).into(2)

  process indexing_custom_genome {
    tag "Indexes supplied reference FASTA file (Samtools)"

    publishDir "$params.outdir/reference", mode: 'copy'

    input:
    file reference from ch_reference

    output:
    file '*.fai'

    script:
    
    """
    samtools faidx ${reference[0]}
    """
    
   }

   process building_dictionary {
    tag "Creating Dictionary file (GATK)"

    publishDir "$params.outdir/reference", mode: 'copy'

    input:
    file reference from ch_reference

    output:
    file '*.dict'

    script:
    
    """
    gatk CreateSequenceDictionary -R ${reference[0]}
    """
    
   }

  }



 //skiping main conditional key

/*
 * result.view { it.trim() }
 */


