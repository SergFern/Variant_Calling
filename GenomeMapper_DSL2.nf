#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""

	================================================================
	V A R I A N T  C A L L E R  - I R Y C I S    v 1.2
	================================================================

    Usage:
    The typical command for running the pipeline is as follows:
    ./nextflow_GenomeMapper.nf [OPTIONS]

    Options:

      --indir [DIR]                           The input directory, all fastq files or csv files in this directory will be processed. (default: "data")
      --outdir [DIR]                          The output directory where the results will be saved (default: "my-results")
      --genome <GRCh37 | GRCh38 | [FILE]>     Reference genome to undergo the maping. Options: GRCh37, GRCh38, [/path/to/reference.fasta] (default: GRCh37)
      --adapters [FILE]                       Adapter file to trimm reads by. (Trimmomatic adapters provided in $baseDir/adapter)
      --region_intervals [BED FILE]           Complete path to specific genomic region in .list format (without chr) to constrict mapping and variant calling. Necessary for Whole Exome Sequencing and Panels. (default: NO_FILE)
      --dbSNP [FILE]                          Automatically provided when selecting GRCh37 or GRCh38, if a custom reference is used and a custom_dbSNP is not provided base recalibration will not be performed. (default: NO_FILE)

      --paired <true | false>                 Execute pipleine in single-end or paired-end mode. If "--paired true" then all fastq files in $params.indir will be processed as samples from the same experiment.
                                              If "--paired false" a csv with the single-end files path and their IDs will be used to identify the fasq files. Options: true, false (default: true)

      --reads [GLOB]                          Glob pattern to identify paired-end reads or the single-end csv file. All data must be compressed. (default: "$baseDir/data/*R{1,2}*.fastq.gz (if paired) "$baseDir/data/*.csv" (if single-end))
                                              CSV format: SampleID, [path/to/read].fastq
      --aln <bwa | bowtie2>                   Aligner chosen to map reads to the reference. Options: bwa, bowtie2 (default: bwa)
      --vc  <freebayes | gatk | varscan>      Variant caller to use for variant calling. Overrideed by --skip_variant_calling Options: gatk, freebayes, varscan (default:gatk)
                                              By default the variant caller will execute in single sample mode. For joint variant calling use jointVariantCalling.nf pipeline.
      --common_id "[STRING]"                  Id by which to identify all samples as coming from the same experiment. Assumed to be leading the file name. (default: first two characters of file name are used as experiment identifier)

      --skip_variant_calling <true | false>   Skips variant calling process entirely, only perform the alignment (default: true)
      --remove_duplicates <true | false>      Remove marked as duplicated reads. Options: true, false (default: false)
      --min_alt_fraction [NUM]                Freebayes specific option, minimumn threshold at which allele frequency is considered real. (default: 0.2)

    """.stripIndent()
}

// Show help message if --help specified
if (params.help){
  helpMessage()
  exit 0
}

if(params.paired){

  //Reads must be read twice, both as a tupple and as an array
   Channel
      .fromFilePairs(params.reads, size: 2, checkIfExists: true, flat:true)
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
      .take( params.dev ? params.number_of_inputs : -1 )
      .ifEmpty {error "File ${params.reads} not parsed properly"}
      .into{ [ch_samples, ch_sampleName, ch_FlowCell_lane] }

   ch_FlowCell_lane.splitFastq(record: true, pe: true, limit: 1).map{it -> new Tuple(it[0].split("_")[0], it[1].readHeader.split(":")[2,3,9].join("."))}.set{ch_RG_ID}

 }
//Fuse Ids and samples to manage SampleId as a tuple together with the samples they identify.
ch_RG_ID.concat(ch_samples_with_id).groupTuple().map{ it -> [[it[0],it[1][0]],it[1][1]] }.set{ ch_samples }

ch_dbSNP = file(params.dbSNP)

def region_interval = params.region_intervals != 'NO_FILE' ? "-L ${params.region_intervals} -ip 100 ":''
def ploidy = params.ploidy != 'n' || params.ploidy == 'y' && params.ploidy.getClass() == java.lang.Integer ? "--ploidy ${params.ploidy} ":''

log.info """\

================================================================
V A R I A N T  C A L L E R  - I R Y C I S    v 1.2
================================================================
genome               : $params.genome
reads                : $params.reads
adapters             : $params.adapters

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

