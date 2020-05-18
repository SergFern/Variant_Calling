

ch_samples_qc = .fromPath(params.reads, checkIfExists: true)

//------------------------------------------------------------FASTQ Quality Check-----------------------------------------------


File f = new File("${params.outdir}/quality_checks");
f.mkdirs();

process fastq_Quality_Check {
  tag "Quality check on raw fastq files", mode: 'copy'

  // publishDir "$params.outdir/quality_checks"

  input:
    set sampleId, samples from ch_samples_qc

  script:
  
  """
  fastqc -t ${params.threads} ${samples[0]} ${samples[1]} -o ${params.working_dir}/${params.outdir}/quality_checks
  """

}

//---------------------------------------------BamQC & MultiQC----------------------------------------

ch_bamqc = .fromPath("${params.outdir}/alignment")

if(!params.skipBamQC){
  process bam_Quality_Check{    
    tag "Quality check on processed bamqc files"
    publishDir "$params.outdir/quality_checks", mode: 'copy'

    input:
      set sampleId, bam_file, bai_file from ch_bamqc

    output:

    file("*.tag") into ch_goMultiQC
    file "*stats"

    //when:
    //!params.skipBamQC

    def gff = params.custom_gff != 'NO_FILE' ? "-gff ${params.custom_gff} ":''

    params.tagBam = true

    script:

      """
      bamqc -bam ${bam_file[0]} -nt ${params.threads} ${gff}-outdir ${params.working_dir}/${params.outdir}/quality_checks --java-mem-size=${params.java_mem}
      echo multiqc_go > multiqc_go.tag
      """
  }
}


process Generating_integral_QC_report{

  tag "Make a general report integrating all QC"
  publishDir "$params.outdir/quality_checks", mode: 'copy'
  //when:
  //params.tagBam && params.tagFastQC

  input:
   file('tag') from ch_goMultiQC

  output:

    file "*report.html"
    file "*data"

  script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc --force ${params.working_dir}/${params.outdir}/quality_checks
    """
}