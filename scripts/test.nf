#!/usr/bin/env nextflow

//Interesting script illustrating the use of arrays in channels.

params.java = '/home/chm2059/chm2059/lib/jre1.8.0_25/bin/java'
params.varscan = '/home/chm2059/chm2059/lib/VarScan.v2.3.9.jar'
params.samtools = '/home/chm2059/chm2059/lib/samtools-1.1/samtools'

bams = Channel.from(
    [
        'HL01',
        file('/home/chm2059/chm2059/HL01/HL01.gatk.bam'),
        file('/home/chm2059/chm2059/HL01/HL01.gatk.bam.bai'),
        file('/home/chm2059/chm2059/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'),
        file('/home/chm2059/chm2059/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.dict'),
        file('/home/chm2059/chm2059/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai')
    ],
    [
        'HL02',
        file('/home/chm2059/chm2059/HL02/HL02.gatk.bam'),
        file('/home/chm2059/chm2059/HL02/HL02.gatk.bam.bai'),
        file('/home/chm2059/chm2059/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'),
        file('/home/chm2059/chm2059/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.dict'),
        file('/home/chm2059/chm2059/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai')
    ]
)


process varscan_indel {

  storeDir "${baseDir}/Varscan/${prefix}"
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=8G -pe smp 2 -l h_rt=6:00:00 -l os=rhel5.4|rhel6.3'
  stageInMode 'copy'

  input:
    set prefix, file(bam_file), file(bam_index_file), file(ref), file(ref_dict), file(ref_fai) from bams

  output:
    set file('*indels.vcf') into varscan_indel_out

  """
  ${params.samtools} mpileup \
    -d10000 -f ${ref} \
    -q 15 \
    ${bam_file} \
    | ${params.java} -Xmx10g -jar ${params.varscan} mpileup2indel \
    --p-value 0.05 \
    --output-vcf > \
    ${prefix}.varscan.indels.vcf
  """
}

