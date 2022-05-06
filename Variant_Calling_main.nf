#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""

	================================================================
	V A R I A N T  C A L L E R  - I R Y C I S  - DS2    v 0.1
	================================================================

    Usage:
    The typical command for running the pipeline is as follows:
    ./Variant_Calling_main.nf --BAM_input_dir [DIR] --genome [FILE]

    Options:

    --BAM_input_dir [DIR]
    --outdir [DIR]                          The output directory where the results will be saved (default: "my-results")
    --genome <GRCh37 | GRCh38 | [FILE]>     Reference genome to undergo the maping. Options: GRCh37, GRCh38, [/path/to/reference.fasta] (default: GRCh37)
    --region_intervals [BED FILE]           Complete path to specific genomic region in .list format (without chr) to constrict mapping and variant calling. Necessary for Whole Exome Sequencing and Panels. (default: NO_FILE)
    --dbSNP [FILE]                          Automatically provided when selecting GRCh37 or GRCh38, if a custom reference is used and a custom_dbSNP is not provided base recalibration will not be performed. (default: NO_FILE)

    --vc  <freebayes | gatk | varscan>      Variant caller to use for variant calling. Overrideed by --skip_variant_calling Options: gatk, freebayes, varscan (default:gatk)
                                            By default the variant caller will execute in single sample mode. For joint variant calling use jointVariantCalling.nf pipeline.
    --min_alt_fraction [NUM]                Freebayes specific option, minimumn threshold at which allele frequency is considered real. (default: 0.2)

""".stripIndent()
}

// Show help message if --help specified
if (params.help){
helpMessage()
exit 0
}

log.info """\

================================================================
V A R I A N T  C A L L E R  - I R Y C I S - DS2    v 0.1
================================================================
BAM_input_dir           : $params.BAM_input_dir
genome reference        : $params.seqRef
variant caller          : $params.vc
min_alt_fraction        : $params.min_alt_fraction
outdir                  : $params.outdir
================================================================
"""

include { Variant_Calling_single } from './modules/Variant_Calling.nf'

workflow {
    //data cannot by symlinks
    data = channel.fromFilePairs("${params.BAM_input_dir}/*.{bam,bai}", checkIfExists: true).take( params.dev ? params.number_of_inputs : -1 )
    Variant_Calling_single(data)
    data.subscribe onNext: { println it }, onComplete: { println 'Done' }

}