#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputFastqs;
    getInputAlignments;
} from "${projectDir}/modules/alignmentPipeline.nf"

include {
    getFastqQualityReports;
    getAlignmentQualityReports;
    getMultiQcFastqReports;
} from "${projectDir}/modules/getReadsQualityReports.nf"

workflow {
    println "\nFASTQ/BAM/CRAM Quality Reports\n"
    if( params.input_ftype.toUpperCase() == "FASTQ" ) {
        println "INPUT FILE TYPE IS FASTQ\n"
        fastq = getInputFastqs()
        fastqReports = getFastqQualityReports(fastq).collect()
        getMultiQcFastqReports(fastqReports)
    }
    else {
        println "INPUT FILE TYPE IS ALIGNMENT [BAM/CRAM]\n"
        alignment = getInputAlignments()
        alignmentReports = getAlignmentQualityReports(alignment).collect()
        getMultiQcFastqReports(alignmentReports)
    }
}

