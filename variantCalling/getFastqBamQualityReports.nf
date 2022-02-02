#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFastq;
    getReadsQualityReports;
    getMultiQcFastqReports;
} from "${projectDir}/modules/computeFastqBamQualityReports.nf"

workflow {
    println "\nFASTQ/BAM Quality Reports begins here\n"
    if( params.inputFileType == "FASTQ" ) {
        println "INPUT FILE TYPE IS FASTQ\n"
        fastq = getInputFastqs()
        fastqReports = getReadsQualityReports(fastq).collect().view()
        //getMultiQcFastqReports(fastqReports)
    }
    else if( params.inputFileType == "BAM" ) {
        println "INPUT FILE TYPE IS BAM\n"
        bam = getInputBams()
        fastqReports = getReadsQualityReports(bam).collect().view()
        //getMultiQcFastqReports(bamReports)
    }
    else { error "\nERROR: You must specify a file type! Options are FASTQ and BAM (case sensitive)\n" }
}

