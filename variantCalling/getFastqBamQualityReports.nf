#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputFastqs;
    getInputBams;
} from "${projectDir}/modules/alignmentPipeline.nf"

include {
    getFastqQualityReports;
    getBamQualityReports;
    getMultiQcFastqReports;
} from "${projectDir}/modules/getReadsQualityReports.nf"

workflow {
    println "\nFASTQ/BAM Quality Reports begins here\n"
    if( params.inputFileType == "FASTQ" ) {
        println "INPUT FILE TYPE IS FASTQ\n"
        fastq = getInputFastqs().view()
        fastqReports = getFastqQualityReports(fastq).collect()
        getMultiQcFastqReports(fastqReports)
    }
    else if( params.inputFileType == "BAM" ) {
        println "INPUT FILE TYPE IS BAM\n"
        bam = getInputBams()
        bamReports = getBamQualityReports(bam).collect()
        getMultiQcFastqReports(bamReports)
    }
    else { error "\nERROR: You must specify a file type! Options are FASTQ and BAM (case sensitive)\n" }
}

