#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFastq;
    getFastqQualityReports
    getMultiQcFastqReports;
} from "${projectDir}/modules/computeFastqBamQualityReports.nf"

workflow {
    println "\nFASTQ/BAM Quality Reports begins here\n"

    bamFastq = getBamFastq().view()
    fastqc_reports = getFastqQualityReports(bamFastq).collect().view()
    //getMultiQcFastqReports(fastqc_reports)
}

