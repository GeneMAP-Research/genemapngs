#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputBams;
    sortBamByName;
    convertBamToFastq;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nBAM2FASTQ workflow begins here\n"
    bam = getInputBams()
    bamSortedByName = sortBamByName(bam)
    fastq = convertBamToFastq(bamSortedByName)
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
