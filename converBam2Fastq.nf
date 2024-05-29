#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputAlignments;
    sortAlignmentByName;
    convertAlignmentToFastq;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nBAM2FASTQ workflow begins here\n"
    bam = getInputAlignments()
    bamSortedByName = sortAlignmentByName(bam)
    fastq = convertAlignmentToFastq(bamSortedByName)
}

workflow.onComplete { println "\nDone! Check results in ${params.output_dir}\n" }
