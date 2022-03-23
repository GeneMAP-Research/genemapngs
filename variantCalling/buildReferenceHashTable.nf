#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    buildRefHashTable
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    fastaRef = params.fastaRef
    hash = buildRefHashTable(fastaRef).view()
}

workflow.onComplete { println "\nDone!\n" }
