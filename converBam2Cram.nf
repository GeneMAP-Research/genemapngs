#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputBams;
    sortBamByName;
    convertBamToCram;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nBAM2CRAM\n"
    bam = getInputBams()
    cram = convertBamToCram(bam)
}

workflow.onComplete { 
    println "\nDone! Check results in ${params.outputDir}\n" 
}
