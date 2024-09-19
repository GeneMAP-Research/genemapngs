#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//nextflow.enable.moduleBinaries = true

include {
    getInputFastqs;
    getSEInputFastqs;
    getAlignmentDir;
    sortAlignmentByName;
    sortCram;
    indexAlignment;
    fixAlignmentTags;
    fixAlignmentMate;
    mergeMultiLaneAlignment;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nMERGE ALIGNMENT FILES\n"
    alignment = getAlignmentDir().view()
    mergeMultiLaneAlignment(alignment).view()
}

workflow.onComplete { 
    println "\nDone! Check results in ${params.output_dir}\n" 
}
