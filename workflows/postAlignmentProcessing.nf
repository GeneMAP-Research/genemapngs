#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//nextflow.enable.moduleBinaries = true

include {
    getAlignment;
    getAlignmentDir;
    markDuplicatesGatk;
    markDuplicates;
    markDupSambam;
    indexAlignment as indexMarkedAlignment;
    recalibrateBaseQualityScores;
    applyBaseQualityRecalibrator;
    markDuplicatesSpark;
    fixAlignmentTags;
    fixAlignmentMate;
    recalibrateBaseQualityScoresSpark;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow BQSR {
    println "\nAlignment workflow begins here\n"

    println "INPUT FILE TYPE IS ALIGNMENT (BAM/CRAM)\n"
    alignment = getAlignmentDir()

    if(!(params.build == 't2t')) {
        recalTable = recalibrateBaseQualityScores(alignment)
        alignment.combine(recalTable, by: 0).set { applyBQSR_input }
        recalibrated = applyBaseQualityRecalibrator(applyBQSR_input)
    }
    else {
        recalTable = recalibrateBaseQualityScores(alignment)
        alignment.combine(recalTable, by: 0).set { applyBQSR_input }
        recalibrated = applyBaseQualityRecalibrator(applyBQSR_input)
    }

}

//workflow.onComplete { 
//    println "\nDone! Check results in ${params.output_dir}\n" 
//}
