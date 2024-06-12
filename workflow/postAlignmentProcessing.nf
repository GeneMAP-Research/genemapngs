#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//nextflow.enable.moduleBinaries = true

include {
    getInputBams;
    sortBamByName;
    convertBamToFastq;
    convertSamToBam;
    convertBamToCram;
    sortBam;
    indexBam;
    indexAndCopyBam;
    markDuplicatesGatk;
    markDuplicates;
    indexBam as indexMarkedBam;
    recalibrateBaseQualityScores;
    applyBaseQualityRecalibrator;
    markDuplicatesSpark;
    fixBamTags;
    fixAlignmentMate;
    recalibrateBaseQualityScoresSpark;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nAlignment workflow begins here\n"
    if( params.input_ftype == "FASTQ" ) {
        error: "INPUT FILE TYPE MUST BE BAM\n"
    }
    else if( params.input_ftype == "BAM" ) {
        println "INPUT FILE TYPE IS BAM\n"
        inputbam = getInputBams()
        bam = fixAlignmentMate(inputbam)
    }
    else { error "\nERROR: You must specify a file type! Options are FASTQ and BAM (case sensitive)\n" }

    //cram = convertBamToCram(fixedsam)

    if(params.buildVersion == 't2t') {
        if(params.sparkMode == false) {
            sortedBam = sortBam(bam)
            //indexedBam = indexBam(sortedBam)
            markedBam = markDuplicates(indexedBam)
            markedIndexedBam = indexMarkedBam(markedBam)
        }
        else {
            markedBam = markDuplicatesSpark(bam)
            fixedBam = fixBamTags(markedBam)
        }
    } 
    else {
        if(params.sparkMode == false) {
            sortedBam = sortBam(bam)
            //indexedBam = indexBam(sortedBam)
            markedBam = markDuplicates(sortedBam)
            markedIndexedBam = indexAndCopyBam(markedBam)
            //cram = convertBamToCram(sortedBam)
            //markedIndexedBam = indexMarkedBam(finalBam)
            recalTable = recalibrateBaseQualityScores(markedIndexedBam)
            markedIndexedBam.combine(recalTable, by: 0).set { applyBQSR_input }
            recalBam = applyBaseQualityRecalibrator(applyBQSR_input)
            cram = convertBamToCram(recalBam)
        }
        else {
            markedBam = markDuplicatesSpark(bam)
            fixedBam = fixBamTags(markedBam)
            recalBam = recalibrateBaseQualityScoresSpark(fixedBam)
            cram = convertBamToCram(recalBam)
        }
    }
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
