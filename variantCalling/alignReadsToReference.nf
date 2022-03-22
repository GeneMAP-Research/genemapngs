#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputFastqs;
    getInputBams;
    sortBamByName;
    convertBamToFastq;
    alignReadsToReference;
    dragenAligner;
    convertSamToBam;
    sortBam;
    indexBam;
    markDuplicates;
    indexBam as indexMarkedBam;
    recalibrateBaseQualityScores;
    applyBaseQualityRecalibrator;
    markDuplicatesSpark;
    fixBamTags;
    recalibrateBaseQualityScoresSpark;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nAlignment workflow begins here\n"
    if( params.inputFileType == "FASTQ" ) {
        println "INPUT FILE TYPE IS FASTQ\n"
        fastq = getInputFastqs()
    }
    else if( params.inputFileType == "BAM" ) {
        println "INPUT FILE TYPE IS BAM\n"
        bam = getInputBams()
        bamSortedByName = sortBamByName(bam)
        fastq = convertBamToFastq(bamSortedByName)
    }
    else { error "\nERROR: You must specify a file type! Options are FASTQ and BAM (case sensitive)\n" }

    if( params.aligner == "DRAGMAP" ) {
        sam = dragenAligner(fastq)
    }
    else {
        sam = alignReadsToReference(fastq)
    }

    bam = convertSamToBam(sam)

    if(params.sparkMode == false) {
        sortedBam = sortBam(bam)
        indexedBam = indexBam(sortedBam)
        markedBam = markDuplicates(indexedBam)
        markedIndexedBam = indexMarkedBam(markedBam)
        recalTable = recalibrateBaseQualityScores(markedIndexedBam)
        markedIndexedBam.combine(recalTable, by: 0).set { applyBQSR_input }
        applyBaseQualityRecalibrator(applyBQSR_input)
    }
    else {
        markedBam = markDuplicatesSpark(bam)
        fixedBam = fixBamTags(markedBam)
        recalBam = recalibrateBaseQualityScoresSpark(fixedBam).view()
    }
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
