#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputFastqs;
    getInputBams;
    sortBamByName;
    convertBamToFastq;
    alignReadsToReference;
    convertSamToBam;
    sortBam;
    indexBam;
    buildBamIndex;
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
        fastq = getInputFastqs().view()
    }
    else if( params.inputFileType == "BAM" ) {
        println "INPUT FILE TYPE IS BAM\n"
        bam = getInputBams().view()
        bamSortedByName = sortBamByName(bam)
        fastq = convertBamToFastq(bamSortedByName)
    }

    sam = alignReadsToReference(fastq)
    bam = convertSamToBam(sam)

    if (params.sparkMode == false) {
        sortedBam = sortBam(bam)
        indexedBam = indexBam(sortedBam)
        markedBam = markDuplicates(indexedBam)
        MarkedIndexedBam = indexMarkedBam(markedBam)
        recalTable = recalibrateBaseQualityScores(MarkedIndexedBam)
        MarkedIndexedBam.combine(recalTable, by: 0).set { applyBQSR_input }
        applyBaseQualityRecalibrator(applyBQSR_input)
    }
    else {
        markedBam = markDuplicatesSpark(bam)
        fixedBam = fixBamTags(markedBam)
        recalBam = recalibrateBaseQualityScoresSpark(fixedBam).view()
    }
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
