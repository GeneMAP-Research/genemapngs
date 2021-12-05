#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getFastq;
    alignReadsToReference;
    convertSamToBam;
    sortBam;
    indexBam;
    markDuplicates;
    indexBam as indexMarkedBam;
    recalibrateBaseQualityScores;
    applyBaseQualityRecalibrator;
    //getFinalBam
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nAlignment begins here\n"
    fastq = getFastq()
    sam = alignReadsToReference(fastq)
    bam = convertSamToBam(sam)
    sortedBam = sortBam(bam)
    indexedBam = indexBam(sortedBam)
    markedBam = markDuplicates(indexedBam)
    MarkedIndexedBam = indexMarkedBam(markedBam)
    recalTable = recalibrateBaseQualityScores(MarkedIndexedBam)
    MarkedIndexedBam.combine(recalTable, by: 0).set { applyBQSR_input }
    applyBaseQualityRecalibrator(applyBQSR_input)
    //getFinalBam(MarkedIndexedBam)
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
