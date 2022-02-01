#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getFastq;
    alignReadsToReference;
    convertSamToBam;
    //sortBam;
    //indexBam;
    buildBamIndex;
    markDuplicates;
    fixBamTags;
    //indexBam as indexMarkedBam;
    recalibrateBaseQualityScores;
    //applyBaseQualityRecalibrator;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nAlignment workflow begins here\n"
    fastq = getFastq()
    sam = alignReadsToReference(fastq)
    bam = convertSamToBam(sam)

    /*--- GATK SPARK PIPELINES ---*/
    //indexedBam = buildBamIndex(bam)
    markedBam = markDuplicates(bam)
    fixedBam = fixBamTags(markedBam)
    recalBam = recalibrateBaseQualityScores(fixedBam).view()

/*
*    sortedBam = sortBam(bam)
*    indexedBam = indexBam(sortedBam)
*    markedBam = markDuplicates(indexedBam)
*    MarkedIndexedBam = indexMarkedBam(markedBam)
*    recalTable = recalibrateBaseQualityScores(MarkedIndexedBam)
*    MarkedIndexedBam.combine(recalTable, by: 0).set { applyBQSR_input }
*    applyBaseQualityRecalibrator(applyBQSR_input)
*/

}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
