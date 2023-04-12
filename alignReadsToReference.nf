#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputFastqs;
    getSEInputFastqs;
    getInputBams;
    sortBamByName;
    convertBamToFastq;
    alignReadsToReference;
    dragenAligner;
    tmapAligner;
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
        if(params.pe == true) {
            println "PAIRED END READS\n"
            fastq = getInputFastqs()
        }
        else {
            println "SINGLE END READS\n"
            fastq = getSEInputFastqs().view()
        }
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
    else if( params.aligner == "TMAP" ) {
        sam = tmapAligner(fastq)
    }
    else {
        sam = alignReadsToReference(fastq)
    }

    bam = convertSamToBam(sam)

    if(params.buildVersion == 't2t') {
        if(params.sparkMode == false) {
            sortedBam = sortBam(bam)
            indexedBam = indexBam(sortedBam)
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

}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
