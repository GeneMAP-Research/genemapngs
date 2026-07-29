#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*********************************************************************************************************
    WORKFLOW RUN MODES:
    - EACH PIPELINE CAN BE CALLED TO RUN A SPECIFIC TASK DEPEINDING ON THE INPUT
      DATA AND DESIRED OUTPUT
        -> QC
        -> TRIM
        -> ALIGN
        -> MERGEALIGN
        -> CALL
        -> FILTER
        -> ANNOTATE

    - USER MAY DESIRE TO PROCESS THEIR RAW DATA IN FASTQ OR ALIGNMENT (BAM/CRAM) FORMAT TO VCF 
      IN ONE RUN IN WHICH CASE THE WORKFLOW WILL START FROM 'ALIGN' TO 'CALL' AND OUTPUT FOR EACH 
      PIPELINE WILL BE SAVED IN USER-SPERCIFIED OUTPUT DIRECTORY.
        -> RAW2VCF

    - USER MAY DESIRE TO PROCESS THEIR RAW DATA IN FASTQ OR ALIGNMENT (BAM/CRAM) FORMAT TO ANNOTATED VCF 
      IN ONE RUN IN WHICH CASE THE WORKFLOW WILL START FROM 'ALIGN' TO 'CALL' AND OUTPUT FOR EACH 
      PIPELINE WILL BE SAVED IN USER-SPERCIFIED OUTPUT DIRECTORY.
        -> RAW2ANNOTATE

    [IMPORTANT]: 
        - THIS WORKFLOW IS SPECIFICALLY DESIGNED TO MAKE 'TRIM' STAND-ALONE TO ALLOW USERS INSPECT THEIR 
          TRIMMED DATA AND PROCEED TO RUN ENTIRE WORKFLOWONLY WHEN CONFIDENT WITH THE FASTQ FILES.

        - WHEN 'RAW2VCF' or 'RAW2ANNOTATE' the workflow will run with default parameters for all 
          but the 'ALIGN' pipeline.



**********************************************************************************************************/

include { 
    jobCompletionMessage;
    checkWkflow;
    checkParams;
    checkFaiIndex;
    checkGatkSeqDict;
    checkBwaIndex;
    checkBwa2Index
} from "${projectDir}/includes/workflowInspection.nf"

include { TEST } from "${projectDir}/workflows/test.nf"
include { QC } from "${projectDir}/workflows/getQualityReports.nf"
include { TRIM } from "${projectDir}/workflows/trimReads.nf"
include { ALIGN } from "${projectDir}/workflows/alignReadsToReference.nf"
include { MERGEALIGN } from "${projectDir}/workflows/mergeMultiLaneAlignments.nf"
include { CALL } from "${projectDir}/workflows/callVariants.nf"
include { FILTER } from "${projectDir}/workflows/filterVariantCalls.nf"
include { ANNOTATE } from "${projectDir}/workflows/annotateVarints.nf"

include { BQSR } from "${projectDir}/workflows/postAlignmentProcessing.nf"

//include { BAM2FASTQ } from "${projectDir}/workflows/converBam2Fastq.nf"
//include { BAM2CRAM } from "${projectDir}/workflows/converBam2Cram.nf"

workflow {

    checkWkflow()

    if(params.wkflow.toUpperCase() == "TEST") {
        TEST()
    }

    if(params.wkflow.toUpperCase() == "QC") {
        checkParams()
        QC()
    }

    if(params.wkflow.toUpperCase() == "TRIM") {
        TRIM()
    }

    if(params.wkflow.toUpperCase() == "ALIGN") {
        ALIGN()
    }

    if(params.wkflow.toUpperCase() == "BQSR") {
        BQSR()
    }

    if(params.wkflow.toUpperCase() == "MERGEALIGN") {
        MERGEALIGN()
    }

    if(params.wkflow.toUpperCase() == "CALL" || params.wkflow.toUpperCase() == "SCALL" || params.wkflow.toUpperCase() == "JCALL") {
        CALL()
    }

    if(params.wkflow.toUpperCase() == "FILTER") {
        FILTER()
    }

    if(params.wkflow.toUpperCase() == "ANNOTATE") {
        ANNOTATE()
    }
}


workflow.onComplete {
    msg = jobCompletionMessage()
    if(params.email == 'NULL') {
        println msg
    } 
    else {
        println msg
        sendMail(
            //from: '[ei-ngs team email]',
            to: params.email,
            subject: "[EI-NGS ${params.wkflow.toUpperCase()} WORKFLOW EXECUTION STATUS]",
            body: msg
        )
    }
}

workflow.onError{
  println "workflow execution stopped with the following message: ${workflow.errorMessage}"
}
