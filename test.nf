#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
  println "\n        ALIGNMENT WORKFLOW: TEST"

  println "pe=${params.pe}                                   #// optional: true, false [dfault: false]   (Whether reads are paired-end or single end)"
  println "aligner=${params.aligner}                         #// optional: BWA, DRAGMAP [default: BWA]"
  println "ftype=${params.input_ftype}                       #// required: FASTQ, BAM, CRAM     (input file type)"
  println "input_dir=${params.input_dir}                     #// required"
  println "output_dir=${params.output_dir}                   #// required"
  println "output_prefix=${params.output_prefix}             #// required"
  println "single_caller=${params.single_caller}             #// options: gatk, deepvariant, dysgu, manta"
  println "exome=${params.exome}                             #// for manta structural variant calling, specify whether WES or WGS"
  println "joint_caller=${params.joint_caller}               #// options: gatk, glnexus"
  println "gvcf_dir=${params.gvcf_dir}                       #// if GVCF files already exist and only joint calling is required"
  println "spark=${params.spark}                             #// Select GATK spark mode to run GATK commands with multi-threading (true) or not (false)"
  println "threads=${params.threads}"
  println "njobs=${params.njobs}"
  println ""

  plink()
	
}



workflow.onComplete { 
  println "Workflow completed at: ${workflow.complete}"
  println "     Execution status: ${ workflow.success ? 'OK' : 'failed'}"
}

workflow.onError{
  println "workflow execution stopped with the following message: ${workflow.errorMessage}"
}

process plink() {

  // directives
  tag "processing ..."
  label 'plink'
  label 'test'
  //debug true
  echo true
  
  script:
    """
    echo "PLINK2 is used for the test as it is light-weight and easily pulled from docker hub"

    plink2 \
      --help \
      --file
    """
}
