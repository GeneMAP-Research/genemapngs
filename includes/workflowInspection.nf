def jobCompletionMessage() {
   return msg = """
          Pipeline execution summary
          ---------------------------
          Completed at: ${workflow.complete}
          Duration    : ${workflow.duration}
          exit status : ${workflow.success ? 'OK' : 'failed' }
          Success     : ${workflow.success}
          workDir     : ${workflow.workDir}
          outputDir   : ${params.output_dir}
          """
          .stripIndent()
}

def checkWkflow() {
   if(params.wkflow == "" || params.wkflow.toUpperCase() == "NULL") {
      wkflowOptions = """
      options:

            test: Run test to see if workfow installed correctly.
              qc: Check FASTQ or Alignment (BAM/CRAM) quality.
            trim: Trim adapters and poor quality bases from reads (input is FASTQ or BAM/CRAM).
           align: Align/map reads to reference and post-alignment processing (input is FASTQ or BAM/CRAM).
      mergealign: Megre Alignment (BAM/CRAM) files.
            call: Perform variant calling (both single and joint sample) in one run.
           scall: Perform only sinlge sample variant calling to generate gVCF files.
           jcall: Perform only joint (multi-sample) variant calling with pre-existing gVCF files.
          filter: Filter variant calls in VCF/BCF files.
        annotate: Annotate variants with ANNOVAR
         raw2vcf: Run entire workflow from 'align' to 'call' (i.e., 'scall' + 'jcall')
    raw2annotate: Run entire workflow from 'align' to 'annotate'

        [NOTE] For 'trim' and 'align', BAM/CRAM input is first converted to FASTQ

      """
      error: println "\nPlease specify a workflow to run using '--wkflow <option>' \n${wkflowOptions}"
   }
}

def checkParams() {
   if(params.wkflow.toUpperCase() == "QC" && params.input_dir == "NULL") {
      usage = """
      Usage: main.nf --wkflow qc <options>
         options:
         --------
         --wgs                : Specify this flag if your data is whole-genome sequence (it runs whole exome - wes - by default)
                                 This is important for resource allocation.
         --ftype              : Input file type; FASTQ, BAM, CRAM [default: FASTQ].
         --input_dir          : (required) Path to FASTQ/BAM/CRAM files.
         --output_dir         : (optional) Results will be saved to parent of input directory ['input_dir/../'].
         --threads            : number of computer cpus to use [default: 4].
         --njobs              : (optional) number of jobs to submit at once [default: 4]
         --help               : print this help message.
      """
      error: println "\nPlease provide all required arguments!\n${usage}" 
   }
}

def checkFaiIndex() {
   if(params.fastaRef == 'NULL') {
      error: "Please select a reference to use by adding in the profile builder. E.g., -profile singularity,slurm,hg19"
   } else {
      channel.fromPath(params.fastaRef + '.fai', checkIfExists: true).ifEmpty(error: ".fai reference index could not be found! Please create one in the reference directory with SAMTOOLS.")
   }
}

def checkGatkSeqDict() {
   channel.fromPath(params.fastaRef.replaceFirst(/fa/,"dict"), checkIfExists: true).ifEmpty(error: "Reference dictionary '.dict' required by GATK could not be found! Please create one in the reference directory.")
}

def checkBwaIndex() {
   channel.fromPath(params.fastaRef + '.sa', checkIfExists: true).ifEmpty(error: "BWA index files could not be found! Please index the reference with BWA.")
   channel.fromPath(params.fastaRef + '.pac', checkIfExists: true).ifEmpty(error: "Please index the reference with BWA")
   channel.fromPath(params.fastaRef + '.ann', checkIfExists: true).ifEmpty(error: "Please index the reference with BWA")
   channel.fromPath(params.fastaRef + '.amb', checkIfExists: true).ifEmpty(error: "Please index the reference with BWA")
}

def checkBwa2Index() {
   channel.fromPath(params.fastaRef + '.0123', checkIfExists: true).ifEmpty(error: "BWA-MEM2 index files could not be found! Please index the reference with BWA-MEM2")
   channel.fromPath(params.fastaRef + '.bwt.2bit.64', checkIfExists: true).ifEmpty(error: "Please index the reference with BWA-MEM2")
}