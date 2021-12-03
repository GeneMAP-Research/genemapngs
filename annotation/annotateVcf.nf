#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getVCF;
    annotateVcfFiles;
} from "${projectDir}/modules/annovarAnnotation.nf"


workflow {
    println "\nANNOVAR annotation starts here\n"
    vcf = getVCF().view()
    annotateVcfFiles(vcf)
}

