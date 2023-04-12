#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getVCF;
    annovarGRCh37;
    annovarGRCh38;
    minimalAnnovarGRCh37;
    minimalAnnovarGRCh38;
} from "${projectDir}/modules/annovarAnnotation.nf"


workflow {
    println "\nANNOVAR annotation starts here\n"

    msg = "\nERROR: You must select a Build Version! Options are hg19 and hg38\n" 

    vcf = getVCF()

    if (params.minimal == true) {
       if (params.buildVersion == "hg19") {
           minimalAnnovarGRCh37(vcf)
       }
       else if (params.buildVersion == "hg38") {
           minimalAnnovarGRCh38(vcf)
       }
       else { error: println msg }
    }
    else {
       if (params.buildVersion == "hg19") {
           annovarGRCh37(vcf)
       }
       else if (params.buildVersion == "hg38") {
           annovarGRCh38(vcf)
       }
       else { error: println msg }
    }
}

