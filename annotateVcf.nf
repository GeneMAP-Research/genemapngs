#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getVCF;
    annovarGRCh37;
    annovarGRCh38;
} from "${projectDir}/modules/annovarAnnotation.nf"


workflow {
    println "\nANNOVAR annotation starts here\n"

    vcf = getVCF()

    if (params.buildVersion == "hg19") {
        annovarGRCh37(vcf)
    }
    else if (params.buildVersion == "hg38") {
        annovarGRCh38(vcf)
    }
    else { error: "\nERROR: You must select a Build Version! Options are hg19 and hg38\n" }
}

