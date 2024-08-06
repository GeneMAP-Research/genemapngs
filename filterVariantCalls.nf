#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getVcf;
    getVcfIndex;
    getThousandGenomesReference;
    vqsrSnp;
    vqsrIndel;
    applyVqsrSnp;
    applyVqsrIndel;
    mergeVCFs;
    getVcfIndex as indexFilteredVcf;
    splitMultiallelicSnvs;
    getVcfIndex as indexMultiallelicSplittedVcf;
    leftnormalizeSnvs;
    filterGatkCalls;
    filterGlnexusCalls;
    getVcfStats;
    getVcfStats as getVcfStats_b;
    getCleanVcf;
    //plotVcfStats;
} from "${projectDir}/modules/vcfQualityMetricsAndFiltering.nf"

workflow {
    println "\nWorkflow starts here\n"

    vcf = getVcf()
    vcf_index = getVcfIndex(vcf)

    if(params.joint_caller.toUpperCase() == "GATK") {
        recalTable_tranches_snp = vqsrSnp(vcf_index)
        recalTable_tranches_indel = vqsrIndel(vcf_index)
        recalibrated_vcf_snp = applyVqsrSnp(recalTable_tranches_snp)
        recalibrated_vcf_indel = applyVqsrIndel(recalTable_tranches_indel)
        recalibrated_vcf_snp.combine(recalibrated_vcf_indel).flatten().set { filter_input }
        filtered = filterGatkCalls(filter_input).collect()
        mergedVcf = mergeVCFs(filtered).view()
        //merged_vcf_index = indexFilteredVcf(mergedVcf)
    } else { 
        vcf.combine(vcf_index).set { vcfstats_input }     
        getVcfStats(vcfstats_input).view()
        mergedVcf = filterGlnexusCalls(vcf, vcf_index).view()
    }

    if(params.left_norm == true) {
        multiallelicSplitted = splitMultiallelicSnvs(mergedVcf)
        multiallelicSplitted
            .map { vcf, index -> vcf }
            .set { multiallelicSplitted_vcf }
        multiallelicSplitted_vcf_index = indexMultiallelicSplittedVcf(multiallelicSplitted_vcf)
        leftnormalized_vcf = leftnormalizeSnvs(multiallelicSplitted)
        vcfstats = getVcfStats_b(leftnormalized_vcf)
        //plotVcfStats( vcfstats )
        getCleanVcf(leftnormalized_vcf).view()
    }
}

workflow.onComplete { println "\nDone filtering VCF file!\n" }
