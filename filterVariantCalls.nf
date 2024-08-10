#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getVcf;
    getVcfIndex;
    splitVcfs;
    getVcfIndex as indexSnpVcf;
    getVcfIndex as indexIndelVcf;
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
        splitVcfs(vcf_index)
            .multiMap { snp, indel -> 
                snp: snp
                indel: indel
            }
            .set { vcf_split }
        indexSnpVcf(vcf_split.snp)
            .set { snp_indexed }
        indexIndelVcf(vcf_split.indel)
            .set { indel_indexed }
        vqsrSnp( snp_indexed )
            .set { recalTable_tranches_snp }
        vqsrIndel( indel_indexed )
            .set { recalTable_tranches_indel }
        applyVqsrSnp(recalTable_tranches_snp).view()
            .set { recalibrated_vcf_snp }
        applyVqsrIndel(recalTable_tranches_indel).view()
            .set { recalibrated_vcf_indel }
        recalibrated_vcf_snp
            .combine(recalibrated_vcf_indel, by:0)
            .view()
            .set { recalibrated }
        mergedVcf = mergeVCFs( recalibrated ).view()
        //filtered = filterGatkCalls(filter_input).collect().view()
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
