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
    filterVCF;
    getVcfStats;
    //plotVcfStats;
} from "${projectDir}/modules/vcfQualityMetricsAndFiltering.nf"

workflow {
    println "\nWorkflow starts here\n"

    vcf = getVcf()
    vcf_index = getVcfIndex(vcf)

    if(params.jointCaller.toUpperCase() == "GATK") {
        recalTable_tranches_snp = vqsrSnp(vcf, vcf_index)
        recalTable_tranches_indel = vqsrIndel(vcf, vcf_index)
        recalibrated_vcf_snp = applyVqsrSnp(vcf, vcf_index, recalTable_tranches_snp)
        recalibrated_vcf_indel = applyVqsrIndel(vcf, vcf_index, recalTable_tranches_indel)
        recalibrated_vcf_snp.combine(recalibrated_vcf_indel).flatten().set { filter_input }
        filtered = filterVCF(filter_input).collect()
        mergedVcf = mergeVCFs(filtered).view()
        //merged_vcf_index = indexFilteredVcf(mergedVcf)
    } else { 
        mergedVcf = vcf
        merged_vcf_index = vcf_index
    }

    multiallelicSplitted_vcf = splitMultiallelicSnvs(mergedVcf, merged_vcf_index)
    //multiallelicSplitted_vcf_index = indexMultiallelicSplittedVcf(multiallelicSplitted_vcf)
    leftnormalized_vcf = leftnormalizeSnvs(multiallelicSplitted_vcf)
    vcfstats = getVcfStats(leftnormalized_vcf)
    //plotVcfStats( vcfstats )
}

workflow.onComplete { println "\nDone filtering VCF file!\n" }
