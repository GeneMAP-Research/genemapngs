#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getVcf;
    getVcfIndex;
    getChromosomes;
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
    transformVcfWithPlink;
    getVcfIndex as getTransformedVcfIndex;
    getVcfStats;
    //plotVcfStats;
    getVcfIntersect;    
} from "${projectDir}/modules/vcfQualityAndFilter.nf"

workflow {
    println "\nWorkflow starts here\n"

    vcf = getVcf()
    vcf_index = getVcfIndex(vcf)
    //chrom = getChromosomes()
    //kgp = getThousandGenomesReference()
    recalTable_tranches_snp = vqsrSnp(vcf, vcf_index)
    recalTable_tranches_indel = vqsrIndel(vcf, vcf_index)
    recalibrated_vcf_snp = applyVqsrSnp(vcf, vcf_index, recalTable_tranches_snp)
    recalibrated_vcf_indel = applyVqsrIndel(vcf, vcf_index, recalTable_tranches_indel)

    recalibrated_vcf_snp.combine(recalibrated_vcf_indel).flatten().set { filter_input }
    filtered = filterVCF(filter_input).collect()
    mergedVcf = mergeVCFs(filtered).view()
    merged_vcf_index = indexFilteredVcf(mergedVcf)
    multiallelicSplitted_vcf = splitMultiallelicSnvs(mergedVcf, merged_vcf_index)
    multiallelicSplitted_vcf_index = indexMultiallelicSplittedVcf(multiallelicSplitted_vcf)
    leftnormalized_vcf = leftnormalizeSnvs(multiallelicSplitted_vcf, multiallelicSplitted_vcf_index)

/*
    transformedVcf = transformVcfWithPlink(leftnormalized_vcf)
    transformedVcf_index = getTransformedVcfIndex(transformedVcf)
    vcfstats = getVcfStats(leftnormalized_vcf)
    //plotVcfStats( vcfstats )

    transformedVcf
        .combine(transformedVcf_index)
        .set { transformedVcf_input }

    chrom
        .combine(transformedVcf_input)
        .map { chr, vcf, vcf_index -> 
            tuple( "${chr}", vcf, vcf_index ) 
        }
        .join(kgp)
        .set { input }

    vcf_intersect = getVcfIntersect( input )
    vcf_intersect
        .collect()
        .set { intersects }

    intersects
        .collect()
        .map { files -> tuple( files) }
        .view()
    //intersect_path = getIntersectPath().collect().view()
    //concatenateIntersectVcf( intersects )
*/
}

