#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getVcf;
    getVcfIndex;
    getChromosomes;
    getThousandGenomesReference;
    vqsr;
    applyVQSR;
    getVcfIndex as indexFilteredVcf;
    splitMultiallelicSnvs;
    getVcfIndex as indexMultiallelicSplittedVcf;
    leftnormalizeSnvs;
    filtereVCF;
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
    chrom = getChromosomes()
    kgp = getThousandGenomesReference()
    recalTable_tranches = vqsr(vcf, vcf_index)
    recalibrated_vcf = applyVQSR(vcf, vcf_index, recalTable_tranches)
    filtered_vcf = filtereVCF(recalibrated_vcf)
    filtered_vcf_index = indexFilteredVcf(filtered_vcf)
    multiallelicSplitted_vcf = splitMultiallelicSnvs(filtered_vcf, filtered_vcf_index)
    multiallelicSplitted_vcf_index = indexMultiallelicSplittedVcf(multiallelicSplitted_vcf)
    leftnormalized_vcf = leftnormalizeSnvs(multiallelicSplitted_vcf, multiallelicSplitted_vcf_index)
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
}

