#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFileSet;
    getCramFileSet;
    getAlignmentFileSet;
    getAlignmentGenomicIntervals;
    haplotypeCaller;
    haplotypeCallerWithIntervals;
    concatGvcfs;
    haplotypeCallerSpark;
    //indexVcf;
    combineGvcfs;
    //combineGvcfsPerInterval;
    getPedFile;
    genotypeGvcfs;
    deepVariantCaller;
    glnexusJointCaller;
    dysguCallSvs;
    indexVcf;
    dysguMergeVcfs;
    getBamChuncks;
    mantaCallSvs;
    mergeMantaCandidateSvCalls;
    mergeMantaCandidateSmallIndelCalls;
    mergeMantaDiploidSvCalls;
    convertBcfToVcf;
    createGenomicsDb;
    createGenomicsDbPerInterval;
    callVariantsFromGenomicsDB;
} from "${projectDir}/modules/variantCallingPipeline.nf"

include {
    indexAlignment;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nVariant calling begins here\n"

/*
    if(params.input_ftype.toUpperCase() == "CRAM") {
        bamFileSet = getCramFileSet()
    }
    else {
        bamFileSet = getBamFileSet()
    }
*/

    bamFileSet = getAlignmentFileSet()

    // STRUCTURAL VARIANT SINGLE SAMPLE CALLERS - DYSGU | MANTA
    if(params.single_caller == 'dysgu') {
        dysguCallSvs(bamFileSet)
            .set { vcf }
        indexVcf(vcf)
            .collect()
            .set { vcfs }
        dysguMergeVcfs(vcfs)
    }
    else if(params.single_caller == 'manta') {
        bamFileSet
            .map { bamName, bamFile, bamIndex -> tuple(bamFile, bamIndex) }
            .collect()
            .set { manta_input }
        getBamChuncks(manta_input)
            .flatten().view()
            .set { bamlist }
        mantaCallSvs(bamlist)
            .collect().view()
            .set { manta_calls }
        mergeMantaDiploidSvCalls(manta_calls).view()
        mergeMantaCandidateSmallIndelCalls(manta_calls).view()
        mergeMantaCandidateSvCalls(manta_calls).view()
    }
    else { // SMALL VARIANTS SINGLE SAMPLE CALLERS - DEEPVARIANT | GATK
        if(params.single_caller == 'deepvariant') {
            gvcf = deepVariantCaller(bamFileSet)
        }
        else { // SINGLE CALLER DEFAULTS TO GATK
            gvcf = haplotypeCaller(bamFileSet)

            //bamFileSet.first().view().set { getIntervalInput }
            //genomicIntervals = getAlignmentGenomicIntervals(getIntervalInput).flatten()
            //hapltotype_caller_input = bamFileSet.combine(genomicIntervals)
            //gvcf = haplotypeCallerWithIntervals(hapltotype_caller_input).collect().view()
        }
    
        gvcfList = gvcf.collect().view()

        // COMBINE GVCFS PER INTERVAL/CHROMOSOME FOR COMPUTATIONAL EFFICIENCY
        channel.from(1..22,'X','Y','MT')
               .collect()
               .flatten()
               .map { chr -> "chr${chr}" }
               .combine(gvcfList.toList())
               .set { per_chrom_gvcf_input }



        // // JOINT CALLERS FOR SMALL VARIANTS - GLNEXUS | GATK 
        // if(params.joint_caller == 'glnexus')  {
        //     bcf = glnexusJointCaller(gvcfList).view()
        //     vcf = convertBcfToVcf(bcf).view()
        // } 
        // else {

        //     /************************************************** 
        //     *          JOINT CALLER DEFAULTS TO GATK          *
        //     *  Use 'CombineGVCFs' with less than 100 samples  *
        //     *        Otherwise, use 'GenomicsDBImport'        *
        //     **************************************************/

        //     sampleCount = gvcfList.size()
        //     if( sampleCount < 100 ) {
        //         combinedGvcf = combineGvcfsPerInterval(per_chrom_gvcf_input)
        //         ped = getPedFile(combinedGvcf)
        //         combinedGvcf.combine(ped).set { join_call_input }
        //         vcf = genotypeGvcfs(join_call_input)
        //     } 
        //     else {
        //         combinedGvcf = createGenomicsDbPerInterval(per_chrom_gvcf_input).view()
        //         vcf = callVariantsFromGenomicsDB(combinedGvcf).view()

        //     //combinedGvcf = combineGvcfsPerInterval(gvcfList)
        //     //getGenomicIntervals(gvcfList)
        //     //combinedGvcf = createGenomicsDb(gvcfList).view()

        //     }
        // }
    }
}

workflow.onComplete { println "\nDone! Check results in ${params.output_dir}\n" }
