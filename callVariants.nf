#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getBamFileSet;
    getCramFileSet;
    haplotypeCaller;
    haplotypeCallerSpark;
    //indexVcf;
    combineGvcfs;
    combineGvcfsPerChromosome;
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
    createGenomicsDbPerChromosome;
    callVariantsFromGenomicsDB;
    getGenomicIntervals;
} from "${projectDir}/modules/variantCallingPipeline.nf"

include {
    indexBam;
} from "${projectDir}/modules/alignmentPipeline.nf"

workflow {
    println "\nVariant calling begins here\n"

    if(params.inputFileType.toUpperCase() == "CRAM") {
        bamFileSet = getCramFileSet()
    }
    else {
        bamFileSet = getBamFileSet()
    }

    // STRUCTURAL VARIANT SINGLE SAMPLE CALLERS - DYSGU | MANTA
    if(params.singleCaller == 'dysgu') {
        dysguCallSvs(bamFileSet)
            .set { vcf }
        indexVcf(vcf)
            .collect().view()
            .set { vcfs }
        dysguMergeVcfs(vcfs).view()
    }
    else if(params.singleCaller == 'manta') {
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
        if(params.singleCaller == 'deepvariant') {
            gvcf = deepVariantCaller(bamFileSet)
        }
        else { // SINGLE CALLER DEFAULTS TO GATK
            if(params.sparkMode == true) {
               gvcf = haplotypeCallerSpark(bamFileSet)
            } else {
               gvcf = haplotypeCaller(bamFileSet)
            }
        }
    
        gvcfList = gvcf.collect().view()

        // COMBINE GVCFS PER INTERVAL/CHROMOSOME FOR COMPUTATIONAL EFFICIENCY
        channel.from(1..22,'X','Y','MT')
               .collect()
               .flatten()
               .map { chr -> "chr${chr}" }
               .combine(gvcfList.toList())
               .set { per_chrom_gvcf_input }


        // JOINT CALLERS FOR SMALL VARIANTS - GLNEXUS | GATK 
        if(params.jointCaller == 'glnexus')  {
            bcf = glnexusJointCaller(gvcfList).view()
            vcf = convertBcfToVcf(bcf).view()
        } 
        else {

            /************************************************** 
            *          JOINT CALLER DEFAULTS TO GATK          *
            *  Use 'CombineGVCFs' with less than 100 samples  *
            *        Otherwise, use 'GenomicsDBImport'        *
            **************************************************/

            sampleCount = gvcfList.size()
            /*if( sampleCount < 100 ) {
                combinedGvcf = combineGvcfsPerChromosome(per_chrom_gvcf_input)
                ped = getPedFile(combinedGvcf)
                combinedGvcf.combine(ped).set { join_call_input }
                vcf = genotypeGvcfs(join_call_input)
            } 
            else {
                combinedGvcf = createGenomicsDbPerChromosome(per_chrom_gvcf_input).view()
                vcf = callVariantsFromGenomicsDB(combinedGvcf).view()

            //combinedGvcf = combineGvcfsPerChromosome(gvcfList)
            //getGenomicIntervals(gvcfList)
            //combinedGvcf = createGenomicsDb(gvcfList).view()

            }*/
        }
    }
}

workflow.onComplete { println "\nDone! Check results in ${params.outputDir}\n" }
