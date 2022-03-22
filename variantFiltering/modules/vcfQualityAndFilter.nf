def getVcf() {
    return channel.fromPath( params.inputDir + params.vcf )
}

def getChromosomes() {
    return channel.of(1..22, 'X')
}

def getIntersectPath() {
    return channel.fromPath( params.outputDir + "chr*", type: 'dir' )
}

def getThousandGenomesReference() {
    return channel.fromFilePairs( params.referenceDir + "chr*1kg.phase3.v5a.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error "\nAn error occurred! Please check that the reference file and its index '.tbi' exist...\n" }
	   	  .map { chr, ref_file -> 
			 tuple( chr.replaceFirst(/chr/,""), ref_file.first(), ref_file.last())
	    	  }
}

process getVcfIndex() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path input_vcf
    output:
        path "*.tbi"
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf}
        """
}

process vqsrSnp() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'vqsr'
    input:
        path(input_vcf)
        path(vcf_index)
    output:
        publishDir path: "${params.outputDir}"
        tuple \
              path("${params.outputPrefix}.snp.recal"), \
              path("${params.outputPrefix}.snp.recal.idx"), \
              path("${params.outputPrefix}.snp.tranches")
    script:
        """
        gatk \
           --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
           VariantRecalibrator \
           -R ${params.fastaRef} \
           -tranche 100.0 -tranche 99.95 -tranche 99.9 \
           -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
           -tranche 95.0 -tranche 94.0 \
           -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           -V ${input_vcf} \
           --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
           --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.snpRef} \
           --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omniRef} \
           --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
           -an QD \
           -an MQ \
           -an MQRankSum \
           -an ReadPosRankSum \
           -an FS \
           -an SOR \
           -mode SNP \
           -O ${params.outputPrefix}.snp.recal \
           --tranches-file ${params.outputPrefix}.snp.tranches \
           --rscript-file ${params.outputPrefix}.snp.plots.R
        """
}

process vqsrIndel() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'vqsr'
    input:
        path(input_vcf)
        path(vcf_index)
    output:
        publishDir path: "${params.outputDir}"
        tuple \
              path("${params.outputPrefix}.indel.recal"), \
              path("${params.outputPrefix}.indel.recal.idx"), \
              path("${params.outputPrefix}.indel.tranches")
    script:
        """
        gatk \
           --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
           VariantRecalibrator \
           -R ${params.fastaRef} \
           -V ${input_vcf} \
           -tranche 100.0 -tranche 99.95 -tranche 99.9 \
           -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
           -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 \
           -tranche 92.0 -tranche 91.0 -tranche 90.0 \
           --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.hapmap} \
           --resource:1000G,known=false,training=true,truth=false,prior=12.0 ${params.millsIndel} \
           --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.omniRef} \
           --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
           -an QD \
           -an DP \
           -an MQRankSum \
           -an ReadPosRankSum \
           -an FS \
           -an SOR \
           -mode INDEL \
           --max-gaussians 4 \
           -O ${params.outputPrefix}.indel.recal \
           --tranches-file ${params.outputPrefix}.indel.tranches \
           --rscript-file ${params.outputPrefix}.indel.plots.R
        """
}

process applyVqsrSnp() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'applyVqsr'
    input:
        path(input_vcf)
        path(vcf_index)
        tuple \
            path(recal_table), \
            path(recal_index), \
            path(tranches)   
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path "${params.outputPrefix}.snp.vqsr.vcf.gz"
    script:
        """
        gatk \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            ApplyVQSR \
            -V ${input_vcf} \
            --recal-file ${recal_table} \
            --mode SNP \
            --truth-sensitivity-filter-level 99.9 \
            --create-output-variant-index true \
            --tranches-file ${tranches} \
            -R ${params.fastaRef} \
            -O ${params.outputPrefix}.snp.vqsr.vcf.gz
        """
}

process applyVqsrIndel() {
    tag "VCF supplied: ${input_vcf}"
    label 'gatk'
    label 'applyVqsr'
    input:
        path(input_vcf)
        path(vcf_index)
        tuple \
            path(recal_table), \
            path(recal_index), \
            path(tranches)   
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path "${params.outputPrefix}.indel.vqsr.vcf.gz"
    script:
        """
        gatk \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            ApplyVQSR \
            -V ${input_vcf} \
            --recal-file ${recal_table} \
            --mode INDEL \
            --truth-sensitivity-filter-level 99.9 \
            --create-output-variant-index true \
            --tranches-file ${tranches} \
            -R ${params.fastaRef} \
            -O ${params.outputPrefix}.indel.vqsr.vcf.gz
        """
}

process mergeVCFs() {
    tag "merging VCFs"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            path(snp_vcf), \
            path(indel_vcf)
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path"${params.outputPrefix}.snp.indel.vqsr.vcf.gz"
    script:
        """
        for i in $snp_vcf $indel_vcf; do bcftools index -ft --threads ${task.cpus} \$i; done

        bcftools \
            concat \
            -a -d all \
            --threads ${task.cpus} \
            -Oz \
            $snp_vcf \
            $indel_vcf \
            -o ${params.outputPrefix}.snp.indel.vqsr.vcf.gz
        """
}

process filterVCF() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path input_vcf
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path "${input_vcf.baseName}.filtered.vcf.gz"
    script:
        """
        bcftools \
            view \
            -i \'FILTER=="PASS"\' \
            --threads ${task.cpus} \
            ${input_vcf} | \
            bcftools \
                view \
                -i \'INFO/DP>=${params.minDP}\' \
                --threads ${task.cpus} | \
                bcftools \
                    view \
                    -i \'GQ>=${params.minGQ}\' \
                    --threads ${task.cpus} \
                    -Oz \
                    -o "${input_vcf.baseName}.filtered.vcf.gz"
        """
}

process splitMultiallelicSnvs() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path(input_vcf)
        path(vcf_index)
    output:
        publishDir path: "${params.outputDir}"
        path "${params.outputPrefix}-tmp.vcf.gz"
    script:
        """
        bcftools \
            norm \
            -m-both \
            --threads ${task.cpus} \
            -Oz \
            -o "${params.outputPrefix}-tmp.vcf.gz" \
            ${input_vcf}
        """
}

process leftnormalizeSnvs() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path(input_vcf)
        path(vcf_index)
    output:
        publishDir path: "${params.outputDir}/annotate/", mode: 'copy'
        path "${params.outputPrefix}-filtered-leftnorm.vcf.gz"
    script:
        """
        bcftools \
            norm \
            -f ${params.fastaRef} \
            --threads ${task.cpus} \
            -Oz \
            -o "${params.outputPrefix}-filtered-leftnorm.vcf.gz"  \
            ${input_vcf}
        """
}

process getVcfStats() {
    tag "VCF supplied: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path input_vcf
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path "${input_vcf.baseName}.vcfstats.txt"
    script:
        """
        bcftools \
            stats \
            -F "${params.fastaRef}" \
            -s - \
            "${input_vcf}" > \
            "${input_vcf.baseName}.vcfstats.txt"
        """
}

process plotVcfStats() {
    tag "Processing ${vcfstat}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path vcfstat
    output:
        publishDir path: "${params.outputDir}${vcfstat.baseName}", mode: 'copy'
        path "./*"
    script:
        """
        mkdir -p "${params.outputDir}${vcfstat.baseName}"
        plot-vcfstats \
           -p "${params.outputDir}${vcfstat.baseName}" \
           "${vcfstat}"
        """
}

process transformVcfWithPlink() {
    tag "PLINK2 EXPORT: ${input_vcf.baseName}-plink.vcf.gz"
    label 'plink2'
    label 'mediumMemory'
    input:
        path input_vcf
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path "${input_vcf.baseName}-plink.vcf.gz"
    script:
        """
        plink2 \
            --vcf ${input_vcf} \
            --aec \
            --export vcf-4.2 bgz id-paste='iid' \
            --threads ${task.cpus} \
            --out "${input_vcf.baseName}-plink"
        """
}

process getVcfIntersect() {
    tag "processing chr${chrom}"
    label 'bcftools'
    label 'mediumMemory'
    cache 'lenient'
    input:
        tuple val(chrom), path(input_vcf), path(vcf_index), path(ref_vcf), path(ref_index)
    output:
        publishDir path: "${params.outputDir}chr${chrom}/", mode: 'copy'
        path "000{2,3}.vcf.gz"
    script:
        """
        mkdir -p "${params.outputDir}/chr${chrom}"
        
        bcftools \
            isec \
            -r ${chrom} \
            --threads ${task.cpus} \
            -p . \
            -Oz \
            -o "chr${chrom}" \
            ${input_vcf} \
            ${ref_vcf}
        """
}

process concatenateIntersectVcf() {
    tag "CONCAT"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path intersects
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path ""
    script:
        """
        bcftools \
            concat \
            --threads ${task.cpus} \
            -Oz \
            -o "${params.outputPrefix}-kgp-intersect.vcf.ga" \
            ${intersects}
        """

}
