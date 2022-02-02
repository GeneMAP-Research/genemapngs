def getBamFiles() {
    return channel.fromFilePairs( params.bamDir + "*.{bam,bai}", size: 2 )
                  .ifEmpty { error "\nERROR: Could not locate a file!\n" }
                  .map { bamName, bamFileset -> tuple(bamName, bamFileset) }
}

process callVariants() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            val(bamName), \
            path(bamFile)
    output:
        publishDir path: "${params.outputDir}/gvcfs/"
        path "${bamName}.gvcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \
            HaplotypeCallerSpark \
            -I ${bamFile[1]} \
            -R ${params.fastaRef} \
            --dbsnp ${params.dbsnp} \
            --native-pair-hmm-threads ${task.cpus} \
            --lenient true \
            -ERC GVCF \
            -OBI 'false' \
            -O "${bamName}.gvcf.gz" \
            -- \
            --spark-master local[${task.cpus}]
        """
}

process combineGvcfs() {
    tag "Combining GVCFs to ${params.outPrefix}.gvcf.gz"
    label 'gatk' 
    label 'variantCaller'
    input:
        path gvcfList
    output:
        publishDir path: "${params.outputDir}/vcf/"
        path "${params.outPrefix}.gvcf.{gz,gz.tbi}"
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

        gatk \
            CombineGVCFs \
            -R ${params.fastaRef} \
            --arguments_file gvcf.list \
            --dbsnp ${params.dbsnp} \
            -O "${params.outPrefix}.gvcf.gz"
        """
}

process getPedFile() {
    tag "Writing PED file"
    label 'bcftools'
    label 'smallMemory'
    echo true
    input:
        tuple \
            path(gvcf), \
            path(gvcfIndex)
    output:
        publishDir path: "${params.outputDir}/ped/"
        path "${params.outPrefix}.tmp.ped"
    script:
        """
        if [ ${params.ped} == NULL ]; then
           bcftools query -l ${gvcf} | \
               awk '{print \$1,\$1,"0","0","-9","-9"}' \
               > ${params.outPrefix}.ped
           ped="${params.outPrefix}.ped"
        else
           ped="${params.ped}"
        fi
        cp \${ped} ${params.outPrefix}.tmp.ped
        """
}

process joinCallVariants() {
    tag "Writing genotypes to ${params.outPrefix}.vcf.gz"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            path(gvcf), \
            path(gvcfIndex), \
            path(ped)
    output:
        publishDir path: "${params.outputDir}/vcf/", mode: 'copy'
        path "${params.outPrefix}.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            GenotypeGVCFs \
            -R ${params.fastaRef} \
            --dbsnp ${params.dbsnp} \
            -ped ${ped} \
            -V ${gvcf} \
            -O "${params.outPrefix}.vcf.gz"
        """
}

/*
process indexVcf() {
    tag "processing ${vcfName}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            val(vcfName), \
            path(vcf)
    output:
        publishDir path: "${params.outputDir}/gvcfs/"
        tuple \
            val(vcfName), \ 
            path("${vcf}"), \
            path("${vcf}.tbi")
    script: 
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${vcf}
        """
}
*/
