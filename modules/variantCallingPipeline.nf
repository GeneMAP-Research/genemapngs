def getBamFile() {
    return channel.fromFilePairs( params.outputDir + '/bam/' + "*.bam", size: 1, flat: true )
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
                  //.map { bamName, bamFile -> tuple(bamName, bamFile) }
}

def getBamIndex() {
    return channel.fromFilePairs( params.outputDir + '/bam/' + "*.bai", size: 1, flat: true )
                  .ifEmpty { error "\nERROR: Could not locate a file! \n" }
                  //.map { bamName, bamIndex -> tuple(bamName, bamIndex) }
}

process callVariants() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/gvcfs/"
        path "${bamName}.g.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            HaplotypeCaller \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            -ERC GVCF \
            -OBI 'false' \
            -O "${bamName}.g.vcf.gz" 
        """
}

//            --java-options "-XX:ConcGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g" \

process callVariantsSpark() {
    tag "processing ${bamName}"
    label 'gatk'
    label 'variantCaller'
    input:
        tuple \
            val(bamName), \
            path(bamFile), \
            path(bamIndex)
    output:
        publishDir path: "${params.outputDir}/gvcfs/"
        path "${bamName}.g.vcf.{gz,gz.tbi}"
    script:
        """
        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            HaplotypeCallerSpark \
            -I ${bamFile} \
            -R ${params.fastaRef} \
            --native-pair-hmm-threads ${task.cpus} \
            -ERC GVCF \
            -OBI 'false' \
            -O "${bamName}.g.vcf.gz" \
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
        path "${params.outPrefix}.g.vcf.{gz,gz.tbi}"
    script:
        """
        for file in ${gvcfList}; do
            if [ \${file##*.} != "tbi" ]; then
                echo "-V \${file}";
            fi
        done > gvcf.list

        gatk \
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            CombineGVCFs \
            -R ${params.fastaRef} \
            --arguments_file gvcf.list \
            --dbsnp ${params.dbsnp} \
            -O "${params.outPrefix}.g.vcf.gz"
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
            --java-options "-XX:ConcGCThreads=${task.cpus} -Xms${task.memory.toGiga()}g -Xmx${task.memory.toGiga()}g -XX:ParallelGCThreads=${task.cpus}" \
            GenotypeGVCFs \
            -R ${params.fastaRef} \
            --dbsnp ${params.dbsnp} \
            -ped ${ped} \
            -V ${gvcf} \
            -O "${params.outPrefix}.vcf.gz"
        """
}

