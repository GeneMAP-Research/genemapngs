def getVCF() {
    return channel.fromPath( params.inputDir + "*.vcf*" )
                  .ifEmpty { println "\nERROR: No VCF file found!\n" }
}

process annotateVcfFiles() {
    tag "processing ${vcfFile}"
    label 'annovar'
    label 'annovarMem'

    input:
        path(vcfFile)
    output:
        publishDir path: "${params.outputDir}/annovarOutput/", mode: 'copy'
        path "${vcfFile.baseName}*multianno.{vcf,txt}"
    script:
        """
        table_annovar.pl \
            ${params.inputDir}/${vcfFile} \
            ${params.annovarDB} \
            -buildver hg19 \
            -out ${vcfFile.baseName} \
            -remove \
            -protocol refGene,knownGene,gwasCatalog,genomicSuperDups,cytoBand,exac03,avsnp150,abraom,dbnsfp33a,regsnpintron,esp6500siv2_all,SAS.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,AMR.sites.2015_08,ALL.sites.2015_08,AFR.sites.2015_08,nci60,clinvar_20210501,refGeneWithVer,revel,mitimpact24,intervar_20180118,hrcr1,gme,gnomad211_exome,gene4denovo201907,dbscsnv11,dbnsfp31a_interpro,cosmic70 \
            -operation g,g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
            -nastring '.' \
            -vcfinput \
            -polish 
        """
}
