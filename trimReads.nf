#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getInputFastqs;
} from "${projectDir}/modules/alignmentPipeline.nf"

include {
    getNpAdapter;
    getT3uAdapter;
    getT2pAdapter;
    getT3pAdapter;
    getT2sAdapter;
    getT3sAdapter;
    readsTrimmer;
} from "${projectDir}/modules/getReadsQualityReports.nf"

workflow {
    println "\nFastq trimming begins here\n"
    fastq = getInputFastqs()

    if(params.adapter.toUpperCase() == "NP") {
        adapter = getNpAdapter()
        fastq
            .combine(adapter)
            .set { trim_input }
    }
    else if(params.adapter.toUpperCase() == "T3U") {
        adapter = getT3uAdapter()
        fastq
            .combine(adapter)
            .set { trim_input }
    }
    else if(params.adapter.toUpperCase() == "T2P") {
        adapter = getT2pAdapter()
        fastq
            .combine(adapter)
            .set { trim_input }
    }
    else if(params.adapter.toUpperCase() == "T3P") {
        adapter = getT3pAdapter()
        fastq
            .combine(adapter)
            .set { trim_input }
    }
    else if(params.adapter.toUpperCase() == "T2S") {
        adapter = getT2sAdapter()
        fastq
            .combine(adapter)
            .set { trim_input }
    }
    else if(params.adapter.toUpperCase() == "T3S") {
        adapter = getT3sAdapter()
        fastq
            .combine(adapter)
            .set { trim_input }
    }
    else {
        error: println """
    Please select and adapter!
    --------------------------
    Options:
        NP   --> NexteraPE-PE.fa
        T3U  --> TruSeq3-PE-2.fa [Illumina universal]
        T2P  --> TruSeq2-PE.fa
        T2S  --> TruSeq2-SE.fa
        T3P  --> TruSeq3-PE.fa
        T3S  --> TruSeq3-SE.fa

    Edit the config file or use '--adapter NP' on the command line
    """
    }

    trimmed_reads = readsTrimmer(trim_input).view()

}

