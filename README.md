READS QUALITY, ALIGNMENT, VARIANT CALLING, VARIANT FILTERING, AND ANNOTATION WORKFLOW
---
This repository contains workflows/pipelines for the following tasks:
- Reads quality assessment: ```getFastqBamQualityReports.nf```
- Reads alignment to reference: ```alignReadsToReference.nf``` 
- Variant calling: ```callVariants.nf``` 
- Variant fltering: ```filterVariantCalls.nf```
- Variant annotation: ```annotateVariants.nf```

It takes as input FASTQ or BAM files which are specified in the ```nextflow.config``` file.

Basic Usage
---

* Get reads quality reports
```
$ nextflow run getFastqBamQualityReports.nf -w /path/to/work-directory/ -profile local
```

* Align reads to reference
```
$ nextflow run alignReadsToReference.nf -w /path/to/work-directory/ -profile chpc
```

* Call variants from aligned BAM files
```
$ nextflow run callVariants.nf -w /path/to/work-directory/ -profile ucthpc
```

IMPORTANT NOTE
----
You must build the indices of all references in the same directory as the references for this workflow to run successfully

- Build BWA index
```
$ bwa index ref
```

- Build GATK index
```
$ gatk CreateSequenceDictionary -R ref
```

- Build .fai index with SAMTOOLS FAIDX
```
$ samtools faidx ref
```
