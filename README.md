READS QUALITY, ALIGNMENT, VARIANT CALLING, VARIANT FILTERING, AND ANNOTATION WORKFLOW
---
This repository contains workflows/pipelines for the following tasks:
- Reads quality assessment: ```getQualityReports.nf```
- Reads alignment to reference: ```alignReadsToReference.nf``` 
- Variant calling: ```callVariants.nf``` 
- Variant fltering: ```filterVariantCalls.nf```
- Variant annotation: ```annotateVariants.nf```

Documentation
---
Go [here]("https://genemap-research.github.io/docs/workflows/ngs/") for usage.


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
