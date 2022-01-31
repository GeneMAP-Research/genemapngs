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
