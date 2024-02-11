#!/usr/bin/env bash
#SBATCH --account=""
#SBATCH --partition=""
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:0
#SBATCH --job-name=""
#SBATCH --mail-user=""
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_%j.log
pwd; hostname; date

module load software/nextflow-23.04.3

cd <run directory>

work="<work directory>"

#rm -rf ${work}

nextflow run wes/getFastqBamQualityReports.nf -w ${work} -profile local -resume --threads 4 --njobs 10

#nextflow run wes/trimReads.nf -w ${work} -profile uctada -resume

#nextflow run wes/alignReadsToReference.nf -w ${work} -profile sadacclong,hg38 -resume --with-report --njobs 40 --threads 11 --sparkMode false

#nextflow wes/converBam2Cram.nf -w ${work} -profile sadacclong,hg38 -resume --njobs 40 --threads 11

#rsync -r <CRAM source directory> <CRAM destination directory>

#nextflow run wes/postAlignmentProcessing.nf -w ${work} -profile sadacclong,hg38 -resume --with-report --njobs 20 --threads 20 --sparkMode false

#nextflow run wes/callVariants.nf -w ${work} -profile sadacclong,hg38 -resume --njobs 40 --threads 11

#nextflow run wes/genotypegvcf.nf -w ${work} -profile sadacclong,hg38 -resume

#nextflow run wes/filterVariantCalls.nf -w ${work} -profile sadacclong,hg38 -resume

