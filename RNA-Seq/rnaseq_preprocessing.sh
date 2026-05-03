#!/bin/bash
#SBATCH --account gendood-preclinomics
#SBATCH --qos bbdefault
#SBATCH --time 96:0:0
#SBATCH --nodes 1
#SBATCH --ntasks 2

date

set -e

module purge; module load bluebear
module load bear-apps/2024a; module load Nextflow/25.04.6

cd /rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/

nextflow run \
    nf-core/rnaseq \
    -resume \
    --input /rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/samplesheet.csv \
    --outdir /rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/results \
    --fasta /rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/results/genome/GRCh38.p14.genome.fa \
    --gtf /rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/results/genome/gencode.v45.annotation.gtf \
    --save_reference \
    -c /rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/nextflow.config \
    -profile bluebear