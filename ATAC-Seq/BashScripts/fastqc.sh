#!/bin/bash
#SBATCH --ntasks=15
#SBATCH --time=80:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2023a
module load FastQC/0.11.9-Java-11

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/

fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B1-1-C_R1_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B1-1-C_R2_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B2-1-KD_R1_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B2-1-KD_R2_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B3-1-C-T_R1_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B3-1-C-T_R2_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B4-1-KD-T_R1_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B4-1-KD-T_R2_001.fastq.gz


fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C1-1-C_R1_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C1-1-C_R2_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C2-1-KD_R1_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C2-1-KD_R2_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C3-1-C-T_R1_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C3-1-C-T_R2_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C4-1-KD-T_R1_001.fastq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/QC/fastqc/ /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C4-1-KD-T_R2_001.fastq.gz