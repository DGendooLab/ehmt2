#!/bin/bash
#SBATCH --ntasks=20
#SBATCH --time=2880:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2023a
module load Trim_Galore/0.6.10-GCCcore-12.3.0

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/


trim_galore --paired /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B1-1-C_R1_001.fastq.gz /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B1-1-C_R2_001.fastq.gz --gzip -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/

trim_galore --paired /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B2-1-KD_R1_001.fastq.gz /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B2-1-KD_R2_001.fastq.gz --gzip -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/

trim_galore --paired /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B3-1-C-T_R1_001.fastq.gz /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B3-1-C-T_R2_001.fastq.gz --gzip -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/

trim_galore --paired /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B4-1-KD-T_R1_001.fastq.gz /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B4-1-KD-T_R2_001.fastq.gz --gzip -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/



trim_galore --paired /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C1-1-C_R1_001.fastq.gz /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C1-1-C_R2_001.fastq.gz --gzip -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/

trim_galore --paired /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C2-1-KD_R1_001.fastq.gz /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C2-1-KD_R2_001.fastq.gz --gzip -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/

trim_galore --paired /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C3-1-C-T_R1_001.fastq.gz /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C3-1-C-T_R2_001.fastq.gz --gzip -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/

trim_galore --paired /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C4-1-KD-T_R1_001.fastq.gz /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C4-1-KD-T_R2_001.fastq.gz --gzip -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/
