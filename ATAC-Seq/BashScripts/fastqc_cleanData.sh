#!/bin/bash
#SBATCH --ntasks=15
#SBATCH --time=120:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2023a
module load FastQC/0.11.9-Java-11

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc

# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/B1-1-C_R1_001_val_1.fq.gz
# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/B1-1-C_R2_001_val_2.fq.gz
# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/B2-1-KD_R1_001_val_1.fq.gz
# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/B2-1-KD_R2_001_val_2.fq.gz
# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/B3-1-C-T_R1_001_val_1.fq.gz
# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/B3-1-C-T_R2_001_val_2.fq.gz
# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/B4-1-KD-T_R1_001_val_1.fq.gz
# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/B4-1-KD-T_R2_001_val_2.fq.gz


# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/C1-1-C_R1_001_val_1.fq.gz
# fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/C1-1-C_R2_001_val_2.fq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/C2-1-KD_R1_001_val_1.fq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/C2-1-KD_R2_001_val_2.fq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/C3-1-C-T_R1_001_val_1.fq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/C3-1-C-T_R2_001_val_2.fq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/C4-1-KD-T_R1_001_val_1.fq.gz
fastqc -t 12 -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/fastqc /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CleanData/C4-1-KD-T_R2_001_val_2.fq.gz