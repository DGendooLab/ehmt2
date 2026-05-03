#!/bin/bash
#SBATCH --ntasks=30
#SBATCH --time=4320:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2023a
module load Bowtie2/2.5.4-GCC-12.3.0

cd /rds/projects/g/gendood-preclinomics/EHMT2/SAM_Bowtie2/

# B1-1-C
bowtie2 -p 25 -x /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index_hg38/Bowtie2_index -1 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B1-1-C_R1_001.fastq.gz -2 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B1-1-C_R2_001.fastq.gz -S /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/B1-1-C.sam

# B2-1-KD
bowtie2 -p 25 -x /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index_hg38/Bowtie2_index -1 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B2-1-KD_R1_001.fastq.gz -2 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B2-1-KD_R2_001.fastq.gz -S /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/B2-1-KD.sam

# B3-1-C-T
bowtie2 -p 25 -x /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index_hg38/Bowtie2_index -1 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B3-1-C-T_R1_001.fastq.gz -2 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B3-1-C-T_R2_001.fastq.gz -S /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/B3-1-C-T.sam

# B4-1-KD-T
bowtie2 -p 25 -x /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index_hg38/Bowtie2_index -1 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B4-1-KD-T_R1_001.fastq.gz -2 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/B4-1-KD-T_R2_001.fastq.gz -S /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/B4-1-KD-T.sam

# C1-1-C
bowtie2 -p 25 -x /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index_hg38/Bowtie2_index -1 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C1-1-C_R1_001.fastq.gz -2 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C1-1-C_R2_001.fastq.gz -S /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/C1-1-C.sam

# C2-1-KD
bowtie2 -p 25 -x /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index_hg38/Bowtie2_index -1 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C2-1-KD_R1_001.fastq.gz -2 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C2-1-KD_R2_001.fastq.gz -S /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/C2-1-KD.sam

# C3-1-C-T
bowtie2 -p 25 -x /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index_hg38/Bowtie2_index -1 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C3-1-C-T_R1_001.fastq.gz -2 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C3-1-C-T_R2_001.fastq.gz -S /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/C3-1-C-T.sam

# C4-1-KD-T
bowtie2 -p 25 -x /rds/projects/g/gendood-preclinomics/EHMT2/References_Index/Bowtie2_index_hg38/Bowtie2_index -1 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C4-1-KD-T_R1_001.fastq.gz -2 /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/00_fastq/C4-1-KD-T_R2_001.fastq.gz -S /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/C4-1-KD-T.sam