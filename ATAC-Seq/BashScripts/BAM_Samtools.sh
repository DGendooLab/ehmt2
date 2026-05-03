#!/bin/bash
#SBATCH --ntasks=30
#SBATCH --time=4320:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2024a
module load SAMtools/1.21-GCC-13.3.0

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/

# B1-1-C
samtools view -@ 25 -b /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/B1-1-C.sam -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B1-1-C.bam

# B2-1-KD
samtools view -@ 25 -b /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/B2-1-KD.sam -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B2-1-KD.bam

# B3-1-C-T
samtools view -@ 25 -b /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/B3-1-C-T.sam -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B3-1-C-T.bam

# B4-1-KD-T
samtools view -@ 25 -b /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/B4-1-KD-T.sam -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/B4-1-KD-T.bam

# C1-1-C
samtools view -@ 25 -b /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/C1-1-C.sam -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C1-1-C.bam

# C2-1-KD
samtools view -@ 25 -b /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/C2-1-KD.sam -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C2-1-KD.bam

# C3-1-C-T
samtools view -@ 25 -b /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/C3-1-C-T.sam -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C3-1-C-T.bam

# C4-1-KD-T
samtools view -@ 25 -b /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/SAM_Bowtie2/C4-1-KD-T.sam -o /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/BAM_Samtools/C4-1-KD-T.bam
