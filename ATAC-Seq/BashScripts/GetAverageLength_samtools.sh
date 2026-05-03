#!/bin/bash
#SBATCH --ntasks=30
#SBATCH --time=2880:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2023a
module load SAMtools/1.21-GCC-12.3.0

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/

samtools stats /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B1-1-C_rmdup.bam | grep "average length"

samtools stats /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B2-1-KD_rmdup.bam | grep "average length"

samtools stats /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B3-1-C-T_rmdup.bam | grep "average length"

samtools stats /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B4-1-KD-T_rmdup.bam | grep "average length"

samtools stats /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C1-1-C_rmdup.bam | grep "average length"

samtools stats /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C2-1-KD_rmdup.bam | grep "average length"

samtools stats /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C3-1-C-T_rmdup.bam | grep "average length"

samtools stats /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C4-1-KD-T_rmdup.bam | grep "average length"
