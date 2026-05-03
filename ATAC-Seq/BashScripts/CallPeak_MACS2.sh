#!/bin/bash
#SBATCH --ntasks=30
#SBATCH --time=4320:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2022a
module load MACS2/2.2.9.1-foss-2022a

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CallPeak_MACS2/

# B2 vs B1

macs2 callpeak -t /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B2-1-KD_rmdup.bam -c /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B1-1-C_rmdup.bam -f BAM -g hs -n B2_vs_B1 --outdir /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CallPeak_MACS2/B2_vs_B1/ --nomodel --shift -130 --extsize 260 -p 0.05

# B4 vs B3
macs2 callpeak -t /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B4-1-KD-T_rmdup.bam -c /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B3-1-C-T_rmdup.bam -f BAM -g hs -n B4_vs_B3 --outdir /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CallPeak_MACS2/B4_vs_B3/ --nomodel --shift -131 --extsize 262 -p 0.05

# C2 vs C1
macs2 callpeak -t /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C2-1-KD_rmdup.bam -c /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C1-1-C_rmdup.bam -f BAM -g hs -n C2_vs_C1 --outdir /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CallPeak_MACS2/C2_vs_C1/ --nomodel --shift -134 --extsize 268 -p 0.05

# C4 vs C3
macs2 callpeak -t /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C4-1-KD-T_rmdup.bam -c /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C3-1-C-T_rmdup.bam -f BAM -g hs -n C4_vs_C3 --outdir /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/CallPeak_MACS2/C4_vs_C3 --nomodel --shift -134 --extsize 268 -p 0.05