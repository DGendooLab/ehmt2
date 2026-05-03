#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --time=1440:0
#SBATCH --qos=bbdefault

set -e
module purge; module load bluebear

module load bear-apps/2022b
module load R/4.3.1-foss-2022b
module load picard/2.27.5-Java-11

cd /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B1-1-C_rmdup.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B1-1-C_rmdup.insertsize -H /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B1-1-C_rmdup.hist.pdf

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B2-1-KD_rmdup.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B2-1-KD_rmdup.insertsize -H /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B2-1-KD_rmdup.hist.pdf

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B3-1-C-T_rmdup.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B3-1-C-T_rmdup.insertsize -H /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B3-1-C-T_rmdup.hist.pdf

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B4-1-KD-T_rmdup.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B4-1-KD-T_rmdup.insertsize -H /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/B4-1-KD-T_rmdup.hist.pdf

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C1-1-C_rmdup.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C1-1-C_rmdup.insertsize -H /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C1-1-C_rmdup.hist.pdf

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C2-1-KD_rmdup.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C2-1-KD_rmdup.insertsize -H /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C2-1-KD_rmdup.hist.pdf

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C3-1-C-T_rmdup.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C3-1-C-T_rmdup.insertsize -H /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C3-1-C-T_rmdup.hist.pdf

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -I /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C4-1-KD-T_rmdup.bam -O /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C4-1-KD-T_rmdup.insertsize -H /rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DeDuplication_PICARD/C4-1-KD-T_rmdup.hist.pdf