# featurecounts.sh

#!/bin/bash

#PBS -l nodes=1:ppn=1,mem=4gb,walltime=4:00:00 

module load  subread/1.5.2

dataDir="scRNA_editing/study_SRA/STAR_2pass"
outDir="scRNA_editing/study_SRA/featureCounts"


featureCounts -Q 0 -T 1 -a /wehisan/bioinf/Bioinformatics/SNPchipdata/haloom/genome/Homo_sapiens.GRCh38.91/Homo_sapiens.GRCh38.91.gtf \
-o $outDir/STAR_MYFILE_fCounts_unstranded  $dataDir/STAR_MYFILEAligned.sortedByCoord.out.bam


