#PBS -l nodes=1:ppn=1,mem=12gb,walltime=6:00:00 

module load picard-tools
module load java

dataDir=/scRNA_editing/study_SRA/STAR_2pass
writeDir=/scRNA_editing/study_SRA/gatk
tmpDir=/scRNA_editing/study_SRA/gatk/tmp_scratch

mkdir $tmpDir/MYFILEtmp

genome_fa=/software/GATK_haplotypeCaller/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

java -Xmx12g -jar /stornext/System/data/apps/picard-tools/picard-tools-2.9.4/picard.jar AddOrReplaceReadGroups  \
INPUT=$dataDir/STAR_MYFILEAligned.sortedByCoord.out.bam \
OUTPUT=$writeDir/MYFILEAligned.sortedByCoord.out_readreplace.bam \
RGID=1 \
RGLB=lib1 \
RGPL=illumina \
RGPU=unit1 \
RGSM=1


java -Xmx12g -jar /stornext/System/data/apps/picard-tools/picard-tools-2.9.4/picard.jar MarkDuplicates \
INPUT=$writeDir/MYFILEAligned.sortedByCoord.out_readreplace.bam \
OUTPUT=$writeDir/MYFILE_dedupped.bam  \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=$writeDir/MYFILE.metrics 

rm $writeDir/MYFILEAligned.sortedByCoord.out_readreplace.bam

module purge 
module load gatk/3.7.0
module load java


java -Xmx12g -jar /stornext/System/data/apps/gatk/gatk-3.7.0/GenomeAnalysisTK.jar -T SplitNCigarReads \
-R $genome_fa \
-I $writeDir/MYFILE_dedupped.bam \
-o $writeDir/MYFILE_split.bam \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

rm $writeDir/MYFILE_dedupped.bam 

java -Xmx12g -jar /stornext/System/data/apps/gatk/gatk-3.7.0/GenomeAnalysisTK.jar -T BaseRecalibrator \
-R $genome_fa \
-I $writeDir/MYFILE_split.bam \
-knownSites /wehisan/bioinf/lab_bahlo/public_datasets/dbSNP/sorted.vcf \
-o $writeDir/MYFILE_recal_data.table

java -Xmx12g -jar /stornext/System/data/apps/gatk/gatk-3.7.0/GenomeAnalysisTK.jar -T PrintReads \
-R $genome_fa \
-I $writeDir/MYFILE_split.bam \
-BQSR $writeDir/MYFILE_recal_data.table \
-o $writeDir/MYFILE_split_bqsr.bam

rm $writeDir/MYFILE_split.bam

java -Xmx12g -jar /stornext/System/data/apps/gatk/gatk-3.7.0/GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $genome_fa \
-I $writeDir/MYFILE_split_bqsr.bam \
-dontUseSoftClippedBases \
-stand_call_conf 0 \
-o $writeDir/MYFILE_var_call.vcf

#NB retain $writeDir/MYFILE_split_bqsr.bam
