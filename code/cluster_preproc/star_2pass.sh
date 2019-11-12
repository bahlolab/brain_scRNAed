
# STAR 2pass

#PBS -l nodes=1:ppn=4,mem=16gb,walltime=1:00:00 

module load samtools
module load STAR

ulimit -v 31654904098

genomeDIR=~/STAR_2pass/genome_sjdbOverhang_47


GTF=~/genome/Homo_sapiens.GRCh38.91/Homo_sapiens.GRCh38.91.gtf

outDir=~/scRNA_editing/study_SRA/STAR_2pass

#r1=~/GSE67835_fastq/MYFILE_1.fastq.gz

r1=~/phs000834_fastq/MYFILE_1.fastq.gz

STAR --runThreadN 4 \
--genomeDir $genomeDIR \
--sjdbGTFfile $GTF \
--sjdbOverhang  47 \
--readFilesIn $r1  \
--outFileNamePrefix $outDir/STAR_MYFILE \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic

samtools index $outDir/STAR_MYFILEAligned.sortedByCoord.out.bam


