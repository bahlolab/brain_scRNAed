# STAR 2pass for ERCC spike-in

#PBS -l nodes=1:ppn=4,mem=16gb,walltime=1:00:00 

module load samtools
module load STAR

ulimit -v 31654904098

genomeDIR=~/genome/ERCC_reference/ERCCgenome_sjdb47

GTF=~genome/ERCC_reference/ERCC92.gtf

outDir=~/scRNA_editing/study_SRA/STAR_2pass/ERCC_ref

#r1=~/GSE67835_fastq/MYFILE_1.fastq.gz

r1=~/phs000834_fastq/MYFILE_1.fastq.gz


STAR --runThreadN 4 \
--genomeDir $genomeDIR \
--sjdbGTFfile $GTF \
--sjdbOverhang  47 \
--readFilesIn $r1  \
--outFileNamePrefix $outDir/ERCC_STAR_MYFILE \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic

samtools index $outDir/ERCC_STAR_MYFILEAligned.sortedByCoord.out.bam


