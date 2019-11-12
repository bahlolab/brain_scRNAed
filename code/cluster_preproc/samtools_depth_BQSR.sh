#PBS -l nodes=1:ppn=1,mem=8gb,walltime=0:30:00

module load samtools

dataDir=/scRNA_editing/study_SRA/gatk
outDir=/scRNA_editing/study_SRA/gatk_samtools_BQSRdepth

samtools depth -b dt_siteStats_TDjoin_sites.bed \
	$dataDir/MYFILE_split_bqsr.bam  > $outDir/MYFILE_40K_sites_BQSRdepth.txt
