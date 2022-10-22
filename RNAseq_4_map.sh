
#!/bin/bash
#PBS -l select=1:ncpus=2:mem=10GB
#PBS -l walltime=72:00:00
#PBS -J 1-4

module load star/2.7.1a

module load htseq/0.6.1
module load samtools/1.2
module load gcc/5.4.0
module load intel-suite/2016.3
module load mpi/intel-5.1
module load anaconda/2.1.0      
module load subread/1.5.2

homeDir="/rds/general/user/mmoradim/home/my_projects/RNAseq_hepatitis_Morfopoulou_et_al/"
inputDir="/rds/general/user/mmoradim/ephemeral/RNAseq_hepatitis_Morfopoulou_et_al/trim_output/"
outputDir="/rds/general/user/mmoradim/ephemeral/RNAseq_hepatitis_Morfopoulou_et_al/map_output"
genomeDir="/rds/general/user/mmoradim/home/my_projects/RNAseq_hepatitis_Morfopoulou_et_al/scripts/1_upstream/4_map_index"

scriptsDir="$homeDir/scripts/1_upstream/5_map"

sampleID_file="$scriptsDir/sampleIDs.txt"
cat $sampleID_file | sed 's/.\///' > sample_list
# num_samples=`cat $sample_list | wc -l | tr -d ' '`

sample=`sed "\${PBS_ARRAY_INDEX}q;d" sample_list`
	shortname=`basename $sample`
	echo $shortname
	inFile1=$inputDir/${sample}"_1P.fq.gz"
	inFile01=$(echo "${inFile1}" | tr -d '\r')

	inFile2=$inputDir/${sample}"_2P.fq.gz"
	inFile02=$(echo "${inFile2}" | tr -d '\r')

    STAR --runThreadN 1 \
		 --genomeDir $genomeDir \
		 --readFilesCommand zcat \
		 --readFilesIn $inFile01 $inFile02 \
         --outFileNamePrefix $outputDir/$shortname"_" \
		 --outSAMtype BAM SortedByCoordinate \
		 --runMode alignReads \
		 --quantMode GeneCounts \
		 --outReadsUnmapped Fastx \
		 --chimSegmentMin 12 \
