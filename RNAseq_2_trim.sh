
#!/bin/bash
#PBS -l select=1:ncpus=2:mem=3GB
#PBS -l walltime=72:00:00
#PBS -J 1-4

module load htseq/0.6.1
module load samtools/1.2
module load gcc/5.4.0
module load intel-suite/2016.3
module load mpi/intel-5.1
module load anaconda/2.1.0      
module load subread/1.5.2

homeDir="/rds/general/user/mmoradim/home/my_projects/RNAseq_hepatitis_Morfopoulou_et_al/"
inputDir="/rds/general/user/mmoradim/ephemeral/RNAseq_hepatitis_Morfopoulou_et_al/input"
outputDir="/rds/general/user/mmoradim/ephemeral/RNAseq_hepatitis_Morfopoulou_et_al/trim_output/"

trimmomaticPath="/rds/general/user/mmoradim/home/tools/RNAseq_pipelines_dpt/star_featurecounts/additional_data/trimmomatic-0.36.jar"

scriptsDir="$homeDir/scripts/1_upstream/2_trim"

sampleID_file="$scriptsDir/sampleIDs.txt"
cat $sampleID_file | sed 's/.\///' > sample_list
# num_samples=`cat $sample_list | wc -l | tr -d ' '`

sample=`sed "\${PBS_ARRAY_INDEX}q;d" sample_list`
	shortname=`basename $sample`
	echo $shortname
	inFile1=$inputDir/${sample}"_R1_001.fastq.gz"
	inFile01=$(echo "${inFile1}" | tr -d '\r')

	inFile2=$inputDir/${sample}"_R2_001.fastq.gz"
	inFile02=$(echo "${inFile2}" | tr -d '\r')
		
	java -jar $trimmomaticPath PE -phred33 -baseout $outputDir/$shortname.fq.gz $inFile01 $inFile02  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
