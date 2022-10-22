
#!/bin/bash
#PBS -l select=1:ncpus=2:mem=3GB
#PBS -l walltime=10:00:00
#PBS -J 1-8

homeDir="/rds/general/user/mmoradim/home/my_projects/RNAseq_hepatitis_Morfopoulou_et_al/"

inputDir="/rds/general/user/mmoradim/ephemeral/RNAseq_hepatitis_Morfopoulou_et_al/input"
outputDir="/rds/general/user/mmoradim/ephemeral/RNAseq_hepatitis_Morfopoulou_et_al/qc1_output/"

fastqcPath="/rds/general/user/mmoradim/home/tools/FastQC/fastqc"

scriptsDir="$homeDir/scripts/1_upstream/1_qc1"

sampleID_file="$scriptsDir/sampleIDs.txt"
cat $sampleID_file | sed 's/.\///' > sample_list
# num_samples=`cat $sample_list | wc -l | tr -d ' '`

sample=`sed "\${PBS_ARRAY_INDEX}q;d" sample_list`
	shortname=`basename $sample`
	echo $shortname
	
	inFile1=$inputDir/${sample}"_001.fastq.gz"
	inFile01=$(echo "${inFile1}" | tr -d '\r')

	$fastqcPath $inFile01 -o $outputDir -d $outputDir 
