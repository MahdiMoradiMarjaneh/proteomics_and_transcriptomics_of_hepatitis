
#!/bin/bash
#PBS -l select=1:ncpus=2:mem=10GB
#PBS -l walltime=72:00:00
#PBS -J 1-4

module load samtools/1.2
module load gcc/5.4.0
module load intel-suite/2016.3
module load anaconda/2.1.0      
module load subread/1.5.2

homeDir="/rds/general/user/mmoradim/home/my_projects/RNAseq_hepatitis_Morfopoulou_et_al/"
inputDir="/rds/general/user/mmoradim/ephemeral/RNAseq_hepatitis_Morfopoulou_et_al/map_output"
outputDir="/rds/general/user/mmoradim/ephemeral/RNAseq_hepatitis_Morfopoulou_et_al/count_output"

annotationFile="/rds/general/user/mmoradim/home/my_projects/RNAseq_hepatitis_Morfopoulou_et_al/scripts/1_upstream/4_map_index/Homo_sapiens.GRCh38.106.gtf"

scriptsDir="$homeDir/scripts/1_upstream/7_count"

sampleID_file="$scriptsDir/sampleIDs.txt"
cat $sampleID_file | sed 's/.\///' > sample_list
# num_samples=`cat $sample_list | wc -l | tr -d ' '`

sample=`sed "\${PBS_ARRAY_INDEX}q;d" sample_list`
	shortname=`basename $sample`
	echo $shortname
	inFile1=$inputDir/${sample}"_Aligned.sortedByCoord.out.bam"
	inFile01=$(echo "${inFile1}" | tr -d '\r')
	
	featureCounts -p -J -C  -a $annotationFile -o $outputDir/$shortname"-gene.fc" -s 2 $inFile01
