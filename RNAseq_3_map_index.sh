module load star/2.7.1a 
cd $PBS_O_WORKDIR 
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ./Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ./Homo_sapiens.GRCh38.106.gtf --sjdbOverhang 75