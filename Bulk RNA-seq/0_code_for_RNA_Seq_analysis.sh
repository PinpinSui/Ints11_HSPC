#Author: Pinpin Sui
#e-Mail: suipinpin0722@163.com
##########################################code for Bulk RNA-seq sequencing
#Step 1 Trimming for raw sequenced data
java -jar trimmomatic-0.38.jar  PE \
sample_R1.fastq sample_R2.fastq \
sample_R1_paired.fq sample_R1_unpaired.fq \
sample_R2_paired.fq sample_R2_unpaired.fq \
ILLUMINACLIP:/path/to/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Step 2 Quality control
mkdir fastQC_result_sample
rawdata_1=sample_R1_paired.fq
rawdata_2=sample_R2_paired.fq
/path/to/fastqc -o fastQC_result_sampple $rawdata_1 $rawdata_2 -f fastq -t 2 

#Step 3 Maaping to the mm10 reference
rawdata_1=sample_R1_paired.fq
rawdata_2=sample_R2_paired.fq

/path/to/STAR --runThreadN 3 \
--readFilesIn $rawdata_1 $rawdata_2 \
--quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix sample_ \
--genomeDir /path/to/STARIndex/ 

#Step 4 Raw counts calculating for each gene
/path/to/htseq-count -f bam -s no sample_Aligned.sortedByCoord.out.bam \
/path/to/mm10_genes.gtf > sample.count.txt


