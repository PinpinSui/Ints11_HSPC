#Author: Pinpin Sui
#e-Mail: suipinpin0722@163.com
#Step 1 Raw data trimming
/path/to/trim_galore -q 25 --phred33 --length 35 -e 0.1 --stringency 2  -o ./02_clean/ $rawdatapath/sample.fastq.gz
 
#Step 2 Quality control
/path/to/fastqc -t 4  ../02_clean/*gz -o clean/
/path/to/multiqc ./*zip

#Step 3 Mapping to the mm10 reference
zcat sample_trimmed.fq.gz  > sample.fastq
bowtie2 -x /path/to/mm10/Bowtie2Index/genome -U sample.fastq  | samtools sort -@ 5 -O bam -o sample.bam
sambamba markdup -r sample.bam  sample.sambamba.rmdup.bam
samtools flagstat sample.bam
samtools flagstat sample.sambamba.rmdup.bam
mv sample.sambamba.rmdup.bam sample.last.bam
samtools index sample.last.bam
bedtools bamtobed -i sample.last.bam  > sample.bed

#Step 4 Peak calling
/path/to/macs2 callpeak -t ../04_align/sample.last.bam -c ../04_align/input.last.bam -f BAM -g mm  -B -q 0.05 -n sample

#Step 5  Experimental evaluation with deeptools
/path/to/bamCoverage --normalizeUsing CPM -b ../04_align/sample_1.last.bam -o sample_1.bw
...
/path/to/bamCoverage --normalizeUsing CPM -b ../04_align/sample_n.last.bam -o sample_n.bw
/path/to/computeMatrix reference-point  --referencePoint TSS  -p 15  \
-b 10000 -a 10000    \
-R /path/to/mm10_hgTables.txt  \
-S ../sample_1.bw ... ../sample_2.bw \
--skipZeros  -o matrix_samples_TSS_10K_final.gz  \
--outFileSortedRegions regions_samples_genes_10K_final.bed

plotHeatmap -m matrix_samples_TSS_10K_final.gz  -out samples_Heatmap_TSS_10K_final.pdf --plotFileFormat pdf  --dpi 720  
plotProfile -m matrix_samples_TSS_10K_final.gz  -out samples_Profile_TSS_10K_final.pdf --plotFileFormat pdf --perGroup --dpi 720


/path/to/computeMatrix scale-regions  -p 15  \
-R /path/to/mm10_hgTables.txt  \
-S ../sample_1.bw ... ../sample_2.bw \
-b 10000 -a 10000  \
--skipZeros -o matrix_samples_body_10K_final.gz

/path/to/plotHeatmap -m matrix_samples_body_10K_final.gz  -out samples_Heatmap_body_10K_final.pdf --plotFileFormat pdf  --dpi 720  
/path/to/plotProfile -m matrix_samples_body_10K_final.gz  -out samples_Profile_body_10K_final.pdf --plotFileFormat pdf --perGroup --dpi 720


#Step 6  Experimental evaluation with ChIPQC
/path/to/Rscript /path/to/ChIPQC.r

