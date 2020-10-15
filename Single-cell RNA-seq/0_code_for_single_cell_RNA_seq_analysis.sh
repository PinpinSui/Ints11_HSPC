#Author Pinpin Sui
#e-Mail suipinpin0722@163.com
#####################################code for data preprocessing

#Step 1 CBL to fastq
/path/to/cellranger mkfastq --run /path/to/CBL_data \
--csv /path/to/cellranger-index.csv --qc \
--output-dir /path/to/output/

#Step 2 Count matrix generating
/path/to/cellranger count --id=KO --fastqs=/path/to/sapmle_prefix \
--sample=sample --transcriptome=/path/to/refdata-cellranger-mm10-3.0.0/ \
--localmem=60 --localcores=12 

#Step 3 Doublets prediction
#Details in 1_Doublets_prediction.py

#Step 4 Clustering and Differential expression analysis with Seurat
#Details in 2_Clustering&Differential_expression_analysis.r

