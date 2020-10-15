## Load libraries
library(ChIPQC)
## Load sample data
samples <- read.csv('/path/to/samplesheet_samples.csv')
chipObj <- ChIPQC(samples, annotation="mm10", chromosomes=NULL) 
ChIPQCreport(chipObj, reportName="ChIP QC report of smples", reportFolder="ChIPQCreport_samples",facet=FALSE)
