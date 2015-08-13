# mrcrna
pipeline for rnaseq analysis at CSC
The pipeline performs spliced alignment using subreads aligner, counting using feature counts module of subreads package.
It then performs differential gene expression analysis using deseq2 package. It also performs differential splicing analysis using limma.

How to run the pipeline

Run the pipeline by executing following code using a jobscript.
         
         Rscript rnapipeline0.1.r file=targets.txt genome=mm9 strandspecific=0 factor1=Group factor2=NULL dir=/csc/rawdata/-----/Unalined/ isPairedEnd=TRUE


targets.txt is  sample information file and needs to be present in working directory.
The file needs to have following columns
sample  Group   InputFile       InputFile2      OutputFile


genome can be either mm9 or hg19
strandspecific (0 for non strand specific or 1 for strand specific)
dir is the directory where raw data is located 

The pipeline was developed on  R/3.2.0

