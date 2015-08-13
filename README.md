# mrcrna
pipeline for rnaseq analysis at CSC

How to run the pipeline
     Run the pipeline by executing following code using a jobscript.
Rscript rnapipeline0.1.r file=targets.txt genome=mm9 strandspecific=0 factor1=Group factor2=NULL dir=/csc/rawdata/-----/Unalined/ isPairedEnd=TRUE

targets.txt is  sample information file and needs to be present in working directory.
genome can be either mm9 or hg19
strandspecific (0 for non strand specific or 1 for strand specific)
