
**Pipeline for rnaseq data analysis at MRC LMS**

The pipeline performs following steps:

 1. Spliced alignment using subreads aligner.
 
 2. Counting using feature counts module of subreads package.
 
 3. Differential gene expression analysis using deseq2 package. 
 
 4. Differential splicing analysis using limma. 
 
 5. Gene ontology and pathway analysis is performed using GoSeq package.
 

**How to run the pipeline**

Run the pipeline by executing following code using a jobscript.
         
         Rscript rnapipeline0.1.r 
         
 Required input files to run the script are,
 

1) targets.r is config file that contains all the required parameters to run the pipeline. Modify it to change the parameters.

2) targets.txt is  sample information file.

Both of these files needs to be present in working directory.

genome can be either mm9 or hg19. Genome index built using rsubreads package needs to be provided. 

strandspecific (0 for non strand specific or 1 for strand specific)

dir is the directory where raw data is located 

The pipeline was tested on  R version 3.3.2 (2016-10-31)

