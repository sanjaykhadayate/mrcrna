# mrcrna
Pipeline for rnaseq data analysis at MRC LMS

The pipeline performs spliced alignment using subreads aligner, counting using feature counts module of subreads package.

It then performs differential gene expression analysis using deseq2 package. 
It also performs differential splicing analysis using limma. Gene ontology and pathway analysis is performed using GoSeq package.

**How to run the pipeline**

Run the pipeline by executing following code using a jobscript.
         
         Rscript rnapipeline0.1.r 

targets.r is config file that contains all the required parameters to run the pipeline. Modify it to change the parameters.

targets.txt is  sample information file and needs to be present in working directory.


genome can be either mm9 or hg19

strandspecific (0 for non strand specific or 1 for strand specific)

dir is the directory where raw data is located 

The pipeline was developed on  R/3.2.0

