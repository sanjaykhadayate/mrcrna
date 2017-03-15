# load libraries
library(Rsubread)
library(limma)
library(edgeR)
library(DESeq2)
library(BiocParallel)
library(GenomicAlignments)
library(rtracklayer)

source("config.r")

# Create output directories
folders<-c("Alignments","Counts","BigWigs","GeneExpression","QC","Splicing","GeneOntology")

for (i in 1:length(folders))  {
  dir.create(paste(getwd(),folders[i], sep="/"))
} 


# read in target file
targets <- readTargets(samplesheet)
targets

#buildindex(basename="mm9",reference="chr1.fa")

index=paste0(resource,genome)

#setup workers
multicoreParam <- MulticoreParam(threads)

rnaaln <- function(x)
{
  read1=paste0(rawdata,targets[x,]$InputFile)
  read2=paste0(rawdata,targets[x,]$InputFile2)
  
  # align reads
  if (isPairedEnd == TRUE)
  {
    subjunc(index=index,readfile1=read1,readfile2=read2,input_format="gzFASTQ",output_format="BAM",output_file=paste0("Alignments/",targets[x,]$OutputFile),nthreads=threads,unique=TRUE,indels=5)
    
  }else {
    subjunc(index=index,readfile1=read1,input_format="gzFASTQ",output_format="BAM",output_file=paste0("Alignments/",targets$OutputFile),nthreads=threads,unique=TRUE,indels=5)
  }
}

getbw <- function(x)
{
  Alignments <- readGAlignments(paste0("Alignments/",targets[x,]$OutputFile))
  Coverage <- coverage(Alignments)
  export.bw(Coverage[1], paste0("BigWigs/",targets[x,]$sample,".bw")) # exporting chr1
  
}

#perform alignments in parallel
bplapply(seq(nrow(targets)),rnaaln)

# Get bigwigs
bplapply(seq(nrow(targets)),getbw)

# Alignment statistics
alignstats<-propmapped(paste0("Alignments/",targets$OutputFile),countFragments=TRUE,properlyPaired=FALSE)
write.table(alignstats,file=paste0("QC/","AlignmentSummary.txt"),sep="\t")

# count numbers of reads mapped to NCBI Refseq genes
featureCountsobj <-featureCounts(files=paste0("Alignments/",targets$OutputFile),annot.inbuilt=genome,nthreads=threads,strandSpecific=strandspecific,isPairedEnd=isPairedEnd)
write.csv(featureCountsobj$counts,file=(paste0("Counts/","Rawcounts.csv")),row.names = T)

if (is.null(factor2))
{
Group<-factor1
design<-formula(~Group)
colData<-cbind(targets$OutputFile,targets$Group)
rownames(colData)<-colData[,1]
colnames(colData)<-c("name","Group")
colData<-data.frame(colData)
dds<-DESeqDataSetFromMatrix(countData= featureCountsobj$counts,colData= targets,design=design)
} else {
  Group1<-factor1
  Group2<-factor2
  design<-formula(~Group1 + Group2)
  colData<-cbind(targets$OutputFile,targets$Group1, targets$Group2)
  rownames(colData)<-colData[,1]
  colnames(colData)<-c("name",Group1,Group2)
  
  colData<-data.frame(colData)
  dds<-DESeqDataSetFromMatrix(countData= featureCountsobj$counts,colData= targets,design=design)
  
}

dds<-DESeq(dds)

res<-results(dds)
resOrdered<-res[order(res$padj),]


write.table(resOrdered,file=paste0("GeneExpression/","DEresult.csv"),sep=",")

# Plot dispersions
pdf(file=paste0("QC/","qc-plots.pdf"))
plotDispEsts(dds, main="Dispersion plot")


# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))
plotPCA(rld,intgroup=grep("^Gr",colnames(colData(dds)),value=T))

dev.off()
# Splicing analysis
countsExons <-featureCounts(files=paste0("Alignments/",targets$OutputFile),annot.inbuilt=genome,nthreads=threads,strandSpecific=strandspecific,isPairedEnd=isPairedEnd,GTF.featureType="exon", GTF.attrType="ID",
                      useMetaFeatures=FALSE, allowMultiOverlap=TRUE)

dgeList <- DGEList(counts=countsExons$counts, genes=countsExons$annotation)

sums <- rowSums(dgeList$counts)
dgeList <- dgeList[sums>10,,keep.lib.sizes=FALSE]
dgeList <- calcNormFactors(dgeList)

if (is.null(factor2))
{
design <- model.matrix(~targets$Group)
} else {
  design <- model.matrix(~targets$Group1 + targets$Group2 )
  
}

voomOut <- voom(dgeList,design,plot=TRUE)
fit <- lmFit(voomOut,design)

spliceOut <- diffSplice(fit, geneid="GeneID")


spliced<-topSplice(spliceOut,coef=2,test="F")

altUsed<-topSplice(spliceOut,coef=2,test="t")

write.table(spliced,file=paste0("Splicing/","AltSplicedgenes.txt"),sep="\t")
write.table(altUsed,file=paste0("Splicing/","DiffUsedExons.txt"),sep="\t")

# functional analysis for DE genes
library(goseq)
anno_version=genome;

newX<- resOrdered[complete.cases(resOrdered$padj),]
degenes<-as.integer(newX$padj<0.05)
names(degenes)<-rownames(newX)
# remove duplicate gene names
degenes<-degenes[match(unique(names(degenes)),names(degenes))]
table(degenes)
# GO and kegg analysis
# (1) fitting the probability weighting function (PWF)
pwf=nullp(degenes,anno_version,'knownGene')
# nullp: probability weighting function
# (2) Using the Wallenius approximation
# change the Keggpath id to name in the goseq output
library(KEGG.db)
xx <- as.list(KEGGPATHID2NAME)

temp <- cbind(names(xx),unlist(xx))
addKeggTogoseq <- function(JX,temp){
  for(l in 1:nrow(JX)){
    if(JX[l,1] %in% temp[,1]){
      JX[l,"term"] <- temp[temp[,1] %in% JX[l,1],2]
      JX[l,"ontology"] <- "KEGG"
    }
  }
  return(JX)
}
functional_analysis=goseq(pwf,anno_version,'knownGene',test.cats=c("GO:BP","GO:MF","KEGG"))
restemp<-addKeggTogoseq(functional_analysis,temp)    # switch Keggpathid to name
write.table(restemp,file=paste0("GeneOntology/","GO_Kegg_Wallenius.txt"),row.names=F,sep="\t")

outnames<-c("countsExons","dds","dgeList","featureCountsobj","fit","res","resOrdered","spliceOut","voomOut","rld")
save(list=outnames,file="analysis.RData")


sessionInfo()

