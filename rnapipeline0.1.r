# load libraries
library(Rsubread)
library(limma)
library(edgeR)
library(DESeq2)


args<- commandArgs()
print(args)

# read in target file
targets <- readTargets(unlist(strsplit(args[6],"="))[2])
targets


dir=unlist(strsplit(args[11],"="))[2]
read1=paste0(dir,targets$InputFile)
read2=paste0(dir,targets$InputFile2)
index=paste0("/csc/skhadaya/resources/",unlist(strsplit(args[7],"="))[2])

# align reads
if (unlist(strsplit(args[12],"="))[2] == TRUE)
{
subjunc(index=index,readfile1=read1,readfile2=read2,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,nthreads=8,tieBreakHamming=TRUE,unique=TRUE,indels=5)

}
else (
subjunc(index=index,readfile1=read1,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,nthreads=8,tieBreakHamming=TRUE,unique=TRUE,indels=5)
)
pdf(file="plots.pdf")


# count numbers of reads mapped to NCBI Refseq genes
fc <-featureCounts(files=targets$OutputFile,annot.inbuilt=unlist(strsplit(args[7],"="))[2],nthreads=8,isPairedEnd=unlist(strsplit(args[12],"="))[2])

#design<-if (unlist(strsplit(args[9],"="))[2] == "") formula(~unlist(strsplit(args[9],"="))[2]) else formula(~unlist(strsplit(args[9],"="))[2]+unlist(strsplit(args[10],"="))[2])
#strsplit(args[9],"=")[2]

#design<-formula(~(unlist(strsplit(args[9],"="))[2]))
design<-formula(~Group)


colData<-cbind(targets$OutputFile,targets$Group)
rownames(colData)<-colData[,1]
colnames(colData)<-c("name","Group")

colData<-data.frame(colData)
dds<-DESeqDataSetFromMatrix(countData= fc$counts,colData= targets,design=design)

dds<-DESeq(dds)

res<-results(dds)
resOrdered<-res[order(res$padj),]

dev.off()
write.table(resOrdered,file="DEgenes.txt",sep="\t")
fc <-featureCounts(files=targets$OutputFile,annot.inbuilt=unlist(strsplit(args[7],"="))[2],nthreads=8,isPairedEnd=unlist(strsplit(args[12],"="))[2],GTF.featureType="exon", GTF.attrType="ID",
useMetaFeatures=FALSE, allowMultiOverlap=TRUE)

dge <- DGEList(counts=fc$counts, genes=fc$annotation)

A <- rowSums(dge$counts)
dge <- dge[A>10,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~targets$Group)
v <- voom(dge,design,plot=TRUE)
 fit <- lmFit(v,design)

ex <- diffSplice(fit, geneid="GeneID")


spliced<-topSplice(ex,coef=2,level="gene")

altUsed<-topSplice(ex,coef=2,level="exon")

write.table(spliced,file="AltSplicedgenes.txt",sep="\t")
write.table(altUsed,file="DiffUsedExons.txt",sep="\t")


sessionInfo()
