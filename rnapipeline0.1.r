# load libraries

library(Rsubread)
library(limma)
library(edgeR)
library(DESeq2)

source("config.r")

#args<- commandArgs()
#print(args)

# read in target file
targets <- readTargets(targets)
targets

#dir=unlist(strsplit(args[11],"="))[2]
read1=paste0(dir,targets$InputFile)
read2=paste0(dir,targets$InputFile2)

# use the igenome information
genomeIndex = list("hg19"="/csc/rawdata/Cscbioinf/bioinfResources/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/SubreadIndex/hg19_index",
                   "mm9"="/csc/rawdata/Cscbioinf/bioinfResources/iGenomes/Mus_musculus/Ensembl/NCBIM37/Sequence/SubreadIndex/mm9_index",
                   "dm3"="/csc/rawdata/Cscbioinf/bioinfResources/iGenomes/Drosophila_melanogaster/Ensembl/BDGP5/Sequence/SubreadIndex/bdgp5_index")
index=genomeIndex[[genome]]

# align reads
if (isPairedEnd == TRUE)
{
subjunc(index=index,readfile1=read1,readfile2=read2,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,nthreads=10,tieBreakHamming=TRUE,unique=TRUE,indels=5)

}else {
subjunc(index=index,readfile1=read1,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,nthreads=10,tieBreakHamming=TRUE,unique=TRUE,indels=5)
}

alignstats<-propmapped(targets$OutputFile,countFragments=TRUE,properlyPaired=FALSE)
write.table(alignstats,file="AlignmentSummary.txt",sep="\t")

# count numbers of reads mapped to iGenome genes
gtfGFiles = list("hg19UCSC"="/csc/rawdata/Cscbioinf/bioinfResources/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf",
              "mm9UCSC"="/csc/rawdata/Cscbioinf/bioinfResources/iGenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf",
              "dm3"="/csc/rawdata/Cscbioinf/bioinfResources/iGenomes/Drosophila_melanogaster/Ensembl/BDGP5/Annotation/Genes/genes.gtf",
              "hg19"="/csc/rawdata/Cscbioinf/bioinfResources/iGenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
              "mm9"="/csc/rawdata/Cscbioinf/bioinfResources/iGenomes/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/genes.gtf")
anno_for_featurecount<-gtfGFiles[[genome]]
fc <-featureCounts(files=targets$OutputFile,annot.ext=anno_for_featurecount,isGTFAnnotationFile=TRUE,
                GTF.featureType="exon",GTF.attrType="gene_id",nthreads=10,
                strandSpecific=strandspecific,isPairedEnd=isPairedEnd)
#fc <-featureCounts(files=targets$OutputFile,annot.inbuilt=genome,nthreads=10,strandSpecific=strandspecific,isPairedEnd=isPairedEnd)

design<-formula(~Group)

colData<-cbind(targets$OutputFile,targets$Group)
rownames(colData)<-colData[,1]
colnames(colData)<-c("name","Group")

colData<-data.frame(colData)
dds<-DESeqDataSetFromMatrix(countData= fc$counts,colData= targets,design=design)

dds<-DESeq(dds)

#####
# Quality control
#####
 # 1. sample distance plot
 # 2. PCA plot
 rlogvalue <- rlog(dds)
 rlogcount <- assay(rlogvalue)
 rlogcount <- rlogcount[!rowSums(rlogcount) == 0,]
 colnames(rlogcount) <-  paste0(colData(dds)$sample)
 library("RColorBrewer")
 mycols <- brewer.pal(8, "Dark2")[1:length(unique(colData(dds)$Group))]
 showcol<-mycols
 names(showcol)<-unique(colData(dds)$Group)
 pcafromdeseq2<-plotPCA(rlogvalue, intgroup="Group",returnData=T) +
             geom_hline(yintercept=0, colour="gray65")+
             geom_vline(xintercept=0, colour="gray65")+
             ggtitle("PCA plot")+
             scale_colour_manual(values=showcol)
 ggsave(pcafromdeseq2,file="qc_PCAplot.pdf")
 
 sampleDists <- as.matrix(dist(t(rlogcount)))
 library(gplots)
 png("qc_heatmap_samples.png", w=1000, h=1000, pointsize=20)
 heatmap.2(as.matrix(sampleDists), key=F, trace="none",
           col=colorpanel(100, "black", "white"),
           ColSideColors=mycols[colData(dds)$Group], RowSideColors=mycols[colData(dds)$Group],
           margin=c(10, 10), main="Sample Distance Matrix")
 dev.off()
 
##########
## DE analysis
##########

pdf(file="plots.pdf")

res<-results(dds)
resOrdered<-res[order(res$padj),]

dev.off()
write.table(resOrdered,file="DEgenes.txt",sep="\t")
#fcexo <-featureCounts(files=targets$OutputFile,annot.inbuilt=genome,nthreads=8,strandSpecific=strandspecific,isPairedEnd=isPairedEnd,GTF.featureType="exon", GTF.attrType="ID",useMetaFeatures=FALSE, allowMultiOverlap=TRUE)
fcexo<-featureCounts(files=targets$OutputFile,annot.ext=anno_for_featurecount,strandSpecific=strandspecific,
                    isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="exon_id",nthreads=8,
                    isPairedEnd=isPairedEnd,useMetaFeatures=FALSE, allowMultiOverlap=TRUE)

dge <- DGEList(counts=fcexo$counts, genes=fcexo$annotation)

A <- rowSums(dge$counts)
dge <- dge[A>10,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~targets$Group)
v <- voom(dge,design,plot=TRUE)
 fit <- lmFit(v,design)

ex <- diffSplice(fit, geneid="GeneID")


spliced<-topSplice(ex,coef=2,test="F")

altUsed<-topSplice(ex,coef=2,test="t")

write.table(spliced,file="AltSplicedgenes.txt",sep="\t")
write.table(altUsed,file="DiffUsedExons.txt",sep="\t")

# functional analysis for DE genes
  library(goseq)
  
  newX<- resOrdered[complete.cases(resOrdered$padj),]
  degenes<-as.integer(newX$padj<0.05)
  names(degenes)<-rownames(newX)
  # remove duplicate gene names
  degenes<-degenes[match(unique(names(degenes)),names(degenes))]
  table(degenes)
  # GO and kegg analysis
  # (1) fitting the probability weighting function (PWF)
    pwf=nullp(degenes,genome,'ensGene')
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
    functional_analysis=goseq(pwf,genome,'ensGene',test.cats=c("GO:BP","GO:MF","KEGG"))
    restemp<-addKeggTogoseq(functional_analysis,temp)    # switch Keggpathid to name
    write.table(restemp,file="GO_Kegg_Wallenius.txt",row.names=F,sep="\t")


sessionInfo()
