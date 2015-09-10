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
subjunc(index=index,readfile1=read1,readfile2=read2,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,nthreads=10,tieBreakHamming=TRUE,unique=TRUE,indels=5)

}else {
subjunc(index=index,readfile1=read1,input_format="gzFASTQ",output_format="BAM",output_file=targets$OutputFile,nthreads=10,tieBreakHamming=TRUE,unique=TRUE,indels=5)
}

alignstats<-propmapped(targets$OutputFile,countFragments=TRUE,properlyPaired=FALSE)
write.table(alignstats,file="AlignmentSummary.txt",sep="\t")

pdf(file="plots.pdf")


# count numbers of reads mapped to NCBI Refseq genes
fc <-featureCounts(files=targets$OutputFile,annot.inbuilt=unlist(strsplit(args[7],"="))[2],nthreads=8,strandSpecific=unlist(strsplit(args[8],"="))[2],isPairedEnd=unlist(strsplit(args[12],"="))[2])

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
fc <-featureCounts(files=targets$OutputFile,annot.inbuilt=unlist(strsplit(args[7],"="))[2],nthreads=8,strandSpecific=unlist(strsplit(args[8],"="))[2],isPairedEnd=unlist(strsplit(args[12],"="))[2],GTF.featureType="exon", GTF.attrType="ID",
useMetaFeatures=FALSE, allowMultiOverlap=TRUE)

dge <- DGEList(counts=fc$counts, genes=fc$annotation)

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
  anno_version=unlist(strsplit(args[7],"="))[2];
  
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
    write.table(restemp,file="GO_Kegg_Wallenius.txt",row.names=F,sep="\t")


sessionInfo()
