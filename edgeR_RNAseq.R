### edgeR DEG analysis

HTSeqCount <- read.delim("./data/HTSeqcounts_RNAseq_Huse.txt",header=T,row.names=1)

### sample informations

sampleInfo <- read.table("./data/sampleInfo.txt",header=T,sep="\t",stringsAsFactors=F)

library(edgeR)

targets <- sampleInfo[match(colnames(HTSeqCount),sampleInfo[,1]),]

Group <- factor(targets$Label)

design <- model.matrix(~0+Group)

colnames(design) <- levels(Group)
rownames(design) <- targets[,1]


## creat Gene DGEList.obj

DGEList.obj <- DGEList(counts=HTSeqCount[,targets$Run],group=Group)

## filtering: keep genes with at least1 cpm in > 3 samples across all groups

isexpr <- rowSums(cpm(DGEList.obj)>=1) >=3

DGEList.obj <- DGEList.obj[isexpr,]

## re-compute library size

DGEList.obj$samples$lib.size <- colSums(DGEList.obj$counts)

## sample normalization

DGEList.obj <- calcNormFactors(DGEList.obj)

## dispersion estimation (common dispersion and tagwise dispersion in one step)

y <- estimateDisp(DGEList.obj,design, tagwise=TRUE, robust=TRUE)

## fit genewise glms

fit <- glmQLFit(y,design)

## DE test by LRT

my.contrasts <- makeContrasts(ATRXpos_vs_ATRXneg = TP53negATRXpos - TP53negATRXneg, levels=design)

qlf <- glmQLFTest(fit,contrast=my.contrasts)

qlf.ranked <- topTags(qlf, n=nrow(qlf$table))

DEGenes <- subset(qlf.ranked$table, FDR < 0.05 & abs(logFC) >= 1)

### calculate log2(RPKM)

exon <- read.table("./data/mm9_exon_merged.bed",header=F,sep="\t",stringsAsFactors=F)

library(splitstackshape)

library(dplyr)

exon_splitGene <- unique(data.frame(cSplit(exon, splitCols="V4", sep=",", direction="long"),stringsAsFactors=FALSE))

exon_splitGene$exon.length <- exon_splitGene$V3-exon_splitGene$V2+1

geneExonLength <- data.frame(summarise(group_by(exon_splitGene,by=V4),gene.length=sum(exon.length)))

colnames(geneExonLength)[1] <- "gene.symbol"

y_RPKM <- rpkm(y, gene.length=geneExonLength[match(rownames(y$counts), geneExonLength[,1]),2], normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)

### visualization ###

library(gplots)

pdf("./heatmap_DEGenes.pdf")

Data <- y_RPKM[rownames(DEGenes),]
Data <- Data[,order(sampleInfo[match(colnames(Data),sampleInfo[,1]),"Label"])]

heatmap.2(Data,col=greenred(256),key=T,symkey=F,density.info="none",trace="none",dendrogram="row",Colv=FALSE,revC=T,labRow=NA,labCol=colnames(Data),cexRow=1,cexCol=1, srtCol=45, keysize=1,scale="row",margins=c(5,10), ColSideColors= c(rep("dark blue",3),rep("dark red",3)))

legend("topright",legend=c("TP53negATRXneg","TP53negATRXpos"),fill=c("dark blue","dark red"),cex=0.8,bty="n")

dev.off()

### volcano plots ###

pdf("./volcanoPlot.pdf", colormodel="cmyk", useDingbats=FALSE)

par(mai=rep(1.5,4))
a <- qlf.ranked$table[,"logFC"]
b <- -log10(qlf.ranked$table[,"PValue"])

smoothScatter(a,b,pch=20,cex=0.5,xlab="log2(fold change)", ylab="-log10(p value)", cex.lab=2.5,cex.axis=2,font=6, xlim=c(-ceiling(max(a)), ceiling(max(a))),mgp = c(3.5, 1, 0))

points(DEGenes[,"logFC"],-log10(DEGenes[,"PValue"]),col="red",pch=20,cex=0.5)

abline(h= sort(-log10(DEGenes[,"PValue"]))[1],col="lightgrey")
abline(v=1,col="lightgrey")
abline(v=-1,col="lightgrey")

dev.off()

### unsupervised analysis ###

## jittered boxplot of Atrx

rpkm_ATRX.matrix$RunID <- sampleInfo[match(rownames(rpkm_ATRX.matrix),sampleInfo[,1]),1]

pdf("jitterplot_rpkm_Atrx.pdf")
ggplot(rpkm_ATRX.matrix,aes(x=Label,y=rpkm,color=RunID)) + geom_boxplot(color="black")+ geom_jitter(position=position_jitter(0.2))
dev.off()
