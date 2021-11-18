### run GSEA for differential genes between clusters or between groups

library(fgsea)

### prepare pre-ranked gene list for GSEA test

### convert gene.symbol to Entrez.ID

library(org.Hs.eg.db)

tmp <- unlist(as.list(org.Hs.egSYMBOL))

EntrezID2Symbol <- data.frame(GeneID=as.integer(names(tmp)),GeneSymbol=tmp,stringsAsFactors=F)

rankedGenes_ATRX_KOvWT_perClusters.list <- lapply(DEgenes_ATRX_KOvWT_perClusters.list, function(x){if(!is.null(x)){
	y = data.frame(x) %>% tibble::rownames_to_column("gene.symbol") %>%
inner_join(EntrezID2Symbol,by=c("gene.symbol"="GeneSymbol"));
ranks = y$avg_logFC;
names(ranks) = y$GeneID;
ranks}
})

### exclude null elements

rankedGenes_ATRX_KOvWT_perClusters.list[unlist(lapply(rankedGenes_ATRX_KOvWT_perClusters.list, is.null))] <- NULL

### get GOBP gene set

library(chipenrich)
library(chipenrich.data)

KEGG <- as.list(geneset.kegg_pathway.hsa@set.gene)

GOBP <- as.list(geneset.GOBP.hsa@set.gene)

### Run GSEA

### DE genes between groups

GSEArs_ATRX_KOvWT_perClusters.list <- lapply(rankedGenes_ATRX_KOvWT_perClusters.list, function(x){
fgsea(pathways = c(KEGG, GOBP), stats = x, minSize=15, maxSize=500, nperm=1000)
})

### DE genes between clusters

GSEArs_clusterMarkers.list <- sapply(as.character(unique(cluster_markers$cluster)), function(x){
y = subset(cluster_markers, cluster==as.integer(x)) %>% inner_join(EntrezID2Symbol, by=c("gene"="GeneSymbol"));
	ranks = y$avg_logFC;
	names(ranks) = y$GeneID;
fgsea(pathways = c(KEGG, GOBP), stats = ranks, minSize=15, maxSize=500, nperm=1000)
}, simplify=F)

### combine results

GSEArs_ATRX_KOvWT_perClusters.df <- bind_rows(GSEArs_ATRX_KOvWT_perClusters.list, .id="cluster") %>% data.frame

GSEArs_clusterMarkers.df <- bind_rows(GSEArs_clusterMarkers.list, .id="cluster") %>% data.frame

### add gene set name

genesetID2Name_KEGG_GOBP <- do.call(rbind,c(as.list(geneset.kegg_pathway.hsa@set.name),as.list(geneset.GOBP.hsa@set.name))) %>% data.frame %>% tibble::rownames_to_column("genesetID") %>% setNames(c("genesetID","geneset.description"))

GSEArs_ATRX_KOvWT_perClusters.df <- left_join(GSEArs_ATRX_KOvWT_perClusters.df, genesetID2Name_KEGG_GOBP, by=c("pathway"="genesetID")) %>% select(1:2,geneset.description,3:9)

GSEArs_clusterMarkers.df <- left_join(GSEArs_clusterMarkers.df, genesetID2Name_KEGG_GOBP, by=c("pathway"="genesetID")) %>% select(1:2,geneset.description,3:9)

### extract sig. results

GSEArs_ATRX_KOvWT_perClusters_padj05.df <- subset(GSEArs_ATRX_KOvWT_perClusters.df, padj < 0.05) %>% group_by(cluster) %>% arrange(padj, desc(abs(NES))) %>% data.frame

GSEArs_clusterMarkers_padj05.df <- subset(GSEArs_clusterMarkers.df, padj < 0.05) %>% group_by(cluster) %>% arrange(padj, desc(abs(NES))) %>% data.frame

### GSEA VlnPlot
pdf("GSEA_enrichPlot.pdf", colormodel="cmyk", useDingbats=FALSE, height=3)

plotEnrichment(GOBP[["GO:0044839"]], rankedGenes_ATRX_KOvWT_perClusters.list[["7"]]) + labs(title="GOBP: cell cycle G2/M phase transition (ATRX KO vs WT in cluster 7)")

plotEnrichment(GOBP[["GO:0044843"]], rankedGenes_ATRX_KOvWT_perClusters.list[["11"]]) + labs(title="GOBP: cell cycle G1/S phase transition (ATRX KO vs WT in cluster 11)")

dev.off()

#####################################################################
### cell cycle scoring

### S phase genes

s.genes <- Seurat::cc.genes$s.genes

### genes in dataset (using integrated data with all genes)

s.genes <- s.genes[s.genes %in% rownames(atrx.integrated.fullFeatures)]

### G2M phase genes

g2m.genes <- Seurat::cc.genes$g2m.genes

### genes in dataset

g2m.genes <- g2m.genes[g2m.genes %in% rownames(atrx.integrated.fullFeatures)]

### calcualte cell cycle score
atrx.integrated.fullFeatures <- CellCycleScoring(object = atrx.integrated.fullFeatures, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

### PCA on cell cycle genes
atrx.integrated.fullFeatures_PCAonCCphases <- RunPCA(atrx.integrated.fullFeatures, features = c(s.genes, g2m.genes))

pdf("PCAplot_byCellCyclePhase.pdf", width=10)

PCAPlot(atrx.integrated.fullFeatures_PCAonCCphases)

PCAPlot(atrx.integrated.fullFeatures_PCAonCCphases, split.by="orig.ident")

dev.off()

### UMAP plot with cell cycle phase annotation

pdf("UMAP_24PCs_cellCyclePhaseAnnot.pdf", width=15)

### all cells together, labeled by clusters

DimPlot(atrx.integrated.fullFeatures, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

### separated groups, labeled by clusters

DimPlot(atrx.integrated.fullFeatures, reduction = "umap", label = TRUE, pt.size = 0.5, split.by="orig.ident")

dev.off()


### summarize the number of cells at each cell cycle phase in each cluster

metaData <- atrx.integrated.fullFeatures@meta.data

cellNum_cellcyclePhase <- metaData %>% group_by(orig.ident, seurat_clusters, Phase) %>% summarise(cellNum.perPhase = n()) %>% data.frame %>% spread(orig.ident, cellNum.perPhase) %>% setNames(c("cluster","Phase","cellNum.ATRXwt","cellNum.ATRXko"))

cellNum_cellcyclePhase %>% mutate(cellPerc.ATRXwt = cellNum.ATRXwt/sum(cellNum.ATRXwt,na.rm=T)*100, cellPerc.ATRXko = cellNum.ATRXko/sum(cellNum.ATRXko,na.rm=T)*100)

cellNum_cellcyclePhase <- cellNum_cellcyclePhase %>% mutate(cellPerc.ATRXwt = cellNum.ATRXwt/sum(cellNum.ATRXwt,na.rm=T)*100, cellPerc.ATRXko = cellNum.ATRXko/sum(cellNum.ATRXko,na.rm=T)*100, cellPerc.diff.KOvWT = cellPerc.ATRXko - cellPerc.ATRXwt)

### bar plots

library(ggplot2)

cellPerc_cellcyclePhase.long <- gather(cellNum_cellcyclePhase[,c(1,2,5,6)], "Type", "cell.perc", 3:4)

cellPerc_cellcyclePhase.long$Type <- ifelse(cellPerc_cellcyclePhase.long$Type=="cellPerc.ATRXwt","ATRX.WT","ATRX.KO")

cellPerc_cellcyclePhase.long$Type <- factor(cellPerc_cellcyclePhase.long$Type,levels=c("ATRX.WT","ATRX.KO"))

cellPerc_cellcyclePhase.long$Phase <- factor(cellPerc_cellcyclePhase.long$Phase,levels=c("G1","S","G2M"))

pdf("percCells_perCellCylePhase_perCluster.pdf", colormodel="cmyk", useDingbats=FALSE, height=5)

ggplot(data= cellPerc_cellcyclePhase.long, aes(x=Type,y=cell.perc, fill=Phase)) + facet_wrap(~cluster, ncol=13) + geom_bar(stat="identity", color="black", position="stack") + scale_fill_manual(values=c("grey","dark blue","dark red")) + theme_bw() + labs(title="Percent of cells at each phase in each cluster", x="Group", y="Percent of cells (%)") + theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, colour="black"), panel.border = element_rect(colour = "grey", fill=NA))

dev.off()

### refine cycling phase annotation

ccScores.df <-  atrx.integrated.fullFeatures@meta.data[,c("S.Score","G2M.Score","Phase","orig.ident")]

ccScores.df$cellID <- rownames(ccScores.df)

### re-define using cycling vs non-cycling cutoff of log(1) (FC = 1)

ccScores.df$Phase_FC1 <- ifelse(ccScores.df$Phase == "G1", "non-cycling", ifelse(ccScores.df$Phase == "S","G1S", as.character(ccScores.df$Phase)))

colnames(ccScores.df)[1] <- "G1S.Score"

ccScores.df <- mutate(ccScores.df, avgCC.Score=log((exp(G1S.Score)+exp(G2M.Score))/2))

metaData <- left_join(Data, ccScores.df)

pdf("/ckoschma/Manuscripts/ATRX/Fig2_UMAP.pdf",height=5,width=10)

### all cells

Idents(atrx.integrated.fullFeatures) <- atrx.integrated.fullFeatures$seurat_clusters
DimPlot(atrx.integrated.fullFeatures, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "orig.ident") + NoLegend() + labs(title="UMAP of all cells")

### cycling cutoff: FC = 1 (ie. phse.score = log(1))

Idents(subset(atrx.integrated.fullFeatures, subset = Phase_FC1!= "non-cycling")) <- atrx.integrated.fullFeatures$seurat_clusters
DimPlot(atrx.integrated.fullFeatures, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "orig.ident") + NoLegend() + labs(title="UMAP of cycling cells (FC > 1)")

ggplot(data= subset(metaData, Phase_FC1!= "non-cycling"), aes(x=UMAP_1, y=UMAP_2, color=avgCC.Score)) + geom_point(size=0.5, alpha=0.5) + scale_colour_gradient(low="blue", high="red") + facet_wrap(~orig.ident) + labs(title="UMAP of cycling cells (FC > 1)") + theme_bw() + theme(plot.title=element_text(size=12,face="bold"), axis.text=element_text(size=16), axis.title=element_text(size=14), strip.text.x=element_text(size=18,face= "bold"), legend.text=element_text(size=16))

### cycling cutoff: FC = 1, annotated by phase

metaData$Phase <- factor(metaData$Phase,levels=c("G1","S","G2M"))

ggplot(data= subset(metaData, Phase_FC1!= "non-cycling"), aes(x=UMAP_1, y=UMAP_2, color=Phase)) + geom_point(size=0.5,alpha=0.5) + scale_colour_manual(values=c("dark blue","dark red")) + facet_wrap(~orig.ident) + labs(title="UMAP of cycling cells (FC > 1)") + theme_bw() + theme(plot.title=element_text(size=12,face="bold"), axis.text=element_text(size=16), axis.title=element_text(size=14), strip.text.x=element_text(size=18,face= "bold"), legend.text=element_text(size=16))

dev.off()

### piechart plot of cycling cell proportions at each condition (WT vs KO)

library(scales)
cellNum_cellcyclePhase_noG1 <- subset(metaData, Phase_FC1!="non-cycling") %>% group_by(orig.ident, seurat_clusters, Phase_FC1) %>% summarise(cellNum.perPhase = n_distinct(cellID)) %>% data.frame %>% spread(orig.ident, cellNum.perPhase) %>% setNames(c("cluster","Phase","cellNum.ATRXwt","cellNum.ATRXko"))

cellNum_cellcyclePhase_noG1 <- cellNum_cellcyclePhase_noG1 %>% mutate(cellPerc.ATRXwt = cellNum.ATRXwt/sum(cellNum.ATRXwt,na.rm=T)*100, cellPerc.ATRXko = cellNum.ATRXko/sum(cellNum.ATRXko,na.rm=T)*100, cellPerc.diff.KOvWT = cellPerc.ATRXko - cellPerc.ATRXwt)

cellPerc_cellcyclePhase_noG1.long <- gather(cellNum_cellcyclePhase_noG1[,c(1,2,5,6)], "Type", "cell.perc", 3:4)

Data <- cellPerc_cellcyclePhase_noG1.long %>% group_by(Type,cluster) %>% summarise(cell.perc=sum(cell.perc,na.rm=T)) %>% data.frame

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

pdf("piechart_cyclingCellProp _perCluster_noG1.pdf", colormodel="cmyk", useDingbats=FALSE, height=5)

ggplot(Data, aes(x="", y=cell.perc, fill=cluster)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + geom_text(aes(label =  paste0(round(cell.perc,0), "%")), position = position_stack(vjust = 0.5)) + facet_wrap(~Type) + scale_fill_manual(values=c("darkolivegreen2","light green","light blue","#55DDE0", "cornflowerblue","darkturquoise","#33658A", "#F6AE2D", "magenta","deeppink2","darkred","#F26419", "#999999")) + blank_theme +   theme(axis.text.x=element_blank())

dev.off()
