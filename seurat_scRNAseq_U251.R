
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(MAST)
library(cowplot)
library(patchwork)

### 1) read in count data: input 10X matrix file, genes file, and barcodes file (aggregated data)

aggr.data <- Read10X(data.dir = "filtered_feature_bc_matrix")

atrx <- CreateSeuratObject(counts = aggr.data, project = "U251_ATRX", min.cells = 3, min.features = 200)

### 2) data preprocessing (rename Idents using sample/group information)

### create sample/group labels

sampleIDs = as.character(read.csv("../cellranger_aggr/aggr.csv")$library_id)

groupLabels <- paste(sampleIDs,c("ATRX_WT","ATRX_KO"), sep=".")

### rename “orig.ident” object using sample/group labels

idents = names(Idents(atrx))

idents <- groupLabels[as.numeric(gsub(".*\\-","",idents))]

Idents(atrx) <- as.factor(idents)

atrx$orig.ident = as.factor(idents)

### 3) data filtering

### tag mitochondrial genes

atrx[["percent.mt"]] <- PercentageFeatureSet(object = atrx, pattern = "^MT-")

### Visualize QC metrics on pre-filtering data

pdf("QCplots_preFiltering.pdf")

VlnPlot(atrx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

dev.off()

### identify outliars (boxplot upper whisker: Q3 + 1.5 *IQR) of number of genes (nFeatures)

nFeature_Outlier <- quantile(atrx$nFeature_RNA,0.75) + 1.5*IQR(atrx$nFeature_RNA)

### remove cells with too many genes and over 25% MT genes

atrx <- subset(x = atrx, subset = nFeature_RNA > 0 & nFeature_RNA < nFeature_Outlier & percent.mt < 25)

### Visualize QC metrics on post-filtering data

pdf("QCplots_postFiltering.pdf")

VlnPlot(atrx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

dev.off()

### 4) batch correction (between different samples/groups)

### split suerat object by conditions (orig.ident)

atrx.list <- SplitObject(atrx, split.by = "orig.ident")

### normalization and find the highly variable features (n=2000, by default)

atrx.list <- lapply(atrx.list, function(x) {
    x = NormalizeData(x)
    FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })


### visualize the top10 highly variable genes of the two separated samples

top10VariableGenes.list <- lapply(atrx.list, function(x){head(VariableFeatures(x),10)})

pdf("top10VariableGenes.pdf")

plots.list <- lapply(atrx.list, function(x){VariableFeaturePlot(x)})
mapply(function(x,y,z){
LabelPoints(plot = x, points = y, repel = TRUE) + labs(title=z)
}, plots.list, top10VariableGenes.list, names(plots.list), SIMPLIFY=F)

dev.off()

### 5) data scaling and visualization by UMAP
### i) Scaling
all.genes <- rownames(atrx)

atrx.integrated.fullFeatures <- IntegrateData(anchorset = atrx.anchors, dims = 1:50, features.to.integrate = all.genes)

DefaultAssay(atrx.integrated.fullFeatures) <- "integrated"

atrx.integrated.fullFeatures <- ScaleData(atrx.integrated.fullFeatures, verbose = TRUE, features=all.genes)

### ii) perform PCA (PCs will be used for downstream tSNE or UMAP clustering)

atrx.integrated.fullFeatures <- RunPCA(atrx.integrated.fullFeatures, npcs = 50, verbose = FALSE)

pdf("PCAvisualization.pdf")

VizDimLoadings(atrx.integrated.fullFeatures, dims = 1:2, reduction = "pca")
DimPlot(atrx.integrated.fullFeatures, reduction = "pca")
DimHeatmap(atrx.integrated.fullFeatures, dims = 1:15, cells = 500, balanced = TRUE)

dev.off()

### iii) Determine the ‘dimensionality’ of the dataset
atrx.integrated.fullFeatures <- JackStraw(atrx.integrated.fullFeatures, dims=50, num.replicate = 100)
atrx.integrated.fullFeatures <- ScoreJackStraw(atrx.integrated.fullFeatures, dims = 1:50)

pdf("JackStraw_Elbow_Plots.pdf")

JackStrawPlot(atrx.integrated.fullFeatures, dims = 1:50)

ElbowPlot(atrx.integrated.fullFeatures, ndims=50)

dev.off()

### identify the optimal number of PCs from Elbow Plot result (< 0.5% decrease in standard deviation as cutoff)

PC.Stdev <- Stdev(object = atrx.integrated.fullFeatures, reduction = "pca")

index <- c()
for (i in 1:(length(PC.Stdev)-1)){
  d1 <- PC.Stdev[i]
  d2 <- PC.Stdev[i+1]
  if (d2 > 0.995*d1) {
    index <- append(index, i+1)
  }
}
pc.dim <- min(index)

### iv) cell clustering

atrx.integrated.fullFeatures <- FindNeighbors(atrx.integrated.fullFeatures, dims = 1:pc.dim)

atrx.integrated.fullFeatures <- FindClusters(atrx.integrated.fullFeatures, resolution = 1)

### v) non-linear dimensional reduction (UMAP/tSNE)

atrx.integrated.fullFeatures <- RunUMAP(atrx.integrated.fullFeatures, reduction = "pca", dims = 1:pc.dim)

### plot UMAP

pdf("UMAP_24PCs.pdf", width=10)

Idents(atrx.integrated.fullFeatures) <- atrx.integrated.fullFeatures$seurat_clusters

DimPlot(atrx.integrated.fullFeatures, reduction = "umap", label = TRUE, pt.size = 0.5, split.by="orig.ident", label.size=6) + NoLegend()

dev.off()

### 6) finding differentially expressed features
cluster_markers <- FindAllMarkers(object = atrx.integrated.fullFeatures, min.pct = 0, logfc.threshold = 0, return.thresh=1 , test.use = "MAST")

### select the significant markers (p_val_adj < 0.05, min.pct >= 0.25, abs(logfc) >= log(1.5))
### logfc is in natral log (ln)
### avg_logFC: Positive values indicate that the feature is more highly expressed in the first group.

cluster_markers_sigUp <- subset(cluster_markers, p_val_adj < 0.05 & avg_logFC  > log(1.5) & (pct.1 >= 0.25 | pct.2 >= 0.25))

### get the frequency of clusters that a gene is present

cluster_markers_sigUp <- cluster_markers_sigUp %>% group_by(gene) %>% mutate(gene.clusterFreq=n_distinct(cluster)) %>% ungroup %>% select(cluster,gene,gene.clusterFreq,1:5) %>% data.frame

### cluster marker plot

### top1 marker genes

pdf("featurePlots_sigUpClusterMarkers_top1.pdf", width=10)

top1_sigUpMarkers <- (cluster_markers_sigUp %>% group_by(cluster) %>% top_n(n=1, wt = avg_logFC) %>% data.frame)$gene

FeaturePlot(atrx.integrated.fullFeatures, features = top1_sigUpMarkers, min.cutoff = "q9", pt.size=0.5, ncol=4)

dev.off()

### top1 marker unique genes (only present in one cluster)

pdf("featurePlots_sigUpClusterMarkers_top1_uniqueCluster.pdf", width=10)

top1_sigUpMarkers_uniqueCluster <- (cluster_markers_sigUp %>% filter(gene.clusterFreq==1) %>% group_by(cluster) %>% top_n(n=1, wt = avg_logFC) %>% data.frame)$gene

FeaturePlot(atrx.integrated.fullFeatures, features = top1_sigUpMarkers_uniqueCluster, min.cutoff = "q9", pt.size=0.5, ncol=4)

dev.off()

### 7) finding DE genes between groups for each individual cluster

### split seurat object by clusters

clusters.list <- SplitObject(atrx.integrated.fullFeatures, split.by = "seurat_clusters")

### run DE gene analysis w/o cutoffs

DEgenes_ATRX_KOvWT_perClusters.list <- lapply(clusters.list, function(x){
Idents(x) = as.factor(x$orig.ident)
if(length(table(Idents(x))) == 2 & all(table(Idents(x)) > 3)){
FindMarkers(x, ident.1 = "981-RR-2.ATRX_KO", ident.2 = "981-RR-1.ATRX_WT", min.pct = 0, logfc.threshold = 0, test.use = "MAST")
}
})

### select the significant DE genes (p_val_adj < 0.05, min.pct >= 0.25, abs(logfc) >= log(1.5))

sigDEgenes_ATRX_KOvWT_clusters.df <- lapply(DEgenes_ATRX_KOvWT_perClusters.list, function(x){
	if(!is.null(x)){
subset(x, p_val_adj < 0.05 & abs(avg_logFC) > log(1.5) & (pct.1 >= 0.25 | pct.2 >= 0.25)) %>% tibble::rownames_to_column("gene")
		}
}) %>% bind_rows(.id= "cluster")
