### R script to run annotatr
### Input:
#1. peak file (PePr output)
#2. genome type (eg. "hg19", "hg38")
#3. plot title
#4. output prefix

myargs <- commandArgs(TRUE)

if(length(myargs) != 4){
    print("usage: peakFile, genome, plotTitle, outPrefix")
    q()
    }


library(annotatr)

### read in data and annotation reference

# extraCols = c(signalValue = 'numeric', pValue = 'numeric', qValue = 'numeric')

# peaks_GR <- read_regions(con=myargs[1], genome=myargs[2], extraCols=extraCols, format="bed")

peaks_GR <- read_regions(con=myargs[1], genome=myargs[2], format="bed")

annots <- builtin_annotations()[grep(myargs[2],builtin_annotations())]

annotations = build_annotations(genome = myargs[2], annotations = annots)

### annotate peaks

peaks_annotated <- annotate_regions(regions=peaks_GR, annotations=annotations,ignore.strand=TRUE,quiet=F)

# randomize regions and do annotation

peaks_random <- randomize_regions(regions=peaks_GR, allow.overlaps=TRUE, per.chromosome=TRUE)

peaks_random_annotated <- annotate_regions(regions=peaks_random,annotations=annotations, ignore.strand=TRUE, quiet=F)

### plot of annotation count with random regions

library(ggplot2)

p <- plot_annotation(annotated_regions = peaks_annotated, annotated_random = peaks_random_annotated, annotation_order = annots, plot_title = myargs[3], x_label = 'Known Gene Annotations',y_label = 'Count')

pdf(paste(myargs[4], "knownGeneRegionAnnot.pdf", sep="_"))
p + theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
dev.off()

### plot regions occurring in pairs of annoations

p <- plot_coannotations(annotated_regions= peaks_annotated, annotation_order=annots, axes_label='Known Gene Annotations',plot_title = myargs[3])

pdf(paste(myargs[4], "knownGeneRegionAnnot_coannot.pdf", sep="_"))
p
dev.off()
