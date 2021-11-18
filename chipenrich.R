### R script to run chipenrich
### Input:
#1. peak file (PePr output)
#2. genome type (eg. "hg19", "hg38")
#3. chipenrich method (eg. "chipenrich","broadenrich")
#4. locus definition (eg. "nearest_tss")
#5. number of cores
#6. output file prefix
#7. output directory


myargs <- commandArgs(TRUE)

if(length(myargs) != 7){
    print("usage: peakFile, genome, method, locusDef, threads, outputPrefix, outputDir")
    q()
    }

peaks <- read.table(myargs[1], header=F,sep="\t",stringsAsFactors=F)

colnames(peaks) <- c("chrom","start","end","score")

library(chipenrich)

result_chipenrich.list <- chipenrich(peaks, genesets=c("GOBP","kegg_pathway"),
genome = myargs[2], method= myargs[3], locusdef=myargs[4], qc_plots=F,
n_cores=myargs[5], out_name= myargs[6], out_path=myargs[7])
