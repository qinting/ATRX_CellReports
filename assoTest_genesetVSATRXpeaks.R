### a function to perform fisher's exact test on the association between ATRX peak and gene set of interest
### 2X2 table: rows are genes associated with ATRX peak and non-ATRX peaks (random regions); columns are the genes in the geneset of interest and genes not in the geneset
### input:
# 1. gene set to be test (a data.frame of genes and pathway/GO "ATRXpeakAnnot_nonDDRgenes_cellReportSup")
# 2. ATRX peaks annotation matrix with "annot.symbol" column
# 3. random region (background) annotaiton matrix with "annot.symbol" column
# 4. matched gene symbol column names (geneColName.geneset, geneColName.peakAnnot, geneColName.bgAnnot)
# 5. genomic regions of interest (promoter, enhancer, CpG, etc.), i.e. restrict the peaks on the regions
### output:
# 1. the 2x2 table
# 2. Fisher's Exact test result
# 3. geneset genes linked to peaks (geneset.peakAnnot)
# 4. geneset genes linked to random regions (geneset.bgAnnot)

genePeakAssoFET <- function(geneset, peakAnnot, bgAnnot, geneColName.geneset, geneColName.peakAnnot, geneColName.bgAnnot, genomicRegion){
  require(dplyr)
  ### get the number of peak-linked genes in geneset and not in geneset
  tmp <- left_join(peakAnnot,geneset,by=setNames(nm=geneColName.peakAnnot,geneColName.geneset))
  geneset.peakAnnot <- tmp[!is.na(tmp$Annot),]
  nonGeneset.peakAnnot <- tmp[is.na(tmp$Annot),]
  genesetNum.peakAnnot <- data.frame(count(group_by(unique(geneset.peakAnnot[,c("annot.symbol","annot.type")]), annot.type),"annot.symbol"))[,-2]
  nonGenesetNum.peakAnnot <- data.frame(count(group_by(unique(nonGeneset.peakAnnot[,c("annot.symbol","annot.type")]), annot.type),"annot.symbol"))[,-2]

  ### get the number of random region (background)-linked genes in geneset and not in geneset
  tmp <- left_join(bgAnnot,geneset,by=setNames(nm=geneColName.peakAnnot,geneColName.geneset))
  geneset.bgAnnot <- tmp[!is.na(tmp$Annot),]
  nonGeneset.bgAnnot <- tmp[is.na(tmp$Annot),]
  genesetNum.bgAnnot <- data.frame(count(group_by(unique(geneset.bgAnnot[,c("annot.symbol","annot.type")]), annot.type),"annot.symbol"))[,-2]
  nonGenesetNum.bgAnnot <- data.frame(count(group_by(unique(nonGeneset.bgAnnot[,c("annot.symbol","annot.type")]), annot.type),"annot.symbol"))[,-2]

  ### generate 2X2 contingent table for gene numbers at the genomic region of interest
  ### merge gene number data.frames
  geneNum.df <- list(genesetNum.peakAnnot,nonGenesetNum.peakAnnot,genesetNum.bgAnnot,nonGenesetNum.bgAnnot) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="annot.type"), .)
  colnames(geneNum.df)[2:5] <- c("genesetNum.peakAnnot","nonGenesetNum.peakAnnot","genesetNum.bgAnnot","nonGenesetNum.bgAnnot")
  geneNum.df[is.na(geneNum.df)] <- 0
  contingentTable <- matrix(as.numeric(geneNum.df[geneNum.df$annot.type==genomicRegion,-1]),2,2,byrow=T, dimnames=list(c("geneAssoPeak","geneAssoBG"),c("geneset","nonGeneset")))

  ### perform FET
  fet.out <- fisher.test(contingentTable)
  return(list("contingentTable"=contingentTable,"fet.output"=fet.out,"geneset.peakAnnot"=geneset.peakAnnot,"geneset.bgAnnot"=geneset.bgAnnot))
}
