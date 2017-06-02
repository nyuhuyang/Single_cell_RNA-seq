#hemberg-lab.github.io/scRNA.seq.course/index.html
#Identification of important genes

#check package
list.of.packages <- "devtools"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("https://bioconductor.org/biocLite.R")
list.of.packages <- c("M3Drop","monocle","DDRTree","pheatmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

library("devtools")
library("DDRTree")
library("pheatmap")
library("M3Drop")
library("monocle")
#### A-1) Setup enviroment and read data
#detect OS and set enviroment
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Single_Cell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Single_Cell");getwd();list.files()}

counts.merged<-read.csv("counts.merged.txt",sep="\t",header =T, row.names = 1)
dim(counts.merged)

#### A-2) Fitting the models
#First we must normalize & QC the dataset. 
#M3Drop contains a built-in function for this which removes cells
#with few detected genes, removes undetected genes, 
#and converts raw counts to CPM.
svg(filename="Renat.Single_Cell.fiter2.svg")
uso_list <- M3Drop::M3DropCleanData(
        counts.merged,
        labels = colnames(counts.merged),
        min_detected_genes = 2000,
        is.counts = TRUE
)
models <- M3Drop::M3DropDropoutModels(uso_list$data)
title(main = "Renat.Single_Cell")
dev.off()

#### A-3) Right outliers
DE_genes <- M3Drop::M3DropDifferentialExpression(
        uso_list$data,
        mt_method = "fdr",
        mt_threshold = 0.1
) 
title(main = "Renat.Single_Cell")

#### A-4) Validation of DE results
#We can also plot the expression levels of these genes
#to check they really are DE genes.
M3Drop::M3DropExpressionHeatmap(
        DE_genes$Gene,
        uso_list$data,
        cell_labels = uso_list$labels
#        key_genes = uso_markers
)
svg(filename="Renat.Single_Cell.svg")
Brennecke_HVG <- M3Drop::BrenneckeGetVariableGenes(
        uso_list$data,
        fdr = 0.01,
        minBiolDisp = 0.5
)


dev.off()
#### B-1) Pseudotime analysis of DE results 
#m3dGenes <- as.character(
#        M3Drop::M3DropDifferentialExpression(counts.merged),
# Selet gene list
m3dGenes<-rownames(DE_genes) # length(rownames(DE_genes))=38
m3dGenes<-Brennecke_HVG # length(Brennecke_HVG) = 307
m3dGenes<-rownames(uso_list$data) # length(rownames(uso_list$data)) = 7942
m3dGenes<-rownames(counts.merged) # length(rownames(counts.merged)) = 26364
write.csv(Brennecke_HVG,"GeneList_307.csv")
write.csv(rownames(uso_list$data),"GeneList_7942.csv")
write.csv(rownames(counts.merged),"GeneList_26364.csv")





#Now run monocle:
d <- counts.merged[which(rownames(counts.merged) %in% m3dGenes), ]
d <- as.matrix(d[!duplicated(rownames(d)), ])

pd <- data.frame(timepoint = colnames(d))
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = m3dGenes)
fd <- new("AnnotatedDataFrame", data=fd)       
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)


dCellData <- newCellDataSet(d,phenoData = pd, featureData = fd, expressionFamily = tobit())
dCellData <- setOrderingFilter(dCellData, which(geneNames %in% m3dGenes))
dCellData <- estimateSizeFactors(dCellData)
dCellDataSet <- reduceDimension(dCellData, pseudo_expr = 1) # takes long time
dCellDataSet <- orderCells(dCellDataSet, reverse = TRUE)
svg(filename="plot_cell_trajectory.svg")
plot_cell_trajectory(dCellDataSet)
dev.off()
