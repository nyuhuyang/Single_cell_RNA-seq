#hemberg-lab.github.io/scRNA.seq.course/index.html
#Identification of important genes

#check package
list.of.packages <- c("devtools","reshape")
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
library("reshape")
#### A-1) Setup enviroment and read data
#detect OS and set enviroment
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Single_Cell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Single_Cell");getwd();list.files()}

count_matrix<-read.csv("counts.merged.txt",sep="\t",header =T, row.names = 1)
dim(count_matrix)

#2.1 create a CellDataSet object
d <- as.matrix(count_matrix)

pd <- data.frame(timepoint = colnames(d))
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = rownames(d))
fd <- new("AnnotatedDataFrame", data=fd)       
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)


HSMM <- newCellDataSet(d,phenoData = pd,
                       featureData = fd,
                       expressionFamily = tobit())

#2.5 Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
print(head(pData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
#length(expressed_genes)=8542


#If you are using RPC values to measure expresion, 
#as we are in this vignette, it's also good to look at the distribution
#of mRNA totals across the cells:
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
qplot(Total_mRNAs, data=pData(HSMM), geom="density") +
        geom_vline(xintercept=lower_bound) +
        geom_vline(xintercept=upper_bound)
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
                     pData(HSMM)$Total_mRNAs < upper_bound]
HSMM <- detectGenes(HSMM, min_expr = 0.1)
# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
# Plot the distribution of the standardized gene expression values.
qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') +
        xlab("Standardized log(FPKM)") +
        ylab("Density")

#'from' cannot be NA, NaN or infinite 

















#### A-2) Fitting the models
#First we must normalize & QC the dataset. 
#M3Drop contains a built-in function for this which removes cells
#with few detected genes, removes undetected genes, 
#and converts raw counts to CPM.
#svg(filename="Renat.Single_Cell.fiter2.svg")
uso_list <- M3Drop::M3DropCleanData(
        count_matrix,
        labels = colnames(count_matrix),
        min_detected_genes = 2000,
        is.counts = TRUE
)
models <- M3Drop::M3DropDropoutModels(uso_list$data)
title(main = "Renat.Single_Cell")
dev.off()

#### B-1) Pseudotime analysis of DE results 
#m3dGenes <- as.character(
#        M3Drop::M3DropDifferentialExpression(count_matrix),
# Selet gene list

m3dGenes<-rownames(uso_list$data) # length(rownames(uso_list$data)) = 7942

write.csv(rownames(uso_list$data),"GeneList_7942.csv")




#Now run monocle:
d <- count_matrix[which(rownames(count_matrix) %in% m3dGenes), ]
d <- as.matrix(d[!duplicated(rownames(d)), ])

pd <- data.frame(timepoint = colnames(d))
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = m3dGenes)
fd <- new("AnnotatedDataFrame", data=fd)       
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)


HSMM <- newCellDataSet(d,phenoData = pd, featureData = fd, expressionFamily = tobit())
HSMM <- setOrderingFilter(HSMM, which(geneNames %in% m3dGenes))
HSMM <- estimateSizeFactors(HSMM)
HSMM <- reduceDimension(HSMM, pseudo_expr = 1) # takes long time
HSMM <- orderCells(HSMM, reverse = TRUE)
#svg(filename="plot_cell_trajectory.svg")
plot_cell_trajectory(HSMM)
dev.off()

#"State" is just Monocle's term for the segment of the tree.
#The function below is handy for identifying the State
#which contains most of the cells from time zero. 
#We can then pass that to orderCells:
GM_state <- function(cds){
        if (length(unique(pData(cds)$State)) > 1){
                T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
                return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
        }else { return (1) }
}
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by="Pseudotime")





#So the root should have cells that expresshigh levels of proliferation markers.
#We can use the jitter plot to pick figure out which state corresponds to rapid 
#proliferation:
        
blast_genes <- row.names(subset(fData(HSMM), gene_short_name %in%  c("CXCL2", "NFKBIA", "FTH1")))
plot_genes_jitter(HSMM[blast_genes,], grouping="State", min_expr=0.1)




#To confirm that the ordering is correct we can select a couple of markers
#of myogenic progress. Plotting these genes demonstrates that ordering looks good:
HSMM<-detectGenes(HSMM,min_expr=10)
HSMM_expressed_genes<-rownames(subset(fData(HSMM), 
                                              num_cells_expressed >= 10))
#length(HSMM_expressed_genes)=333
HSMM_filtered<- HSMM[HSMM_expressed_genes,]

my_genes<-rownames(subset(fData(HSMM_filtered),gene_short_name %in% 
                         c("CXCL2", "NFKBIA", "FTH1")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="State")


#5.2 Finding genes that distinguish cell type or state


        