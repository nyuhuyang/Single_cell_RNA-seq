#Identification of important genes
list.of.packages <- c("devtools","dplyr","pheatmap","VGAM", "irlba",
                      "matrixStats", "igraph", "combinat", "fastICA",
                      "grid", "reshape2", "plyr", "parallel", "methods")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("https://bioconductor.org/biocLite.R")
list.of.packages <- c("GenomeInfoDb","DESeq2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

#check package
library("DESeq2")
#1. Setup enviroment and read data
#detect OS and set enviroment
if (Sys.info()[['sysname']]=="Darwin"){
    setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Renat_Single_Cell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
    setwd("C:/Users/User/Dropbox/Public/Olivier/R/Renat_Single_Cell");getwd();list.files()}


#7. Compare the gene expression between cluster using DESeq
# read gene list and HSMM S4 data
HSMM_ordering_genes <-read.csv("HSMM_ordering_genes.csv")
HSMM_ordering_genes<-as.character(unlist(HSMM_ordering_genes))
HSMM_epi_3<-readRDS("HSMM_epi_3")

# create sample table for DESeq
sample_epi.table <- pData(HSMM_epi_3)
head(sample_epi.table) 
# Creat count matrix
count_matrix<-read.csv("counts.merged.txt",sep="\t",header =T, row.names = 1)
count_matrix_epi <-count_matrix[,(colnames(count_matrix)  %in% rownames(sample_epi.table))]

all(rownames(sample_epi.table)  == colnames(count_matrix_epi)) #Test the consistant of epi sample name

## Creating a DESeqDataSet object
#create ddsHTSeq
ddsHTSeq <- DESeqDataSetFromMatrix( countData = count_matrix_epi,
                                    colData = sample_epi.table,
                                    design= ~ Cluster)
ddsHTSeq
### Normalization for sequencing depth
####Pre-filtering and switch factor
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) >1, ] #removing rows with 1 reads


cts <- counts(ddsHTSeq)
geoMeans <- apply(cts, 1, function(row) exp(sum(log(row[row != 0]))/length(row)))
dds <- estimateSizeFactors(ddsHTSeq, geoMeans=geoMeans)
#use above command if below error occur
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

#run DESeq
dds$Cluster <- factor(dds$Cluster, levels=c("1","2","3"))
#ddsHTSeq$Cluster <- droplevels(ddsHTSeq$Cluster)
dds <- DESeq(dds)
saveRDS(dds,file="dds") #save dds for other applications, like shiny

#Extract results
#For testing at a different threshold, we provide the `alpha` to *results*,
#so that the mean filtering is optimal for our new FDR threshold.

res_cluster_2vs1 <- results(dds,contrast = c("Cluster","2","1"), alpha=0.1, lfcThreshold=1)
res_cluster_3vs2 <- results(dds,contrast = c("Cluster","3","2"), alpha=0.1, lfcThreshold=1) # takes long time, don't know why
res_cluster_3vs1 <- results(dds,contrast = c("Cluster","3","1"), alpha=0.1, lfcThreshold=1)


# Examining results tables and summary
table(res_cluster_2vs1$padj < 0.1)
table(res_cluster_3vs2$padj < 0.1)
table(res_cluster_3vs1$padj < 0.1)
summary(res_cluster_2vs1)
summary(res_cluster_3vs2)
summary(res_cluster_3vs1)


#Exporting results to CSV files

res_cluster_2vs1 <- res_cluster_2vs1[order(res_cluster_2vs1$padj),]
res_cluster_3vs2 <- res_cluster_3vs2[order(res_cluster_3vs2$padj),]
res_cluster_3vs1 <- res_cluster_3vs1[order(res_cluster_3vs1$padj),]
res_cluster_2vs1 <- res_cluster_2vs1[complete.cases(res_cluster_2vs1),] #remove NA
res_cluster_3vs2 <- res_cluster_3vs2[complete.cases(res_cluster_3vs2),] #remove NA
res_cluster_3vs1 <- res_cluster_3vs1[complete.cases(res_cluster_3vs1),] #remove NA

write.csv(as.data.frame(res_cluster_2vs1),file="DESeq_result_cluster_2vs1.csv")
write.csv(as.data.frame(res_cluster_3vs2),file="DESeq_result_cluster_3vs2.csv")
write.csv(as.data.frame(res_cluster_3vs1),file="DESeq_result_cluster_3vs1.csv")

#RUN GESA on https://david.ncifcrf.gov/ by uploading g ene names


### Visualizing results

#The MA-plot provides a global view of the differential genes, 
#with the log2 fold change on the y-axis over the mean of normalized counts:
par(mfrow=c(1,3),cex=1.5) 
plotMA(res_cluster_2vs1, ylim=c(-5,5),main ="cluster2 genes/ cluster1 genes")
plotMA(res_cluster_3vs2, ylim=c(-5,5),main ="cluster3 genes/ cluster2 genes",yaxt='n',ylab="")
plotMA(res_cluster_3vs1, ylim=c(-5,5),main ="cluster3 genes/ cluster1 genes",yaxt='n',ylab="")

#Examine the counts for the top gene, sorting by p-value:
res_cluster_2vs1 <-read.csv("DESeq_result_cluster_2vs1.csv",row.names = 1)
res_cluster_3vs2 <-read.csv("DESeq_result_cluster_3vs2.csv",row.names = 1)
res_cluster_3vs1 <-read.csv("DESeq_result_cluster_3vs1.csv",row.names = 1)
dds<-readRDS("dds")
par(mfrow=c(1,1))
par(oma=c(2,2,2,2))
d<-data.frame(x= str(0), y= integer(0)) # create a empty data.frame

for( i in 1:9){
    d1 <- plotCounts(dds, gene=rownames(res_cluster_2vs1)[i], intgroup="Cluster",returnData=T)
    d1$gene_short_name<-rep(rownames(res_cluster_2vs1)[i],nrow(d1)) # add one more column for facet_warp
    d <- rbind(d,d1)
}
ggplot(d, aes(x=Cluster, y=count)) + 
        facet_wrap(~gene_short_name)+
        geom_point(position=position_jitter(w=0.2,h=0),size=3)+
        ggtitle("Finding genes that disinguish Clusters")+
        theme(plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
              axis.text=element_text(size=20),
              strip.text = element_text(size=25),
              axis.title=element_text(size=40, colour = "black"))+
        coord_trans(y = "log10")

d <- plotCounts(dds, gene=which(rownames(dds)=="SCGB1A1"), intgroup="Cluster",
                returnData=T)
ggplot(d, aes(x=Cluster, y=count)) +
    geom_point(position=position_jitter(w=0.2,h=0),size=3)+
    labs(title = (rownames(dds)[which(rownames(dds)=="SCGB1A1")]))+
    theme(axis.text=element_text(size=40),
          axis.title=element_text(size=40, colour = "black"),
          plot.title = element_text(hjust = 0.5,size=40,face="bold"))


#Make normalized counts plots for the top 9 genes:
par(mfrow=c(3,3))
for (i in 1:9)  {
    plotCounts(dds, gene=i, intgroup="Cluster")
}
    ggplot(d, aes(x=Cluster, y=count)) +
        geom_point(position=position_jitter(w=0.2,h=0),size=3)+
        labs(title = (rownames(dds)[order(res_cluster_2vs1$padj)[i]]))+
        theme(axis.text=element_text(size=40),
              axis.title=element_text(size=40, colour = "black"),
              plot.title = element_text(hjust = 0.5,size=40,face="bold"))
}
#A more sophisticated plot of counts:

library(ggplot2)
data <- plotCounts(dds, gene=which(rownames(res_cancer)=="SPI1"), intgroup=c("Condition","SPI1"), returnData=TRUE)
ggplot(data, aes(x=Condition, y=count, col=SPI1,size = SPI1))+
    geom_point(position=position_jitter(width=.1,height=0))+
    labs(title = "SPI1 expression level in Healthy vs Cancer group")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=40)) + 
    scale_y_log10()

ggplot(data, aes(x=SPI1, y=count, col=Condition,size = SPI1))+
    geom_point(position=position_jitter(width=.1,height=0))+
    labs(title = "SPI1 expression level in unmutated vs mutated group")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),text = element_text(size=40)) + 
    scale_y_log10()

#A sorted results table:

res_cancer_Sort <- res_cancer[order(res_cancer$padj),]
head(res_cancer_Sort)

res_SPI1_Sort <- res_SPI1[order(res_SPI1$padj),]
head(res_SPI1_Sort)
#A heatmap of the top genes:

library(pheatmap)

topgenes <- head(rownames(res_cancer_Sort),100)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Condition","SPI1")])
pheatmap(mat, annotation_col=df, 
         fontsize=15,fontsize_row=7,fontsize_col=13,cex=1.05,
         main ="Top 100 differentially expressed genes between Healthy and Cancer")

topgenes <- head(rownames(res_SPI1_Sort),100)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("Condition","SPI1")])
pheatmap(mat, annotation_col=df, 
         fontsize=15,fontsize_row=7,fontsize_col=13,cex=1.05,
         main ="Top 100 differentially expressed genes between mutated SPI1 and unmutated SPI1")
