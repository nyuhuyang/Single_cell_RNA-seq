#Identification of important genes
list.of.packages <- c("devtools","dplyr","pheatmap","VGAM", "irlba",
                      "matrixStats", "igraph", "combinat", "fastICA",
                      "grid", "reshape2", "plyr", "parallel", "methods")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#check package
library("devtools")
library("DDRTree")
library("pheatmap")
library("M3Drop")
library("monocle")
library("reshape")

#### A-1) Setup enviroment and read data
#detect OS and set enviroment
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Renant_Single_Cell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Renant_Single_Cell");getwd();list.files()}

count_matrix<-read.csv("counts.merged.txt",sep="\t",header =T, row.names = 1)
dim(count_matrix)

#2.1 create a CellDataSet object
d <- as.matrix(count_matrix)

pd <- data.frame(timepoint = colnames(d))
rownames(pd)<-colnames(d)
fd <- data.frame(gene_short_name = rownames(d))
rownames(fd) <-rownames(d)
pd <- new("AnnotatedDataFrame", data=pd)
fd <- new("AnnotatedDataFrame", data=fd)       
geneNames <- rownames(d)


HSMM <- newCellDataSet(d,phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=negbinomial.size())

#2.5 Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
print(head(pData(HSMM)))
HSMM<-estimateSizeFactors(HSMM, locfunc=genefilter::shorth ) # if Size_Factor=NA, there are too many zeroes, 
print(head(pData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes)


#If you are using RPC values to measure expresion, 
#as we are in this vignette, it's also good to look at the distribution
#of mRNA totals across the cells:
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
                     pData(HSMM)$Total_mRNAs < upper_bound]
# Log-transform each value in the expression matrix.
#\L <- log(exprs(HSMM[expressed_genes,])+1)
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
#\melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
# Plot the distribution of the standardized gene expression values.
#\qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') +
#\        xlab("Standardized log(FPKM)") +
#\        ylab("Density")

#3 Classifying and counting cells of different types
#3.1 Classifying cells with CellTypeHierarchy

id_ <-function(gene){
    #return gene id
    return(row.names(subset(fData(HSMM), gene_short_name == gene)))
}
cth <- newCellTypeHierarchy()
#Epithelial
#KRT19 - all (majority) of epithelial; + negative for non-epithelial (DCN, PECAM1, LAPTM5, CD68)----------
cth <- addCellType(cth, "KRT19+ epithelial Cells", classify_func=function(x) {
    x[id_("KRT19"),] >= 1 &
    x[id_("DCN"),] == 0 & x[id_("PECAM1"),] == 0 & x[id_("LAPTM5"),] == 0 & x[id_("CD68"),] == 0 })
#Non-epithelial structural
#Stromal+endothelial - GPX3 + negative for leukocyte markers (PTPRC, LAPTM5, SRGN)------------
#Stromal/fibroblasts - DCN, COL6A1, TIMP3, PDGFRA
#Endothelial - VWF, EMCN, PECAM1, CDH5
cth <- addCellType(cth, "Stromal+\nEndothelial+\nFibroblasts Cells", classify_func=function(x) {
    x[id_("GPX3"),] >= 1 & x[id_("PTPRC"),] == 0 & x[id_("LAPTM5"),] == 0 & x[id_("SRGN"),] == 0  | #Stromal+endothelial
    x[id_("DCN"),] >= 1 & x[id_("COL6A1"),] >= 1 & x[id_("TIMP3"),] >= 1 & x[id_("PDGFRA"),] >= 1}) #Stromal/fibroblasts
 #subscript out of bounds
#  x[id_("VWF"),] >= 1 & x[id_("EMCN"),] >= 1 & x[id_("PECAM1"),] >= 1 & x[id_("CDH5"),] >= 1   }) #Endothelial


#Immune
#Leukocytes (all)  - PTPRC, LAPTM5, SRGN----------
#Macrophages - CD68, MARCO, LYZ-------
#T cells - CD3G
cth <- addCellType(cth, "Immune Cells", classify_func=function(x) {
    (x[id_("PTPRC"),] >= 1 | x[id_("LAPTM5"),] >= 1 | x[id_("SRGN"),] >= 1)| #Leukocytes (all)
    (x[id_("CD68"),] >= 1 & x[id_("MARCO"),] >= 1 & x[id_("LYZ"),] >= 1) }) #Macrophages
#subscript out of bounds
#    x[id_("CD3G"),] >= 1   })           #T cells

suppressWarnings(HSMM <- classifyCells(HSMM, cth, 0.1))

table(pData(HSMM)$CellType)

#3.2 Unsupervised cell clustering
#We can filter genes based on average expression level, and we can additionally select genes
#that are unusually variable across cells. These genes tend to be highly informative about cell state.
#estimateDispersions only makes sense right now for negbinomial and negbinomial.size expression families.
HSMM<-estimateDispersions(HSMM)
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
#3.3 Semi-supervised cell clustering with known marker genes
#Before we just picked genes that were highly expressed and highly variable.
#Now, we'll pick genes that co-vary with our markers.
#In a sense, we will be building a large list of genes to use as markers,
#so that even if a cell does not have KRT19, it might be recognizable as a Epithelial cells based on other genes.

#It then removes all the "Unknown"and "Ambiguous" functions before identifying genes that are differentially expressed between the types.

marker_diff <- markerDiffTable(HSMM[expressed_genes,], 
                               cth,
                               residualModelFormulaStr="~num_genes_expressed", # obtain 0 genes using ~CellType, there is no cluster yet
                               cores=1) #Long cth will trigure "subscript out of bounds". It takes long time.

#The function markerDiffTable then returns a data frame of test results, 
#and you can use this to pick the genes you want to use for clustering.
#Often it's best to pick the top 10 or 20 genes that are most specific for each cell type.
#This ensures that the clustering genes aren't dominated by markers for one cell type.
#You generally want a balanced panel of markers for each type if possible.
#Monocle provides a handy function for ranking genes by how restricted their expression is for each type.

candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)

#Below shows top three marker genes for each cell type
selectTopMarkers(marker_spec, 3)

#To cluster the cells, we'll choose the top 500 markers for each of these cell types:
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM_3 <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM_3)
HSMM_3 <- reduceDimension(HSMM_3, max_components=2, num_dim = 6, reduction_method = 'tSNE',
                        residualModelFormulaStr="~num_genes_expressed", verbose = T)
HSMM_3 <- clusterCells(HSMM_3, num_clusters=10) # increase num_clusters from 3 to 10 for cell trajectories analysis
plot_cell_clusters(HSMM_3, 1, 2, color="CellType")


#3.4 Imputing cell type
#we've reduce the number of "contaminating" Endo Cells in the Epithelial cluster, and vice versa.
#what about the "Unknown" cells?
#When a cluster is composed of more than a certain percentage (in this case, 10%) of a certain type,
#all the cells in the cluster are set to that type.
#If a cluster is composed of more than one cell type, the whole thing is marked "Ambiguous".
#If there's no cell type thats above the threshold, the cluster is marked "Unknown"
suppressWarnings(HSMM_3 <- clusterCells(HSMM_3,
                     num_clusters=10,
                     frequency_thresh=0.45, #Original frequency_thresh=0.1 {0.4,0.52}
                     cell_type_hierarchy=cth))
table(pData(HSMM_3)$CellType)  

#To control the order and color, change the factor levels to (Epi,Endo, Immune)
pData(HSMM_3)$CellType <-factor(pData(HSMM_3)$CellType,levels=c("KRT19+ epithelial Cells",
                                                                "Stromal+\nEndothelial+\nFibroblasts Cells",
                                                                "Immune Cells"))

#Pie

pie <- ggplot(pData(HSMM_3), aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)

pie + coord_polar(theta = "y") +theme(axis.title.x=element_blank(), 
                                      axis.title.y=element_blank(),
                                      legend.key.size = unit(2, "cm"),
                                      text = element_text(size=40),
                                      legend.text = element_text(size=30))


#All_Cell_Cluster: 
plot_cell_clusters(HSMM_3, 1, 2, color="CellType") +
    theme(text = element_text(size=40),
          legend.justification = 'right', 
          legend.position=c(0.75,0.82),
          legend.key.height = grid::unit(1, "in"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))

#Each_Cell_Cluster
plot_cell_clusters(HSMM_3, 1, 2, color="CellType") +facet_wrap(~CellType)+ 
    theme(text = element_text(size=40),
          legend.justification = 'right', 
          legend.position=c(0.7,0.7),
          legend.key.width = grid::unit(1, "in"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))

#Marker genes level in three Cell types 
plot_cell_clusters(HSMM_3, 1,2, color="CellType", markers=c("SRGN","KRT19","GPX3"))+
    theme(text = element_text(size=35),
          legend.text = element_text(size=30))+
    guides(colour = guide_legend(override.aes = list(size=10)))

#Marker genes level in Epithelial Cells 
plot_cell_clusters(HSMM_3, 1,2, color="CellType", markers=c("KRT19","KRT5","MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1"))+
    theme(text = element_text(size=50),
          legend.justification = 'right', 
          legend.position=c(1.0,-0.1),
          legend.key.height = grid::unit(0.5, "in"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))


#Finally, we subset the CellDataSet object to create HSMM_epi, 
#which includes only myoblasts. We'll use this in the
HSMM_epi <- HSMM_3[,pData(HSMM_3)$CellType == "KRT19+ epithelial Cells"]
HSMM_epi <- estimateDispersions(HSMM_epi)

#4 Constructing single cell trajectories
#4.3 Unsupervised ordering
table(pData(HSMM_epi)$Cluster)
#what if we don't have time series data?
#Below are two methods to select genes that require no knowledge of the design of the experiment at all.

#4.3.1 Selecting genes with high dispersion across cells

#4.3.2 Selecting genes based on PCA loading

#4.4 Unsupervised feature selection based on density peak clustering
#To use dpFeature, we firrst select superset of feature genes as genes expressed in at least 5% of all the cells.
HSMM_epi_3 <- detectGenes(HSMM_epi, min_expr=0.1)
fData(HSMM_epi_3)$use_for_ordering <- fData(HSMM_epi_3)$num_cells_expressed > 0.05 * ncol(HSMM_epi_3)
plot_ordering_genes(HSMM_epi_3)

#We will then run reduceDimension with t-SNE as the reduction method on those top PCs and project them further
#down to two dimensions.
HSMM_epi_3 <- reduceDimension(HSMM_epi_3, max_components=2, norm_method = 'log', num_dim = 6,
                             reduction_method = 'tSNE', verbose = T)
#Then we can run density peak clustering to identify the clusters on the 2-D t-SNE space.
HSMM_epi_3 <- clusterCells(HSMM_epi_3, verbose = F)

#To control the order and color, change the cluster levels from 1,2,3 to 1,3,2
pData(HSMM_epi_3)$Cluster <-factor(pData(HSMM_epi_3)$Cluster,levels=c("1","3","2"))

#After the clustering, we can check the clustering results.

plot_cell_clusters(HSMM_epi_3, 1,2, color='Cluster', markers=c("KRT19","KRT5","MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1"))+
    theme(text = element_text(size=50),
          legend.justification = 'right', 
          legend.position=c(1.0,-0.15),
          legend.key.height = grid::unit(0.5, "in"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))

#We also provide the decision plot for users to check the p; b for each cell and decide the threshold for defining the cell clusters.
#skip plot_rho_delta(HSMM_epi_3, rho_threshold = 2, delta_threshold = 4 )+
#skip     theme(text = element_text(size=50))

#skip HSMM_epi_3 <- clusterCells(HSMM_epi_3,
#skip                          rho_threshold = 2,
#skip                          delta_threshold = 4,
#skip                          skip_rho_sigma = T,
#skip                          verbose = F)



#After we confirm the clustering makes sense, 
#we can then perform differential gene expression test as a way to extract the genes that distinguish them.
HSMM_expressed_genes <- row.names(subset(fData(HSMM_epi_3), num_cells_expressed >= 10))
clustering_DEG_genes <- differentialGeneTest(HSMM_epi_3[HSMM_expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 1) #takes long time
#We will then select the top 1000 significant genes as the ordering genes.
clustering_DEG_genes<-clustering_DEG_genes[order(clustering_DEG_genes$qval),]
write.csv(clustering_DEG_genes,"clustering_DEG_genes_cluster.csv")

HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
write.csv(HSMM_ordering_genes,"HSMM_ordering_genes.csv",row.names = F)
plot_ordering_genes(HSMM_epi_3)
HSMM_epi_3 <- setOrderingFilter(HSMM_epi_3, ordering_genes = HSMM_ordering_genes)
HSMM_epi_3 <- reduceDimension(HSMM_epi_3)
HSMM_epi_3 <- orderCells(HSMM_epi_3)
GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$Pseudotime)[,"0"] # replace $Hours with $Pseudotime
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    }else { return (1) }
}
HSMM_epi_3 <- orderCells(HSMM_epi_3, root_state=GM_state(HSMM_epi_3))
plot_cell_trajectory(HSMM_epi_3, color_by="State")+
    theme(text = element_text(size=40))+
    guides(colour = guide_legend(override.aes = list(size=10)))

plot_cell_trajectory(HSMM_epi_3, color_by="Cluster")+
    theme(text = element_text(size=40))+
    guides(colour = guide_legend(override.aes = list(size=10)))

plot_cell_trajectory(HSMM_epi_3, color_by="Pseudotime")+
    theme(text = element_text(size=40))+
    guides(colour = guide_legend(override.aes = list(size=10)))

plot_cell_trajectory(HSMM_epi_3, color_by="Cluster") + facet_wrap(~State, nrow=1)

#We can check the final state results as following:
plot_cell_clusters(HSMM_epi_3, 1,2, color="Cluster", markers=HSMM_ordering_genes[1:15])+
    theme(text = element_text(size=30),
          legend.key.height = grid::unit(0.75, "in"),
          legend.justification = 'right', 
          legend.position=c(1.0,-0.15))+
    guides(colour = guide_legend(override.aes = list(size=10)))


plot_genes_branched_pseudotime(HSMM_epi_3[HSMM_ordering_genes[1:15],],
#                               branch_point=1,
                               color_by="Cluster",
                               ncol=4)+
    theme(text = element_text(size=30),
          legend.key.height = grid::unit(0.75, "in"),
          legend.justification = 'right', 
          legend.position=c(0.9,-0.05))+
    guides(colour = guide_legend(override.aes = list(size=10)))


#4.5 Semi-supervised ordering with known marker genes
#5 Differential expression analysis
#5.1 Basic differential analysis
#5.2 Finding genes that distinguish cell type or state

#to_be_tested <- row.names(subset(fData(HSMM_epi_3),
#                                 gene_short_name %in% c("SCGB3A2","CXCL8","HSPA5",
#                                                        "SCGB3A1","FOS","RRAD","ANXA1",
#                                                        "EGR1","ID1","AHR","TP53BP2",
#                                                        "JUN","NEDD4L","STK17A","EMP2")))
#diff_test_res <- differentialGeneTest(HSMM_epi_3[to_be_tested,],
#                                      fullModelFormulaStr="~Cluster")
#skip diff_test_res <- subset(diff_test_res, qval < 0.1)
#diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_jitter(HSMM_epi_3[HSMM_ordering_genes[1:15],], grouping="Cluster", color_by="Cluster",
                  nrow=4, ncol=NULL, plot_trend=TRUE)+
    theme(text = element_text(size=30),
          legend.key.height = grid::unit(0.75, "in"),
          legend.justification = 'right', 
          legend.position=c(0.9,0.1))+
    guides(colour = guide_legend(override.aes = list(size=2)))



#5.3 Finding genes that change as a function of pseudotime
marker_genes <- row.names(subset(fData(HSMM_epi_3),
          gene_short_name %in% c("SCGB1A1","SCGB3A2","CXCL8","HSPA5",
                                            "SCGB3A1","FOS","RRAD","ANXA1",
                                         "EGR1","ID1","AHR","TP53BP2",
                                         "JUN","NEDD4L","STK17A","EMP2")))
#5.4 Clustering genes by pseudotemporal expression pattern
diff_test_res <- differentialGeneTest(HSMM_epi_3[marker_genes,],
                                      fullModelFormulaStr="~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1)) #Don't remove this
plot_pseudotime_heatmap(HSMM_epi_3[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap= T)
plot_genes_branched_pseudotime(HSMM_epi_3[Epi_subtype_gene,],
                               #                               branch_point=1,
                               color_by="Cluster",
                               ncol=2)+
    theme(text = element_text(size=30),legend.position="none")
#5.5 Multi-factorial differential expression analysis
#In below example, Monocle tests three genes for differential expression between
#epi end immune, while subtracting the effect of Hours

#6 Analyzing branches in single-cell trajectories

plot_cell_trajectory(HSMM_epi_3, color_by="Pseudotime")

#BEAM takes as input a CellDataSet that's been ordered with orderCells and the name of a branch point in
#the trajectory. It returns a table of significance scores for each gene. 
#Genes that score significant are said to be branch-dependent in their expression.

HSMM_epi_3_res <- BEAM(HSMM_epi_3, branch_point=1, cores = 1) #It takes long long time
BEAM_res<-HSMM_epi_3_res
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
#You can visualize changes for all the genes that are significantly branch dependent using a special type of heatmap.
plot_genes_branched_heatmap(HSMM_epi_3[row.names(subset(BEAM_res, qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 2,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
#We can plot a couple of these genes, such as Pdpn and Sftpb
HSMM_epi_genes <- row.names(subset(fData(HSMM_epi_3), gene_short_name %in% c("HSPA6", "BAG3")))
plot_genes_branched_pseudotime(HSMM_epi_3[HSMM_epi_genes,],
                               branch_point=1,
                               color_by="Pseudotime",
                               ncol=1)
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

rownames(sample_epi.table) == colnames(count_matrix_epi) #Test the consistant of epi sample name
## Creating a DESeqDataSet object
library(DESeq2)
#create ddsHTSeq
ddsHTSeq <- DESeqDataSetFromMatrix( countData = count_matrix_epi,
                                    colData = sample_epi.table,
                                    design= ~ Cluster)
ddsHTSeq
### Normalization for sequencing depth
####Pre-filtering and switch factor
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) >0, ] #removing rows with 0 reads
ddsHTSeq$Cluster <- factor(ddsHTSeq$Cluster, levels=c("1","2","3")) 

cts <- counts(ddsHTSeq)
geoMeans <- apply(cts, 1, function(row) exp(sum(log(row[row != 0]))/length(row)))
dds <- estimateSizeFactors(ddsHTSeq, geoMeans=geoMeans)
#use above command if below error occur
#Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
#every gene contains at least one zero, cannot compute log geometric means

#run DESeq

dds <- DESeq(dds)


#Extract results
#For testing at a different threshold, we provide the `alpha` to *results*,
#so that the mean filtering is optimal for our new FDR threshold.

res_cluster_2vs1 <- results(dds,contrast = c("Cluster","2","1"), alpha=0.5, lfcThreshold=1)
res_cluster_3vs1 <- results(dds,contrast = c("Cluster","3","1"), alpha=0.5, lfcThreshold=1)
res_cluster_3vs2 <- results(dds,contrast = c("Cluster","3","2"), alpha=0.5, lfcThreshold=1) # takes long time, don't know why
res2 <- results(dds, alpha=0.5,lfcThreshold=1)

# Examining results tables and summary
table(res_cluster_2vs1$padj < 0.1)
table(res_cluster_3vs1$padj < 0.1)
table(res_cluster_3vs2$padj < 0.1)
table(res2$padj < 0.1)
summary(res_cluster_2vs1)
summary(res_cluster_3vs1)
summary(res_cluster_3vs2)
summary(res2)


#Exporting results to CSV files

resO5rdered <- res2[order(res2$padj),]
resO5rdered <- resO5rdered[complete.cases(resO5rdered),] #remove NA
res05 <- resO5rdered[resO5rdered$padj<0.05,]
write.csv(as.data.frame(res05),file="DESeq_result_cluster123.csv")
#RUN GESA on https://david.ncifcrf.gov/ by uploading g ene names


### Visualizing results

#The MA-plot provides a global view of the differential genes, 
#with the log2 fold change on the y-axis over the mean of normalized counts:
par(mfrow=c(2,2))
plotMA(res_cluster_2vs1, ylim=c(-5,5))
plotMA(res_cluster_3vs1, ylim=c(-5,5))
plotMA(res_cluster_3vs2, ylim=c(-5,5))
plotMA(res2, ylim=c(-5,5))


#A p-value histogram:

hist(res_cluster_2vs1$pvalue[res_cluster_2vs1$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")


#Examine the counts for the top gene, sorting by p-value:

par(mfrow=c(1,1))
par(oma=c(2,2,2,2))
d <- plotCounts(dds, gene=which.min(res_cluster_2vs1$padj), intgroup="Cluster",
                returnData=T)
ggplot(d, aes(x=Cluster, y=count)) +
        geom_point(position=position_jitter(w=0.2,h=0),size=3)+
        labs(title = (rownames(dds)[gene]))+
        theme(axis.text=element_text(size=40),
              axis.title=element_text(size=40, colour = "black"),
              plot.title = element_text(hjust = 0.5,size=40,face="bold"))

d <- plotCounts(dds, gene=which(rownames(dds)=="FOS"), intgroup="Cluster",
                returnData=T)
ggplot(d, aes(x=Cluster, y=count)) +
    geom_point(position=position_jitter(w=0.2,h=0),size=3)+
    labs(title = (rownames(dds)[which(rownames(dds)=="FOS")]))+
    theme(axis.text=element_text(size=40),
          axis.title=element_text(size=40, colour = "black"),
          plot.title = element_text(hjust = 0.5,size=40,face="bold"))


#Make normalized counts plots for the top 9 genes:
par(mfrow=c(3,3))
for (i in 1:9)  {
    d <- plotCounts(dds, gene=order(res_cluster_1vs2$padj)[4], intgroup="Cluster",
                    returnData=T)
    ggplot(d, aes(x=Cluster, y=count)) +
        geom_point(position=position_jitter(w=0.2,h=0),size=3)+
        labs(title = (rownames(dds)[gene]))+
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