#Identification of important genes
########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages <- c("VGAM", "irlba", "matrixStats", "igraph", 
                           "combinat", "fastICA","grid", "gridExtra",
                           "reshape2", "plyr", "parallel", "methods")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("DDRTree", "monocle")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

libraries(list.of.bio.packages,list.of.cran.packages)
memory.limit(size=65000)


# #####################################################################
# 
#  2.a Load and normalize data
#  https://www.bioconductor.org/packages/devel/bioc/vignettes/monocle/inst/doc/monocle-vignette.pdf
#
# ####################################################################
#====== A-1) Setup enviroment and read data (Required)============
#detect OS and set enviroment
# load the pipeline data by specifying a pipestance path
if (Sys.info()[['sysname']]=="Darwin"){
        WD <- "/Users/yah2014/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        WD <- "C:/Users/User/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}

#count_matrix<-read.csv(paste0(WD,"/Cell_KB-1/outs/sample123.csv"),row.names = 1)
#dim(count_matrix)

#2.1=======create a CellDataSet object (Required) from CellRanger GeneBCMatrix====
#http://cole-trapnell-lab.github.io/monocle-release/docs/
gbm_cds <- newCellDataSet(exprs(gbm),
                          phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                          featureData = new("AnnotatedDataFrame", data = fData(gbm)),
                          lowerDetectionLimit=0.5,
                          expressionFamily=negbinomial.size())
colnames(fData(gbm_cds))[2] <-"gene_short_name"

#2.2=================Filtering low-quality cells(Recommended)==================
gbm_cds <- detectGenes(gbm_cds, min_expr = 0.1)
print(head(fData(gbm_cds)))
print(head(pData(gbm_cds)))
gbm_cds<-estimateSizeFactors(gbm_cds, locfunc=genefilter::shorth ) # if Size_Factor=NA, there are too many zeroes, 
print(head(pData(gbm_cds)))
expressed_genes <- row.names(subset(fData(gbm_cds), num_cells_expressed >= 10))
length(expressed_genes)
gbm_cds<-gbm_cds[expressed_genes,]
#

#If you are using RPC values to measure expresion, 
#as we are in this vignette, it's also good to look at the distribution
#of mRNA totals across the cells:
pData(gbm_cds)$Total_mRNAs <- Matrix::colSums(exprs(gbm_cds))
gbm_cds <- gbm_cds[,pData(gbm_cds)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(gbm_cds)$Total_mRNAs)) + 2*sd(log10(pData(gbm_cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(gbm_cds)$Total_mRNAs)) - 2*sd(log10(pData(gbm_cds)$Total_mRNAs)))
gbm_cds <- gbm_cds[,pData(gbm_cds)$Total_mRNAs > lower_bound &
                     pData(gbm_cds)$Total_mRNAs < upper_bound]
qplot(Total_mRNAs, data=pData(gbm_cds), geom="density") +
        geom_vline(xintercept=lower_bound) +
        geom_vline(xintercept=upper_bound)

# 2.3==============Quanlity control (Recommended)====================================
# Log-transform each value in the expression matrix.
L <- log(exprs(gbm_cds)+1)
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily"
melted_dens_df <- melt(t(scale(t(L))))
# Plot the distribution of the standardized gene expression values.
qplot(value, geom="density", data=melted_dens_df) + stat_function(fun = dnorm, size=0.5, color='red') +
        xlab("Standardized log(FPKM)") +
        ylab("Density") #Slow.. 
#===================================================================================
# #####################################################################
# 
#  3. Classifying and counting cells of different types
#  https://www.bioconductor.org/packages/devel/bioc/vignettes/monocle/inst/doc/monocle-vignette.pdf
#
# ####################################################################

#
#3.1 Classifying cells by type (Recommended)=============================================
id_ <-function(gene){
    #return gene id
    return(row.names(subset(fData(gbm_cds), gene_short_name == gene)))
}
cth <- newCellTypeHierarchy()
#Epithelial
#---(Recommended)KRT19 - all (majority) of epithelial; + negative for non-epithelial (DCN, PECAM1, LAPTM5, CD68)----------
cth <- addCellType(cth, "Epithelial Cells", classify_func=function(x) {
    x[id_("KRT19"),] >= 1 &  x[id_("DCN"),] == 0 & x[id_("LAPTM5"),] == 0 & x[id_("CD68"),] == 0 })


#Non-epithelial structural
#---(Recommended)Stromal+endothelial+fibroblasts - GPX3 + negative for leukocyte markers (PTPRC, LAPTM5, SRGN)------------
cth <- addCellType(cth, "Stromal Cells+\nEndothelial Cells+\nFibroblasts", classify_func=function(x) {
        x[id_("GPX3"),] >= 1 & x[id_("PTPRC"),] == 0 & x[id_("LAPTM5"),] == 0 & x[id_("SRGN"),] == 0})

#---Stromal/fibroblasts - DCN, COL6A1, TIMP3, PDGFRA-------
cth <- addCellType(cth, "Stromal+\nfibroblasts Cells", classify_func=function(x) {
    x[id_("DCN"),] >= 1 & x[id_("COL6A1"),] >= 1 & x[id_("TIMP3"),] >= 1 & x[id_("PDGFRA"),
        parent_cell_type_name = "Stromal Cells+\nEndothelial Cells+\nFibroblasts"] >= 1}) #Stromal/fibroblasts
#---(Recommended)fibroblasts - ANPEP-------
cth <- addCellType(cth, "Fibroblasts", classify_func=function(x) {
        x[id_("ANPEP"),] >= 1 },
        parent_cell_type_name = "Stromal Cells+\nEndothelial Cells+\nFibroblasts") #fibroblasts

#---Stromal/ - MMRN2,CD248------- #Expression of stromal cell markers in distinct compartments of human skin cancers
cth <- addCellType(cth, "Stromal+ Cells", classify_func=function(x) {
        x[id_("MMRN2"),] >= 1 & x[id_("CD248"),] >= 1},
        parent_cell_type_name = "Stromal Cells+\nEndothelial Cells+\nFibroblasts") #Stromal

#----Endothelial - VWF, EMCN, PECAM1, CDH5-------
cth <- addCellType(cth, "Endothelial Cells", classify_func=function(x) {
  x[id_("VWF"),] >= 1 & x[id_("EMCN"),] >= 1 & x[id_("CDH5"),] >= 1},
        parent_cell_type_name = "Stromal+\nEndothelial Cells") #Endothelial

#Immune http://www.abcam.com/primary-antibodies/immune-cell-markers-poster
#---(Recommended)Leukocytes (all)  - PTPRC, LAPTM5, SRGN----------
cth <- addCellType(cth, "Leukocytes", classify_func=function(x) {
        (x[id_("PTPRC"),] >= 1 | x[id_("LAPTM5"),] >= 1 | x[id_("SRGN"),] >= 1)}) #Leukocytes (all)
                
#---Macrophages - CD68, MARCO, LYZ-------
cth <- addCellType(cth, "Macrophages", classify_func=function(x) {
        (x[id_("CD68"),] >= 1 & x[id_("MARCO"),] >= 1 & x[id_("LYZ"),] >= 1) },
        parent_cell_type_name = "Leukocytes") #Macrophages
#---(Recommended)T cells - CD3D-------
cth <- addCellType(cth, "T cells", classify_func=function(x) {
    x[id_("CD3D"),] >= 1 },
    parent_cell_type_name = "Leukocytes") #T cells
#---(Recommended)Monocytes - CD14-------
cth <- addCellType(cth, "Monocytes", classify_func=function(x) {
        x[id_("CD14"),] >= 1 },
        parent_cell_type_name = "Leukocytes") #Monocytes
#---Adipocytes---------
cth <- addCellType(cth, "Beige Adipocytes", classify_func=function(x) {
        x[id_("TNFRSF9"),] >= 1 & x[id_("TMEM2"),] >= 1}) #Beige Adipocytes
#===========================================================================

HSMM <- classifyCells(gbm_cds, cth, 0.1)
#Error in if (type_res[cell_name] == TRUE) next_nodes <- c(next_nodes,  : 
#missing value where TRUE/FALSE needed
#check id_("") exist or not
table(pData(HSMM)$CellType)

# generate Pie

pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)

pie + coord_polar(theta = "y") +theme(axis.title.x=element_blank(), 
                                      axis.title.y=element_blank(),
                                      legend.key.size = unit(1.5, "cm"),
                                      text = element_text(size=30),
                                      legend.text = element_text(size=20))

#generate table
cellType.table <- read.delim("CellType2.txt") #previously as
ss <- tableGrob(cellType.table) #library(gridExtra)
grid.arrange(ss)

#3.2 Unsupervised cell clustering(Alternative)========================
#Clustering cells without marker genes

#We can filter genes based on average expression level, and we can additionally select genes
#that are unusually variable across cells. These genes tend to be highly informative about cell state.
#estimateDispersions only makes sense right now for negbinomial and negbinomial.size expression families.
HSMM<-estimateDispersions(HSMM)
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)
#Now we're ready to try clustering the cells:
HSMM_1 <- reduceDimension(HSMM, max_components=2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM_1 <- clusterCells(HSMM_1, num_clusters=9)  # increase num_clusters for cell trajectories analysis
plot_cell_clusters(HSMM_1, 1, 2, color="CellType") +
        theme(text = element_text(size=20),
              legend.justification = 'right', 
              legend.position=c(0.20,0.78))+ 
        guides(colour = guide_legend(override.aes = list(size=9)))
#here can test sets of markers
plot_cell_clusters(HSMM_1, 1,2, color="CellType", markers=c("CD14","ANPEP","CD3G","KRT19"))+
        theme(text = element_text(size=30),
              legend.position="right",
              legend.key.width = unit(2, "cm"))+ 
        guides(colour = guide_legend(override.aes = list(size=10)))
#Monocle allows us to subtract the effects of "uninteresting" sources of variation
#to reduce their impact on the clustering.
HSMM_1 <- clusterCells(HSMM_1, num_clusters=2)
#All_Cell_Cluster 
plot_cell_clusters(HSMM_1, 1, 2, color="CellType") +
        theme(text = element_text(size=40),
              legend.position="none")+ 
        guides(colour = guide_legend(override.aes = list(size=10)))+
        facet_wrap(~Cluster)



#3.3 Clustering cells using marker genes(Recommended)========================
#Semi-supervised cell clustering with known marker genes
#Before we just picked genes that were highly expressed and highly variable.
#Now, we'll pick genes that co-vary with our markers.
#In a sense, we will be building a large list of genes to use as markers,
#so that even if a cell does not have KRT19, it might be recognizable as a Epithelial cells based on other genes.

#It then removes all the "Unknown"and "Ambiguous" functions before identifying genes that are differentially expressed between the types.

marker_diff <- markerDiffTable(HSMM_1[expressed_genes,], 
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
HSMM_2 <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM_2)
HSMM_2 <- reduceDimension(HSMM_2, max_components=2, num_dim = 6, reduction_method = 'tSNE',
                        residualModelFormulaStr="~num_genes_expressed", verbose = T)
HSMM_2 <- clusterCells(HSMM_2, num_clusters=10) # increase num_clusters from 3 to 10 for cell trajectories analysis
plot_cell_clusters(HSMM_2, 1, 2, color="CellType")


#3.4 Imputing cell type
#we've reduce the number of "contaminating" Cells.
#what about the "Unknown" cells?
#When a cluster is composed of more than a certain percentage (in this case, 10%) of a certain type,
#all the cells in the cluster are set to that type.
#If a cluster is composed of more than one cell type, the whole thing is marked "Ambiguous".
#If there's no cell type thats above the threshold, the cluster is marked "Unknown"
suppressWarnings(HSMM_2 <- clusterCells(HSMM_1,
                     num_clusters=10,
                     frequency_thresh=0.4, #Original frequency_thresh=0.1 {0.4,0.52}
                     cell_type_hierarchy=cth))
table(pData(HSMM_2)$CellType)

#To control the order and color, change the factor levels to (Epi,Endo, Immune)
#???pData(HSMM_2)$CellType <-factor(pData(HSMM_2)$CellType,levels=c("Fibroblasts",
#                                                                "Stromal Cells+\nEndothelial Cells",
#                                                                "Myeloid cells",
#                                                                "Monocytes"))

#Pie

pie <- ggplot(pData(HSMM_2), aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)

pie + coord_polar(theta = "y") +theme(axis.title.x=element_blank(), 
                                      axis.title.y=element_blank(),
                                      legend.key.size = unit(2, "cm"),
                                      text = element_text(size=40),
                                      legend.text = element_text(size=30))


#All_Cell_Cluster: 
plot_cell_clusters(HSMM_2, 1, 2, color="CellType") +
    theme(text = element_text(size=40),
          legend.justification = 'right', 
          legend.position=c(0.28,0.86),
          legend.key.height = grid::unit(0.8, "in"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))

#Each_Cell_Cluster
plot_cell_clusters(HSMM_2, 1, 2, color="CellType") +facet_wrap(~CellType)+ 
    theme(text = element_text(size=40),
          legend.justification = 'right', 
          legend.position=c(0.35,0.2),
          legend.key.width = grid::unit(0.8, "in"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))

#Marker genes level in two Cell types 
plot_cell_clusters(HSMM_2, 1,2, color="CellType", markers=c("SRGN","DCN"))+
    theme(text = element_text(size=35),
          legend.text = element_text(size=30))+
    guides(colour = guide_legend(override.aes = list(size=10)))



#Finally, we subset the CellDataSet object to create HSMM_immue, 
#which includes only myoblasts. We'll use this in the
HSMM_epi <- HSMM_2[,pData(HSMM_2)$CellType == "KRT19+ epithelial Cells"]
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
HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
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
saveRDS(HSMM_epi_3,file="HSMM_epi_3") #save HSMM_epi_3 for other applications, like shiny
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
