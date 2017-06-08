#hemberg-lab.github.io/scRNA.seq.course/index.html
#Identification of important genes

#check package
list.of.packages <- c("devtools","reshape","dplyr","pheatmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("https://bioconductor.org/biocLite.R")
list.of.packages <- c("M3Drop","monocle","DDRTree","DESeq2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


library("devtools")
library("DDRTree")
library("pheatmap")
library("M3Drop")
library("monocle")
library("reshape")
library("DESeq2")

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
#\print(head(fData(HSMM)))
#\print(head(pData(HSMM)))
HSMM<-estimateSizeFactors(HSMM, locfunc=genefilter::shorth ) # if Size_Factor=NA, there are too many zeroes, 
#\print(head(pData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
#length(expressed_genes)


#If you are using RPC values to measure expresion, 
#as we are in this vignette, it's also good to look at the distribution
#of mRNA totals across the cells:
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
#\qplot(Total_mRNAs, data=pData(HSMM), geom="density") +
#\        geom_vline(xintercept=lower_bound) +
#\        geom_vline(xintercept=upper_bound)
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
                     pData(HSMM)$Total_mRNAs < upper_bound]
HSMM <- detectGenes(HSMM, min_expr = 0.1)
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
#KRT5 - rare population of basal stem cells
cth <- addCellType(cth, "KRT5 rare population of basal stem cells", classify_func=function(x) {
    x[id_("KRT5"),] >= 1 })
#MUC1 - majority of differentiated cells
cth <- addCellType(cth, "MUC1 majority of differentiated cells", classify_func=function(x) {
    x[id_("MUC1"),] >= 1 })
#SCGB3A2 - early secretory cell differentiation
cth <- addCellType(cth, "SCGB3A2 early secretory cell differentiation", classify_func=function(x) {
    x[id_("SCGB3A2"),] >= 1 })
#SCGB1A1 - non-mucous club secretory cells
cth <- addCellType(cth, "SCGB1A1 non-mucous club secretory cells", classify_func=function(x) {
    x[id_("SCGB1A1"),] >= 1 })
#SCGB3A1 - non-mucous club secretory cells
cth <- addCellType(cth, "SCGB3A1 non-mucous club secretory cells", classify_func=function(x) {
    x[id_("SCGB3A1"),] >= 1 })
#SFTPB - secretory cells (majority of cells)
cth <- addCellType(cth, "SFTPB secretory cells (majority of cells)", classify_func=function(x) {
    x[id_("SFTPB"),] >= 1 })
#FOXJ1 - ciliated cells (could be rare in this region)
cth <- addCellType(cth, "FOXJ1 ciliated cells", classify_func=function(x) {
        x[id_("FOXJ1"),] >= 1 })

#CDH1 - Cadherin for all epithelial)
cth <- addCellType(cth, "Cadherin+ epithelial Cells", classify_func=function(x) {
    x[id_("CDH1"),] >= 1 })


#Non-epithelial structural
#Stromal+endothelial - GPX3 + negative for leukocyte markers (PTPRC, LAPTM5, SRGN)------------
#Stromal/fibroblasts - DCN, COL6A1, TIMP3, PDGFRA
#Endothelial - VWF, EMCN, PECAM1, CDH5
cth <- addCellType(cth, "Stromal + endothelial + fibroblasts Cells", classify_func=function(x) {
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
pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)

pie + coord_polar(theta = "y") +theme(axis.title.x=element_blank(), 
                                      axis.title.y=element_blank(),
                                      text = element_text(size=25))

#3.2 Unsupervised cell clustering
#estimateDispersions only makes sense right now for negbinomial and negbinomial.size expression families.
HSMM<-estimateDispersions(HSMM)
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
#plot_ordering_genes(HSMM)
# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
#plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',


HSMM_1 <- reduceDimension(HSMM, max_components=2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = T)
HSMM_1 <- clusterCells(HSMM_1,num_clusters=2)
## Distance cutoff calculated to 1.243471
plot_cell_clusters(HSMM_1, 1,2, color="CellType", markers=c("GPX3","KRT19","SRGN"))+
    theme(text = element_text(size=30),
          legend.position="right",
          legend.key.width = unit(2, "cm"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))

plot_cell_clusters(HSMM_1, 1,2, color="CellType", markers=c("KRT19","KRT5","MUC1","SCGB3A2","SCGB1A1","SCGB3A1","SFTPB","FOXJ1"))+
    theme(text = element_text(size=30),
          legend.position="right",
          legend.key.width = unit(2, "cm"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))

plot_cell_clusters(HSMM_1, 1, 2, color="num_genes_expressed")

#Monocle allows us to subtract the effects of \uninteresting" sources of variation
#to reduce their impact on the clustering.
HSMM_2 <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE',
#                        residualModelFormulaStr="~CellType + num_genes_expressed",
                        verbose = T)
HSMM_2 <- clusterCells(HSMM_2, num_clusters=3)
plot_cell_clusters(HSMM_2, 1, 2, color="CellType")

#Now that we've accounted for some unwanted sources of variation, 
#we're ready to take another crack at classifying the cells by unsupervised clustering:

HSMM_3 <- clusterCells(HSMM_2, num_clusters=4)
plot_cell_clusters(HSMM_3, 1, 2, color="Cluster") + facet_wrap(~CellType)+
    theme(text = element_text(size=25),
          legend.justification = 'right', 
          legend.position=c(0.8,0.3),
          legend.key.width = unit(2, "cm"))+ 
    guides(colour = guide_legend(override.aes = list(size=10)))


#3.3 Semi-supervised cell clustering with known marker genes
marker_diff <- markerDiffTable(HSMM[expressed_genes,], 
                               cth,
#                               residualModelFormulaStr="~num_genes_expressed", # fail using ~CellType
                               cores=1) #Long cth will trigure "subscript out of bounds". It takes long time.
 
candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)

#Below shows top three marker genes for each cell type
selectTopMarkers(marker_spec, 3)

#To cluster the cells, we'll choose the top 500 markers for each of these cell types:
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',
HSMM_4 <- reduceDimension(HSMM, max_components=2, num_dim = 5, reduction_method = 'tSNE',
                        residualModelFormulaStr="~num_genes_expressed", verbose = T)
#HSMM_4 <- clusterCells(HSMM_4, num_clusters=3)
#plot_cell_clusters(HSMM_4, 1, 2, color="CellType")


#3.4 Imputing cell type
HSMM_4 <- clusterCells(HSMM_4,
                     num_clusters=5,
                     frequency_thresh=0.1,
                     cell_type_hierarchy=cth)

plot_cell_clusters(HSMM_4, 1, 2, color="CellType", markers = c("KRT19", "GPX3","SRGN"))
pie <- ggplot(pData(HSMM_4), aes(x = factor(1), fill = factor(CellType))) +
    geom_bar(width = 1)
pie + coord_polar(theta = "y") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
#Finally, we subset the CellDataSet object to create HSMM_epi, 
#which includes only myoblasts. We'll use this in the
HSMM_epi <- HSMM_1[,pData(HSMM_1)$CellType == "KRT19+ epithelial Cells"]
HSMM_epi <- estimateDispersions(HSMM_epi)


#4 Constructing single cell trajectories
#4.3 Unsupervised ordering
diff_test_res <- differentialGeneTest(HSMM_epi[expressed_genes,],fullModelFormulaStr="~num_genes_expressed")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

#4.3.1 Selecting genes with high dispersion across cells
disp_table <- dispersionTable(HSMM_epi)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.5 & 
                        dispersion_empirical >= 1 * dispersion_fit)$gene_id

HSMM_epi <- setOrderingFilter(HSMM_epi, ordering_genes)
plot_ordering_genes(HSMM_epi)
HSMM_epi <- reduceDimension(HSMM_epi, max_components=2)
HSMM_epi <- orderCells(HSMM_epi)
plot_cell_trajectory(HSMM_epi, color_by="Pseudotime")
plot_cell_trajectory(HSMM_epi, color_by="State")
GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$Pseudotime)[,"0"] # replace $Hours with $Pseudotime
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    }else { return (1) }
}
HSMM_epi <- orderCells(HSMM_epi, root_state=GM_state(HSMM_epi))
plot_cell_trajectory(HSMM_epi, color_by="Pseudotime")
plot_cell_trajectory(HSMM_epi, color_by="State") + facet_wrap(~State, nrow=1)
#We can use the jitter plot to pick figure out which state corresponds to rapid proliferation:
blast_genes <- row.names(subset(fData(HSMM_epi), gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
plot_genes_jitter(HSMM_epi[blast_genes,], grouping="State", min_expr=0.1)

#To confirm that the ordering is correct we can select a couple of markers of myogenic progress. 
#Plotting these genes demonstrates that ordering looks good:
HSMM_expressed_genes <- row.names(subset(fData(HSMM_epi), num_cells_expressed >= 10))
HSMM_filtered <- HSMM_epi[HSMM_expressed_genes,]

my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="State") # replace Hours with State

#4.3.2 Selecting genes based on PCA loading
HSMM_epi <- HSMM_epi[HSMM_expressed_genes,]
exprs_filtered <- t(t(exprs(HSMM_epi)/pData(HSMM_epi)$Size_Factor))
#nz_genes <- which(exprs_filtered !== 0, arr.ind = T)
exprs_filtered <- log(exprs_filtered + 1)
# Calculate the variance across genes without converting to a dense
# matrix:
expression_means <- Matrix::rowMeans(exprs_filtered)
expression_vars <- Matrix::rowMeans((exprs_filtered - expression_means)^2)
# Filter out genes that are constant across all cells:
genes_to_keep <- expression_vars > 0
exprs_filtered <- exprs_filtered[genes_to_keep,]
expression_means <- expression_means[genes_to_keep]
expression_vars <- expression_vars[genes_to_keep]
# Here's how to take the top PCA loading genes, but using
# sparseMatrix operations the whole time, using irlba. Note
# that the v matrix from irlba is the loading matrix
set.seed(0)
irlba_pca_res <- irlba(t(exprs_filtered),
                       nu=0,
                       center=expression_means,
                       scale=sqrt(expression_vars),
                       right_only=TRUE)$v
row.names(irlba_pca_res) <- row.names(exprs_filtered)
# Here, we will just
# take the top 200 genes from components 2 and 3.
# Component 1 usually is driven by technical noise.
# We could also use a more principled approach,
# similar to what dpFeature does below
PC2_genes <- names(sort(abs(irlba_pca_res[, 2]), decreasing = T))[1:200]
PC3_genes <- names(sort(abs(irlba_pca_res[, 3]), decreasing = T))[1:200]
ordering_genes <- union(PC2_genes, PC3_genes)
#Using these to order the cells as above yields the following trajectory:
HSMM_epi <- setOrderingFilter(HSMM_epi, ordering_genes)
HSMM_epi <- reduceDimension(HSMM_epi, max_components=2)
HSMM_epi <- orderCells(HSMM_epi)
HSMM_epi <- orderCells(HSMM_epi, root_state=GM_state(HSMM_epi))
plot_cell_trajectory(HSMM_epi, color_by="State")

#4.4 Unsupervised feature selection based on density peak clustering
#To use dpFeature, we firrst select superset of feature genes as genes expressed in at least 5% of all the cells.
HSMM_epi <- detectGenes(HSMM_epi, min_expr=0.1)
fData(HSMM_epi)$use_for_ordering <- fData(HSMM_epi)$num_cells_expressed > 0.05 * ncol(HSMM_epi)

plot_pc_variance_explained(HSMM_epi, return_all = F) #look at the plot and decide how many di-
#mensions you need. It is determined by a huge drop of variance at that dimension. pass that num-
#ber to num_dim in the next function.

#We will then run reduceDimension with t-SNE as the reduction method on those top PCs and project them further
#down to two dimensions.
HSMM_epi <- reduceDimension(HSMM_epi, max_components=2, norm_method = 'log', num_dim = 3,
                            reduction_method = 'tSNE', verbose = T)

HSMM_epi <- clusterCells(HSMM_epi, verbose = F)
plot_cell_clusters(HSMM_epi, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_epi, color_by = 'as.factor(State)')  # replace Hours with State
plot_rho_delta(HSMM_epi, rho_threshold = 2, delta_threshold = 4 )

HSMM_epi <- clusterCells(HSMM_epi,
                         rho_threshold = 2,
                         delta_threshold = 4,
                         skip_rho_sigma = T,
                         verbose = F)

#We can check the final clustering results as following:
plot_cell_clusters(HSMM_epi, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_epi, color_by = 'as.factor(State)') # replace Hours with State
#After we confirm the clustering makes sense, 
#we can then perform differential gene expression test as a way to extract the genes that distinguish them.
clustering_DEG_genes <- differentialGeneTest(HSMM_epi[HSMM_expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 1) #takes long time
#We will then select the top 1000 significant genes as the ordering genes.
HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
HSMM_epi <- setOrderingFilter(HSMM_epi, ordering_genes = HSMM_ordering_genes)
HSMM_epi <- reduceDimension(HSMM_epi)
HSMM_epi <- orderCells(HSMM_epi)
HSMM_epi <- orderCells(HSMM_epi, root_state=GM_state(HSMM_epi))
plot_cell_trajectory(HSMM_epi, color_by="State") # replace Hours with State

#4.5 Semi-supervised ordering with known marker genes

cth <- newCellTypeHierarchy()

#KRT19 - all (majority) of epithelial; + negative for non-epithelial (DCN, PECAM1, LAPTM5, CD68)----------
cth <- addCellType(cth, "KRT19 all(majority) of epithelial", classify_func=function(x) {
    x[id_("KRT19"),] >= 1 & x[id_("DCN"),] == 0 & x[id_("PECAM1"),] == 0 & x[id_("LAPTM5"),] == 0 & x[id_("CD68"),] == 0 })
#KRT5 - rare population of basal stem cells
cth <- addCellType(cth, "KRT5 rare population of basal stem cells", classify_func=function(x) {
    x[id_("KRT5"),] >= 1 })
#MUC1 - majority of differentiated cells
cth <- addCellType(cth, "MUC1_majority of differentiated cells", classify_func=function(x) {
    x[id_("MUC1"),] >= 1 })
#SCGB3A2 - early secretory cell differentiation
cth <- addCellType(cth, "SCGB3A2 early secretory cell differentiation", classify_func=function(x) {
    x[id_("SCGB3A2"),] >= 1 })
#SCGB1A1 - non-mucous club secretory cells
cth <- addCellType(cth, "SCGB1A1 non-mucous club secretory cells", classify_func=function(x) {
    x[id_("SCGB1A1"),] >= 1 })
#SCGB3A1 - non-mucous club secretory cells
cth <- addCellType(cth, "SCGB3A1 non-mucous club secretory cells", classify_func=function(x) {
    x[id_("SCGB3A1"),] >= 1 })
#SFTPB - secretory cells (majority of cells)
cth <- addCellType(cth, "SFTPB secretory cells (majority of cells)", classify_func=function(x) {
    x[id_("SFTPB"),] >= 1 })
#FOXJ1 - ciliated cells (could be rare in this region)
cth <- addCellType(cth, "FOXJ1 ciliated cells", classify_func=function(x) {
    x[id_("FOXJ1"),] >= 1 })


HSMM_epi <- classifyCells(HSMM_epi, cth)
#Now we select the set of genes that co-vary (in either direction) with these two \bellweather" genes:
marker_diff <- markerDiffTable(HSMM_epi[HSMM_expressed_genes,],
                               cth,
                               cores=1) #takes long time
#semisup_clustering_genes <- row.names(subset(marker_diff, qval < 0.05))
semisup_clustering_genes <- row.names(marker_diff)[order(marker_diff$qval)][1:500]

#Using the top 1000 genes for ordering produces a trajectory that's highly similar to the one we obtained with
#unsupervised methods, but it's a little \cleaner".
HSMM_epi <- setOrderingFilter(HSMM_epi, semisup_clustering_genes)
#plot_ordering_genes(HSMM_epi)
HSMM_epi <- reduceDimension(HSMM_epi, max_components=2)
HSMM_epi <- orderCells(HSMM_epi)
HSMM_epi <- orderCells(HSMM_epi, root_state=GM_state(HSMM_epi))
plot_cell_trajectory(HSMM_epi, color_by="CellType") + theme(legend.position="right")
#To confirm that the ordering is correct, we can select a couple of markers of myogenic progress. 
#In this experiment, one of the branches corresponds to cells that successfully fuse to form myotubes,
#and the other to those that fail to fully differentiate. 
#We'll exclude the latter for now, but you can learn more about tools for dealing with branched trajectories in section 6.

HSMM_filtered <- HSMM_epi[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("KRT19", "SCGB3A2", "FOXJ1")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_branched_pseudotime(cds_subset,
                               branch_point=1,
                               color_by="State",  # replace Hours with State
                               ncol=1)


