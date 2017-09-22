########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages <- c("easypackages","Seurat","dplyr",
                           "Matrix")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c()
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
library(easypackages)
libraries(list.of.bio.packages,list.of.cran.packages)

# #####################################################################
# 
#  1. Seurat Clustering 
#      http://satijalab.org/seurat/pbmc3k_tutorial.html
# ####################################################################
#======1.0 Setup enviroment and read data (Required)============
# detect OS and set enviroment
# load the pipeline data by specifying a pipestance path
if (Sys.info()[['sysname']]=="Darwin"){
        WD <- "/Users/yah2014/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        WD <- "C:/Users/User/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}

#======1.1 Setup the Seurat Object (Required)=========================
# Load the Kidney dataset
#Kidney.data <- Read10X(data.dir = "./Cell_KB-1/outs/filtered_gene_bc_matrices/hg19/")
Kidney.data <- Read10X(data.dir = "./Cell_KB-1/outs/raw_gene_bc_matrices/hg19/") #include all
# Examine the memory savings between regular and sparse matrices
#dense.size <- object.size(x = as.matrix(x = Kidney.data))
#dense.size
#sparse.size <- object.size(x = Kidney.data)
#sparse.size

# Initialize the Seurat object with the raw (non-normalized data).
# Keep all genes expressed in >= 3 cells (~0.1% of the data). 
# Keep all cells with at least 200 detected genes
Kidney <- CreateSeuratObject(raw.data = Kidney.data,
                           min.cells = 3,
                           min.genes = 200,
                           project = "10X_Kidney")
str(Kidney@raw.data) #Don't use summary()
#======1.2 Standard pre-processing workflow =====
#======1.2.1 QC and selecting cells for further analysis====
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", 
                   x = rownames(x = Kidney@data),
                   value = TRUE)
percent.mito <- Matrix::colSums(Kidney@raw.data[mito.genes, ])/Matrix::colSums(Kidney@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
Kidney <- AddMetaData(object = Kidney, 
                    metadata = percent.mito, 
                    col.name = "percent.mito")
VlnPlot(object = Kidney, 
        features.plot = c("nGene", "nUMI", "percent.mito"), 
        nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = Kidney, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = Kidney, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
Kidney <- FilterCells(object = Kidney, 
                      subset.names = c("nGene", "percent.mito"), 
                      low.thresholds = c(200, -Inf), 
                      high.thresholds = c(9000, 0.5)) #2500-> 50000, 0.05->0.5

str(Kidney@data) #Don't use summary()
#=====1.3 Normalizing the data========
Kidney <- NormalizeData(object = Kidney, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)


#====1.4 Detection of variable genes across the single cells===
# The parameters here identify ~2,000 variable genes, and represent typical 
# parameter settings for UMI data that is normalized to a total of 1e4 molecules.
par(mfrow = c(1, 1))
Kidney <- FindVariableGenes(object = Kidney, 
                          mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, 
                          x.high.cutoff = 3, 
                          y.cutoff = 0.5)
length(x = Kidney@var.genes)


#=====1.5 Scaling the data and removing unwanted sources of variation=====
Kidney <- ScaleData(object = Kidney, 
                  vars.to.regress = c("nUMI", "percent.mito")) #takes some time

#=====1.6 Perform linear dimensional reduction========
Kidney <- RunPCA(object = Kidney, 
               pc.genes = Kidney@var.genes, 
               do.print = TRUE, 
               pcs.print = 1:5, 
               genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = Kidney, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = Kidney, pcs.use = 1:2)
PCAPlot(object = Kidney, dim.1 = 1, dim.2 = 2)
# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
Kidney <- ProjectPCA(object = Kidney, do.print = FALSE)

PCHeatmap(object = Kidney, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = Kidney, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

#======1.7 Determine statistically significant principal components=====
# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
Kidney <- JackStraw(object = Kidney, num.replicate = 100, do.print = FALSE) #It takes long time
JackStrawPlot(object = Kidney, PCs = 1:12)

PCElbowPlot(object = Kidney)

#======1.8 Cluster the cells ======
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
Kidney <- FindClusters(object = Kidney, 
                     reduction.type = "pca", 
                     dims.use = 1:10, 
                     resolution = 0.8,  #orginal 0.6,Further subdivisions 0.8
                     print.output = 0, 
                     save.SNN = TRUE)
PrintFindClustersParams(object = Kidney) 
# While we do provide function-specific printing functions, the more general
# function to print calculation parameters is PrintCalcParams().

#======1.9 Run Non-linear dimensional reduction (tSNE)========
Kidney <- RunTSNE(object = Kidney, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = Kidney)

#======1.10 Finding differentially expressed genes (cluster biomarkers)===
# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = Kidney, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = Kidney, ident.1 = 5, ident.2 = c(0, 3), 
                                min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
Kidney.markers <- FindAllMarkers(object = Kidney, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25) #It takes long time
Kidney.markers %>% group_by(cluster) %>% top_n(2, avg_diff)

# Four tests for differential expression
# ROC test ("roc"),
# t-test ("t"),
# LRT test based on zero-inflated data ("bimod", default),
# LRT test based on tobit-censoring models ("tobit") 
# The ROC test returns the "classification power" 
# for any individual marker (ranging from 0 - random, to 1 - perfect).

cluster1.markers <- FindMarkers(object = Kidney, ident.1 = 0, thresh.use = 0.25, 
                                test.use = "roc", only.pos = TRUE)
VlnPlot(object = Kidney, features.plot = c("MS4A1", "CD79A"))

# you can plot raw UMI counts as well
VlnPlot(object = Kidney, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)

# find marker gene from proteinatlas.org
# 1.type cell type (eg. epithelial)
# 2.download tsv
# 3.open tsv with excel
# 4.copy the first 16 gene names
# 5.paste in word, replace ^p with ","

FeaturePlot(object = Kidney, 
            features.plot = c("GLS","GPX3","MME","SLC22A2",
                              "CD3E","IL7R","CD8A","FCER1A",
                              "MS4A7","LYZ","GNLY","GZMB",
                              "KLRD1","NKG7","FGF9","SFRP1"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

FeaturePlot(object = Kidney, 
            features.plot = c("CD14","KDR","FLT1","ACVRL1",
                              "TEK","SELP","VWF","CD93",
                              "CD79A","MS4A1","FGF1","NET1",
                              "FGF9","DCN"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

# DoHeatmap generates an expression heatmap for given cells and genes.
# In this case, we are plotting the top 20 markers 
# (or all markers if less than 20) for each cluster.
top10 <- Kidney.markers %>% group_by(cluster) %>% top_n(3, avg_diff)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = Kidney, genes.use = top10$gene,
          slim.col.label = TRUE, remove.key = TRUE)

#======1.11 Assigning cell type identity to clusters========
current.cluster.ids <- 1:16

new.cluster.ids <- c("0) CD4 T cells",
                     "1) kidney cells",
                     "2) kidney cells",
                     "3) kidney cells",
                     "4) CD8 T cells",
                     "5) kidney cells",
                     "6) Dendritic Cells\n & Monocytes",
                     "7) NK cells &\n CD8 T cells",
                     "8) Fibroblasts",
                     "9) Monocytes",
                     "10) Monocytes",
                     "11) Endothelial &\n Megakaryocytes",
                     "12) B cells",
                     "13) Fibroblast &\n Epithelial",
                     "14) Unkown(FGF9 positive)",
                     "15) Stromal Cells",
                     "16) Unkown")
table(Kidney@ident)
Kidney@ident <- plyr::mapvalues(x = Kidney@ident,
                              from = current.cluster.ids,
                              to = new.cluster.ids)
TSNEPlot(object = Kidney, no.legend = TRUE, do.label = TRUE,
         label.size = 6)
cluster16.markers <- FindMarkers(object = Kidney, ident.1 = 16, 
                                min.pct = 0.25)
cluster16 <- cluster16.markers[order(cluster16.markers$p_val),] %>% head(100) %>% rownames()

FeaturePlot(object = Kidney, 
            features.plot = cluster16[1:16], 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
#-----1.12 Further subdivisions within cell types--------
# First lets stash our identities for later
Kidney <- StashIdent(object = Kidney, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the
# SNN, and can simply put: Kidney <- FindClusters(Kidney,resolution = 0.8)
Kidney <- FindClusters(object = Kidney, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.8, print.output = FALSE)
# Demonstration of how to plot two tSNE plots side by side, and how to color
# points based on different criteria
plot1 <- TSNEPlot(object = Kidney, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
plot2 <- TSNEPlot(object = Kidney, do.return = TRUE, group.by = "ClusterNames_0.6", 
                  no.legend = TRUE, do.label = TRUE)
plot_grid(plot1, plot2)


# Find discriminating markers
tcell.markers <- FindMarkers(object = Kidney, ident.1 = 0, ident.2 = 1)

# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we
# can see that CCR7 is upregulated in C0, strongly indicating that we can
# differentiate memory from naive CD4 cells.  cols.use demarcates the color
# palette from low to high expression
FeaturePlot(object = Kidney, features.plot = c("S100A4", "CCR7"), cols.use = c("green", 
                                                                             "blue"))
Kidney <- SetAllIdent(object = Kidney, id = "ClusterNames_0.7")
save(Kidney, file = "./datasets/Kidney_20170920.Rda")




