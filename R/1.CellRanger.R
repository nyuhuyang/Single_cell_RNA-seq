#on linux server
#cd ~/share
#wget --no-check-certificate --content-disposition --ask-password --user=youremail@your.institution.edu -i links.txt
#tar -xvfz
#cd ~/share/Cell_KB-1/outs
#cellranger mat2csv filtered_gene_bc_matrices_h5.h5 sample123.csv
#rsync -rav yah2014@pascal.med.cornell.edu:/home/yah2014/share/Cell_KB-1/outs/sample123.csv ./

########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################
#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-1.1.0.R")
library(cellrangerRkit)

memory.limit(size=25000)

# #####################################################################
# 
#  1 Load and normalize data
#  http://cf.10xgenomics.com/supp/cell-exp/cellrangerrkit-PBMC-vignette-knitr-1.1.0.pdf
#
# ####################################################################

# =========1.1 load the pipeline data by specifying a pipestance path (Required)====
if (Sys.info()[['sysname']]=="Darwin"){
        WD <- "/Users/yah2014/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        WD <- "C:/Users/User/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}

# ----------1.2 load gbm (Alternative)---------------------------------------------------
gbm <- load_cellranger_matrix(paste0(WD,"/Cell_KB-1"))
analysis_results <- load_cellranger_analysis_results(paste0(WD,"/Cell_KB-1"))

#gbm is an object that stores the barcode filtered gene, 
#expression matrix and metadata, such as gene symbols and barcode ID
dim(exprs(gbm)) # expression matrix
fData(gbm) # data frame of genes
pData(gbm) # data frame of cell barcodes
saveRDS(gbm,file="gbm")

# =========1.1 Load gbm (Recommended)=========================================================
gbm <- readRDS("gbm")


#=====QC (Recommended) ======================================
anyNA(exprs(gbm))
dim(gbm)
counts_per_cell<- colSums(exprs(gbm)) # mean count per cell
genes_per_cell <- apply(exprs(gbm), 2, function(c)sum(c!=0)) # mean gene per cell
mean(counts_per_cell)
median(genes_per_cell)
qplot(counts_per_cell, data=pData(gbm), geom="density")+coord_cartesian(xlim = c(0, 30000))
qplot(genes_per_cell, data=pData(gbm), geom="density")+scale_x_continuous(limits = c(0, 10000))
qplot(counts_per_cell,genes_per_cell, data=pData(gbm))
#============================================================================

#The variable analysis_results contains pre-computed results for 
#principle component analysis (PCA) dimensional reduction, 
#t-SNE (t-Distributed Stochastic Neighbor Embedding) projection.

#access the t-SNE projection and plot the cells colored by UMI counts
tsne_proj <- analysis_results$tsne
visualize_umi_counts(gbm,tsne_proj[c("TSNE.1","TSNE.2")],
                     limits=c(3,4),marker_size=0.05)
#Figure 1: t-SNE projection where each cell is colored by log10 of UMI counts. 
#Color scale represents log10 of UMI counts

#extending reading for UMI:
#http://www.nature.com/nmeth/journal/v11/n2/full/nmeth.2772.html
#http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1778.html
#https://hemberg-lab.github.io/scRNA.seq.course/unique-molecular-identifiers-umis.html

#filter unexpressed genes, normalize the UMI counts for each barcode,
#and use the log-transformed gene-barcode matrix.
use_genes <- get_nonzero_genes(gbm)
gbm_bcnorm <- normalize_barcode_sums_to_median(gbm[use_genes,])
gbm_log <- log_gene_bc_matrix(gbm_bcnorm,base=10)
print(dim(gbm_log))

visualize_umi_counts(gbm_log,tsne_proj[c("TSNE.1","TSNE.2")],
                     limits=c(3,4),marker_size=0.05)

# #####################################################################
#  
# 1. 2 Visualizing signatures of known gene markers
#
# ####################################################################

genes <- c("KRT19", "GPX3", "CD3G", "ANPEP","PTPRC","CD14")
tsne_proj <- analysis_results$tsne
visualize_gene_markers(gbm_log,genes,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(0,1.5))

#Figure 2: t-SNE projection where each cell is colored by normalized expression
#of the marker in the cell. Color scale represents the normalized expression of the marker.


# #####################################################################
#  
# 1. 3 Unbiased analysis using clustering results
#
# ####################################################################
#=======3.1 Visualizing clustering results==========================
n_clu <- 2:10
x=2
km_res <- analysis_results$kmeans # load pre-computed kmeans results
clu_res <- sapply(n_clu, function(x)
        km_res[[paste("kmeans",x,"clusters",sep="_")]]$Cluster)
colnames(clu_res) <- sapply(n_clu, function(x) paste("kmeans",x,sep="."))
visualize_clusters(clu_res,tsne_proj[c("TSNE.1","TSNE.2")])

#Figure 3: t-SNE projection where each cell is colored by cluster
#assigned by kmeans clustering. ID on the right represents the cluster ID.

#===========3.2 Analyzing cluster specific genes===================

#we focus on the k-means clustering result above with 5 clusters
example_K <- 9 # number of clusters (use "Set3" for brewer.pal below if example_K > 8)
example_col <- rev(brewer.pal(example_K,"Set3")) # customize plotting colors
cluster_result <- analysis_results$kmeans[[paste("kmeans",example_K,"clusters",sep="_")]]
visualize_clusters(cluster_result$Cluster,tsne_proj[c("TSNE.1","TSNE.2")])#,colour=example_col)

#Figure 4: t-SNE projection where each cell is colored by
#cluster assigned by kmeans clustering. ID on the right represents the cluster ID.


#The function prioritize_top_genes identifies markers that
#are up-regulated in particular clusters of cells.

# sort the cells by the cluster labels
cells_to_plot <- order_cell_by_clusters(gbm, cluster_result$Cluster)
# order the genes from most up-regulated to most down-regulated in each cluster
prioritized_genes <- prioritize_top_genes(gbm, cluster_result$Cluster, "sseq", min_mean=0.5)

#In this case, we output all the top 10 gene symbols for the 5 clusters to file.
#output_folder <-"/path_to_your_local_folder/pbmc_data_public/pbmc3k/gene_sets"
write_cluster_specific_genes(prioritized_genes, WD, n_genes=20)



#use the prioritized genes to plot a heat-map 
#where the top three most up-regulated genes in each cluster are displayed.

# create values and axis annotations for pheatmap
gbm_pheatmap(log_gene_bc_matrix(gbm), prioritized_genes, cells_to_plot,
             n_genes=10, colour=example_col, limits=c(-1,2))

#Figure 5: Heatmap of scaled expression of top 3 genes (row) in each cell (column). 
#The horizontal and vertical bars around the heat map represents cluster ID
#assigned by kmeans clustering. Color scale represents scaled expression of the gene.

cell_composition(cluster_result$Cluster,
                 anno=c("monocytes","T cells","NK cells","megakaryocytes","B cells",
                        "X","Y","Z","ZZ"))
sessionInfo()

