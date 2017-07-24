#https://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/vignette.html

#Introduction to scater: Single-cell analysis toolkit for expression in R
########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages <- c("cowplot")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("scater")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

libraries(list.of.bio.packages,list.of.cran.packages)
memory.limit(size=65000)

########################################################################
#
# 1.5 Scater
#
########################################################################

#====== 1.5.0 Setup enviroment and read data (Required)============
#detect OS and set enviroment
# load the pipeline data by specifying a pipestance path
if (Sys.info()[['sysname']]=="Darwin"){
        WD <- "/Users/yah2014/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        WD <- "C:/Users/User/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}


##------1.5.1 quickstart (Alternative)--------------------------
gbm <- readRDS("gbm")
## form an SCESet object
pd <- new("AnnotatedDataFrame", data = pData(gbm))
rownames(pd) <- pd$barcode
sce <- newSCESet(countData = exprs(gbm), phenoData = pd)

## filter no exprs
keep_feature <- rowSums(exprs(sce) > 0) > 0
sce <- sce[keep_feature,]

## quick start, eval=FALSE
sce <- calculateQCMetrics(sce, feature_controls = 1:40)


scater_gui(sce) #quick-start-gui
##-----------------------------------------------------------------------------


## ====1.5.2 sceset make sceset counts only (Recommended)==========================
gbm <- readRDS("gbm")
pd <- new("AnnotatedDataFrame", data = pData(gbm))
rownames(pd) <- pd$barcode
gene_df <- data.frame(Gene = fData(gbm)$symbol)
rownames(gene_df) <- rownames(exprs(gbm))
fd <- new("AnnotatedDataFrame", data = gene_df)
sce <- newSCESet(countData = exprs(gbm), phenoData = pd,
                            featureData = fd)
sce
dim(sce)
## ==========================================================================

## ----1.5.3 sceset make sceset exprs only (Alternative, practise)-----------------------------------------
example2 <- newSCESet(exprsData = log2(calculateCPM(sce) + 1))
pData(example2) ## data frame with 0 columns and 1883 rows
fData(example2) ## data frame with 0 columns and 32738 rows

## compare following three functions
counts(example2)[1:10, 1:3]
#counts(object): returns the matrix of read counts. 
#If no counts are defined for the object, 
#then the counts matrix slot is simpy NULL.

exprs(example2)[1:10, 1:3]
#exprs(object): returns the matrix of feature expression values. 
#Typically these should be log2(counts-per-million) values or 
#log2(reads-per-kilobase-per-million-mapped), appropriately normalised.

get_exprs(sce, "counts")[1:10, 1:3]
# access any assay data from the object with the get_exprs function.

## add some count data to SCESet if necessary
set_exprs(example2, "counts") <- counts(sce)
counts(example2) <- get_exprs(sce, "counts")
example2
counts(example2)[1:10, 1:3]

## sceset demo replacement
gene_df <- data.frame(Gene = rownames(exprs(gbm)))
rownames(gene_df) <- gene_df$Gene
fd <- new("AnnotatedDataFrame", data = gene_df)
## replace featureData
fData(sce) <- fd
## replace phenotype data
pData(sce) <- pd
## replace expression data to be used
exprs(sce) <- log2(calculateCPM(sce) + 1)


## ----1.5.4  Plots of expression values (Alternative, practise)-----------------
## plot sceset blocking
#plot shows the cumulative proportion across the top 500 features
plot(sce, block1 = "Mutation_Status", block2 = "Treatment",
     colour_by = "Cell_Cycle", nfeatures = 300, exprs_values = "counts")

## plot expression
# plotexpression function plot xpression values for a subset of genes or features
plotExpression(sce, rownames(sce)[1:6],
               x = "Mutation_Status", exprs_values = "exprs", colour = "Treatment")

## plot expression theme bw
plotExpression(sce, rownames(sce)[7:12],
               x = "Mutation_Status", exprs_values = "counts", colour = "Cell_Cycle",
               show_median = TRUE, show_violin = FALSE,  xlab = "Mutation Status",
               log = TRUE)
##--------------------------------------------------------------------

########################################################################
#
# 1.6 Scater Quality Control
#
########################################################################
#The scater package puts a focus on aiding with quality control (QC) 
#and pre-processing of single-cell RNA-seq data before further downstream analysis.

#We see QC as consisting of three distinct steps:
        
# 1.QC and filtering of features (genes)
# 2.QC and filtering of cells
# 3.QC of experimental variables

## 1.6.1 calc QCmetrics
sce <- calculateQCMetrics(sce, feature_controls = 1:20)
varLabels(sce)  #list-pata-qc
names(fData(sce)) #list-fdata-qc
#=============================================

## calc-qc-metrics multi feature controls(Alternative)-------------------------
#More than one set of feature controls can be defined if desired.
sce <- calculateQCMetrics(
        sce, feature_controls = list(controls1 = 1:20, controls2 = 500:1000),
        cell_controls = list(set_1 = 1:5, set_2 = 31:40))
varLabels(sce)


## 1.6.2 QC and filtering of features(Recommended)======================================
#filtering out unwanted features, with very low overall expression

# plot-qc-expression, fig.height=7.5, fig.width=8.5
keep_feature <- rowSums(counts(sce) > 0) > 4
sce <- sce[keep_feature,]
## Plot QC
plotQC(sce, type = "highest-expression", exprs_values = "counts")

## ----plot-qc-expression-cell-controls(Alternative)---------------------------------
p1 <- plotQC(sce[, !sce$is_cell_control],
             type = "highest-expression")
p2 <- plotQC(sce[, sce$is_cell_control],
             type = "highest-expression")
multiplot(p1, p2, cols = 2)

## plot-qc-exprs-freq-vs-mean-default(Recommended)=====
plotQC(sce, type = "exprs-freq-vs-mean")

## -plot-qc-exprs-mean-vs-freq-defined-feature-set(Alternative)------
feature_set_1 <- fData(sce)$is_feature_control_controls1
plotQC(sce, type = "exprs-freq-vs-mean", feature_set = feature_set_1)

## plot-fdata(Recommended)===
plotFeatureData(sce, aes(x = n_cells_exprs, y = pct_total_counts))

## 1.6.3 QC and filtering of cells ======================================
## plot-pdata
plotPhenoData(sce, aes(x = total_counts, y = total_features,
                                  colour = Mutation_Status))

## ----plot-pdata-cont-col, fig.show = TRUE----------------------------------
plotPhenoData(sce, aes(x = Mutation_Status, y = total_features,
                                  colour = log10_total_counts))

## ----plot-pdata-col-gene-exprs---------------------------------------------
plotPhenoData(sce, aes(x = total_counts, y = total_features,
                                  colour = Gene_1000))

## ----plot-pdatacol-gene-exprs-2, fig.show = FALSE--------------------------
plotPhenoData(sce, aes(x = pct_counts_feature_controls,
                                  y = total_features, colour = Gene_0500))

## ----plot-pdatacol-gene-exprs-3, fig.show = FALSE--------------------------
plotPhenoData(sce, aes(x = pct_counts_feature_controls,
                                  y = pct_counts_top_50_features,
                                  colour = Gene_0001))

## ----plot-pdata-pct-exprs-controls-----------------------------------------
plotPhenoData(sce, aes(x = total_features,
                                  y = pct_counts_feature_controls,
                                  colour = Mutation_Status)) +
        theme(legend.position = "top") +
        stat_smooth(method = "lm", se = FALSE, size = 2, fullrange = TRUE)

## ----plot-pca-default------------------------------------------------------
plotPCA(sce)

## ----plot-pca-cpm, eval=FALSE----------------------------------------------
#  plotPCA(sce, exprs_values = "cpm")

## ----plot-pca-feature-controls, fig.show = FALSE---------------------------
plotPCA(sce, feature_set = fData(sce)$is_feature_control)

## ----plot-pca-4comp-colby-shapeby, fig.height=5.5--------------------------
plotPCA(sce, ncomponents = 4, colour_by = "Treatment",
        shape_by = "Mutation_Status")

## ----plot-pca-4comp-colby-sizeby-exprs, fig.height=5.5---------------------
plotPCA(sce, colour_by = "Gene_0001", size_by = "Gene_1000")

## ----plot-tsne-1comp-colby-sizeby-exprs, fig.height=5.5--------------------
plotTSNE(sce, colour_by = "Gene_0001", size_by = "Gene_1000")

## ----plot-difmap-1comp-colby-sizeby-exprs, fig.height=5.5------------------
plotDiffusionMap(sce, colour_by = "Gene_0001", size_by = "Gene_1000")

## ----plot-pca-4comp-colby-shapeby-save-pcs, fig.show = FALSE---------------
sce <- plotPCA(sce, ncomponents = 4,
                          colour_by = "Treatment", shape_by = "Mutation_Status",
                          return_SCESet = TRUE, theme_size = 12)
head(reducedDimension(sce))

## ----plot-reduceddim-4comp-colby-shapeby, fig.show=FALSE-------------------
plotReducedDim(sce, ncomponents = 4, colour_by = "Treatment",
               shape_by = "Mutation_Status")

## ----plot-reduceddim-4comp-colby-sizeby-exprs, fig.show = FALSE------------
plotReducedDim(sce, ncomponents = 4, colour_by = "Gene_1000",
               size_by = "Gene_0500")

## ----plot-pca-outlier------------------------------------------------------
sce <- plotPCA(sce, pca_data_input = "pdata", 
                          detect_outliers = TRUE, return_SCESet = TRUE)


## ----plot-qc-expl-variables-all, warning=FALSE-----------------------------
plotQC(sce, type = "expl")

## ----plot-qc-expl-variables-select-variables-------------------------------
plotQC(sce, type = "expl",
       variables = c("total_features", "total_counts", "Mutation_Status", "Treatment",
                     "Cell_Cycle"))

## ----plot-qc-pairs-pc------------------------------------------------------
plotQC(sce, type = "expl", method = "pairs", theme_size = 6)

## ----plot-qc-find-pcs-pcs-vs-vars, fig.width=8, fig.height=7---------------
p1 <- plotQC(sce, type = "find-pcs", variable = "total_features",
             plot_type = "pcs-vs-vars")
p2 <- plotQC(sce, type = "find-pcs", variable = "Cell_Cycle",
             plot_type = "pcs-vs-vars")
multiplot(p1, p2, cols = 2)

## ----plot-qc-find-pcs-pairs, fig.width=10, fig.height=7--------------------
plotQC(sce, type = "find-pcs", variable = "total_features",
       plot_type = "pairs-pcs")

## ----plot-qc-find-pcs-pairs-2, fig.show=FALSE------------------------------
plotQC(sce, type = "find-pcs", variable = "Cell_Cycle",
       plot_type = "pairs-pcs")

## ----cell-pairwise-distance-matrices-euclidean, eval=TRUE------------------
cell_dist <- dist(t(exprs(sce)))
cellPairwiseDistances(sce) <- cell_dist
plotMDS(sce)

## ----cell-pairwise-distance-matrices-canberra, eval=TRUE, fig.show = FALSE----
cell_dist <- dist(t(counts(sce)), method = "canberra")
cellPairwiseDistances(sce) <- cell_dist
plotMDS(sce, colour_by = "Mutation_Status")

## ----feature-pairwise-distance-matrices, eval=FALSE------------------------
#  feature_dist <- dist(exprs(sce))
#  featurePairwiseDistances(sce) <- feature_dist
#  limma::plotMDS(featDist(sce))

## ----kallisto-demo-kallisto-test-data, eval=FALSE--------------------------
#  ################################################################################
#  ### Tests and Examples
#  
#  # Example if in the kallisto/test directory
#  setwd("/home/davis/kallisto/test")
#  kallisto_log <- runKallisto("targets.txt", "transcripts.idx", single_end=FALSE,
#              output_prefix="output", verbose=TRUE, n_bootstrap_samples=10)
#  
#  sce_test <- readKallistoResults(kallisto_log, read_h5=TRUE)
#  sce_test

## ----kallisto-cell-cycle-example, eval=FALSE-------------------------------
#  setwd("/home/davis/021_Cell_Cycle/data/fastq")
#  system("wc -l targets.txt")
#  ave_frag_len <- mean(c(855, 860, 810, 760, 600, 690, 770, 690))
#  
#  kallisto_test <- runKallisto("targets.txt",
#                               "Mus_musculus.GRCm38.rel79.cdna.all.ERCC.idx",
#                               output_prefix="kallisto_output_Mmus", n_cores=12,
#                               fragment_length=ave_frag_len, verbose=TRUE)
#  sce_kall_mmus <- readKallistoResults(kallisto_test, read_h5=TRUE)
#  sce_kall_mmus
#  
#  sce_kall_mmus <- readKallistoResults(kallisto_test)
#  
#  sce_kall_mmus <- getBMFeatureAnnos(sce_kall_mmus)
#  
#  head(fData(sce_kall_mmus))
#  head(pData(sce_kall_mmus))
#  sce_kall_mmus[["start_time"]]
#  
#  counts(sce_kall_mmus)[sample(nrow(sce_kall_mmus), size=15), 1:6]
#  
#  ## Summarise expression at the gene level
#  sce_kall_mmus_gene <- summariseExprsAcrossFeatures(
#      sce_kall_mmus, exprs_values="tpm", summarise_by="feature_id")
#  
#  datatable(fData(sce_kall_mmus_gene))
#  
#  sce_kall_mmus_gene <- getBMFeatureAnnos(
#      sce_kall_mmus_gene, filters="ensembl_gene_id",
#      attributes=c("ensembl_gene_id", "mgi_symbol", "chromosome_name",
#                   "gene_biotype", "start_position", "end_position",
#                   "percentage_gc_content", "description"),
#      feature_symbol="mgi_symbol", feature_id="ensembl_gene_id",
#      biomart="ensembl", dataset="mmusculus_gene_ensembl")
#  
#  datatable(fData(sce_kall_mmus_gene))
#  
#  ## Add gene symbols to featureNames to make them more intuitive
#  new_feature_names <- featureNames(sce_kall_mmus_gene)
#  notna_mgi_symb <- !is.na(fData(sce_kall_mmus_gene)$mgi_symbol)
#  new_feature_names[notna_mgi_symb] <- fData(sce_kall_mmus_gene)$mgi_symbol[notna_mgi_symb]
#  notna_ens_gid <- !is.na(fData(sce_kall_mmus_gene)$ensembl_gene_id)
#  new_feature_names[notna_ens_gid] <- paste(new_feature_names[notna_ens_gid],
#            fData(sce_kall_mmus_gene)$ensembl_gene_id[notna_ens_gid], sep="_")
#  sum(duplicated(new_feature_names))
#  featureNames(sce_kall_mmus_gene) <- new_feature_names
#  head(featureNames(sce_kall_mmus_gene))
#  tail(featureNames(sce_kall_mmus_gene))
#  sum(duplicated(fData(sce_kall_mmus_gene)$mgi_symbol))
#  
