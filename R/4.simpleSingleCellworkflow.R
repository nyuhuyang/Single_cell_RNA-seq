########################################################################
#
# A step-by-step workflow for low-level analysis of 
# single-cell RNA-seq data with Bioconductor
#
########################################################################
#https://www.bioconductor.org/help/workflows/simpleSingleCell/
#https://f1000research.com/articles/5-2122/v1

########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages <- c("cowplot")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("http://bioconductor.org/workflows.R")
workflowInstall("simpleSingleCell")

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("minqa", "nloptr", "lme4", "SparseM", 
                           "MatrixModels", "modeltools", "prettyunits", 
                           "pbkrtest", "quantreg", "flexmix", "prabclus", 
                           "diptest", "trimcluster", "progress", "lmtest", 
                           "xts", "Rsamtools", "GenomicAlignments", "e1071",
                           "pls", "car", "cvTools", "fpc", "GGally", 
                           "kernlab", "mclust", "sROC", "sp", "vcd", 
                           "laeken", "TTR", "rtracklayer", "bookdown", 
                           "sgeostat", "robCompositions", "FNN", "VIM", 
                           "smoother", "scatterplot3d", "R.oo", "R.methodsS3",
                           "zoo", "statmod", "graph", "GenomicFeatures", 
                           "BiocStyle", "mvoutlier", "destiny", "R.utils", 
                           "readxl", "scran", "RBGL","scran", 
                           "TxDb.Mmusculus.UCSC.mm10.ensGene")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
library(easypackages)
libraries(list.of.cran.packages,list.of.bio.packages)

########################################################################
#
#  A.  Count loading
# 
# ######################################################################
#detect OS and set enviroment
# load the pipeline data by specifying a pipestance path
if (Sys.info()[['sysname']]=="Darwin"){
        WD <- "/Users/yah2014/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        WD <- "C:/Users/User/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}

# 10x Chromium Single Cell RNA-seq data was reviously saved as S4 object "gbm"
gbm <- readRDS("gbm")
#gbm <-gbm[1:2000,1:40] #pilot test
pd <- new("AnnotatedDataFrame", data = pData(gbm))
rownames(pd) <- pd$barcode
gene_df <- data.frame(Gene = fData(gbm)$symbol)
rownames(gene_df) <- rownames(exprs(gbm))
fd <- new("AnnotatedDataFrame", data = gene_df)
sce <- newSCESet(countData = exprs(gbm), phenoData = pd,
                        featureData = fd)
sce
dim(sce)
#We identify the rows corresponding to ERCC spike-ins and mitochondrial genes. 
is.spike <- grepl("^ERCC", fData(sce)$Gene)
is.mito <- grepl("^mt-", fData(sce)$Gene)
#calculate quality control metrics such as the total number of counts 
#or the proportion of counts in mitochondrial genes or spike-in transcripts.
sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
head(pData(sce),1)
head(fData(sce),1)
library(scran)  #?
setSpike(sce) <- "ERCC" #?

########################################################################
#
#  B. Quality control on the cells
# 
# ######################################################################
# library size(counts per cell) =  total sum of counts across all features.
# expressed features (genes per cell) = the number of features with non-zero counts
#  Cells with small library sizes and few expressed genes are like to be of poor quality.
par(mfrow=c(1,2))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
qplot(counts_per_cell, data=sce$total_counts/1e6, geom="density")
qplot(counts_per_cell, data=pData(gbm), geom="density")+coord_cartesian(xlim = c(0, 30000))



