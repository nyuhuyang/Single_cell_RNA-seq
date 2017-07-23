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

#source("http://bioconductor.org/workflows.R")
#workflowInstall("simpleSingleCell")

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("scater","org.Mm.eg.db")
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
is.spike <- grepl("^ERCC", fData(sce)$Gene) # External RNA Controls Consortium
is.mito <- grepl("^mt-", fData(sce)$Gene)
#calculate quality control metrics such as the total number of counts 
#or the proportion of counts in mitochondrial genes or spike-in transcripts.
sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
head(pData(sce),1)
head(fData(sce),1)
setSpike(sce) <- "ERCC"

########################################################################
#
#  B. Quality control on the cells
# 
# ######################################################################
# library size(counts per cell) =  total sum of counts across all features.
# expressed features (genes per cell) = the number of features with non-zero counts
#  Cells with small library sizes and few expressed genes are like to be of poor quality.
par(mfrow=c(1,2))
hist(sce1$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=30, col="grey80", ylab="Number of cells")
hist(sce1$total_features, xlab="Number of expressed genes", main="", 
     breaks=30, col="grey80", ylab="Number of cells")
#remove cells with log-library sizes that are more than 3 median absolute deviations (MADs)
#below the median log-library size.
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
table(libsize.drop)

feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
table(feature.drop)

sce <- sce[,!(libsize.drop | feature.drop )]#subsetting by column


#Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data quality
dim(sce)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           Remaining=ncol(sce))


#-------QC on mt- & ERCC (Alternative)-----------------------------------------
par(mfrow=c(1,2))
#High mitochondrial proportion indicates poor-quality of cells.
hist(sce$pct_counts_feature_controls_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
#High proportion spike-in RNA indicates loss of endogenous RNA in low-quality cells.
hist(sce$pct_counts_feature_controls_ERCC, xlab="ERCC proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

#use the MAD-based definition of outliers to remove putative low-quality cells from the dataset
mito.drop <- isOutlier(sce$pct_counts_feature_controls_Mt, nmads=3, type="higher")
spike.drop <- isOutlier(sce$pct_counts_feature_controls_ERCC, nmads=3, type="higher")
table(mito.drop)
table(spike.drop)

sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
#Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data quality
dim(sce)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(sce))


#QC by PCA (Alternative)--------------------------------------------------
fontsize <- theme(axis.text=element_text(size=12), 
                  axis.title=element_text(size=16),
                  plot.title = element_text(size=22))
plotPCA(sce, pca_data_input="pdata") + fontsize +
        ggtitle("PCA based on the quality metrics for each cell")

########################################################################
#
#  C. Classification of cell cycle phase
# 
# ######################################################################

anno <- select(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
ensembl <- anno$ENSEMBL[match(rownames(sce), anno$SYMBOL)]
assignments <- cyclone(sce, mm.pairs, gene.names=ensembl)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)













