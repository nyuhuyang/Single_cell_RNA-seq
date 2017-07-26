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
list.of.bio.packages <- c("scater","org.Hs.eg.db","scran")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
library(easypackages)
libraries(list.of.cran.packages,list.of.bio.packages)

########################################################################
#
#  1.5.1  Count loading (Required)
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
#exprs(gbm) are counts, exprs(sce) are fpkm(maybe)

########################################################################
#
#  B. Quality control on the cells (Recommended)
# 
# ######################################################################
# library size(counts per cell) =  total sum of counts across all features.
# expressed features (genes per cell) = the number of features with non-zero counts
#  Cells with small library sizes and few expressed genes are like to be of poor quality.
par(mfrow=c(1,3))


hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=30, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=30, col="grey80", ylab="Number of cells")
# remove cells with log-library sizes that are more than 3 median absolute deviations (MADs)
# below the median log-library size.
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
table(libsize.drop)

feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
table(feature.drop)

sce <- sce[,!(libsize.drop | feature.drop )]#subsetting by column
hist(sce$total_features, xlab="Number of expressed genes", main="nmads=3", 
     breaks=30, col="grey80", ylab="Number of cells")



# Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data quality
dim(sce)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           Remaining=ncol(sce))

mean(sce$total_counts)
median(sce$total_features)

#-------QC on mt- & ERCC (Alternative)--

#-------QC on mt- & ERCC (Alternative)-----------------------------------------
par(mfrow=c(1,2))
# High mitochondrial proportion indicates poor-quality of cells.
hist(sce$pct_counts_feature_controls_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
# High proportion spike-in RNA indicates loss of endogenous RNA in low-quality cells.
hist(sce$pct_counts_feature_controls_ERCC, xlab="ERCC proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

# use the MAD-based definition of outliers to remove putative low-quality cells from the dataset
mito.drop <- isOutlier(sce$pct_counts_feature_controls_Mt, nmads=3, type="higher")
spike.drop <- isOutlier(sce$pct_counts_feature_controls_ERCC, nmads=3, type="higher")
table(mito.drop)
table(spike.drop)

sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
# Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data quality
dim(sce)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(sce))


#------QC by PCA (Alternative)--------------------------------------------------
fontsize <- theme(axis.text=element_text(size=12), 
                  axis.title=element_text(size=16),
                  plot.title = element_text(size=22))
plotPCA(sce, pca_data_input="pdata") + fontsize +
        ggtitle("PCA based on the quality metrics for each cell")

########################################################################
#
#  C. Classification of cell cycle phase (Alternative)
# 
# ######################################################################
hm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
anno <- select(org.Hs.eg.db, keys=as.character(fData(sce)$Gene), keytype="SYMBOL", column="ENSEMBL")
# try to construct anno by sce
# anno <- data.frame(SYMBOL= fData(sce)$Gene, ENSEMBL =rownames(fData(sce)) )
ensembl <- anno$ENSEMBL[match(fData(sce)$Gene, anno$SYMBOL)]
assignments <- cyclone(exprs(sce), hm.pairs, gene.names=ensembl) # takes long time
par(mfrow=c(1,1))
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
# G1 phase <- (G1>0.5 & G1> G2/M)
sce_G1 <- sce[,assignments$phases=="G1"]
# G2/M phase <- (G2/M >0.5 & G2/M >G1)
# s phase <- (G1 <= 0.5 & G2/M <= 0.5)

########################################################################
#
#  D. Filtering out low-abundance genes (recommended)
# 
# ######################################################################
#low-abundance genes = an average count below a filter threshold of 1
ave.counts <- calcAverage(sce)
keep <- ave.counts >= 0.01
table(keep)
par(mfrow=c(1,1))
hist(log10(ave.counts), 
     breaks=100, main="Before gene filtering", col="grey80",
     xlab=expression(Log[10]~"average count for each gene"))
abline(v=log10(1), col="blue", lwd=2, lty=2)
abline(v=log10(0.01), col="red", lwd=2, lty=2)

#In general, we prefer the mean-based filter as it tends to be less aggressive.
sce <- sce[keep,]
sce <- calculateQCMetrics(sce)

mean(sce$total_counts)
median(sce$total_features)
saveRDS(sce,"sce")
#-----------highly expressed genes (Alternative)--------------------------------
#should generally be dominated by constitutively expressed transcripts
#such as those for ribosomal or mitochondrial proteins. 
fontsize <- theme(axis.text=element_text(size=12), 
                  axis.title=element_text(size=16),
                  plot.title = element_text(size=22))
plotQC(sce, type = "highest-expression", n=50) + fontsize #long long time

#-----select non-zero counts genes at least n cells (Alternative)-------------
numcells <- nexprs(sce, byrow=TRUE)
alt.keep <- numcells >= 10
sum(alt.keep)
smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"), 
              ylab="Number of expressing cells")
is.ercc <- isSpike(sce, type="ERCC")
points(log10(ave.counts[is.ercc]), numcells[is.ercc], col="red", pch=16, cex=0.5)


#---Checking for important technical factors(Alternative)-----
plotExplanatoryVariables(sce, variables=c("counts_feature_controls_ERCC", 
                                          "log10_counts_feature_controls_ERCC")) + fontsize
#Error in if (any(matrixStats::rowVars(exprs_mat) == 0)) stop("Some features have zero variance.
#Please filter out features with zero variance (e.g. all zeros).") : 
#missing value where TRUE/FALSE needed

########################################################################
#
#  E. Normalization of cell-specific biases (Alternative)
# 
# ######################################################################
## Using the deconvolution method to deal with zero counts
sce <- computeSumFactors(sce, sizes=seq(20, 80, 5))
summary(sizeFactors(sce))
par(mfrow=c(1,2))
plot(sizeFactors(sce), colSums(exprs(sce))/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor",
     main="Before normalization")
## Computing separate size factors for spike-in transcripts
#sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)
#Warning message:
#        In .local(x, ...) : zero spike-in counts during spike-in normalization

#Applying the size factors to normalize gene expression
sce.norm <- normalize(sce)
plot(sizeFactors(sce), colSums(exprs(sce.norm))/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor",
     main="After normalization")


########################################################################
#
#  F. Identifying highly variable genes from the normalized log-expression
#  (Alternative)
#
# ######################################################################
var.fit <- trendVar(sce.norm, trend="loess", use.spikes=FALSE, span=0.2)
#Error in out[lower] <- x[lower] * (left.val/left.edge) : 
#NAs are not allowed in subscripted assignments
fitkeep<-names(na.omit(var.fit$mean)) #remove NaN
sce.fitkeep <- sce.norm[fitkeep,]
dim(sce.fitkeep)
var.fit <- trendVar(sce.fitkeep, trend="loess", use.spikes=FALSE, span=0.2)
var.out <- decomposeVar(sce.fitkeep, var.fit)
par(mfrow=c(1,1))
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
cur.spike <- isSpike(sce.fitkeep)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

# HVG = significantly greater than zero at a (FDR) of 5%
# HVG = biological component greater than or equal to 0.5
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
nrow(hvg.out)
head(hvg.out)
plotExpression(sce.fitkeep, rownames(hvg.out)) + fontsize


               