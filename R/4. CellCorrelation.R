########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages <- c("easypackages","stringr","pheatmap","dplyr",
                           "RColorBrewer","stringi")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("oligo","affy", "impute","GEOquery",
                          "monocle","geneplotter")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
library(easypackages)
libraries(list.of.bio.packages,list.of.cran.packages)


# #####################################################################
# 
#  4. Run correlation test of different cell types
#       http://www.ebi.ac.uk
#       https://f1000research.com/articles/5-1384/v1
# ####################################################################
#====== Setup enviroment and read data (Required)============
#detect OS and set enviroment
# load the pipeline data by specifying a pipestance path
if (Sys.info()[['sysname']]=="Darwin"){
        WD <- "/Users/yah2014/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        WD <- "C:/Users/User/Dropbox/Public/Olivier/R/Muthukumar_Single_Cell"
        setwd(WD);getwd();list.files()}

# ----------4.1.0 load Affrimatrix data (Alternative)---------------
eList <- getGEO("GSE3239") #it takes some time
eList
names(eList)
eData_tissue <- eList[[1]] #extract expressionset
eData_tissue
names(pData(eData_tissue))#there is usually a lot of unnecessary stuff here
pData(eData_tissue)$title #enough cell information from title
sampleNames(eData_tissue) <-pData(eData_tissue)$title
#saveRDS(eData_tissue,"eData_tissue")
par(mfrow=c(1,1))
boxplot(eData_tissue)
#normData <- rma(eData_tissue) #  unable to find an inherited method for function 'rma' for signature '"ExpressionSet"'
#normData


#--------clean up the phenotype information(Alternative)-------------
filename <- sampleNames(eData_tissue)
pData(eData_tissue)$filename <- filename
sampleNames <- gsub(" ", ".", filename) #replace space to .
sampleNames <- sub(".Cel.", "", sampleNames) #delete "cell" and "cells"
sampleNames <- sub(".Pooled", "", sampleNames) #delete "Pooled"
sampleNames <- sub("^Human.", "", sampleNames) #delete "Pooled"
sampleNames(eData_tissue) <- sampleNames

sampleNames <- sub(".\\d$", "", sampleNames) #delete .number
sampleNames <- str_extract(sampleNames,'\\w+$') #Extract last word in string
sampleNames <- sub("Ectoepithelial$", "Epithelial", sampleNames) #delete "Ecto"
sampleNames <- sub("Endothelial$", "Endothelials", sampleNames) 
sampleNames <- sub("Muscle$", "Muscles", sampleNames) 
sampleNames <- sub("Stromal$", "Stromals", sampleNames)
table(sampleNames)
pData(eData_tissue)$group <- sampleNames
pData(eData_tissue)$group
saveRDS(eData_tissue,"eData_tissue")
#=======4.1.1 saveRDS(eData_tissue,"eData_tissue")(Required)=========================
eData_tissue <-readRDS("eData_tissue")

#---------4.1.2 Quality control-----------------
exp_raw <- log2(exprs(eData_tissue)+1)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     title = pData(eData_tissue)$title,
                     group = pData(eData_tissue)$group)

(qplot(PC1, PC2, data = dataGG, color = group,
       main = "PCA plot of the raw data (log-transformed)", size = I(2),
       asp = 1.0, geom = "text",
       label = group)
        + scale_colour_brewer(palette = "Set2"))
boxplot(exp_raw, target = "core",
        main = "Boxplots of log2-intensities for the raw data")
dists <- as.matrix(dist(t(exp_raw), method = "manhattan"))
colnames(dists) <- NULL
diag(dists) <- NA
rownames(dists) <-  pData(eData_tissue)$title
hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)

pheatmap(dists, col = rev(hmcol), clustering_distance_rows = "Manhattan",
         clustering_distance_cols = "Manhattan")

#--------------test spearman correlation(Alternative)---------------------
exp_raw <- log2(exprs(eData_tissue)+1)
c <- cor(exp_raw, method="spearman")
diag(c) <- NA
colnames(c) <- pData(eData_tissue)$group
pheatmap(c,
         col = rev(hmcol),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         fontsize_row = 9,
         fontsize_col = 9,
         fontsize =18,
         main="Test spearman correlation: pure tissue cell types")

#------4.2 Filtering based on correlation(Alternative)--------------
colnames(eData_tissue)
select_samples <- c("Cardiac.Stromal.Cells.2",
                    "Cardiac.Stromal.Cells.3",
                    "Cervical.Ectoepithelial.1",
                    "Cervical.Ectoepithelial.2",
                    "Lung.Fibroblasts.1",
                    "Lung.Fibroblasts.2",
                    "Prostate.Smooth.Muscle.1",
                    "Prostate.Smooth.Muscle.3",
                    "Umbilical.Artery.Endothelial.Cell.3",
                    "Umbilical.Artery.Endothelial.Cell.4"
                    )
selected_samples <- sampleNames(eData_tissue) %in% select_samples
table(selected_samples)
selected_tissue <- eData_tissue[,selected_samples]
sampleNames(selected_tissue)
pData(selected_tissue)$group
selected_tissue_exprs <- log2(exprs(selected_tissue)+1)

#---------4.2.1 test spearman correlation(recommend)---------------------
c <- cor(exp_raw, method="spearman")
diag(c) <- NA
colnames(c) <- pData(selected_tissue)$group
pheatmap(c,
         col = rev(hmcol),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         fontsize_row = 9,
         fontsize_col = 9,
         fontsize =18,
         main="Test spearman correlation: pure tissue cell types")


#======= 4.3 Make Gene.symbol_probe_ensembl annotation(recommend)=================
#http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.htmlannotation(eData_tissue)
platf <- getGEO(annotation(eData_tissue), AnnotGPL=TRUE)
anot <- data.frame(attr(dataTable(platf), "table"))
Gene_Probe <- anot[,c("Gene.symbol","ID")]
Gene_Probe<- Gene_Probe[Gene_Probe$Gene.symbol!="",] #filter empty field
head(Gene_Probe,60)
dim(Gene_Probe)
Gene_Probe_test<- as.data.frame(str_split_fixed(Gene_Probe$Gene.symbol, "///",3)) #split to 3 columns
head(Gene_Probe_test,60)
Gene_Probe$Gene.symbol <-Gene_Probe_test$V1 #replace orignal Gene.symbol
head(Gene_Probe,60)

#===4.3.1 open GSE64655annotation.csv
GSE64655_annotation <-read.csv("GSE64655_annotation.csv",fill = TRUE)
head(GSE64655_annotation)
dim(GSE64655_annotation)

# inner_join to Make Gene.symbol_probe_ensembl 
Gene.symbol_probe_ensembl <- inner_join(Gene_Probe, GSE64655_annotation,by = "Gene.symbol")
head(Gene.symbol_probe_ensembl)
dim(Gene.symbol_probe_ensembl)

# inner_join to make tissue_exprs(Recommend)=========================
tissue_exprs <- log2(exprs(eData_tissue)+1)
class(tissue_exprs)
tissue_exprs <- as.data.frame(tissue_exprs)
head(tissue_exprs)
tissue_exprs$ID <- rownames(tissue_exprs)
head(tissue_exprs)
tissue_exprs_ensembl <-inner_join(tissue_exprs, 
                                           Gene.symbol_probe_ensembl,
                                           by = "ID")
head(tissue_exprs_ensembl)
#remove duplicate
dup <- duplicated(tissue_exprs_ensembl$Gene.ID)
tissue_exprs_ensembl <- tissue_exprs_ensembl[!dup,]
rownames(tissue_exprs_ensembl) <-tissue_exprs_ensembl$Gene.ID
tissue_exprs_ensembl<- tissue_exprs_ensembl[,sampleNames(eData_tissue)]
head(tissue_exprs_ensembl)
class(tissue_exprs_ensembl) #ready for join

# inner_join to make selected_tissue_exprs(Alternative)----------------------------
class(selected_tissue_exprs)
selected_tissue_exprs <- as.data.frame(selected_tissue_exprs)
head(selected_tissue_exprs)
selected_tissue_exprs$ID <- rownames(selected_tissue_exprs)
head(selected_tissue_exprs)
selected_tissue_exprs_ensembl <-inner_join(selected_tissue_exprs, 
                                           Gene.symbol_probe_ensembl,
                                           by = "ID")
head(selected_tissue_exprs_ensembl)
#remove duplicate
dup <- duplicated(selected_tissue_exprs_ensembl$Gene.ID)
selected_tissue_exprs_ensembl <- selected_tissue_exprs_ensembl[!dup,]
rownames(selected_tissue_exprs_ensembl) <-selected_tissue_exprs_ensembl$Gene.ID
selected_tissue_exprs_ensembl<- selected_tissue_exprs_ensembl[,select_samples]
head(selected_tissue_exprs_ensembl)
class(selected_tissue_exprs_ensembl) #ready for join

#----4.3.2 find ensembl ID for multi_gene(test)-----------------------------------------
Gene_Probe <- Gene_Probe[grepl("///", Gene_Probe$Gene.symbol),"Gene.symbol"]
multi_genes<- as.data.frame(str_split_fixed(multi_genes, "///",3)) #split to 3 columns
head(multi_genes)
#run one
colnames(multi_genes)[1] <- "Gene.symbol"
head(multi_genes)
multi_gene_test <- inner_join(multi_genes, GSE64655_annotation,by = "Gene.symbol")
head(multi_gene_test)

#run two
colnames(multi_genes)[1] <- "V1"
colnames(multi_genes)[2] <- "Gene.symbol"
head(multi_genes)
multi_gene_test <- inner_join(multi_genes, GSE64655_annotation,by = "Gene.symbol")
head(multi_gene_test)

#run two
colnames(multi_genes)[2] <- "V2"
colnames(multi_genes)[3] <- "Gene.symbol"
head(multi_genes)
multi_gene_test <- inner_join(multi_genes, GSE64655_annotation,by = "Gene.symbol") 
head(multi_gene_test)




# ===========4.5 load Immune cell RNA-seq (Required)=======================
#use excel to remove header
##replace space with _, especially for column names
ImmuneCell_raw <- read.table("GSE64655_Normalized_transcript_expression_in_human_immune_cells.txt")
rownames(ImmuneCell_raw) <-ImmuneCell_raw[,1]
colnames(ImmuneCell_raw) <-sapply(ImmuneCell_raw[1,],as.character,drop=FALSE) #as.character a list
ImmuneCell_raw <- ImmuneCell_raw[-1,]
ImmuneCell_raw <- ImmuneCell_raw[,-1]
head(ImmuneCell_raw[,1:4])

# subset the d0 samples
day.zero <- grepl("d0",colnames(ImmuneCell_raw))
colnames(ImmuneCell_raw[,day.zero])
ImmuneCell_exprs <- ImmuneCell_raw[,day.zero]

HD3_name <-sub("_d0$", "", colnames(ImmuneCell_exprs))  #remove "_d0" at end
HD3_name
colnames(ImmuneCell_exprs) <- sub("^.*_", "", HD3_name) #remove "HD3*_" at begin. ?regex
head(ImmuneCell_exprs[,1:4])
dim(ImmuneCell_exprs)
#========replace column name(Recommend)==================
col_names <- colnames(ImmuneCell_exprs)
col_names
col_names <-sub("^PBMC", "PBMC", col_names)
col_names <-sub("^mDC", "dendritic.cells", col_names)
col_names <-sub("^Mono", "monocytes", col_names)
col_names <-sub("^Neut", "neutrophils", col_names)
col_names <-sub("^B", "B.cells", col_names)
col_names <-sub("^T", "T.cells", col_names)
col_names <-sub("^NK", "natural.killer", col_names)
col_names
colnames(ImmuneCell_exprs) <- col_names
head(ImmuneCell_exprs[,1:4])
class(ImmuneCell_exprs)  #ready for join

#-----test spearman correlation for pure imune cell(Alternative)----------------
# convert data.frame to matrix
is.matrix(ImmuneCell_exprs)
ImmuneCell_exprs.char <- as.matrix(ImmuneCell_exprs)
head(ImmuneCell_exprs.char[,1:4])
class(ImmuneCell_exprs.char[1,1])

ImmuneCell_exprs.matrix <- apply(ImmuneCell_exprs.char,2,as.numeric)
dim(ImmuneCell_exprs.matrix)
head(ImmuneCell_exprs.matrix[,1:4])
class(ImmuneCell_exprs.matrix)
rownames(ImmuneCell_exprs.matrix) <- rownames(ImmuneCell_exprs)
colnames(ImmuneCell_exprs.matrix) <- colnames(ImmuneCell_exprs)


c <- cor(ImmuneCell_exprs.matrix, method="spearman")
diag(c) <- NA
colnames(c) <- NULL
hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)
pheatmap(c, col = rev(hmcol), clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         main ="test spearman correlation for pure immune cell type")


# ===========4.6 Load SingleCell data (Required)=======================
HSMM <-readRDS("HSMM")
table(pData(HSMM)$CellType)
# filter HSMM 
table(pData(HSMM)$Total_mRNAs >5500)
HSMM <-HSMM[,pData(HSMM)$Total_mRNAs >5500]
# order HSMM by celltype
HSMM <-HSMM[,order(pData(HSMM)$CellType)]
# rename samples
sampleNames(HSMM) <-make.unique(as.character(pData(HSMM)$CellType))
sampleNames(HSMM)
head(exprs(HSMM)[,1:6])
# insert Ensembl column to the 1st
SingleCell_exprs <- as.data.frame(exprs(HSMM))
head(SingleCell_exprs[,1:5])
class(SingleCell_exprs) #ready for join
# ===========4.7 test spearman correlation for SingleCell vs immune cells vs tissue cell (Required)=======================


#==========4.7.1 insert Ensembl column to the 1st
tissue  <- tissue_exprs_ensembl
ImmuneCell <- ImmuneCell_exprs
SingleCell <- SingleCell_exprs

tissue$Ensembl <- rownames(tissue) # add Ensembl 
tissue <- tissue[,c("Ensembl",colnames(tissue_exprs_ensembl))] 
head(tissue[,1:5])

ImmuneCell$Ensembl <- rownames(ImmuneCell) #add Ensembl  
ImmuneCell <- ImmuneCell[,c("Ensembl",colnames(ImmuneCell_exprs))] 
head(ImmuneCell)

SingleCell$Ensembl <- rownames(SingleCell) #add Ensembl  
SingleCell <- SingleCell[,c("Ensembl",colnames(SingleCell_exprs))] 
head(SingleCell[,1:5])

#=======merge expression profile===============
tissue_Immune <- inner_join(tissue, ImmuneCell,by = "Ensembl") 
dim(tissue_Immune)
head(tissue_Immune[,1:5])

tissue_Immune_Single <- inner_join(tissue_Immune, SingleCell,by = "Ensembl") 
dim(tissue_Immune_Single)
head(tissue_Immune_Single[,1:5])

rownames(tissue_Immune_Single) <-tissue_Immune_Single$Ensembl
tissue_Immune_Single <- tissue_Immune_Single[,-1]
head(tissue_Immune_Single[,1:5])

# convert data.frame to matrix==================
is.matrix(tissue_Immune_Single)
tissue_Immune_Single.char <- as.matrix(tissue_Immune_Single)
head(tissue_Immune_Single.char[,1:4])
class(tissue_Immune_Single.char[1,1])

tissue_Immune_Single.matrix <- apply(tissue_Immune_Single.char,2,as.numeric)
dim(tissue_Immune_Single.matrix)
head(tissue_Immune_Single.matrix[,1:4])
class(tissue_Immune_Single.matrix)
rownames(tissue_Immune_Single.matrix) <- rownames(tissue_Immune_Single)
colnames(tissue_Immune_Single.matrix) <- colnames(tissue_Immune_Single)
is.matrix(tissue_Immune_Single.matrix)


c <- cor(tissue_Immune_Single.matrix, method="spearman")
diag(c) <- NA
colnames(c) <- NULL
rownames(c) <- NULL
hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)

pheatmap(c,
         col = rev(hmcol),
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         #gaps_row=gaps_row$gap,
         #         gaps_col=gaps_col$gap,
         #         cellheight = 7,
         #         cellwidth = 10,
         #         border_color=NA,
         fontsize_row = 8.5,
         fontsize_col = 8.5,
         fontsize =15,
         main="Test spearman correlation: Single Cell vs pure cell types"
         #         filename = "TEST_12cat.png",
         #annotation_row = annotdf,
         #annotation_colors = mycolors
)

# pheatmap small size-----------
#select columns
n_col <-ncol(SingleCell_exprs):ncol(tissue_Immune_Single)
colnames(c)[n_col]
#select rows
n_row <-1:ncol(exprs(HSMM))
rownames(c)[n_row]

pheatmap(c[n_row,n_col],
         col = rev(hmcol),
         clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan",
         #gaps_row=gaps_row$gap,
         #         gaps_col=gaps_col$gap,
         #         cellheight = 7,
         #         cellwidth = 10,
         #         border_color=NA,
         fontsize_row = 8.5,
         fontsize_col = 15,
         fontsize =15,
         main="Test spearman correlation: Single Cell vs pure cell types"
         #         filename = "TEST_12cat.png",
         #annotation_row = annotdf,
         #annotation_colors = mycolors
)
