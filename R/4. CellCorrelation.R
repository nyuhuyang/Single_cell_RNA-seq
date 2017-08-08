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
sampleNames <- sub("Epithelial$", "Epithelials", sampleNames)
sampleNames <- sub("Endothelial$", "Endothelials", sampleNames) 
sampleNames <- sub("Muscle$", "Muscles", sampleNames) 
sampleNames <- sub("Stromal$", "Stromals", sampleNames)
table(sampleNames)
pData(eData_tissue)$group <- sampleNames
pData(eData_tissue)$group
#----------regroup eData_tissue in alphabetic order(Alternative)-----
eData_Endo <- eData_tissue[,pData(eData_tissue)$group=="Endothelials"]
eData_Epit <- eData_tissue[,pData(eData_tissue)$group=="Epithelials"]
eData_Fibr <- eData_tissue[,pData(eData_tissue)$group=="Fibroblasts"]
eData_Musc <- eData_tissue[,pData(eData_tissue)$group=="Muscles"]
eData_Stro <- eData_tissue[,pData(eData_tissue)$group=="Stromals"]

#prepare files for the new expressoinset
eData_names <- list(names=c(eData_Endo,eData_Epit,eData_Fibr,
                            eData_Musc,eData_Stro))
pData(eData_names$names[[1]])$group
group <- lapply(eData_names$names,function(x) pData(x)$group)
group <- unlist(group)

pData(eData_names$names[[1]])$title
title <- lapply(eData_names$names,function(x) pData(x)$title)
title <- unlist(title)

sampleNames(eData_names$names[[1]])
sampleNames <- lapply(eData_names$names,function(x) sampleNames(x))
sampleNames <- unlist(sampleNames)
sampleNames <- as.character(sampleNames)

head(exprs(eData_names$names[[1]]))
All.exprs <- lapply(eData_names$names,function(x) exprs(x))
All.exprs <- do.call(cbind,All.exprs)
class(All.exprs)
head(All.exprs)

#----Construct the new ExpressoinSet(Alternative)------------------------
samples <- data.frame(title = title,group = group)
rownames(samples) <-sampleNames

pd <- new("AnnotatedDataFrame",data=samples)
head(pd)

annotation(eData_tissue)
tissue_new <- new("ExpressionSet", exprs = All.exprs, 
                  phenoData = pd,
                  annotation = annotation(eData_tissue))
tissue_new
saveRDS(tissue_new,"tissue_new")
#=======4.1.1 saveRDS(eData_tissue,"eData_tissue")(Required)=========================
#eData_tissue <-readRDS("eData_tissue")
eData_tissue <-readRDS("tissue_new")

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
write.csv(ImmuneCell_exprs,"ImmuneCell_exprs.csv")
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
# order HSMM by celltype
HSMM <-HSMM[,order(pData(HSMM)$CellType)]
# rename samples
sampleNames(HSMM) <-make.unique(as.character(pData(HSMM)$CellType))
sampleNames(HSMM)
head(exprs(HSMM)[,1:6])
# prepare expression data.frame
SingleCell.all_exprs <- as.data.frame(exprs(HSMM))
head(SingleCell.all_exprs[,1:5])
class(SingleCell.all_exprs) #ready for join

#-----for part of single cell--------------
# filter HSMM --------------
table(pData(HSMM)$Total_mRNAs >5500)
HSMM <-HSMM[,pData(HSMM)$Total_mRNAs >5500]
# order HSMM by celltype
HSMM <-HSMM[,order(pData(HSMM)$CellType)]
# rename samples
sampleNames(HSMM) <-make.unique(as.character(pData(HSMM)$CellType))
sampleNames(HSMM)
head(exprs(HSMM)[,1:6])
# prepare expression data.frame
SingleCell_exprs <- as.data.frame(exprs(HSMM))
head(SingleCell_exprs[,1:5])
class(SingleCell_exprs) #ready for join


# ===========4.7 test spearman correlation for SingleCell vs immune cells vs tissue cell (Required)=======================


#==========4.7.1 insert Ensembl column to the 1st==========
tissue  <- tissue_exprs_ensembl
ImmuneCell <- ImmuneCell_exprs
SingleCell <- SingleCell_exprs
SingleCell.all <- SingleCell.all_exprs

tissue$Ensembl <- rownames(tissue) # add Ensembl 
tissue <- tissue[,c("Ensembl",colnames(tissue_exprs_ensembl))] 
head(tissue[,1:5])

ImmuneCell$Ensembl <- rownames(ImmuneCell) #add Ensembl  
ImmuneCell <- ImmuneCell[,c("Ensembl",colnames(ImmuneCell_exprs))] 
head(ImmuneCell)
#-----for part of single cell--------------
SingleCell$Ensembl <- rownames(SingleCell) #add Ensembl  
SingleCell <- SingleCell[,c("Ensembl",colnames(SingleCell_exprs))] 
head(SingleCell[,1:5])

#-----for all single cell--------------
SingleCell.all$Ensembl <- rownames(SingleCell.all) #add Ensembl  
SingleCell.all <- SingleCell.all[,c("Ensembl",colnames(SingleCell.all_exprs))] 
head(SingleCell.all[,1:5])

#=======4.7.2 merge expression profile===============
tissue_Immune <- inner_join(tissue, ImmuneCell,by = "Ensembl") 
dim(tissue_Immune)
head(tissue_Immune,5)

#-----for part of single cell--------------
tissue_Immune_Single <- inner_join(tissue_Immune, SingleCell,by = "Ensembl") 
dim(tissue_Immune_Single)
head(tissue_Immune_Single[,1:5])

rownames(tissue_Immune_Single) <-tissue_Immune_Single$Ensembl
tissue_Immune_Single <- tissue_Immune_Single[,-1]
head(tissue_Immune_Single[,1:5])

#-----for all single cell--------------
tissue_Immune_Single.all <- inner_join(tissue_Immune, SingleCell.all,by = "Ensembl") 
dim(tissue_Immune_Single.all)
head(tissue_Immune_Single.all[,1:5])

rownames(tissue_Immune_Single.all) <-tissue_Immune_Single.all$Ensembl
tissue_Immune_Single.all <- tissue_Immune_Single.all[,-1]
head(tissue_Immune_Single.all[,1:5])

# convert data.frame to matrix==================
#-----for all single cell--------------
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

saveRDS(tissue_Immune_Single.matrix,"tissue_Immune_Single")

#-----for all single cell--------------
is.matrix(tissue_Immune_Single.all)
tissue_Immune_Single.all.char <- as.matrix(tissue_Immune_Single.all)
head(tissue_Immune_Single.all.char[,1:4])
class(tissue_Immune_Single.all.char[1,1])

tissue_Immune_Single.all.matrix <- apply(tissue_Immune_Single.all.char,2,as.numeric)
dim(tissue_Immune_Single.all.matrix)
head(tissue_Immune_Single.all.matrix[,1:4])
class(tissue_Immune_Single.all.matrix)
rownames(tissue_Immune_Single.all.matrix) <- rownames(tissue_Immune_Single.all)
colnames(tissue_Immune_Single.all.matrix) <- colnames(tissue_Immune_Single.all)
is.matrix(tissue_Immune_Single.all.matrix)

saveRDS(tissue_Immune_Single.all.matrix,"tissue_Immune_Single.all")
#====4.7.3 spearman correlation(required)===============
#-----for part of single cell--------------
tissue_Immune_Single.matrix <- readRDS("tissue_Immune_Single")
c <- cor(tissue_Immune_Single.matrix, method="spearman")
#-----for all single cell--------------
tissue_Immune_Single.matrix <- readRDS("tissue_Immune_Single.all")
c <- cor(tissue_Immune_Single.all.matrix, method="spearman")

#=====4.8 make grouped pheatmap=================================

tissue_new <-readRDS("tissue_new")
ImmuneCell_exprs <- read.csv("ImmuneCell_exprs.csv",row.names = 1)
#========make annotdf==========

category <- c(as.character(pData(tissue_new)$group),
              rep("Immunes",ncol(ImmuneCell_exprs)),
              rep("scRNA-seq",ncol(SingleCell_exprs)))

Annotdf <- data.frame(row.names = colnames(tissue_Immune_Single.matrix), 
                      category = category )  
dim(Annotdf)
head(Annotdf)

#color
newCols <- colorRampPalette(grDevices::rainbow(length(unique(Annotdf$category))))
mycolors <- newCols(length(unique(Annotdf$category)))
names(mycolors) <- unique(Annotdf$category)
mycolors <- list(category = mycolors)

#=========make gap_row============
gaps_row <- as.data.frame(table(Annotdf))
gaps_row
gaps_row <- gaps_row[c(1,2,3,5,7,4,6),] #switch sequence
gaps_row$gap <-NA #insert empty column
gaps_row$gap[1] <-gaps_row$Freq[1]
for(i in 2:length(gaps_row$gap)){
        gaps_row$gap[i] <- gaps_row$gap[i-1] + gaps_row$Freq[i]}
gaps_row


#==========make pheatmap============
diag(c) <- NA
hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)
pheatmap(c,
         col = rev(hmcol),
         cluster_rows = FALSE,
         cluster_cols = FALSE,         
         #clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         gaps_row=gaps_row$gap,
         #         gaps_col=gaps$col,
         #         cellheight = 7,
         #         cellwidth = 10,
         #border_color=NA,
         fontsize_row = 8,
         #fontsize_col = 8,
         fontsize =18,
         show_colnames = F,
         main="Spearman coefficient: Single Cell vs pure cell types",
         #         filename = "TEST_12cat.png",
         annotation_row = Annotdf,
         annotation_colors = mycolors
)

# pheatmap small size=======
#select columns
n_col <-(ncol(tissue_new)+ncol(ImmuneCell_exprs)+1):ncol(tissue_Immune_Single)
colnames(c)[n_col]
#select rows
tissue_row <-1:(ncol(tissue_new))
rownames(c)[tissue_row]

Immune_row <-(ncol(tissue_new)+1):(ncol(tissue_new)+ncol(ImmuneCell)-1)
rownames(c)[Immune_row]

tissue_Immune_row <-1:(ncol(tissue_new)+ncol(ImmuneCell)-1)
rownames(c)[tissue_Immune_row]

pheatmap(c[tissue_Immune_row,n_col],
#         scale="column",
         col = rev(hmcol),
         cluster_rows = FALSE,
         cluster_cols = FALSE,         
         #clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         gaps_row=gaps_row$gap[1:5],
         #         cellheight = 7,
         #         cellwidth = 10,
         #border_color=NA,
         fontsize_row = 8,
         fontsize_col = 8,
         fontsize =15,
         #show_colnames = T,
         main="Spearman coefficient: Single Cell vs pure cell types",
         #         filename = "TEST_12cat.png",
         annotation_row = Annotdf,
         annotation_colors = mycolors
)

#=====4.9 cell cluster=================================
#replace row cell names
ImmuneCell.names <- colnames(ImmuneCell_exprs)
ImmuneCell.names <- sub(".\\d$", "", ImmuneCell.names) #delete .number

#-----for part of single cell--------------
category.2 <- c(as.character(pData(tissue_new)$group),
              ImmuneCell.names,
              rep("scRNA-seq",ncol(SingleCell_exprs)))

rownames(c) <-category.2
#Find the maximum position for each col of a matrix
Max.coe <-max.col(t(c[tissue_Immune_row,n_col]))
rownames(c[tissue_Immune_row,n_col])[Max.coe]
cell.cluster <- data.frame(Monocle=colnames(c[tissue_Immune_row,n_col]),
                           Spearman=rownames(c[tissue_Immune_row,n_col])[Max.coe])
head(cell.cluster)
ss <- gridExtra::tableGrob(cell.cluster[1:28,])
gridExtra::grid.arrange(ss) #make table
ss <- gridExtra::tableGrob(cell.cluster[29:56,])
gridExtra::grid.arrange(ss) #make table
ss <- gridExtra::tableGrob(cell.cluster[57:67,])
gridExtra::grid.arrange(ss) #make table

#-----for all single cell--------------
category.3 <- c(as.character(pData(tissue_new)$group),
                ImmuneCell.names,
                rep("scRNA-seq",ncol(SingleCell.all_exprs)))

rownames(c) <-category.3
#select columns
single.all_col <-(ncol(tissue_new)+ncol(ImmuneCell_exprs)+1):ncol(tissue_Immune_Single.all)
#Find the maximum position for each col of a matrix
Max.coe <-max.col(t(c[tissue_Immune_row,single.all_col]))
rownames(c[tissue_Immune_row,n_col])[Max.coe]
cell.cluster <- data.frame(Monocle=as.character(pData(HSMM)$CellType),
                           Spearman=rownames(c[tissue_Immune_row,single.all_col])[Max.coe])
dim(cell.cluster)
head(cell.cluster)
tail(cell.cluster)
table(cell.cluster$Monocle)
table(cell.cluster$Spearman)
write.csv(cell.cluster,"cell.cluster.csv")
