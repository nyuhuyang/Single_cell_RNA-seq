#Identification of important genes
########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages <- c("pheatmap","dplyr")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("monocle","geneplotter")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
library(easypackages)
libraries(list.of.bio.packages,list.of.cran.packages)


# #####################################################################
# 
#  3.5 Classifying and counting cells of different types
#       http://www.proteinatlas.org
#  https://bioinfo.uth.edu/scrnaseqdb/index.php
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

HSMM <-readRDS("HSMM")
table(pData(HSMM)$CellType)
# filter HSMM 
table(pData(HSMM)$Total_mRNAs >5300)
HSMM <-HSMM[,pData(HSMM)$Total_mRNAs >5300]
# order HSMM by celltype
HSMM <-HSMM[,order(pData(HSMM)$CellType)]
# rename samples
sampleNames(HSMM) <-make.unique(as.character(pData(HSMM)$CellType))

head(pData(HSMM)$CellType)
head(exprs(HSMM)[,1:3])

#======3.5.1 read files from proteinatlas.org (Required)============
#search kidney tissue specific cell type, download TAB files
#if tab files can't be downloaded, chose genes in frist two pages manually
#creat TAB file with excel, including delete blank rows
tab.files <-list.files(pattern = "\\.tab$") #list all files with a specified extension
tab.files
n <- length(tab.files)
n
celltype <- gsub('.tab', '',tab.files)  #"B_Cells.tab" => "B_Cells"
celltype
select_columns <- c("Gene","Ensembl","RNA.tissue.category")

tabfileslist <- list() #creat a empty list

#create tab files list
for(i in 1:n){
        if(i==1){nrow=0;print(paste0("nrow = ",nrow))} #only for the first run
        tab <- tab.files[i] #get file name
        tabfileslist[[tab]] <- read.delim(tab) #store variable to list
        #select 100 genes and 4 columns
        tabfileslist[[tab]] <- tabfileslist[[tab]][1:20,select_columns]
        #add priority number: prioties the best marker
        tabfileslist[[tab]]$priority <- rownames(tabfileslist[[tab]])
        #filter out genes with RNA.TS score not greater than 0, remove RNA.Ts=NA
#        tabfileslist[[tab]] <-subset(tabfileslist[[tab]],
#                                     !is.na(tabfileslist[[tab]]$RNA.TS))
        #add CellType to category
        tabfileslist[[tab]]$CellType <-rep(celltype[i],
                                               length(tabfileslist[[tab]]$Gene))
        nrow <-nrow + length(tabfileslist[[tab]]$CellType)
        print(dim(tabfileslist[[tab]]))
        print(paste0("nrow = ",nrow))
}

#convert list to data frame

Gene_Ensembl_celltype <- do.call(rbind.data.frame, tabfileslist) 
Gene_Ensembl_celltype <- droplevels(Gene_Ensembl_celltype)

#Granulocyte and monocytes do have 100 marker genes
Gene_Ensembl_celltype <- Gene_Ensembl_celltype[complete.cases(Gene_Ensembl_celltype),]
dim(Gene_Ensembl_celltype)
head(Gene_Ensembl_celltype[,1:5])

#make duplicated gene name unique, by adding .number
#/   Gene<- make.unique(as.character(Gene_Ensembl_celltype$Gene), sep = "_")
#remove duplicated gene name with lower priority number 
dup <- duplicated(Gene_Ensembl_celltype$Gene)
dup_gene <-droplevels(unique(Gene_Ensembl_celltype[dup,"Gene"]))
dup_gene_priority <- Gene_Ensembl_celltype[Gene_Ensembl_celltype$Gene %in% dup_gene,
                                  c("Gene","priority")]
dup_gene_priority <- dup_gene_priority[order(as.numeric(dup_gene_priority$priority)),
                                             c("Gene","priority")]
dup_gene_priority <- dup_gene_priority[order(dup_gene_priority$Gene),
                                       c("Gene","priority")]
dim(dup_gene_priority)
dup_gene_priority <- dup_gene_priority[duplicated(dup_gene_priority$Gene),
                                       c("Gene","priority")] #doesn't mark the first duplicates
dim(dup_gene_priority)
unique.gene <- !(rownames(Gene_Ensembl_celltype) %in% rownames(dup_gene_priority))
Gene_Ensembl_celltype <- Gene_Ensembl_celltype[unique.gene,]
dim(Gene_Ensembl_celltype)
head(Gene_Ensembl_celltype[,1:5])

#=====3.5.2 pheatmap with certain celltype===============
"Make matrix as 'Gene_exprs' dim=(170,50)"
#"nrow = 170"
#table(pData(HSMM)$CellType=="Ambiguous")=50
# First make matrix as "Gene_Ensembl" dim=(155,1)
Gene_Ensembl <-data.frame(row.names = Gene_Ensembl_celltype$Gene,
                          Ensembl = Gene_Ensembl_celltype$Ensembl,
                          priority = Gene_Ensembl_celltype$priority)
head(Gene_Ensembl)
dim(Gene_Ensembl)

# Then make matrix "Ensembl_exprs" dim=(44,51)

Ensembl_exprs <-subset(exprs(HSMM),rownames(exprs(HSMM)) %in% Gene_Ensembl$Ensembl)
Ensembl_exprs <- as.data.frame(Ensembl_exprs)
Ensembl_exprs$Ensembl <- rownames(Ensembl_exprs)
Ensembl_exprs <- Ensembl_exprs[,c("Ensembl",colnames(exprs(HSMM)))]
dim(Ensembl_exprs)
head(Ensembl_exprs[,1:3])

# make matrix Gene_Ensembl_exprs dim=(155,51) then Gene_exprs dim=(~,50)

Gene_Ensembl_exprs <- dplyr::left_join(Gene_Ensembl, Ensembl_exprs,by = "Ensembl") #didn't use inner_join, need to name the rows
#Warning message:
#Column `Ensembl` joining factor and character vector, coercing into character vector 
all(Gene_Ensembl_exprs$Ensembl == Gene_Ensembl$Ensembl)
rownames(Gene_Ensembl_exprs) <- rownames(Gene_Ensembl)
dim(Gene_Ensembl_exprs)
head(Gene_Ensembl_exprs[,1:3])
Gene_Ensembl_exprs<-Gene_Ensembl_exprs[complete.cases(Gene_Ensembl_exprs),]
Gene_exprs <- Gene_Ensembl_exprs[,-c(1:2)]
dim(Gene_exprs)
head(Gene_exprs[,1:3])

# make annotdf
#Gene_Ensembl_celltype$Gene <- Gene # replace Gene with Gene.number
tail(Gene_Ensembl_celltype)
Annotdf <- data.frame(row.names = Gene_Ensembl_celltype$Gene, 
                      name = Gene_Ensembl_celltype$Gene, #have to add one more column to keep as dataframe
                      category = Gene_Ensembl_celltype$CellType )  
Annotdf <- Annotdf[rownames(Annotdf) %in% rownames(Gene_exprs),] #remove genes with NA exprs
annotdf <- data.frame(row.names = rownames(Annotdf), 
                      category =Annotdf$category) #have to create it twice?
                      
dim(annotdf)
head(annotdf)

#===="pheatmap with Ambiguous celltype"=====
#color
newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotdf$category))))
mycolors <- newCols(length(unique(annotdf$category)))
names(mycolors) <- unique(annotdf$category)
mycolors <- list(category = mycolors)

#calculate gap row===
gaps_row <- as.data.frame(table(annotdf))
gaps_row
gaps_row$gap <-NA #insert empty column
gaps_row$gap[1] <-gaps_row$Freq[1]
for(i in 2:length(gaps_row$gap)){
        gaps_row$gap[i] <- gaps_row$gap[i-1] + gaps_row$Freq[i]}
gaps_row

#calculate gap row===
gaps_col <- as.data.frame(table(pData(HSMM)$CellType))
gaps_col
gaps_col <- gaps_col[gaps_col$Freq!=0,]
gaps_col$gap <-NA #insert empty column
gaps_col$gap[1] <-gaps_col$Freq[1]
for(i in 2:length(gaps_col$gap)){
        gaps_col$gap[i] <- gaps_col$gap[i-1] + gaps_col$Freq[i]}
gaps_col

pheatmap(Gene_exprs,
#         scale="row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row=gaps_row$gap,
#         gaps_col=gaps_col$gap,
#         cellheight = 7,
#         cellwidth = 10,
#         border_color=NA,
         fontsize_row = 5,
         main="Grouped Cells with marker genes",
         #         filename = "TEST_12cat.png",
         annotation_row = annotdf,
         annotation_colors = mycolors
)

