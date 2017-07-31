#Identification of important genes
########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages <- c()
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

HSMM_Amb <- HSMM[,pData(HSMM)$CellType == "Ambiguous"]
HSMM_Epi <- HSMM[,pData(HSMM)$CellType == "Epithelial Cells"]
HSMM_Fib <- HSMM[,pData(HSMM)$CellType == "Fibroblasts"]
HSMM_Gra <- HSMM[,pData(HSMM)$CellType == "Granulocyte"]
HSMM_Leu <- HSMM[,pData(HSMM)$CellType == "Leukocytes"]
HSMM_Mon <- HSMM[,pData(HSMM)$CellType == "Monocytes"]
HSMM_Str <- HSMM[,pData(HSMM)$CellType == "Stromal Cells+\nEndothelial Cells+\nFibroblasts"]
HSMM_T <- HSMM[,pData(HSMM)$CellType == "T cells"]
HSMM_Unk <- HSMM[,pData(HSMM)$CellType == "Unknown"]

head(pData(HSMM_Amb)$CellType)
head(exprs(HSMM_Amb)[,1:3])

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
select_columns <- c("Gene","Ensembl","RNA.tissue.category","RNA.TS")

tabfileslist <- list() #creat a empty list

#create tab files list
for(i in 1:n){
        if(i==1){nrow=0;print(paste0("nrow = ",nrow))} #only for the first run
        tab <- tab.files[i] #get file name
        tabfileslist[[tab]] <- read.delim(tab) #store variable to list
        #select 100 genes and 4 columns
        tabfileslist[[tab]] <- tabfileslist[[tab]][1:100,select_columns]
        #filter out genes with RNA.TS score not greater than 0
        tabfileslist[[tab]] <-subset(tabfileslist[[tab]],
                                     tabfileslist[[tab]][,"RNA.TS"]>0)
        #add CellType to category
        tabfileslist[[tab]][,"CellType"] <-rep(celltype[i],
                                               length(tabfileslist[[tab]][,"Gene"]))
        nrow <-nrow + length(tabfileslist[[tab]][,"CellType"])
        print(dim(tabfileslist[[tab]]))
        print(paste0("nrow = ",nrow))
}

#=====3.5.2 pheatmap with grouped celltype===============
myexprs <- do.call(rbind.data.frame, tabfileslist) #unlist to data frame

#make duplicated gene name unique
Gene<- make.unique(as.character(myexprs$Gene), sep = ".")

Gene_Ensembl <-data.frame(row.names = Gene,
                          Ensembl = myexprs$Ensembl)

exprs.Amb <-subset(exprs(HSMM_Amb),rownames(exprs(HSMM_Amb)) %in% myexprs$Ensembl)
dim(exprs.Amb)
myexprs <- merge.data.frame(Gene_Ensembl,exprs.Amb,
                            by.x =  Gene_Ensembl$Ensembl)

annotdf <- data.frame(row.names = Gene, 
                      category = myexprs$CellType )  



#-------------------
mymat <- matrix(rexp(720, rate=.1), ncol=12)
colnames(mymat) <- c(rep("treatment_1", 3), rep("treatment_2", 3), rep("treatment_3", 3), rep("treatment_4", 3))
rownames(mymat) <- paste("gene", 1:dim(mymat)[1], sep="_")

annotdf <- data.frame(row.names = rownames(mymat), 
                      category = rep(paste0("Category_", seq(12)), each=5) )  

newCols <- colorRampPalette(grDevices::rainbow(length(unique(annotdf$category))))
mycolors <- newCols(length(unique(annotdf$category)))
names(mycolors) <- unique(annotdf$category)
mycolors <- list(category = mycolors)

pheatmap(mymat,
         scale="row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
         gaps_col=c(3,6,9),
         cellheight = 6,
         cellwidth = 20,
         border_color=NA,
         fontsize_row = 6,
         main="Genes grouped by categories",
         #         filename = "TEST_12cat.png",
         annotation_row = annotdf,
         annotation_colors = mycolors
)

