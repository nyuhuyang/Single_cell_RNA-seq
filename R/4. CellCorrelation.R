########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages <- c("easypackages","stringr","pheatmap","dplyr")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("oligo","AnnotationDbi", "impute","GEOquery",
                          "monocle","geneplotter")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
library(easypackages)
libraries(list.of.bio.packages,list.of.cran.packages)


# #####################################################################
# 
#  4. Run correlation test of different cell types
#       http://www.ebi.ac.uk
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

# ===========4.1 load processed data (Recommend)=======================
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

eList <- getGEOSuppFiles("GSE3239")
eList
list.files("GSE3239")
untar("GSE3239/GSE3239_RAW.tar", exdir = "GSE3239/CEL")
list.files("GSE3239/CEL")
celfiles <- list.files("GSE3239/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)
rawData
head(exprs(rawData))
head(pData(rawData))

# clean up the phenotype information
filename <- sampleNames(eData_tissue)
pData(eData_tissue)$filename <- filename
sampleNames <- gsub(" ", ".", filename) #replace space to .
sampleNames(eData_tissue) <- sampleNames
sampleNames <- sub(".\\d$", "", sampleNames) #delete .number
sampleNames <- sub(" Cel.", "", sampleNames) #delete "cell" and "cells"
sampleNames <- sub(" Pooled", "", sampleNames) #delete "Pooled"
sampleNames <- sub('^.* ([[:alnum:]]+)$', '\\1', sampleNames) #keep last word
sampleNames <- sub("Ectoepithelial$", "Epithelial", sampleNames) #delete "Ecto"
sampleNames <- sub("Endothelial$", "Endothelials", sampleNames) 
sampleNames <- sub("Muscle$", "Muscles", sampleNames) 
sampleNames <- sub("Stromal$", "Stromals", sampleNames)
table(sampleNames)
pData(eData_tissue)$group <- sampleNames
pData(eData_tissue)$group
saveRDS(eData_tissue,"eData_tissue")
eData_tissue <-readRDS('eData_tissue')
