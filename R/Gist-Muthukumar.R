#generate table in R markdown
cellType.table <- read.delim("CellType.txt")
kable(cellType.table, captioin ="Cell Type")


#generate table in R
x <- data.frame(row.names=paste("Name",1:10))
x[,1] <- 1:10
x[,2] <- sample(1:100,10)
x[,3] <- sample(LETTERS[1:26],10)
colnames(x) <- c("Value 1", "Value 2", "Label")

#Plot your table with table Grob in the library(gridExtra)
ss <- tableGrob(x)

#Make a scatterplot of your data
k <- ggplot(x,aes(x=x$"Value 1",y=x$"Value 2")) + 
        geom_point()

#Arrange them as you want with grid.arrange
grid.arrange(ss)

#2017-07-14
#try NK cells and other myloide cells

#CD3G(-),NCAM1(+),KLRD1 (+),NCR1(+)	, CD3-, CD56+, CD94+, NKp46+
        
#assaign each cluster to cell type

#2017-07-18
#FlowJo pineer in single cell level clustering using machine learning
#k-means clustering, tsne 
#http://www.the-scientist.com/?articles.view/articleNo/49228/title/Picking-Out-Patterns/

#https://followthedata.wordpress.com/2015/12/21/list-of-deep-learning-implementations-in-biology/
ADAGE #Analysis using Denoising Autoencoders for Gene Expression 
#https://github.com/greenelab/adage

