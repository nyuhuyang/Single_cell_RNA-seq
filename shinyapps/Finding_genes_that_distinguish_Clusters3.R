source("https://bioconductor.org/biocLite.R")
list.of.packages <- c("GenomeInfoDb","DESeq2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

library("shiny")
library("DESeq2")
library("ggplot2")





dds<-readRDS("dds")

d<-data.frame(x= str(0), y= integer(0)) # create a empty data.frame
#Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : 
#invalid graphics state-----use below command


server <-function(input, output) {
    output$plot <- renderPlot({


        if(input$Group== "Cluster2 and Cluster1") {
            res_cluster_2vs1 <-read.csv("DESeq_result_cluster_2vs1.csv",row.names = 1)
            for( i in input$range[1]:(input$range[2])){
            d1 <- plotCounts(dds, gene=rownames(res_cluster_2vs1)[i], intgroup="Cluster",returnData=T)
            d1$gene_short_name<-rep(rownames(res_cluster_2vs1)[i],nrow(d1)) # add one more column for facet_warp
            d <- rbind(d,d1)}
        }
        
        if(input$Group== "Cluster3 and Cluster2") {
            res_cluster_3vs2 <-read.csv("DESeq_result_cluster_3vs2.csv",row.names = 1)
            for( i in input$range[1]:(input$range[2])){
                d1 <- plotCounts(dds, gene=rownames(res_cluster_3vs2)[i], intgroup="Cluster",returnData=T)
                d1$gene_short_name<-rep(rownames(res_cluster_3vs2)[i],nrow(d1)) # add one more column for facet_warp
                d <- rbind(d,d1)}
        }
        
        if(input$Group== "Cluster3 and Cluster1") {
            res_cluster_3vs1 <-read.csv("DESeq_result_cluster_3vs1.csv",row.names = 1)
            for( i in input$range[1]:(input$range[2])){
                d1 <- plotCounts(dds, gene=rownames(res_cluster_3vs1)[i], intgroup="Cluster",returnData=T)
                d1$gene_short_name<-rep(rownames(res_cluster_3vs1)[i],nrow(d1)) # add one more column for facet_warp
                d <- rbind(d,d1)}
        }
        p =  ggplot(d, aes(x=Cluster, y=count)) 

        
        p<- p+facet_wrap(~gene_short_name)+
            geom_point(position=position_jitter(w=0.2,h=0),size=3)+
            ggtitle(isolate(input$title))+
            theme(plot.title = element_text(hjust = 0.5,size = 40, face = "bold"),
                  axis.text=element_text(size=20),
                  strip.text = element_text(size=25),
                  axis.title=element_text(size=40, colour = "black"))+
            coord_trans(y = "log10")
        

            print(p)
    }, height=1000)
    
}
ui <- fluidPage(
    titlePanel("Monocle parameters"),
    
    textInput(inputId = "title", 
              label = "Write a title",
              value = "Finding differentially expressed genes"),
    sidebarPanel(
        
        helpText("Create ggplot with gene expression information from DESeq:"),
        
        sliderInput(inputId = "range", 
                    label = "Range of interest:", 
                    value = c(1,9), step=9, min = 1, max = 450),
        selectInput('Group', 'DESeq between:', c("Cluster2 and Cluster1",
                                                 "Cluster3 and Cluster2",
                                                 "Cluster3 and Cluster1"))
#        actionButton(inputId = "go", 
#                     label = "Update")
    ),
    
    mainPanel(
        plotOutput('plot')
    )
)
shinyApp(ui = ui, server = server)
