library("shiny")
library("devtools")
library("DDRTree")
library("pheatmap")
library("monocle")
library("reshape")
library("rdrop2")

HSMM_ordering_genes <-read.csv("HSMM_ordering_genes.csv")
HSMM_ordering_genes<-as.character(unlist(HSMM_ordering_genes))
HSMM_epi_3<-readRDS("HSMM_epi_3")


server <-function(input, output) {

    output$plot <- renderPlot({
        p <- plot_genes_jitter(HSMM_epi_3[HSMM_ordering_genes[input$num:(input$num+15)],],
                          grouping=input$group, color_by=input$group,
                          nrow=4, ncol=NULL, plot_trend=TRUE)+
            ggtitle(input$title)+
            theme(plot.title = element_text(hjust = 0.5),
                  text = element_text(size=30),
                  legend.position="none")
        
        print(p)
        
    }, height=1000)
    
}
ui <- fluidPage(
    dataset<-diamonds,
    titlePanel("Monocle parameters"),
    
    textInput(inputId = "title", 
              label = "Write a title",
              value = "Finding genes that disinguish Clusters"),
    sidebarPanel(
        
        sliderInput(inputId = "num", 
                    label = "Select a group of 16 genes to analysis. \n
                    Choose the postion for the 1st gene in the group.", 
                    value = 1, step=16, round=0, min = 1, max = 977),
        
        selectInput('group', 'Group by', c("Cluster","State","Pseudotime"))
    ),
    
    mainPanel(
        plotOutput('plot')
    )
)
shinyApp(ui = ui, server = server)
