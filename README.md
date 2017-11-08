# Analysis of single cell RNA-seq data

List of software packages for single-cell RNA-seq analysis
List of workflow and tutorials for single-cell RNA-seq analysis

## Software packages

- [Monocle](http://cole-trapnell-lab.github.io/monocle-release/) - [R] - Differential expression and time-series analysis for single-cell RNA-Seq.
- [scater](http://bioconductor.org/packages/release/bioc/html/scater.html) - [R] - Scater places an emphasis on tools for quality control, visualisation and pre-processing of data before further downstream analysis, filling a useful niche between raw RNA-sequencing count or transcripts-per-million data and more focused downstream modelling tools such as monocle, scLVM, SCDE, edgeR, limma and so on.
- [Seurat](http://satijalab.org/seurat/) - [R] - It contains easy-to-use implementations of commonly used analytical techniques, including the identification of highly variable genes, dimensionality reduction (PCA, ICA, t-SNE), standard unsupervised clustering algorithms (density clustering, hierarchical clustering, k-means), and the discovery of differentially expressed genes and markers.



## Tutorials and workflows
 
- [Aaron Lun's Single Cell workflow on Bioconductor](http://bioconductor.org/help/workflows/simpleSingleCell/) - [R] - This article describes a computational workflow for basic analysis of scRNA-seq data using software packages from the open-source Bioconductor project.

- [Cole Trapnell's Single Cell workflow](http://cole-trapnell-lab.github.io/monocle-release/docs/) - [R] - Monocle introduced the strategy of ordering single cells in pseudotime, placing them along a trajectory corresponding to a biological process such as cell differentiation by taking advantage of individual cell's asynchronous progression of those processes.
example:
Reproducing the original pseudo-time clustering plot in Shiny :https://nyhuyang.shinyapps.io/single_cell_clusters/
Distribution based on gene expression levels (y-axis) in Shiny: https://nyhuyang.shinyapps.io/single_cell_expression/

- [Hemberg lab's Analysis of single cell RNA-seq data](http://hemberg-lab.github.io/scRNA.seq.course/) - [R] - This course we will discuss some of the questions that can be addressed using scRNA-seq as well as the available computational and statistical methods avialable. The course is taught through the University of Cambridge Bioinformatics training unit
