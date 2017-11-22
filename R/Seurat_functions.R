# subset by protocol ================

FilterCellsBy <- function(object,name){
        cells.name <- grepl(name,unlist(object@data@Dimnames))
        cells.name <- unlist(object@data@Dimnames)[cells.name]
        object.name <-SubsetData(object = object,
                                 cells.use =cells.name)
        print(table(object.name@meta.data$protocol)) #test
        return(object.name)
}

# rename ident back to 0,1,2,3...
RenameIdentBack <- function(object){
        old.ident.ids <- levels(object@ident)
        new.cluster.ids <- 1:length(old.ident.ids)-1
        object@ident <- plyr::mapvalues(x = object@ident,
                                        from = old.ident.ids,
                                        to = new.cluster.ids)
        return(object)
}

# rename ident by gene expression
RenameIdentbyGene <- function(object,label.genes =NULL,accept.low = 1){
        
        if(is.null(label.genes)){print("At least one label gene and one gene only")}
        
        label.genes <- capitalize(tolower(label.genes))
        label.genes <- label.genes[label.genes %in% object@raw.data@Dimnames[[1]]]
        
        #split seurat object by Gene expression level
        object.label.genes <- SubsetData(object = object,
                                         subset.name = label.genes,
                                         accept.low = accept.low)
        
        object.nolabel.genes <- SubsetData(object = object,
                                           subset.name = label.genes,
                                           accept.high = accept.low-0.0001)
        # rename two objects with ident 0 and 1, respectively
        old.ident.ids <- levels(object.label.genes@ident)
        new.cluster.ids <- rep(0,nlevels(object.label.genes@ident))
        object.label.genes@ident <- plyr::mapvalues(x = object.label.genes@ident,
                                                    from = old.ident.ids,
                                                    to = new.cluster.ids)
        
        old.ident.ids <- levels(object.nolabel.genes@ident)
        new.cluster.ids <- rep(1,nlevels(object.nolabel.genes@ident))
        object.nolabel.genes@ident <- plyr::mapvalues(x = object.nolabel.genes@ident,
                                                      from = old.ident.ids,
                                                      to = new.cluster.ids)         
        # Merge two Seurat objects
        object1 <- MergeSeurat(object.label.genes, object.nolabel.genes,
                               names.field = 1,
                               do.normalize =FALSE,
                               add.cell.id1= 0,add.cell.id2=1)
        return(object1)
}

# Differential Expression Testing with average nUMI,among different Ident====
DiffTTestbyIdent <- function(object){
        
        cluster.markers <- NULL
        avg_UMI <- NULL
        DiffTTestResults <- NULL
        
        for(i in 1:nlevels(object@ident)){
                # find markers for every cluster compared to all remaining cells
                cluster.markers[[i]] <- FindMarkers(object = object,
                                                    ident.1 = i-1, #careful with cluster 0
                                                    test.use = "bimod",
                                                    thresh.use = -Inf,
                                                    min.pct = -Inf,
                                                    min.cells = -Inf) 
                # resluts = "p_val",  "avg_diff", "pct.1", "pct.2"
                
                cluster.markers[[i]]$p_val <- p.adjust(cluster.markers[[i]]$p_val,"BH")
                colnames(cluster.markers[[i]])[1] <- "adj.p_val"
                # resluts = "adj.p_val",  "avg_diff", "pct.1", "pct.2"
                
                avg_UMI[[i]] <-rowMeans(as.matrix(x = object@data[, WhichCells(object = object,
                                                                               ident = i-1)]))
                avg_UMI[[i]] <-data.frame(avg_UMI = avg_UMI[[i]])
                DiffTTestResults[[i]] <- cbind(cluster.markers[[i]],
                                               avg_UMI[[i]][,"avg_UMI"][match(rownames(cluster.markers[[i]]), 
                                                                              rownames(avg_UMI[[i]]))])
                DiffTTestResults[[i]] <- DiffTTestResults[[i]][,c(1,2,5)]
                colnames(DiffTTestResults[[i]])[3] <- "avg_UMI"
                # resluts = "adj.p_val",  "avg_diff", "avg_UMI"
        }
        for(i in 1:nlevels(object@ident)){
                print(dim(DiffTTestResults[[i]]))
        }
        return(DiffTTestResults)
}


# Calculate Log2fold change
Log2fold <- function(object){
        #Limma's "Log(FC)" = mean(log2(Group1)) - mean(log2(Group2))
        # = Sum(log2(Group1))/n1 - Sum(log(all other groups))/(N-n1)     
        sum_UMI <- data.frame(matrix(0,nrow = object@data@Dim[1], # number of genes
                                     ncol = nlevels(object@ident))) # number of idents
        nCells <- data.frame(matrix(0,nrow = 1,ncol = nlevels(object@ident)))
        
        for(i in 1:nlevels(object@ident)){
                exprs <- as.matrix(x = object@data[, WhichCells(object = object,
                                                                ident = i-1)])
                sum_UMI[,i] <- rowSums(log2(exprs+1))
                nCells[1,i] <- ncol(exprs)
                print(head(as.data.frame(sum_UMI[,i]),3));print(nCells[1,i])
        }
        total.cells <- object@data@Dim[2]
        #       Log2FC = Sum(log2(Group1))/n1 - Sum(log(all other groups))/(N-n1) 
        Log2FC <- NULL
        for(i in 1:nlevels(object@ident)){
                Log2FC[[i]] = data.frame(Log2FC= sum_UMI[,i]/nCells[1,i]- 
                                                 rowSums(sum_UMI[,-i])/(total.cells-nCells[1,i]))
                rownames(Log2FC[[i]]) <-rownames(exprs)
        }
        print(lapply(Log2FC,head,3))
        return(Log2FC)
}
# make DiffTTestResults into one data.frame and rename columns
CbindRename <- function(x = x, y= NULL){
        x0 <- x[[1]]
        y0 <- y[[1]]
        name <- c(names(x0),names(y0))
        Rownames <- rownames(x0)
        x0 <- cbind(x0,y0[match(rownames(x0),rownames(y0)),])
        for(i in 2:length(x)){
                x0 <- cbind(x0,x[[i]][match(rownames(x0),
                                            rownames(x[[i]])),])
                x0 <- cbind(x0,y[[i]][match(rownames(x0),
                                            rownames(y[[i]])),])
        }
        
        name <- paste0(rep(paste0(name,".C"),length(x)), 
                       rep(0:(length(x)-1),each = length(name)))
        colnames(x0) <- name
        rownames(x0) <- Rownames
        print(head(x0[1:6]))
        return(x0)
}
