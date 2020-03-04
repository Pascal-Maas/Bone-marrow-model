setwd("Data/")

library(Matrix)
library(Matrix.utils)
library(VennDiagram)
library(corrplot)

#' @description This function will combine matrices from files into a single matrix.
#' This is done by an inner join to only selected the genes that are shared between 
#' the matrices.
get_combined_matrices <- function(){
    data_files <- list.files()
    M <- Matrix()
    sapply(data_files, function(file){
        m <- head(get(load(file)), 100) # head for testing
        if (nrow(M) <= 1){
            M <<- m 
        } else {
            # Inner join to only select the shared genes
            M <<- join.Matrix(M, m, rownames(M), rownames(m), all.x = F, all.y = F) 
        }
    })
    return(M)
}

#' @description This function will create a plot where the differences between the datasets
#' are plotted as a fraction. This is done in order to determine if there are samples that contain less 
#' genes than others and if they should be removed in order to retain the number of genes.
correlation_plot_genes <- function(omit = c()){
    if (length(omit) > 0){
        data_files <- list.files()[-omit]
    } else {
        data_files <- list.files()
    }
    genes <- c()
    list_of_genes <- lapply(data_files, function(file){
        g <- rownames(get(load(file)))
        if (length(genes) == 0){
            genes <<- g
        } else {
            genes <<- intersect(genes, g)
        }
        return(g)
    })
    M <- as.matrix(sapply(list_of_genes, function(vec1){
        return(sapply(list_of_genes, function(vec2){
            return(length(setdiff(vec1, vec2)) / min(length(vec1), length(vec2)))
        }))
    }))
    corrplot(M, method = "square", title = paste("Shared genes:", length(genes)), mar=c(10,0,10,0))
    return(list("proportion" = M, "genes" = list_of_genes))
}


#l <- get_combined_matrices()

par(mfrow=c(1,2))
correlation_plot_genes()
obj <- correlation_plot_genes(omit = c(5, 17, 20))

m1 <- head(get(load('SRA779509_SRS3805255.sparse.RData')), 20)
colnames(m1) <- paste0("sample_1_", colnames(m1))

m2 <- head(get(load('SRA779509_SRS3805247.sparse.RData')), 20)

# Columns are cells (tags for mRNA)
# Rows are genes

names <- strsplit(rownames(m1), "_") 
a <- unlist(lapply(names, function(s){ # get Only ENGSG00....
    return(s[2])
}))

summary(m1[which(startsWith(a, "ENSG00000010610") == T),]) #CD4
