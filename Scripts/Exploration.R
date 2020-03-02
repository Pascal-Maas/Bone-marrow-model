setwd("Data/")

get_list_of_matrices <- function(){
    data_files <- list.files()
    return(lapply(data_files, function(file){
        return(t(get(load(file))))
    }))
}

l <- get_list_of_matrices()


m1 <- get(load('SRA779509_SRS3805255.sparse.RData'))
typeof(m1)

m2 <- get(load('SRA779509_SRS3805247.sparse.RData'))
dplyr::bind_rows(m1, m2)
rbind2(m1, m2, fill=T)

# Columns are cells (tags for mRNA)
# Rows are genes

names <- strsplit(rownames(m1), "_") 
a <- unlist(lapply(names, function(s){ # get Only ENGSG00....
    return(s[2])
}))

summary(m1[which(startsWith(a, "ENSG00000010610") == T),]) #CD4

# TODO:
# Transpose so cells are horizontal, genes are vertical
# Determine the labels of each cell (by overexpression of CD4, CD8 etc..? Or maybe other way?)
# Find a way to combine cells from different patients (averaging maybe? Concatenating?)
# Explore which genes cause most variation (PCA)
# Perform feature selection
# Determine 
# Make model
# Combine into pipeline


# dummy data
set.seed(3344)
A = Matrix(matrix(rbinom(16, 2, 0.2), 4))
colnames(A)=letters[1:4]
B = Matrix(matrix(rbinom(9, 2, 0.2), 3))
colnames(B) = letters[3:5]

# finding what's missing
misA = colnames(B)[!colnames(B) %in% colnames(A)]
misB = colnames(A)[!colnames(A) %in% colnames(B)]

misA
misB

misAl = as.vector(numeric(length(misA)), "list")
names(misAl) = misA
misBl = as.vector(numeric(length(misB)), "list")
names(misBl) = misB

## adding missing columns to initial matrices
An = do.call(cbind, c(A, misAl))
Bn = do.call(cbind, c(B, misBl))[,colnames(An)]

# final bind
rbind(An, Bn)

