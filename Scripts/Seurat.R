library(dplyr)
library(Seurat)

# Load the PBMC dataset
input <- "Data/SRA779509_SRS3805245.sparse.RData"
pbmc.data <- get(load(input))
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


# Calculates the percentage of counts coming from Mitochondria genes (starting with MT), per cell
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Shows the QC data when loading SeuratObject. 
# nCount = total RNA molecules
# nFeature = Number of unique RNA molecules
# percent.mt = Percentage mitochondria (as previously calculated)
head(pbmc@meta.data, 5)

# Shows the vignette plot of these metdata. The red area shows the spread
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# These plots show the relationship between the metadata using a scatterplot
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") # Weak negative correlation
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") # Strong positive correlation
CombinePlots(plots = list(plot1, plot2))

# Create a subset of the original set. Only take:
# cells with > 200 unique RNA molecules, (leaked out of broken membrane)
# cells with < 2500 unique RNA molecules, done to filter doubles (will result in too many unique RNAs)
# cells with < 5 percent mitochondria RNA (leaked out of broken membrane)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 

# Scale the counts by 10.000, then normalize using a log transformation
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Standardize the data by calculating the Z-score ((val - avg) / std dev). 
# Outliers are clipped by sqrt(N), where N = amount of cells
# Calculate the variance of each gene across all cells. This measure is done to control for the mean expression.
# Return the top 2.000 genes with the highest standardized variance
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot(plot2)
CombinePlots(plots = list(plot1, plot2))

# Scale the data for PCA, but only select the genes that had a high variance.
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Performs a PCA on the scaled data. It will show the first 5 dimensions
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca") # Plot a regular pca plot of the first 2 components
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") # Show the loadings of the PCA. 
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) # Heatmap of up versus down regulated

# Permutation of Data, sampling 1% to create PCA's. 
# Compares the PCA scores for the sampled genes to the observed PCA to calculate statistical significance 
pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)

# Get p-values per Principle Component
pbmc[["pca"]]@jackstraw@overall.p.values

# Plot the P-values
JackStrawPlot(pbmc, dims = 15:40)

# Alternatively, a Elbowplot can be made
ElbowPlot(pbmc, ndims = 40)


# Select PCA data for modelling:
embeddings <- pbmc[["pca"]]@cell.embeddings # cells with PCs
head(embeddings[,1:5]) # Replace '1:8' with 1:N' where N = number of PCs wanted.

# TODO:
# Combine several matrices together of multiple patients
# Find Golden truth of cells, so an X matrix and y vector can be created for modelling.
# Create a Markdown document that will create plots made. 
