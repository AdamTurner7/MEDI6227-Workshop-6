###MEDI6227 Workshop 6
# Load R libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
# The readRDS function is used to load in .rds objects into your R environment
pbmc <- readRDS(file = "pbmc3k_QCB_Session5(1).rds")
# This function calculates the percentage ribosomal content of the RNA for each cell barcode,
# similar to calculating percentage mitochrondial content in QCB Session 5. This new variable
# is saved into the pbmc@meta.data slot.
pbmc@meta.data$pt.Ribosomal <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")  #Find out about regular expression use in R to understand the utility of '^RP[LS]'

####QC METRICS
# The code below calculates the largest gene and its percentage and will be stashed into
# pbmc@meta.data.
pbmc.nomalat <- pbmc[rownames(pbmc) != "MALAT1", ]
pbmc.nomalat$largest_count <- apply(pbmc.nomalat@assays$RNA@counts, 2, max)
pbmc.nomalat$largest_index <- apply(pbmc.nomalat@assays$RNA@counts, 2, which.max)
pbmc.nomalat$largest_gene <- rownames(pbmc.nomalat)[pbmc.nomalat$largest_index]
pbmc.nomalat$percent.Largest.Gene <- 100 * pbmc.nomalat$largest_count/pbmc.nomalat$nCount_RNA
pbmc@meta.data$largest_gene <- pbmc.nomalat$largest_gene
pbmc@meta.data$percent.Largest.Gene <- pbmc.nomalat$percent.Largest.Gene

# Clean up and deleted pbmc.nomalat
rm(pbmc.nomalat)

# Visualize the new QC metrics you have created today
VlnPlot(pbmc, features = c("pt.Ribosomal", "percent.Largest.Gene"))

###Q13 Run the code chunk below and record the gene symbol and full gene name of the top largest gene in the l.gene.sorted data frame here:
# Explotation of the largest gene and its percentages across cell barcodes
l.genes <- as.data.frame(cbind(as.numeric(pbmc@meta.data$percent.Largest.Gene), pbmc@meta.data$largest_gene))
l.genes.sorted <- l.genes[order(l.genes$V1, decreasing = T), ]
# TMSB4X thymosin beta 4 X-linked

####LINEAR REDUCTION
###RUN PCA USING SEURAT FUNCTION ON SCALED DATA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
####Q14- Explore the pmbc Seurat object and record the locations of all of the pca data in the R chunk below (Hint - use “str(pbmc)” and there are several!)
pbmc@reductions[["pca"]]@cell.embeddings
pbmc@reductions[["pca"]]@feature.loadings
pbmc@reductions[["pca"]]@feature.loadings.projected
pbmc@reductions[["pca"]]@assay.used
pbmc@reductions[["pca"]]@global
pbmc@reductions[["pca"]]@stdev

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

####Q15 Do you recognise any of the genes associated with the first five PCs? Record your thoughts below:
#CST3 DC cell marker
#MALAT1 - should have been filtered out
#HLA - HLA's should be present on professional APC's (all cells will express HLA1)


# Examine and visualize PCA results a few different ways
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

####Q16 QCB workshop Question 16 - What is this plot showing you? (Hint: Type in: ??VizDimLoadings into the Console and hit return to see the R help page) Record your thoughts below:
??VizDimLoadings
#Visualize top genes associated with reduction components - the top genes in 2 dimensions (dims = 1:2) that contain the highest amount of variance in PC_2 7 genes produce significantly more variance??

# Examine and visualize PCA results a few different ways. DimPlot produces a plot of the
# cell.embeddings for the first two PCs (PC1 and PC2).
DimPlot(pbmc, reduction = "pca")

###QCB workshop Question 17 - Do you remember making a plot similar to this in QCB Session 5?
##was used to visualise the distribution of data for normalisation techniques e.g.percent mt and feature count
## FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# We can use the group.by option to colour by any other metadata column, today will will use
# largest_gene. We can also add labels to the plot. Finally we can add a call to the
# NoLegend() function to suppress the automatic colour legend which is drawn.

DimPlot(pbmc, reduction = "pca", group.by = "largest_gene", label = TRUE, label.size = 3, repel = T) +
  NoLegend()

DimPlot(pbmc, reduction = "pca", label = TRUE, label.size = 3, repel = T) + NoLegend()

####Q18a Do you find the relationship between the cell.embeddings and names of the most highly expressed gene per cell interesting? What do you think the green cells are? (Hint: Google search ‘Ferritin’ (FTL, FTH1) AND ‘immune cells’)
#The largest values for cell embeddings to not necessarily link to largest genes
#Ferretin upregulated during infection - green cells may represent population of active (myeloid) immune cells?

# Now examine the next two PCs (PC_3 and PC_4)
DimPlot(pbmc, reduction = "pca", group.by = "largest_gene", label = TRUE, label.size = 2, repel = T,
        dims = 3:4) + NoLegend()
#platelets captured in PC4 not all of the useful information is captured in the first two principal components. 

###Q18b - Copy the RunPCA code above into the R chunk below and change the PCs (try dim = 5:6 next). What do you see?
DimPlot(pbmc, reduction = "pca", group.by = "largest_gene", label = TRUE, label.size = 2, repel = T,
        dims = 5:6) + NoLegend()
#variation is is reduced no longer in distinct clusters lost information on certain cell types e.g. platelets and B/T cells

#Examine variance explained by the PCA analysis, by ranking stdev values for each PC calculated by RunPCA
# Extract the variance (stdev column) captured by each PC (PCs column), examine and plot the
# data
PCA_stdev <- as.data.frame(cbind(colnames(pbmc@reductions$pca@feature.loadings), pbmc@reductions$pca@stdev))
colnames(PCA_stdev) <- c("PCs", "stdev")

plot(rownames(PCA_stdev), PCA_stdev$stdev)
###Q19 Describe the relationship between the variables in the PCA-stdev plot above:
#majority of the variation is captured in PC1:PC5 
#PCs X axis, variation Y axis
#low variance still important for detecting clusters that contribute less varaibility 

#scree plot in seurat
x <- 1/mean(1/as.numeric(PCA_stdev$stdev))  # Calculates the harmonic mean of stdev for the first 50 PCs. We are using this to highlight the 'elbow' in the plot, what other value could you use?
ElbowPlot(pbmc, ndims = 50) + geom_hline(yintercept = x, color = "grey")


###Visualization using Feature Heatmaps ordered according to their PCA scores

# Explore the first PC
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# Explore the first 15 PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Q20 - How many cross-window shapes do you see? How many PCs will you select? Do you think the cell type abundance in the dataset impacts these heatmaps in any way?
#Select first 5 PCs variation seems minimal after that

#perform a resampling test (The Jackstraw procedure) 

# NOTE: This process can take a long time for bigger datasets. More approximate techniques
# such as those implemented in ElbowPlot() can be used to reduce computation time.
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# Plot the results of the JackStraw analysis
JackStrawPlot(pbmc, dims = 1:15)

###Q21- How many PCs and which PCs are significant in the above JackStrawPlot?
#13PCs are significant (1:13)

###Q22 - Using all three methods above, did you reach a consensus in the number of PCs to select for downstream analysis? Which method was most useful and why?
#could take an average 5+5+13/3 = 8
#Jackstraw may be most useful as the others are open to interpretation, Jackstraw has an emperical basis and gives P values

#in this exercise we use 10PCs
###CLUSTERING
# PC_1 to PC_1 used (dims = 1:10) Resolution value set at maximum of 1.2
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 1.2)
#smaller the resolution (e.g. less than 1) bigger the cell clusters

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#Q23 - How many cell clusters do you find with resolution = 1.2?
#11?

# PC_1 to PC_1 used (dims = 1:10) Resolution value set at minimum of 0.4
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.4)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

###Q24 - How many cell clusters do you find with resolution = 0.4? Which is the best resolution to use and when do you know to stop?
#9 this seems appropriate as that is the number of cell types in the PBMC sample

#Run non-linear dimensional reduction (UMAP/tSNE)

# Run UMAP analysis using PC_1 to PC_10 - same as your clustering above. Today we are not
# going to use the tSNE approach
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

###Q25 - Describe what you see in the UMAP above? Which clusters are more closely related to each other? Which cluster has the smallest cell number?
#9 cell clusters some more closely related to eachother0,2,4,6::5,1,7::3::8
#cluster 8 has the smallest number

saveRDS(pbmc, file = "pbmc3k_QCB_Session6.rds")

###FIND DIFFERENTIAL EXPRESSED FEATURES
