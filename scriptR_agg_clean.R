# HUMAN BRAIN COVID ANALYSIS # single nuclei RNA sequencing 10x

library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(tidyr)
library(patchwork)

getwd()
setwd("/.../Analysis_4samplesAggregated")

### read aggregate file ###
agg <- Read10X(data.dir = "/Volumes/bf-cluster2/asimatso/HumanBrain/AGG/outs/count/filtered_feature_bc_matrix")
agg <- CreateSeuratObject(counts = agg, project = "HumanBrain", min.cells = 3, min.features = 200)
table(agg$orig.ident)  #19763 nuclei
agg -> brain

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")

# Get batches based on cell names
sample <- sapply(colnames(GetAssayData(object = brain, slot = "counts")),
                 FUN=function(x){substr(x,18,19)})
head(sample)
tail(sample)

sample <- as.numeric(as.character(sample))
names(sample) <- colnames(GetAssayData(object = brain, slot = "counts"))
class(sample)
head(sample)

brain <- AddMetaData(brain, sample, "sample")
head(brain)
tail(brain)
table(brain$orig.ident) #19763 nuclei before filtering



############# QC metrics #############
VlnPlot(object = brain, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(brain, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

brain_filtered <- subset(brain, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
VlnPlot(object = brain_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

head(brain_filtered@meta.data)
table(brain_filtered$orig.ident) 

##### remove more mithocondrial genes ####
hist(brain_filtered$percent.mt)
VlnPlot(object = brain_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

brain_filtered_2 <- subset(brain, subset =percent.mt < 3)
table(brain_filtered_2$orig.ident) #16260
table(brain_filtered_2$sample) #1120 4417 3318 7405 
VlnPlot(object = brain_filtered_2, features = c("percent.mt"), split.by = "sample")

brain_filtered_2 <- NormalizeData(brain_filtered_2, normalization.method = "LogNormalize", scale.factor = 10000)
brain_filtered_2 <- FindVariableFeatures(brain_filtered_2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(brain_filtered_2), 10)
all.genes <- rownames(brain_filtered_2)
brain_filtered_2 <- ScaleData(brain_filtered_2, features = all.genes)
brain_filtered_2 <- RunPCA(brain_filtered_2, features = VariableFeatures(object = brain_filtered_2))
DimHeatmap(brain_filtered_2, dims = 5:15, cells = 500, balanced = TRUE)

brain_filtered_2 <- FindNeighbors(brain_filtered_2, dims = 1:12)
brain_filtered_2 <- FindClusters(brain_filtered_2, resolution = 0.4)
brain_filtered_2 <- RunUMAP(brain_filtered_2, dims = 1:13)

install.packages("RColorBrewer")
library(RColorBrewer)
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 8))
DimPlot(brain_filtered_2, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8, cols = mycolors)
DimPlot(brain_filtered_2, reduction = "umap", split.by =  "sample", group.by = "sample")

saveRDS(brain_filtered_2, file = "/Volumes/Alba Simats/snRNAseq - HumanBrainCOVID/Ranalysis/HumanBrainCOVID/Analysis_4samplesAggregated/brain_filtered_2.rds")

# add grouping
new.grouping <- c("group")
brain_filtered_2[[new.grouping]] <- new.grouping
colnames(brain_filtered_2@meta.data)

brain_filtered_2$group[brain_filtered_2$sample == "1"] <- "COVID"
brain_filtered_2$group[brain_filtered_2$sample == "2"] <- "CONTROL"
brain_filtered_2$group[brain_filtered_2$sample == "3"] <- "COVID"
brain_filtered_2$group[brain_filtered_2$sample == "4"] <- "CONTROL"
head(brain_filtered_2@meta.data, 10)
tail(brain_filtered_2@meta.data, 10)

DimPlot(brain_filtered_2, reduction = "umap", split.by =  "group", group.by = "group")

#find key markers per cluster
Brainmarkers_2 <- FindAllMarkers(object = brain_filtered_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Brainmarkers_2 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table((Brainmarkers_2 %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)), 'markers_by_cluster_2.tsv', sep='\t')

FeaturePlot(brain_filtered_2, features = c("RBFOX3", "SLIT3", "CX3CR1"))
FeaturePlot(brain_filtered_2, features = c("SNAP25", "PLP1","LHFPL3", "AQP4", "CHI3L1", "DNAH9","CLDN5", "CD74", "CSF1R"), cols = c( "#E6E6E6", "#D93622", "#A60021"))
DotPlot(brain_filtered_2,features = c("SNAP25", "PLP1","LHFPL3", "AQP4", "CHI3L1", "DNAH9","CLDN5", "CD74", "CSF1R" ), group.by = "seurat_clusters")
DotPlot(brain_filtered_2,features = c( "AQP4","CLDN5", "DNAH9","CSF1R", "SNAP25", "PLP1","LHFPL3", "CD8A","CHI3L1"), group.by = "label", cols = c( "#E6E6E6", "#A60021"))


# Add labels to clusters #
new.label <- c("label")
brain_filtered_2[[new.label]] <- new.label
colnames(brain_filtered_2@meta.data)

brain_filtered_2$label[brain_filtered_2$seurat_clusters == "0"] <- "Oligodendrocytes-1"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "1"] <- "Oligodendrocytes-2"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "2"] <- "Microglia"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "3"] <- "Astrocytes"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "4"] <- "Oligodendrocytes-3"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "5"] <- "OPC"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "6"] <- "Neurons-1"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "7"] <- "Oligodendrocytes-4"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "8"] <- "Ependymal cells"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "9"] <- "Neurons-2"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "10"] <- "VSM"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "11"] <- "Endothelial cells"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "12"] <- "Neurons-3"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "13"] <- "T cells"
brain_filtered_2$label[brain_filtered_2$seurat_clusters == "14"] <- "Neurons-4"

DimPlot(brain_filtered_2, reduction = "umap", group.by = "label", cols = mycolors, label = TRUE)


# group cells per celltype #
new.label <- c("label_2")
brain_filtered_2[[new.label]] <- new.label
colnames(brain_filtered_2@meta.data)
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "0"] <- "Oligodendrocytes"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "1"] <- "Oligodendrocytes"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "2"] <- "Microglia"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "3"] <- "Astrocytes"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "4"] <- "Oligodendrocytes"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "5"] <- "OPC"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "6"] <- "Neurons"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "7"] <- "Oligodendrocytes"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "8"] <- "Ependymal cells"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "9"] <- "Neurons"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "10"] <- "VSM"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "11"] <- "Endothelial cells"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "12"] <- "Neurons"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "13"] <- "T cells"
brain_filtered_2$label_2[brain_filtered_2$seurat_clusters == "14"] <- "Neurons"

DimPlot(brain_filtered_2, reduction = "umap", group.by = "label_2", cols = mycolors, label = TRUE)
DotPlot(brain_filtered_2,features = c( "AQP4","CLDN5", "DNAH9","CSF1R", "SNAP25", "PLP1","LHFPL3", "CD8A","CHI3L1"), group.by = "label_2", cols = c( "#E6E6E6", "#A60021"))


# DEG in all cell types #
Idents(object = brain_filtered_2) <- "label_2"
DimPlot(brain_filtered_2, reduction = "umap", cols = mycolors, label = TRUE)

DEGoligo <- FindMarkers(object = brain_filtered_2, ident.1 = "COVID", group.by = "group", subset.ident = "Oligodendrocytes", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(DEGoligo, 'DEGoligo.tsv', sep='\t')

DEGastro <- FindMarkers(object = brain_filtered_2, ident.1 = "COVID", group.by = "group", subset.ident = "Astrocytes", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(DEGastro, 'DEGastro.tsv', sep='\t')

DEGendo <- FindMarkers(object = brain_filtered_2, ident.1 = "COVID", group.by = "group", subset.ident = "Endothelial cells", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(DEGendo, 'DEGendo.tsv', sep='\t')

DEGneurons <- FindMarkers(object = brain_filtered_2, ident.1 = "COVID", group.by = "group", subset.ident = "Neurons", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(DEGneurons, 'DEGneurons.tsv', sep='\t')

DEGopc <- FindMarkers(object = brain_filtered_2, ident.1 = "COVID", group.by = "group", subset.ident = "OPC", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(DEGopc, 'DEGopc.tsv', sep='\t')


###### Filter microglia only ########
Idents(brain_filtered_2) <- "seurat_clusters"
microglia <- subset(x = brain_filtered_2, idents = c("2"))
DimPlot(object = microglia, split.by = "group", cols= mycolors)

Idents(microglia_2) = "group"
microglia_2.markers.COV <- FindMarkers(microglia_2, ident.1 = "COVID", min.pct = 0.25)
head(microglia_2.markers.COV, n = 5)
write.table(microglia_2.markers.COV, 'markers_microglia_2.tsv', sep='\t')

library(EnhancedVolcano)
EnhancedVolcano(microglia_2.markers.COV,
                lab = rownames(microglia_2.markers.COV),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.30,
                pointSize = 3.0,
                labSize = 6.0,
                title = "microglia_2: COVID vs. CONTROL",
                xlim = c(-2.2,2.2),
                ylim = c(0, 30),
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 1,
                legendPosition = 'down',
                legendLabSize = 12,
                legendIconSize = 4.0,
                arrowheads = FALSE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps =20, 
                selectLab = mito.genes
                # selectLab = c("P2RY12", "CX3CR1")
                )


# recluster microglia #
microglia_2 <- NormalizeData(microglia_2, normalization.method = "LogNormalize", scale.factor = 10000)
microglia_2 <- FindVariableFeatures(microglia_2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(microglia_2), 10)
all.genes <- rownames(microglia_2)
microglia_2 <- ScaleData(microglia_2, features = all.genes)
microglia_2 <- RunPCA(microglia_2, features = VariableFeatures(object = brain_filtered_2))
DimHeatmap(microglia_2, dims = 1:10, cells = 500, balanced = TRUE)

microglia_2 <- FindNeighbors(microglia_2, dims = 1:4)
microglia_2 <- FindClusters(microglia_2, resolution = 0.4) # finally this one!!!

microglia_2 <- RunUMAP(microglia_2, dims = 1:4)
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 8))
DimPlot(microglia_2, reduction = "umap", label = FALSE, pt.size = 1, label.size = 10, cols = mycolors, split.by = "sample")
DimPlot(microglia_2, reduction = "umap", split.by =  "group", group.by = "seurat_clusters",cols = mycolors)

#all clusters markers
Brainmarkers.microglia <- FindAllMarkers(object = microglia_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Brainmarkers.microglia %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table((Brainmarkers.microglia %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)), 'markers_by_cluster_microglia.re0.4.tsv', sep='\t')

# IDENTIFY EACH SUBCLUSTER OF MICROGLIA RECLUSTERED # res 0.4
table(microglia_2$sample)
table(microglia_2$group)
FeaturePlot(microglia_2, features = c("P2RY12", "CX3CR1"))
FeaturePlot(microglia_2, features = c("APOE", "SPP1", "TGFBR1", "C1QC", "C1QA", "IL1B"))
DimPlot(microglia_2, reduction = "umap", split.by =  "group", group.by = "seurat_clusters",cols = mycolors)
DimPlot(microglia_2, reduction = "umap", group.by = "group", cols = mycolors)
DimPlot(microglia_2, reduction = "umap", group.by = "seurat_clusters", cols = mycolors)

head(microglia_2)
Idents(microglia_2) <- "RNA_snn_res.0.2" #resolution 0.2
DimPlot(microglia_2, reduction = "umap", label = FALSE, pt.size = 1, label.size = 10, cols = mycolors, group.by = "RNA_snn_res.0.2")
DimPlot(microglia_2, reduction = "umap", label = FALSE, pt.size = 1, label.size = 10, cols = mycolors, group.by = "RNA_snn_res.0.2", split.by = "group")


# Which mitocondrial genes? #microglia c2 from original, without reclustering
mito.genes <- grep(pattern = "^MT-", x = rownames(x =brain_filtered_2), value = TRUE)
mito.genes
VlnPlot(microglia_2, features = mito.genes)
RidgePlot(microglia_2, features = mito.genes)

EnhancedVolcano(microglia_2.markers.COV,
                lab = rownames(microglia_2.markers.COV),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.30,
                pointSize = 3.0,
                labSize = 6.0,
                title = "microglia_2: COVID vs. CONTROL",
                xlim = c(-2,2),
                ylim = c(0, 30),
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                arrowheads = FALSE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps =10, 
                #selectLab = mito.genes
)


# mito.genes in other cell types: 

# oligodendrocytes
Idents(brain_filtered_2) <- "seurat_clusters"
oligo <- subset(x = brain_filtered_2, idents = c("0", "1", "4", "7"))
DimPlot(object = oligo, group.by = "seurat_clusters", cols= mycolors, label = TRUE)

Idents(oligo) <- "group"
oligo.markers.COV <- FindMarkers(oligo, ident.1 = "COVID", min.pct = 0.25)
head(oligo.markers.COV, n = 5)
write.table(oligo.markers.COV, 'markers_oligo.tsv', sep='\t')

# neurons
Idents(brain_filtered_2) <- "seurat_clusters"
neurons <- subset(x = brain_filtered_2, idents = c("12", "6", "9", "14"))
DimPlot(object = neurons, group.by = "seurat_clusters", cols= mycolors, label = TRUE)

Idents(neurons) <- "group"
neurons.markers.COV <- FindMarkers(neurons, ident.1 = "COVID", min.pct = 0.25)
head(neurons.markers.COV, n = 5)
write.table(neurons.markers.COV, 'markers_neurons.tsv', sep='\t')

# astrocytes
Idents(brain_filtered_2) <- "seurat_clusters"
astrocytes <- subset(x = brain_filtered_2, idents = c("3"))
DimPlot(object = astrocytes, group.by = "seurat_clusters", cols= mycolors, label = TRUE)

Idents(astrocytes) <- "group"
astrocytes.markers.COV <- FindMarkers(astrocytes, ident.1 = "COVID", min.pct = 0.25)
head(astrocytes.markers.COV, n = 5)
write.table(astrocytes.markers.COV, 'markers_astrocytes.tsv', sep='\t')

# endothelial
Idents(brain_filtered_2) <- "seurat_clusters"
EC <- subset(x = brain_filtered_2, idents = c("11"))
DimPlot(object = EC, group.by = "seurat_clusters", cols= mycolors, label = TRUE)

Idents(EC) <- "group"
EC.markers.COV <- FindMarkers(EC, ident.1 = "COVID", min.pct = 0.25)
head(EC.markers.COV, n = 5)
write.table(EC.markers.COV, 'markers_EC.tsv', sep='\t')

# T cells
Idents(brain_filtered_2) <- "seurat_clusters"
Tcells <- subset(x = brain_filtered_2, idents = c("13"))
DimPlot(object = Tcells, group.by = "seurat_clusters", cols= mycolors, label = TRUE)

Idents(Tcells) <- "group"
Tcells.markers.COV <- FindMarkers(Tcells, ident.1 = "COVID", min.pct = 0.25)
head(Tcells.markers.COV, n = 5)
write.table(Tcells.markers.COV, 'markers_Tcells.tsv', sep='\t')

# VSM NOT IN COVID
# EpC NOT IN COVID






########## CELLCHAT ##########
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
devtools::install_github("sqjin/CellChat", force = TRUE)
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# only COVID group #
Idents(object = brain_filtered_2) <-  "group"
DimPlot(object = brain_filtered_2)
brain_COVID <- subset(x = brain_filtered_2, idents = "COVID")
DimPlot(object = brain_COVID, group.by = "group", cols= mycolors, split.by = "group")
table(brain_COVID$group) 

# only CONTROL group #
Idents(object = brain_filtered_2) <-  "group"
DimPlot(object = brain_filtered_2)
brain_CONTROL <- subset(x = brain_filtered_2, idents = "CONTROL")
DimPlot(object = brain_CONTROL, group.by = "group", cols= mycolors, split.by = "group")
table(brain_CONTROL$group) 


## Part I: Data input & processing and initialization of CellChat object
head(brain_filtered_2)
cellchat <- createCellChat(brain_filtered_2, assay = "RNA", meta = NULL, group.by = "label_2")
cellchat <- addMeta(cellchat, meta = brain_filtered_2@meta.data)
cellchat <- setIdent(cellchat, ident.use = "label_2") # set "label" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

head(brain_CONTROL)
cellchat <- createCellChat(brain_CONTROL, assay = "RNA", meta = NULL, group.by = "label_2")
cellchat <- addMeta(cellchat, meta = brain_CONTROL@meta.data)
cellchat <- setIdent(cellchat, ident.use = "label_2") # set "label" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

head(brain_COVID)
cellchat <- createCellChat(brain_COVID, assay = "RNA", meta = NULL, group.by = "label_2")
cellchat <- addMeta(cellchat, meta = brain_COVID@meta.data)
cellchat <- setIdent(cellchat, ident.use = "label_2") # set "label" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#remove EpC and VSM in CONTROL
head(brain_CONTROL)
table(brain_CONTROL$label_2)
Idents(brain_CONTROL) <- "label_2"
brain_CONTROL_reduced <- subset(brain_CONTROL, idents = c("Ependymal cells", "VSM"), invert = TRUE)
table(brain_CONTROL_reduced$label_2)
cellchat <- createCellChat(brain_CONTROL_reduced, assay = "RNA", meta = NULL, group.by = "label_2")
cellchat <- addMeta(cellchat, meta = brain_CONTROL_reduced@meta.data)
cellchat <- setIdent(cellchat, ident.use = "label_2") # set "label" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#remove EpC and VSM in ALL (_reduced)
head(brain_filtered_2)
table(brain_filtered_2$label_2)
table(brain_CONTROL$label_2)
table(brain_COVID$label_2)
Idents(brain_filtered_2) <- "label_2"
brain_ALL_reduced <- subset(brain_filtered_2, idents = c("Ependymal cells", "VSM"), invert = TRUE)
table(brain_ALL_reduced$label_2)
cellchat <- createCellChat(brain_ALL_reduced, assay = "RNA", meta = NULL, group.by = "label_2")
cellchat <- addMeta(cellchat, meta = brain_ALL_reduced@meta.data)
cellchat <- setIdent(cellchat, ident.use = "label_2") # set "label" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#from now on, either with COVID or CONTROL samples, or ALL together
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


## Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)# Compute the communication probability and infer cellular communication network
computeAveExpr(cellchat, features = c("CX3C"), type =  "truncatedMean", trim = 0.1) # The function computeAveExpr can help to check the average expression of signaling genes of interest, 
cellchat <- filterCommunication(cellchat, min.cells = 2)
df.net <- subsetCommunication(cellchat) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
cellchat <- computeCommunProbPathway(cellchat) # Infer the cell-cell communication at a signaling pathway level
cellchat <- aggregateNet(cellchat) # Calculate the aggregated cell-cell communication network

# save cellchat files 
saveRDS(cellchat, file = "cellchat_COVID.rds")
saveRDS(cellchat, file = "cellchat_CONTROL_reduced.rds")
saveRDS(cellchat, file = "cellchat_ALL_reduced.rds")

# load cellchat files
cellchat.CONTROL_reduced <- readRDS("cellchat_CONTROL_reduced.rds")
cellchat.ALL_reduced <- readRDS("cellchat_ALL_reduced.rds")
cellchat.COVID <- readRDS("/Volumes/Alba Simats/snRNAseq - HumanBrainCOVID/Ranalysis/HumanBrainCOVID/Analysis_4samplesAggregated/cellchat_COVID.rds")

cellchat <- cellchat.CONTROL
cellchat <- cellchat.ALL
cellchat <- cellchat.COVID

#We can also visualize the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


## Part III: Visualization of cell-cell communication network
pathways.show.all <- cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

# Chord diagram
cellchat.COVID@netP$pathways
cellchat.CONTROL@netP$pathways
cellchat.ALL@netP$pathways
cellchat@netP$pathways

pathways.show <- c("SPP1") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

pathways.show <- c("CX3C") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
par(mfrow = c(1,1), xpd=TRUE)
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(1:11), lab.cex = 0.5,legend.pos.y = 30) #source 4 is microglia
netVisual_chord_gene(cellchat, sources.use = c(1:9), targets.use = c(1:9), signaling = c("CX3C"),legend.pos.x = 8)
netVisual_chord_cell(cellchat.ALL_reduced, sources.use = c(4), targets.use = c(3), signaling = c("CX3C"),legend.pos.x = 8)



##Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Hierarchy plot
vertex.receiver = seq(1,9) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver, group = "group")
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", group = "group")




######  COMPARISON ###### 
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
library(CellChat)
library(patchwork)
data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)

levels(cellchat.CONTROL_reduced@idents)
levels(cellchat.COVID@idents)
object.list <- list(CONTROL = cellchat.CONTROL_reduced, COVID = cellchat.COVID)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

### Comparing the communication probability between two datasets for each L-R pair and each pair of cell groups
# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat, sources.use = c(1,2,3,4,5,6,7), targets.use = c(1,2,3,4,5,6,7))
df.net

#subset to specific signaling
df.net.cx3cr1 <- subsetCommunication(cellchat, signaling = c("CX3C"))

# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")

# Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight") 
# where red (or blue ) colored edges represent increased (or decreased ) signaling in the second dataset compared to the first one


### Signaling ligand-receptor pairs based on the differential gene expression analysis 
pos.dataset = "COVID"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "COVID",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "CONTROL",ligand.logFC = -0.1, receptor.logFC = -0.1)


# Circle plot
pathways.show <- c("CX3C") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

# save cellchat files 
saveRDS(cellchat, file = "cellchat_MERGED.rds")
saveRDS(cellchat, file = "cellchat_MERGED.reduced.rds")
saveRDS(cellchat, file = "cellchat_MERGED.noEpCVSM.rds")

# load cellchat files
cellchat.MERGED <- readRDS("cellchat_MERGED.rds")
cellchat.MERGED.reduced <- readRDS("cellchat_MERGED.reduced.rds")
cellchat.MERGED.reduced <- readRDS("/Volumes/AS/snRNAseq - HumanBrainCOVID/Ranalysis/HumanBrainCOVID/Analysis_4samplesAggregated/cellchat_MERGED.noEpCVSM.rds")







