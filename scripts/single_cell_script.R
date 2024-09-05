library(GEOquery)
library(Matrix)
library(scCATCH)
library(readr)
#
# extract clinical info of 8 patients included in the study (https://pubmed.ncbi.nlm.nih.gov/33691089/)
dataset=getGEO("GSE166992")
clinical=pData(dataset[[1]])
table(clinical$`disease status:ch1`)
data.frame(clinical$title, clinical$`disease status:ch1`)
# run this once only; make sure to create an empty data folder before running this. 
untar("~/Desktop/comp_analysis/single_cell/GSE166992_RAW.tar", exdir="data")

# processing of the exp data to generate expression matrices
samples = paste("IVAR",c(2:7, 9, 10), sep="")
samples
outdir = "~/Desktop/comp_analysis/single_cell/data/";
data.10x = list(); # empty list to hold the exp matrices
data.10x[[1]] <- Read10X(data.dir = "~/Desktop/comp_analysis/single_cell/data/IVAR2");
data.10x[[2]] <- Read10X(data.dir = "~/Desktop/comp_analysis/single_cell/data/IVAR3");
data.10x[[3]] <- Read10X(data.dir = "~/Desktop/comp_analysis/single_cell/data/IVAR4");
data.10x[[4]] <- Read10X(data.dir = "~/Desktop/comp_analysis/single_cell/data/IVAR5");
data.10x[[5]] <- Read10X(data.dir = "~/Desktop/comp_analysis/single_cell/data/IVAR6");
data.10x[[6]] <- Read10X(data.dir = "~/Desktop/comp_analysis/single_cell/data/IVAR7");
data.10x[[7]] <- Read10X(data.dir = "~/Desktop/comp_analysis/single_cell/data/IVAR9");
data.10x[[8]] <- Read10X(data.dir = "~/Desktop/comp_analysis/single_cell/data/IVAR10");
scrna.list = list()
for (i in 1:length(data.10x)) {
  scrna.list[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=100, min.features=700, project=samples[i]);
  scrna.list[[i]][["DataSet"]] = samples[i];
}
table(scrna$DataSet)
class_labels=c(rep("Healthy", 12925), rep("COVID_19", 4872), rep("Healthy", 4789), rep("COVID_19",sum(3867,5559, 3919, 858)), rep("Healthy",10374))
scrna$class_labels=class_labels

# merging 8 expression measurements
scrna <- merge(x=scrna.list[[1]], y=c(scrna.list[[2]],scrna.list[[3]], scrna.list[[4]], scrna.list[[5]], scrna.list[[6]], scrna.list[[7]],scrna.list[[8]]), add.cell.ids = samples, project="COVID");
rm(scrna.list) # remove this to save space
rm(data.10x) # same
# percent mito
mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna), value = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts'));
scrna[['percent.mito']] <- percent.mito*100;
scrna$percent.mito

# percent ribo
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna), value = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts'));
scrna[['percent.ribo']] <- percent.ribo*100;

# dataset QC plots
par(mfrow=c(2,2))
vln=VlnPlot(object = scrna, features = c("percent.mito", "percent.ribo"), pt.size=0, ncol = 2, group.by="class_labels");
plot(vln);
vln <- VlnPlot(object = scrna, features = c("nCount_RNA", "nFeature_RNA"), pt.size=0, group.by="class_labels", y.max=25000)
plot(vln)

# cell-level filtering
scrna_filtered <- subset(scrna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 5)
scrna_filtered
# gene-level filtering
counts <- GetAssayData(object =scrna_filtered, slot = "counts")

# output a logical matrix specifying whether a gene has counts > 0 for a given cell
#nonzero <- counts > 0
# sum all TRUE values and returns TRUE if more than 10 TRUE values per gene
#keep_genes <- Matrix::rowSums(nonzero) >= 10
keep_genes <- Matrix::rowSums(counts > 0) >= 10

# keeping genes expressed in > 10 cells
filtered_counts <- counts[keep_genes, ]
scrna_filtered <- CreateSeuratObject(filtered_counts, meta.data = scrna_filtered@meta.data)
scrna_filtered
scrna

# find highly variable genes
scrna_filtered <- FindVariableFeatures(scrna_filtered, selection.method = "vst", nfeatures = 2000)

# identify 10 most highly variable genes
top10 <- head(VariableFeatures(scrna_filtered), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scrna_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(scrna_filtered)
scrna_filtered <- ScaleData(scrna_filtered, features = all.genes)

# PCA
scrna_filtered <- RunPCA(scrna_filtered, features = VariableFeatures(object = scrna_filtered))
VizDimLoadings(scrna_filtered, dims = 1:2, reduction = "pca")

# determine number of PCs
ElbowPlot(object = scrna_filtered, 
          ndims = 50)

# clustering
scrna_filtered <- FindNeighbors(scrna_filtered, dims = 1:12)
scrna_filtered <- FindClusters(scrna_filtered, resolution = c(0.6))
scrna_filtered <- RunUMAP(scrna_filtered, dims = 1:12)
DimPlot(scrna_filtered, reduction = "umap", label = TRUE)

#cell type annotation using scCATCH
clu_markers <- scCATCH::findmarkergenes(object = scrna_filtered,
                                        species = 'Human',
                                        cluster = 'All',
                                        match_CellMatch = TRUE,
                                        cancer = NULL,
                                        tissue = c('Blood','Peripheral blood'),
                                        cell_min_pct = 0.25,
                                        logfc = 0.25,
                                        pvalue = 0.05)

#final list of cell types using scCATCH
clu_ann <- scCATCH(clu_markers$clu_markers,
                   species = 'Human',
                   tissue = c('Blood','Peripheral blood'))

clu_ann

#run singleR
ref_1=HumanPrimaryCellAtlasData()
pred <- SingleR(test=GetAssayData(scrna_filtered), ref=ref_1, labels=ref_1$label.main)
table(pred$labels,scrna_filtered$class_labels)

# assign consensus cluster identities from both methods
new.cluster.ids <- c("T cell", "Naive T cell", "T cell", "Memory T cell", "Monocyte", "T cell", "Macrophages", "Plasmacytoid dendritic cell",
                     "B cell", "Monocyte", "CD8+ T cell", "T cell", "Basophil", "B cell", "Platelet", "B cell", "Plasmacytoid dendritic cell", "B cell")
length(new.cluster.ids)
scrna_copy=scrna_filtered 
names(new.cluster.ids) <- levels(scrna_copy)
scrna_copy <- RenameIdents(scrna_copy, new.cluster.ids)
DimPlot(scrna_copy, reduction = "umap", pt.size = 0.5)

# calculate cluster-wise proportion of covid vs healthy patients
prop=as.data.frame(prop.table(table(scrna_filtered$seurat_clusters, scrna_filtered$class_labels)))
colnames(prop)=c("Cluster","Label","Proportion")
head(prop)
for(i in 0:17){
  vec=prop[which(prop$Cluster == i),"Proportion"]
  print(sum(vec))
}

# get list of markers for all clusters
markers <- FindAllMarkers(object = scrna_filtered, 
                          logfc.threshold = 0.25)  
