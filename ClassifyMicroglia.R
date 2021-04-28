## Classify microglia based on 

require(Seurat)
require(readxl)

celltypes10x = read.csv('/home/jovyan/HB_ZIK/celltypes10x.csv')

celltypesSS = read.csv('/home/jovyan/HB_ZIK/celltypesSS.csv')

glioblastoma_10x = LoadH5Seurat('/home/jovyan/data/jhubData/HB_ZIK/HB_ZIK/zikaGlioblastomas_10X_SeuratObject.h5Seurat')

glioblastoma_SS = LoadH5Seurat('/home/jovyan/data/jhubData/HB_ZIK/HB_ZIK/zikaGlioblastomas_SS_SeuratObject.h5Seurat')

load('/home/jovyan/HB_ZIK/microgliaData/metadata-ctrl-gbm.RData')
metadata2 = df

load('/home/jovyan/HB_ZIK/microgliaData/unnormalized_counts_ctrl_gbm.RData')
counts2 = counts_unnorm

genes2 = rownames(counts2)
genes2_split = unlist(lapply(1:length(genes2), function(x) strsplit(genes2[x], split = '_')[[1]][1]))

genes_joint = genes2_split
genes_10x = rownames(glioblastoma_10x@assays$RNA@counts)
genes_SS = rownames(glioblastoma_SS@assays$RNA@counts)

all_genes = intersect(genes_joint, genes_10x)

counts2 = counts2[genes2_split %in% all_genes,]
genes2 = rownames(counts2)
genes2_split = unlist(lapply(1:length(genes2), function(x) strsplit(genes2[x], split = '_')[[1]][1]))
counts2 = counts2[match(all_genes, genes2_split),]

glioblastoma_SS = glioblastoma_SS[genes_SS %in% all_genes,]
genes_SS = rownames(glioblastoma_SS@assays$RNA@counts)
glioblastoma_SS = glioblastoma_SS[match(all_genes, genes_SS),]

glioblastoma_10x = glioblastoma_10x[genes_10x %in% all_genes,]
genes_10x = rownames(glioblastoma_10x@assays$RNA@counts)
glioblastoma_10x = glioblastoma_10x[match(all_genes, genes_10x),]

counts = counts2
rownames(counts) = all_genes

microglia = CreateSeuratObject(counts = counts)
microglia$Cluster = metadata2$Cluster

glioblastoma_10x = glioblastoma_10x[,glioblastoma_10x$SampleName %in% celltypes10x$Cellname[celltypes10x$Classification == 'Microglia']]
glioblastoma_SS = glioblastoma_SS[,glioblastoma_SS$SampleName %in% celltypesSS$Cellname[celltypesSS$Classification == 'Microglia']]

# Now classify 10x and SS microglia based on reference:

microglia = NormalizeData(microglia, normalization.method = "LogNormalize", scale.factor = 10000)
microglia = FindVariableFeatures(microglia, selection.method = "vst", nfeatures = 2000)
microglia <- ScaleData(microglia, verbose = FALSE, features = all_genes)
n_dimensions = 30
microglia <- RunPCA(microglia, npcs = n_dimensions, verbose = FALSE)
microglia <- RunUMAP(microglia, reduction = "pca", dims = 1:n_dimensions)
DimPlot(microglia, reduction = "umap", group.by = "Cluster")

glioblastoma_10x = FindVariableFeatures(glioblastoma_10x, selection.method = "vst", nfeatures = 2000)
glioblastoma_10x <- ScaleData(glioblastoma_10x, verbose = FALSE, features = all_genes)
n_dimensions = 30
glioblastoma_10x <- RunPCA(glioblastoma_10x, npcs = n_dimensions, verbose = FALSE)
glioblastoma_10x <- RunUMAP(glioblastoma_10x, reduction = "pca", dims = 1:n_dimensions)
DimPlot(glioblastoma_10x, reduction = "umap", group.by = "Celltype")

glioblastoma_SS = FindVariableFeatures(glioblastoma_SS, selection.method = "vst", nfeatures = 2000)
glioblastoma_SS <- ScaleData(glioblastoma_SS, verbose = FALSE, features = all_genes)
n_dimensions = 30
glioblastoma_SS <- RunPCA(glioblastoma_SS, npcs = n_dimensions, verbose = FALSE)
glioblastoma_SS <- RunUMAP(glioblastoma_SS, reduction = "pca", dims = 1:n_dimensions)
DimPlot(glioblastoma_SS, reduction = "umap", group.by = "Celltype")

features = VariableFeatures(microglia)
microglia.anchors <- FindTransferAnchors(reference = microglia, query = glioblastoma_10x,
                                         dims = 1:n_dimensions, features = features)
predictions <- TransferData(anchorset = microglia.anchors, refdata = as.character(microglia$Cluster), dims = 1:n_dimensions)
glioblastoma_10x <- AddMetaData(glioblastoma_10x, metadata = predictions)
glioblastoma_10x$MicrogliaCluster = predictions$predicted.id
Idents(glioblastoma_10x) = glioblastoma_10x$MicrogliaCluster

microglia.anchors <- FindTransferAnchors(reference = microglia, query = glioblastoma_SS,
                                         dims = 1:n_dimensions, features = features)
predictions <- TransferData(anchorset = microglia.anchors, refdata = as.character(microglia$Cluster), dims = 1:n_dimensions)
glioblastoma_SS <- AddMetaData(glioblastoma_SS, metadata = predictions)
glioblastoma_SS$MicrogliaCluster = predictions$predicted.id
Idents(glioblastoma_SS) = glioblastoma_SS$MicrogliaCluster

cluster_markers = c('LYVE1', 'MRC1', 'ITPR2', 'IPCEF1', 'MARCKS', 'P2RY13',
                    'IL6ST', 'ACY3', 'LPAR5', 'GATM', 'P2RY12', 'ITM2C', 'CX3CR1',
                    'SELPLG', 'DHRS9', 'A2M', 'HTRA1', 'BIN1', 'IDH2', 'PTAFR',
                    'ALDH2', 'HLA-DRB5', 'CYFIP1', 'TXNIP', 'IFI27', 'MX1', 'SCIN',
                    'FKBP5', 'CPM', 'IFITM3', 'IFI6', 'EIF2AK2', 'IFI44L', 'IFI44',
                    'IFIT3', 'PLAC8', 'CLEC5A', 'LDHA', 'PADI2', 'SLC11A1', 'MSR1',
                    'CEBPD', 'CXCR4', 'ACSL1', 'CD163', 'FCGBP', 'MAFB', 'APOE',
                    'APOC1', 'WIPF3', 'CPVL', 'NAP1L1', 'VSIG4', 'GLDN', 'DDIT4',
                    'SSPP1', 'CACNA1A', 'SERPINE1', 'FTL', 'LPL', 'GPNMB', 'TPT1',
                    'HIF1A', 'NAMPT', 'ELL2', 'VEGFA', 'MT2A', 'C5AR1', 'PLAUR',
                    'HBEGF', 'TFRC', 'IL1RN', 'MCL1', 'SGK1', 'CD83', 'NR4A2',
                    'SRGN', 'NR4A3', 'PDK4', 'CDKN1A', 'OTUD1', 'B3GNT5', 'CH25H',
                    'IER3', 'CCL2', 'TNF', 'CCL8', 'PTGS2', 'EGR1', 'BTG2', 'IL1B',
                    'EGR3', 'EGR2', 'IL1A', 'CCL4L2', 'CCL3L3', 'CCL4')

cluster_markers1 = c('CX3CR1', 'TMEM119', 'CSF1R', 'P2RY12', 'P2RY13', 'SELPLG',
                     'MARCKS', 'TREM2', 'APOE', 'HLA-DRA', 'CD74', 'IFI44L',
                     'CCL2', 'CCL4', 'IFNB', 'IL1B')

print(DoHeatmap(glioblastoma_10x, features = cluster_markers1,
                slot = 'scale.data'))

pdf('10X_MicrogliaClusters.pdf', width = 20, height = 10)
print(DoHeatmap(glioblastoma_10x, features = cluster_markers1,
                slot = 'scale.data'))
dev.off()

pdf('SS_MicrogliaClusters.pdf', width = 20, height = 10)
print(DoHeatmap(glioblastoma_SS, features = cluster_markers1,
                slot = 'scale.data'))
dev.off()

SaveH5Seurat(glioblastoma_SS, filename = "/home/jovyan/data/jhubData/HB_ZIK/HB_ZIK/zikaGlioblastomas_SS_microglia.h5Seurat")
Convert("/home/jovyan/data/jhubData/HB_ZIK/HB_ZIK/zikaGlioblastomas_SS_microglia.h5Seurat", dest = "h5ad")

SaveH5Seurat(glioblastoma_10x, filename = "/home/jovyan/data/jhubData/HB_ZIK/HB_ZIK/zikaGlioblastomas_10x_microglia.h5Seurat")
Convert("/home/jovyan/data/jhubData/HB_ZIK/HB_ZIK/zikaGlioblastomas_10x_microglia.h5Seurat", dest = "h5ad")


