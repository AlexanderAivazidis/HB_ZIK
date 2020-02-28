library(Seurat)
library(Matrix)
library(EnsDb.Hsapiens.v75)
library(rhdf5)
library(dplyr)

setwd('/home/jovyan/HB_ZIK/')

setwd('/home/jovyan/HB_ZIK/')
glioblastoma_counts = read.delim('../data/ZikaGlioblastomas/tic-527/study5953-tic527-star-fc-genecounts.txt', row.names = 1)
manifest = read.delim('../HB_ZIK/5953stdy_manifest_14517_170919_HB_ZIK_046_048_050_051.csv', sep = ',', skip = 8)
manifest = manifest[manifest$SUPPLIER.SAMPLE.NAME != "",]
glioblastoma_metadata = data.frame(matrix('', dim(glioblastoma_counts)[2],4))
colnames(glioblastoma_metadata) = c('SampleName', 'Technology', 'ZikaExposure', 'Patient')
glioblastoma_metadata[,1] = unlist(lapply(colnames(glioblastoma_counts), function(x) substring(x,2)))
glioblastoma_metadata[,2] = rep('SmartSeq', dim(glioblastoma_counts)[2])
glioblastoma_metadata[,3] = rep('TRUE', dim(glioblastoma_counts)[2])
glioblastoma_metadata[,4] = substring(manifest$SUPPLIER.SAMPLE.NAME[match(glioblastoma_metadata$SampleName, manifest$SANGER.SAMPLE.ID)], 3, 4)
patients = c('42','42','43', '43', '45', '45', '46', '46') # (info from sample tracker)
geneNames = as.matrix(read.table(paste('../data/ZikaGlioblastomas/cellranger302_count_32771_5953STDY855119', as.character(1), '_GRCh38-3_0_0_premrna/filtered_feature_bc_matrix/filtered_feature_bc_matrix/features.tsv', sep = '')))

geneOccurence = table(geneNames[,2])
for (gene in names(geneOccurence)){
  if (geneOccurence[gene] > 1){
    geneNames[which(geneNames[,2] == gene),2] = paste(geneNames[which(geneNames[,2] == gene),2], '(', geneNames[which(geneNames[,2] == gene),1], ')', sep = '')
  }
}
rownames(glioblastoma_counts) = geneNames[match(rownames(glioblastoma_counts), geneNames[,1]),2]
glioblastoma_counts = glioblastoma_counts[geneNames[,2],]

for (i in 1:8){
  print(i)
  # prefiltered count_matrix:
  data_subset <- as.matrix(Read10X_h5(filename = paste('../data/HB_ZIK/HB_ZIK/cellranger302_count_32771_5953STDY855119', as.character(i), '_GRCh38-3_0_0_premrna/output_filtered.h5', sep = ''), use.names = TRUE))
  sampleNames = rownames(data_subset)
  glioblastoma_counts = cbind(glioblastoma_counts, data_subset)
  metadata_subset = data.frame(matrix('', dim(data_subset)[2],4))
  colnames(metadata_subset) = c('SampleName', 'Technology', 'ZikaExposure', 'Patient')
  metadata_subset[,1] = colnames(data_subset)
  metadata_subset[,2] = rep('10X', dim(data_subset)[2])
  metadata_subset[,3] = rep('FALSE', dim(data_subset)[2])
  metadata_subset[,4] = rep(patients[i], dim(data_subset)[2])
  glioblastoma_metadata = rbind(glioblastoma_metadata, metadata_subset)
}
mitogenes <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")$gene_id
mitogenes = geneNames[,2][match(mitogenes, geneNames[,1])]
mitogenes = mitogenes[!is.na(mitogenes)] 
percent.mt = colSums(glioblastoma_counts[rownames(glioblastoma_counts) %in% mitogenes,])/colSums(glioblastoma_counts)
Glioblastoma <- CreateSeuratObject(glioblastoma_counts, project = 'HB_ZIK', min.cells = 0, min.features = 0)
Glioblastoma$SampleName = colnames(glioblastoma_counts)
Glioblastoma$Technology = glioblastoma_metadata$Technology
Glioblastoma$ZikaExposure = glioblastoma_metadata$ZikaExposure
Glioblastoma$Patient = glioblastoma_metadata$Patient
Glioblastoma$percent.mt = percent.mt
saveRDS(Glioblastoma, file = "../data/ZikaGlioblastomas/zikaGlioblastomas_SeuratObject.rds")
Glioblastoma = readRDS("../data/ZikaGlioblastomas/zikaGlioblastomas_SeuratObject.rds")

Glioblastoma.list <- SplitObject(Glioblastoma, split.by = 'Technology')
i = 1
VlnPlot(Glioblastoma.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'Patient')
plot1 <- FeatureScatter(Glioblastoma.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Glioblastoma.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
Glioblastoma.list[[i]] <- subset(Glioblastoma.list[[i]], subset = nFeature_RNA > 2500 & nFeature_RNA < 12000 & nCount_RNA < 2*10^6 & percent.mt < 0.1)
VlnPlot(Glioblastoma.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'Patient')
i = 2
VlnPlot(Glioblastoma.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'Patient')
plot1 <- FeatureScatter(Glioblastoma.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Glioblastoma.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
Glioblastoma.list[[i]] <- subset(Glioblastoma.list[[i]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA > 0 & nCount_RNA < 4*10^5 & percent.mt < 0.1)
VlnPlot(Glioblastoma.list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'Patient')

for (i in 1:2){
  Glioblastoma.list[[i]] <- NormalizeData(Glioblastoma.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  Glioblastoma.list[[i]] <- FindVariableFeatures(Glioblastoma.list[[i]], selection.method = "vst", nfeatures = 200)
  top25 <- head(VariableFeatures(Glioblastoma.list[[i]]), 25)
  plot1 <- VariableFeaturePlot(Glioblastoma.list[[i]])
  plot2 <- LabelPoints(plot = plot1, points = top25, repel = TRUE)
  CombinePlots(plots = list(plot1, plot2))
  all.genes <- rownames(Glioblastoma.list[[i]])
  Glioblastoma.list[[i]] <- ScaleData(Glioblastoma.list[[i]], features = all.genes)
}

for (i in 1:2){
  Glioblastoma.list[[i]] <- RunPCA(Glioblastoma.list[[i]], features = VariableFeatures(object = Glioblastoma.list[[i]]))
  print(Glioblastoma.list[[i]][["pca"]], dims = 1:5, nfeatures = 5)
  VizDimLoadings(Glioblastoma.list[[i]], dims = 1:2, reduction = "pca")
  DimPlot(Glioblastoma.list[[i]], reduction = "pca")
  DimHeatmap(Glioblastoma.list[[i]], dims = 1, cells = round(dim(Glioblastoma.list[[i]])[2])/2, balanced = TRUE)
  DimHeatmap(Glioblastoma.list[[i]], dims = 1:15, cells = round(dim(Glioblastoma.list[[i]])[2])/2, balanced = TRUE)
}

for (i in 1:2){
  Glioblastoma.list[[i]] <- JackStraw(Glioblastoma.list[[i]], num.replicate = 100)
  Glioblastoma.list[[i]] <- ScoreJackStraw(Glioblastoma.list[[i]], dims = 1:20)
  JackStrawPlot(Glioblastoma.list[[i]], dims = 1:20)
  ElbowPlot(Glioblastoma.list[[i]], ndims = 40)
}





