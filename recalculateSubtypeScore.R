library(Seurat)
library(Matrix)
library(EnsDb.Hsapiens.v75)
library(rhdf5)
library(dplyr)
library(ggplot2)

setwd('/home/jovyan/HB_ZIK/')

gsc_10x = readRDS("../data/jhubData/HB_ZIK/HB_ZIK/zikaGlioblastomas_10X_SeuratObject.rds")
gsc_SS = readRDS("../data/jhubData/HB_ZIK/HB_ZIK/zikaGlioblastomas_SS_SeuratObject.rds")

gsc_sample_names_SS = as.character(read.delim(file = 'gsc_samplenames_SS.csv', sep = ',', row.names = 1, stringsAsFactors = F)[,1])
gsc_sample_names_10x = as.character(read.delim(file = 'gsc_samplenames.csv', sep = ',', row.names = 1, stringsAsFactors = F)[,1])

gsc_10x = gsc_10x[,gsc_sample_names_10x]
gsc_SS = gsc_SS[,gsc_sample_names_SS]

# 13799 cells


# Load Neftel modules:
options(stringsAsFactors = FALSE)
modules = read.delim('/home/jovyan/HB_ZIK/1-s2.0-S0092867419306877-mmc2.csv', skip = 4, sep = ',')
MESlike = c(modules['MES2'][modules['MES2'] != ''], modules['MES1'][modules['MES1'] != ''])
AClike = modules['AC'][modules['AC'] != '']
OPClike = modules['OPC'][modules['OPC'] != '']
NPClike = c(modules['NPC2'][modules['NPC2'] != ''], modules['NPC1'][modules['NPC1'] != ''])
# Again to be consistent with the Neftel paper note which genes have more than 4 average log2 counts:
normCounts = exp(as.matrix(GetAssayData(object = gsc_SS, slot = "data")))
average = log2(rowMeans(normCounts))
keepGenes = average > 4
rm(average, normCounts)
# Define relative expression:
meanExpression = rowMeans(as.matrix(GetAssayData(object = gsc_SS, slot = "data")))[keepGenes]
relativeExpression = as.matrix(GetAssayData(object = gsc_SS, slot = "data"))[keepGenes,]-meanExpression 
# Bin genes into 30 groups based on expression:
b = seq(from = min(meanExpression)-0.0001, to = max(meanExpression)+0.0001, length.out = 30)
bins = .bincode(meanExpression, b, TRUE)
names(bins) = names(meanExpression)
# Subset signatures to include only available genes:
MESlike = MESlike[MESlike %in% names(meanExpression)]
AClike = AClike[AClike %in% names(meanExpression)]
OPClike = OPClike[OPClike %in% names(meanExpression)]
NPClike = NPClike[NPClike %in% names(meanExpression)]
# Extract control gene set for each signature:
MESlike_bins = bins[MESlike]
AClike_bins = bins[AClike]
OPClike_bins = bins[OPClike]
NPClike_bins = bins[NPClike]
bins_reduced = bins[!names(bins) %in% c(MESlike, AClike, OPClike, NPClike)]
MESlike_control = unname(unlist(lapply(MESlike_bins, function(b) sample(names(bins_reduced[bins_reduced == b]), 100, replace = T))))
AClike_control = unname(unlist(lapply(AClike_bins, function(b) sample(names(bins_reduced[bins_reduced == b]), 100, replace = T))))
OPClike_control = unname(unlist(lapply(OPClike_bins, function(b) sample(names(bins_reduced[bins_reduced == b]), 100, replace = T))))
NPClike_control = unname(unlist(lapply(NPClike_bins, function(b) sample(names(bins_reduced[bins_reduced == b]), 100, replace = T))))
# Calculate score:
MESlike_score = colMeans(relativeExpression[MESlike,]) - colMeans(relativeExpression[MESlike_control,]) 
AClike_score = colMeans(relativeExpression[AClike,]) - colMeans(relativeExpression[AClike_control,]) 
OPClike_score = colMeans(relativeExpression[OPClike,]) - colMeans(relativeExpression[OPClike_control,]) 
NPClike_score = colMeans(relativeExpression[NPClike,]) - colMeans(relativeExpression[NPClike_control,]) 
write.csv(as.data.frame(cbind(MESlike_score, AClike_score, OPClike_score, NPClike_score)), file = 'gsc_SS_scores.csv')

# Load Neftel modules:
options(stringsAsFactors = FALSE)
modules = read.delim('/home/jovyan/HB_ZIK/1-s2.0-S0092867419306877-mmc2.csv', skip = 4, sep = ',')
MESlike = c(modules['MES2'][modules['MES2'] != ''], modules['MES1'][modules['MES1'] != ''])
AClike = modules['AC'][modules['AC'] != '']
OPClike = modules['OPC'][modules['OPC'] != '']
NPClike = c(modules['NPC2'][modules['NPC2'] != ''], modules['NPC1'][modules['NPC1'] != ''])
# Again to be consistent with the Neftel paper note which genes have more than 4 average log2 counts:
normCounts = exp(as.matrix(GetAssayData(object = gsc_10x, slot = "data")))
average = log2(rowMeans(normCounts))
keepGenes = average > 4
rm(average, normCounts)
# Define relative expression:
meanExpression = rowMeans(as.matrix(GetAssayData(object = gsc_10x, slot = "data")))[keepGenes]
relativeExpression = as.matrix(GetAssayData(object = gsc_10x, slot = "data"))[keepGenes,]-meanExpression 
# Bin genes into 30 groups based on expression:
b = seq(from = min(meanExpression)-0.0001, to = max(meanExpression)+0.0001, length.out = 30)
bins = .bincode(meanExpression, b, TRUE)
names(bins) = names(meanExpression)
# Subset signatures to include only available genes:
MESlike = MESlike[MESlike %in% names(meanExpression)]
AClike = AClike[AClike %in% names(meanExpression)]
OPClike = OPClike[OPClike %in% names(meanExpression)]
NPClike = NPClike[NPClike %in% names(meanExpression)]
# Extract control gene set for each signature:
MESlike_bins = bins[MESlike]
AClike_bins = bins[AClike]
OPClike_bins = bins[OPClike]
NPClike_bins = bins[NPClike]
bins_reduced = bins[!names(bins) %in% c(MESlike, AClike, OPClike, NPClike)]
MESlike_control = unname(unlist(lapply(MESlike_bins, function(b) sample(names(bins_reduced[bins_reduced == b]), 100, replace = T))))
AClike_control = unname(unlist(lapply(AClike_bins, function(b) sample(names(bins_reduced[bins_reduced == b]), 100, replace = T))))
OPClike_control = unname(unlist(lapply(OPClike_bins, function(b) sample(names(bins_reduced[bins_reduced == b]), 100, replace = T))))
NPClike_control = unname(unlist(lapply(NPClike_bins, function(b) sample(names(bins_reduced[bins_reduced == b]), 100, replace = T))))
# Calculate score:
MESlike_score = colMeans(relativeExpression[MESlike,]) - colMeans(relativeExpression[MESlike_control,]) 
AClike_score = colMeans(relativeExpression[AClike,]) - colMeans(relativeExpression[AClike_control,]) 
OPClike_score = colMeans(relativeExpression[OPClike,]) - colMeans(relativeExpression[OPClike_control,]) 
NPClike_score = colMeans(relativeExpression[NPClike,]) - colMeans(relativeExpression[NPClike_control,]) 
write.csv(as.data.frame(cbind(MESlike_score, AClike_score, OPClike_score, NPClike_score)), file = 'gsc_10x_scores.csv')



