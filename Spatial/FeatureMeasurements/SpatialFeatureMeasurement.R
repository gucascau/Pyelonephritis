
##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@nationwidechildrens.org                         
# Copyright (c) 2023 Kidney and Urology Tract Center
# Nationwide Children's Hospital
# Description: 
#   Basic spatistics 
#    - Cluster Feature plots (inlcuding, cutting infected areas 2,5,10)
#    - Cluster propostion under each condition
#    - Gene expression under each condition (Bubble Plots)
#    - Progeny values at each condition
#    - 
################################################################

## we updated Seurat version
install.packages('Seurat')
install.packages('ggplot2')
#remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(tidyverse)
library(SeuratObject)
#install.packages("harmony")
library(harmony)


library(SPOTlight)
library(scater)
library(scran)
library(SingleCellExperiment)

#install.packages("SeuratDisk")
library(SeuratDisk)
library(tidyverse)
library(cowplot)

library(ggplot2)
library(reshape2)
library(Hmisc)
library(stats)
library(corrplot)
library(RColorBrewer)
## Deconvolve
#install.packages("L1pack")
#install.packages("reshape")
library(L1pack)
library(reshape)
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)
# if (!("xbioc" %in% rownames(inst))) {
#   remotes::install_github("renozao/xbioc", dependencies = FALSE)
# }
# if (!("SCDC" %in% rownames(inst))) {
#   remotes::install_github("meichendong/SCDC", dependencies = FALSE)
# }

suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))

### alllow large memory
library(usethis) 
#usethis::edit_r_environ()

outdir = "/Users/XXW004/Documents/Projects/RuizRosado/SpatialTranscriptomic/Results/Integration_v04142023/FeatureStatisics/"
dir.create(outdir, showWarnings = F)

# to list available datasets in SeuratData you can run AvailableData()

# first we dowload the dataset
# if (!("stxBrain.SeuratData" %in% rownames(InstalledData()))) {
#   InstallData("stxBrain")
# }


# # now we can list what datasets we have downloaded
# InstalledData()
# 
# brain1 <- LoadData("stxBrain", type = "anterior1")
# 
# head(brain2[[]])
# brain2 <- LoadData("stxBrain", type = "posterior1")
# head(brain2[[]])
# brain <- merge(brain1, brain2)
# head(brain2[[]])
# 

setwd("/Users/XXW004/Documents/Projects/RuizRosado/SpatialTranscriptomic/Results/Integration_v04142023/FeatureStatisics/")
#######################################################
# 1. laoding 10X Genomics Visium dataset
#######################################################
#Read in 10X Image and Spatial Data
#aa<-Read10X_Image(image.dir="Pyelonephritis_31/spatial/", 
#                 image.name = "", filter.matrix=TRUE)




### due to the integration take a long time, we ignore the previous and rea accordingly
Pyelonephritis.integrated<-readRDS(file = "../../PreProcess/V03222023_PreliminaryResults/PreliminaryIntegratedSptatil_Pyelonephritis.RDS")
## save the integrated 
## change the order:
?SpatialFeaturePlot

levels(factor(Pyelonephritis.integrated$orig.ident))

Pyelonephritis.integrated$orig.ident<-ordered(factor(Pyelonephritis.integrated$orig.ident), levels=c("PyeloD0_1","PyeloD1_1","PyeloD1_2","PyeloD3_1","PyeloD3_2","PyeloD5_1","PyeloD5_2","PyeloD7_1","PyeloD28_1","PyeloD28_2","PyeloD56_1"))
levels(factor(Pyelonephritis.integrated$orig.ident))

pdf("SpatialFeaturePlots_Pyeloneprhitis_integrated_NOFeature_v0323.pdf", height = 4, width = 40)
PSpatialFeature_NO <-SpatialFeaturePlot(Pyelonephritis.integrated, features = NULL) + theme(legend.position = "top")
PSpatialFeature_NO

dev.off()

pdf("SpatialFeaturePlots_Pyeloneprhitis_integrated_NOFeature_2_v0323.pdf", height = 4, width = 40)
PSpatialFeature_NO <-SpatialFeaturePlot(Pyelonephritis.integrated, features = NULL, stroke = 0) + theme(legend.position = "top")
PSpatialFeature_NO

dev.off()

#Pyelonephritis.integrated<- RunHarmony(Pyelonephritis.integrated, group.by.vars = "orig.ident")
Pyelonephritis.integrated <- RunPCA(Pyelonephritis.integrated, verbose = FALSE)
Pyelonephritis.integrated <- FindNeighbors(Pyelonephritis.integrated, dims = 1:15)
Pyelonephritis.integrated <- FindClusters(Pyelonephritis.integrated, verbose = FALSE, resolution = 0.3)
Pyelonephritis.integrated <- RunUMAP(Pyelonephritis.integrated, dims = 1:15)

SpatialDimPlot(Pyelonephritis.integrated)
#### Run harmony to remove the 
# Pyelonephritis.integrated_harmony<- RunHarmony(Pyelonephritis.integrated, group.by.vars = "orig.ident")
# Pyelonephritis.integrated_harmony <- RunPCA(Pyelonephritis.integrated_harmony, verbose = FALSE)
# Pyelonephritis.integrated_harmony <- FindNeighbors(Pyelonephritis.integrated_harmony, dims = 1:15)
# Pyelonephritis.integrated_harmony <- FindClusters(Pyelonephritis.integrated_harmony, verbose = FALSE, resolution = 0.4)
# Pyelonephritis.integrated_harmony <- RunUMAP(Pyelonephritis.integrated_harmony,reduction = "harmony", dims = 1:15)
# 

## find all markers

Idents(Pyelonephritis.integrated)

#######################################################
##  2. basical details
#######################################################
#2.1 dimension of the activae assay here is the number of features (genes) by samples (spots)
# dim(P31)
# dim(P32)
# dim(P33)
# dim(P34)
# dim(P41)
# dim(P42)
# dim(P43)
# dim(Pyelonephritis)
### 2.1 dimension of the active assay # the number of feature (genes) by spots (spots)
dim(x=Pyelonephritis.integrated)

nrow(x=Pyelonephritis.integrated) ## the number of features
ncol(x=Pyelonephritis.integrated) ### the number of samples

### 2.2 feature levels and sample name
## check how many spots each spatial
table(Pyelonephritis.integrated$orig.ident)
# calculate the feature (gene) and sample names
head(x=rownames(Pyelonephritis.integrated),n=5)
tail(x=rownames(Pyelonephritis.integrated),n=5)
## sample name 
tail(x=colnames(Pyelonephritis.integrated),n=5)

## 2.3. Sample level metadata
class(Pyelonephritis.integrated[[]])
colnames(Pyelonephritis.integrated[[]])
head(Pyelonephritis.integrated@meta.data)

## nFeature spatial; the number of unique gene in each sample

sum(Pyelonephritis.integrated$nFeature_Spatial == colSums(Pyelonephritis.integrated@assays$SCT@counts >0))


## consist of one or more assay objects as individual represetation of single-cell expression data,
## Each Assay stores multiple slots, including raw (counts), normalized (data) and scaled data (scaled.data) as well as a vector of variable features (var.features) and feature-level metadata (meta.features).
#Pyelonephritis.integrated@assays$Spatial@meta.features

## check the mitochondrial gene list
#grep (pattern = "^hb", x=(rownames(P31)), ignore.case=TRUE, value = TRUE)

# 2.3 sample level metadata
class(P31[[]])
colnames(P31[[]])
colnames(Pyelonephritis.integrated[[]])
tail(rownames(Pyelonephritis.integrated[[]]))
# 2.4 the number of unique genes in each sample
#colSums(P31@assays$Spatial@counts)
sum(P31$nFeature_Spatial == colSums(P31@assays$Spatial@counts>0))
# 2.5 the total numebr of detected molecules in each sample
sum(P31$nCount_Spatial == colSums(P31@assays$Spatial@counts))

## 2.4 objects (e.g. Assaty) together with feature level metadata
# a vector of names of assiacted objects 
names(x=P31)

# Spatial details
Pyelonephritis.integrated[['Spatial']] # equivalent to P31@assays$Spatial
# slice detais for image key
Pyelonephritis.integrated[['
                           ']] # equivalent to P31@assays$slice1

Pyelonephritis.integrated[['Spatial']]

# each seurat object consists of one or more assay objects (basis unit of seurat ) as individual representation of single-cell expression data
# example of the assay objuect include RNA, ADT in CITE-seq, or Spatial. Each assay store multiple slots, including raw (counts), normalized (data) and scaled data as well as 
# a vector of varible features (var.features) and feature-level metadata
Pyelonephritis.integrated@assays$Spatial ## 
Pyelonephritis.integrated[['Spatial']]@counts[5:10,1:3]
Pyelonephritis.integrated[['Spatial']]@counts[5:10,1:3]
# Feature level metadata is associated with each individual assay



## quality control

#####################################################################################################################
### 3. Generated a matrix for the control datasets to evaluate the number of spatial counts and number of features
#####################################################################################################################

#Then we run dimensionality reduction and clustering as before.
DefaultAssay(Pyelonephritis.integrated) <- "SCT"

### 
Pyelonephritis.integrated[[]]

# quality control

Pyelonephritis.integrated<-PercentageFeatureSet(Pyelonephritis.integrated, "^mt-", col.name = "percent_mito")
Pyelonephritis.integrated<- PercentageFeatureSet(Pyelonephritis.integrated, "^Hb.*-", col.name = "percent_hb")
Pyelonephritis.integrated<- PercentageFeatureSet(Pyelonephritis.integrated, "^Rp.*", col.name = "percent_rRNA")

# add the group details and change the order of the groups


Pyelonephritis.integrated[["group"]] = Pyelonephritis.integrated$orig.ident

levels(factor(Pyelonephritis.integrated[[]]$group))

## measure how many spots 
# create a list
TotalSpot<- data.frame(name=row.names(table(Pyelonephritis.integrated$group)), value=as.numeric(table(Pyelonephritis.integrated$group)))

sum(Pyelonephritis.integrated[[]]$nFeature_Spatial)/nrow(Pyelonephritis.integrated[[]])/ (sum(Pyelonephritis.integrated[[]]$nCount_Spatial)/sum(Pyelonephritis.integrated[[]]$nFeature_Spatial))

## change the order the rowname

TotalSpot$name<-ordered(factor(TotalSpot$name), levels=c("PyeloD0_1","PyeloD1_1","PyeloD1_2","PyeloD3_1","PyeloD3_2","PyeloD5_1","PyeloD5_2","PyeloD7_1","PyeloD28_1","PyeloD28_2","PyeloD56_1"))

TotalSpot
###use the default colors
library(scales)
###plot total number of spots
cols = gg_color_hue(11)
## plot the spots number distibution
pdf("SpotsPlots_Pyeloneprhitis_integrated_v0323.pdf", height = 3, width = 6)

#show_col(hue_pal()(16))

# These two are equivalent; by default scale_fill_hue() is used
# These two are equivalent; by default scale_colour_hue() is used
pspots<- ggplot(data=TotalSpot, aes(x=name, y=value, fill=name)) +geom_bar(stat="identity") +theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1))+ scale_fill_hue() + ylim(0, 3000)

#+scale_color_brewer(palette="Spectral")
print (pspots)

dev.off()

### plot the number of counts, numbe of features in the spatial 

pdf("FeaturePlots_Pyeloneprhitis_integrated_v0323.pdf", height = 4, width = 20)
PVlnFeature<- VlnPlot(Pyelonephritis.integrated, features = c("nCount_Spatial", "nFeature_Spatial", "percent_rRNA",
                                                              "percent_hb"), group.by ="group", pt.size = 0, ncol = 4) + NoLegend()

PVlnFeature
dev.off()

pdf("FeaturePlots_Pyeloneprhitis_integrated_withSpots_v0323.pdf", height = 4, width = 20)
PVlnFeature<- VlnPlot(Pyelonephritis.integrated, features = c("nCount_Spatial", "nFeature_Spatial", "percent_rRNA",
                                                              "percent_hb"), group.by ="group", pt.size = 0, ncol = 4) + NoLegend()

PVlnFeature + pspots
dev.off()

pdf("SpatialFeaturePlots_Pyeloneprhitis_integrated_SpatialFeatureCounts_v0323.pdf", height = 12, width = 40)
PSpatialFeature <-SpatialFeaturePlot(Pyelonephritis.integrated, features = c("nCount_Spatial", "nFeature_Spatial", "percent_rRNA","percent_hb")) 
# ggplot2::scale_fill_continuous(limits = c(0.0,1.0), breaks = c(0.0, 1.0))
## the following can scale the legend size
#++ ggplot2::scale_fill_continuous(limits = c(0.0,1.0),, breaks = c(0.0, 0.5, 1.0))
PSpatialFeature

dev.off()

?SpatialFeaturePlot


################################################################
# Fig S. draw the biomarkers expression of acute kidney injury 
################################################################

## here we scale the date to use the z-score for plotting the genes in multiple slides
# 
# eg on PBMC small dataset
library(magrittr)

## we here scale the data into z-score 0 -1
Pyelonephritis.integrated.scaled.mat <-
  sapply(rownames(Pyelonephritis.integrated) %>% seq,
         function (i) {
           scales::rescale(GetAssayData(Pyelonephritis.integrated, assay = "SCT", slot = "data")[rownames(Pyelonephritis.integrated)[i],], to = c(0,1)) }
  ) %>% t

# check expression ranges
quantile(Pyelonephritis.integrated.scaled.mat, c(0.01, 1))
# should output:
#1% 100%
#0    1

# add genes names
dimnames(Pyelonephritis.integrated.scaled.mat)[[1]] <- rownames(Pyelonephritis.integrated.scaled.mat)
# set new scaled matrix to the slot "scale.data"
Pyelonephritis.integrated.scaled.mat <- SetAssayData(Pyelonephritis.integrated, slot = "scale.data", Pyelonephritis.integrated.scaled.mat)

# now plot genes as needed.. using slot = "scale.data"

## several biomarkers: 

pdf("SpatialFeaturePlots_Pyeloneprhitis_KidneyInjuredBiomarker_v0323.pdf", height = 12, width = 40)

?SpatialFeaturePlot
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Lcn2","Cst3","Havcr1"),  ncol = 11, crop = TRUE) 

#+ ggplot2::scale_fill_continuous(limits = c(0,6), breaks = c(0,3,6)) +scale_colour_gradient(low = "white", high = "red")
# 
# SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Cst3"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
# 
# SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Havcr1"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)

?SpatialFeaturePlot
dev.off()

Pyelonephritis.integrated$group
VlnPlot(object = Pyelonephritis.integrated,features = c("Lcn2","Cst3","Havcr1"),  ncol = 1, split.by = "group",pt.size = 0)
#DotPlot(Pyelonephritis.integrated,features = c("Lcn2","Cst3","Havcr1"),  split.by = "group", cols = "RdBu")
################################################################
# Fig S. the cell-type cluster proportion across different time points
################################################################
Pyelonephritis.integrated@meta.data$seurat_clusters
### add the cluster proportion within all the datasets
table(Pyelonephritis.integrated@meta.data$group)
## 4.3 add the proportion of each cluster and show on the annotated clusters from all the datasets
length(Pyelonephritis.integrated@meta.data$seurat_clusters)
df<- data.frame(clu=names(table(Pyelonephritis.integrated@meta.data$seurat_clusters)), totalnumber=sprintf("%1.0f",table(Pyelonephritis.integrated@meta.data$seurat_clusters)), percentage= sprintf("%1.2f", 100*table(Pyelonephritis.integrated@meta.data$seurat_clusters)/length(Pyelonephritis.integrated@meta.data$seurat_clusters)))
names(table(Pyelonephritis.integrated@meta.data$seurat_clusters))
table(Pyelonephritis.integrated@meta.data$seurat_clusters)
df

## 5.4 add the proportion of each cluster on each slide

## we paste the group and cluster number 
paste0(Pyelonephritis.integrated@meta.data$group,":",Pyelonephritis.integrated@meta.data$seurat_clusters)

Pyelonephritis.integrated@meta.data$GroupCluster = paste0(Pyelonephritis.integrated@meta.data$group,":",Pyelonephritis.integrated@meta.data$seurat_clusters)
table(Pyelonephritis.integrated@meta.data$GroupCluster)

names(table(Pyelonephritis.integrated@meta.data$group))
sprintf("%1.0f",table(Pyelonephritis.integrated@meta.data$group))

table(Pyelonephritis.integrated@meta.data$group)

Pyelonephritis.integrated@meta.data$GroupCluster
dfSliceTotalnumber
# add the total number of each condition
dfSliceTotalnumber<- data.frame(Condition= names(table(Pyelonephritis.integrated@meta.data$group)), TotalSpotNumber = sprintf("%1.0f",table(Pyelonephritis.integrated@meta.data$group)))
# add the total number of each condition into meta data
Pyelonephritis.integrated@meta.data$TotalclusterNumber<-dfSliceTotalnumber[match(Pyelonephritis.integrated@meta.data$group,dfSliceTotalnumber$Condition ),2]
# add the number of each cluster number under different condition 
dfSliceClusternumber<- data.frame(ConditionCluster= names(table(Pyelonephritis.integrated@meta.data$GroupCluster)), TotalSpotNumberCluster = sprintf("%1.0f",table(Pyelonephritis.integrated@meta.data$GroupCluster)))

Pyelonephritis.integrated@meta.data$TotalclusterPerCondiiton <- dfSliceClusternumber[match(Pyelonephritis.integrated@meta.data$GroupCluster,dfSliceClusternumber$ConditionCluster),2]

# measure the proportion of each cluster under different condition
Pyelonephritis.integrated@meta.data$TotalclusterNumber

Pyelonephritis.integrated@meta.data$ClusterProportion<- 100*(as.numeric(Pyelonephritis.integrated@meta.data$TotalclusterPerCondiiton)/as.numeric(Pyelonephritis.integrated@meta.data$TotalclusterNumber))

## generate a matrix for each group

Pyelonephritis.integrated@meta.data %>% select(c("group","seurat_clusters","TotalclusterNumber",
                                                 "TotalclusterPerCondiiton","ClusterProportion"))

pdf("Frequency_split_bycondition_V0322.pdf", width = 15, height = 4)
ggplot(Pyelonephritis.integrated@meta.data, aes(fill=Pyelonephritis.integrated@meta.data$group, y=as.numeric(Pyelonephritis.integrated@meta.data$ClusterProportion), x=Pyelonephritis.integrated@meta.data$seurat_clusters)) + 
  geom_bar(position="dodge", stat="identity")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                  axis.text.x = element_text(vjust = 1, hjust=1),
                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits=c(0,40))+ xlab("Cluster number") +
  ylab ("Cluster frequency in each condition") + scale_fill_discrete(name = "Infected date")
dev.off()


pdf("Frequency_split_bycluster_V0322.pdf", width = 15, height = 4)
ggplot(Pyelonephritis.integrated@meta.data, aes(fill=Pyelonephritis.integrated@meta.data$seurat_clusters, y=as.numeric(Pyelonephritis.integrated@meta.data$ClusterProportion), x=Pyelonephritis.integrated@meta.data$group)) + 
  geom_bar(position="dodge", stat="identity")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                  axis.text.x = element_text(vjust = 1, hjust=1),
                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits=c(0,40)) 
dev.off()


## we only select the cluster 2, 5 and 10
Pyelonephritis.integrated@meta.data
Pyelonephritis.integrated.infected<-Pyelonephritis.integrated@meta.data %>% filter(`seurat_clusters` %in% c(2,5,10))
pdf("Frequency_split_bycondition_infectedarea_V0322.pdf", width = 6, height = 4)
ggplot(Pyelonephritis.integrated.infected, aes(fill=Pyelonephritis.integrated.infected$group, y=as.numeric(Pyelonephritis.integrated.infected$ClusterProportion), x=Pyelonephritis.integrated.infected$seurat_clusters)) + 
  geom_bar(position="dodge", stat="identity")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                  axis.text.x = element_text(vjust = 1, hjust=1),
                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits=c(0,40)) +
  xlab("Cluster number") +
  ylab ("Cluster frequency in each condition") + scale_fill_discrete(name = "Infected date")
dev.off()




strsplit(Pyelonephritis.integrated@meta.data$GroupCluster,":") 
Pyelonephritis.integrated@meta.data %>% filter (`group` == str_split(Pyelonephritis.integrated@meta.data$GroupCluster,":", simplify =TRUE) [,1])
(length(Pyelonephritis.integrated@meta.data %>% filter (`group` == str_split(Pyelonephritis.integrated@meta.data$GroupCluster,":", simplify =TRUE) [,1])))
### attention:
dfSlice <- data.frame(Date = names(table(Pyelonephritis.integrated@meta.data$GroupCluster)), totalSpotnumber=sprintf("%1.0f",table(Pyelonephritis.integrated@meta.data$GroupCluster)), 
                      Cluster= sprintf("%1.2f", 100*table(Pyelonephritis.integrated@meta.data$GroupCluster)/(Pyelonephritis.integrated@meta.data %>% filter (`group` == unlist(str_split(Pyelonephritis.integrated@meta.data$GroupCluster,":", simplify =FALSE))[1]))))

Pyelonephritis.integrated@meta.data$TotalSpotsCondition<- 
  
  names(table(Pyelonephritis.integrated@meta.data$group))
sprintf("%1.0f",table(Pyelonephritis.integrated@meta.data$group))

Pyelonephritis.integrated@meta.data$TotalclusterPercent <- df[match(Pyelonephritis.integrated@meta.data$seurat_clusters,df$clu),3]
Pyelonephritis.integrated@meta.data$TotalclusterNumber <-df[match(Pyelonephritis.integrated@meta.data$seurat_clusters,df$clu),2]

### add the cluster proportion within differnt time points
Pyelonephritis.integrated@meta.data$group
## add the proportion of each cluster at different time points
Pyelonephritis.integrated@meta.data$TimepointCluster<-paste0(Pyelonephritis.integrated@meta.data$group,":",Pyelonephritis.integrated@meta.data$seurat_clusters)
#length(Pyelonephritis.integrated@meta.data %>% filter (`group` == strsplit(Pyelonephritis.integrated@meta.data$TimepointCluster,":") [[1]][1]) %>% pull(group))

## make a data frame to store the different proportion
table(Pyelonephritis.integrated@meta.data$TimepointCluster)
TimeDf <- data.frame(Timecluster=names(table(Pyelonephritis.integrated@meta.data$TimepointCluster)), percentage= sprintf("%1.2f", 100*table(Pyelonephritis.integrated@meta.data$TimepointCluster)/length(Pyelonephritis.integrated@meta.data %>% filter (`group` == strsplit(Pyelonephritis.integrated@meta.data$TimepointCluster,":") [[1]][1]) %>% pull(group))))
TimeDf

## we then added the TimeClusterPercentage
Pyelonephritis.integrated@meta.data$TimeClusterPercentage <- TimeDf[match(Pyelonephritis.integrated@meta.data$TimepointCluster, TimeDf$Timecluster),2]
Pyelonephritis.integrated@meta.data

Idents(Pyelonephritis.integrated)

pdf("Frequency_split_cluster_V0322.pdf", width = 15, height = 4)
ggplot(Pyelonephritis.integrated@meta.data, aes(fill=Pyelonephritis.integrated@meta.data$group, y=as.numeric(Pyelonephritis.integrated@meta.data$TimeClusterPercentage), x=Idents(Pyelonephritis.integrated))) + 
  geom_bar(position="dodge", stat="identity")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                  axis.text.x = element_text(vjust = 1, hjust=1),
                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits=c(0,40)) 
dev.off()

pdf("FeaturePlots_Pyeloneprhitis_UmapDimPlot_Withlabels_V0322.pdf", width = 5, height = 4)
DimPlot(Pyelonephritis.integrated, label = TRUE,reduction = "umap",  pt.size=0.1)

dev.off()

#df[match(Pyelonephritis.integrated@meta.data$seurat_clusters,df$clu),3]
## set up the default to SCT

## enable joint dimensional reduction and clustering on the undelyiing RNA seq expression data
# DefaultAssay(Pyelonephritis.integrated) <- "SCT"
# 
# st.list = list(D0 = P0, D1_1=P31, D1_2=P32, D3_1=P33, D3_2=P34, D5_1=P41,D5_2=P42,D7_1=P11, D7_2=P12, D28=P28, D56=P43 )
# VariableFeatures(Pyelonephritis.integrated) <-c(VariableFeatures(P0), VariableFeatures(P31),VariableFeatures(P32),VariableFeatures(P33),VariableFeatures(P34),VariableFeatures(P41),VariableFeatures(P42), VariableFeatures(P43),VariableFeatures(P11),VariableFeatures(P12),VariableFeatures(P28))
# #VariableFeatures(Pyelonephritis.integrated) 
