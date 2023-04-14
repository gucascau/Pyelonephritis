##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@nationwidechildrens.org                         
# Copyright (c) 2023 Kidney and Urology Tract Center
# Nationwide Children's Hospital
# Description: 
#   Quality measurements of spatial trancriptomics for Pyelonephritis
#    - The number of spots per sample
#    - The numbe of UMI (per spot)
#    - The number of genes per spot
################################################################



library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(tidyverse)
library(SeuratObject)
install.packages("harmony")
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
install.packages("L1pack")
install.packages("reshape")
library(L1pack)
library(reshape)
if (!("xbioc" %in% rownames(inst))) {
  remotes::install_github("renozao/xbioc", dependencies = FALSE)
}
if (!("SCDC" %in% rownames(inst))) {
  remotes::install_github("meichendong/SCDC", dependencies = FALSE)
}

suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))

outdir = "/Users/XXW004/Documents/Projects/RuizRosado/SpatialTranscriptomic/Results/PreProcess/V03222023_PreliminaryResults/"
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


setwd("/Users/XXW004/Documents/Projects/RuizRosado/SpatialTranscriptomic/Results/PreProcess/V03222023_PreliminaryResults/")
#######################################################
# 1. laoding 10X Genomics Visium dataset
#######################################################
#Read in 10X Image and Spatial Data
#aa<-Read10X_Image(image.dir="Pyelonephritis_31/spatial/", 
 #                 image.name = "", filter.matrix=TRUE)





P0<-Load10X_Spatial(data.dir="../Datasets_V04282022/ControlSample/", 
                    filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                    slice="slice1",  filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P0[['orig.ident']] <-"PyeloD0_1"

P31<-Load10X_Spatial(data.dir="../WithSlideID_V0111/Pyelonephritis_31/", 
                            filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                            slice="slice2",  filter.matrix=TRUE, to.upper=FALSE, image=NULL)

P31[['orig.ident']] <-"PyeloD1_1"

#bb<-Read10X_Image(image.dir="acute_chronic_pyelo/seurat/",
#                 image.name = "tissue_lowres_image.png", filter.matrix=TRUE)
P32<-Load10X_Spatial(data.dir="../WithSlideID_V0111/Pyelonephritis_32", 
                                     filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                                     slice="slice3", filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P32[['orig.ident']] <-"PyeloD1_2"
#dd<-Read10X_Image(image.dir="chronic/seurat/",
 #                 image.name = "tissue_lowres_image.png", filter.matrix=TRUE)
P33<-Load10X_Spatial(data.dir="../WithSlideID_V0111/Pyelonephritis_33", 
                               filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                               slice="slice4", filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P33[['orig.ident']] <-"PyeloD3_1"
#cc<-Read10X_Image(image.dir="control/seurat/",
 #                 image.name = "tissue_lowres_image.png", filter.matrix=TRUE)
P34<-Load10X_Spatial(data.dir="../WithSlideID_V0111/Pyelonephritis_34", 
                               filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                               slice="slice5", filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P34[['orig.ident']] <-"PyeloD3_2"
#aa<-Read10X_Image(image.dir="Pyelonephritis_31/spatial/", 
#                 image.name = "", filter.matrix=TRUE)
P41<-Load10X_Spatial(data.dir="../WithSlideID_V0111/Pyelonephritis_41/", 
                     filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                     slice="slice6", filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P41[['orig.ident']] <-"PyeloD5_1"
#bb<-Read10X_Image(image.dir="acute_chronic_pyelo/seurat/",
#                 image.name = "tissue_lowres_image.png", filter.matrix=TRUE)
P42<-Load10X_Spatial(data.dir="../WithSlideID_V0111/Pyelonephritis_42", 
                     filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                     slice="slice7", filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P42[['orig.ident']] <-"PyeloD5_2"


P11<-Load10X_Spatial(data.dir="../Datasets_V04282022/AcutePyelo/", 
                    filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                    slice="slice8",  filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P11[['orig.ident']] <-"PyeloD7_1"

P12<-Load10X_Spatial(data.dir="../Datasets_V04282022/AcutePyeloneph/", 
                    filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                    slice="slice9",  filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P12[['orig.ident']] <-"PyeloD7_2"

P28<-Load10X_Spatial(data.dir="../Datasets_V04282022/Chronicpyelo/", 
                     filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                     slice="slice10",  filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P28[['orig.ident']] <-"PyeloD28_1"

#dd<-Read10X_Image(image.dir="chronic/seurat/",
#                 image.name = "tissue_lowres_image.png", filter.matrix=TRUE)
P43<-Load10X_Spatial(data.dir="../WithSlideID_V0111/Pyelonephritis_43", 
                     filename="filtered_feature_bc_matrix.h5", assay="Spatial", 
                     slice="slice11", filter.matrix=TRUE, to.upper=FALSE, image=NULL)
P43[['orig.ident']] <-"PyeloD56_1"



# spatial normalization
P0 <- SCTransform(P0, assay = "Spatial", verbose = FALSE)

P31 <- SCTransform(P31, assay = "Spatial", verbose = FALSE)
#P31 <- NormalizeData(P31, verbose = FALSE, assay = "Spatial")

P32 <- SCTransform(P32, assay = "Spatial", verbose = FALSE)
#P32 <- NormalizeData(P32, verbose = FALSE, assay = "Spatial")

P33 <- SCTransform(P33, assay = "Spatial", verbose = FALSE)
#P33 <- NormalizeData(P33, verbose = FALSE, assay = "Spatial")

P34 <- SCTransform(P34, assay = "Spatial", verbose = FALSE)
#P34 <- NormalizeData(P34, verbose = FALSE, assay = "Spatial")

P41 <- SCTransform(P41, assay = "Spatial", verbose = FALSE)
#P41 <- NormalizeData(P41, verbose = FALSE, assay = "Spatial")

P42 <- SCTransform(P42, assay = "Spatial", verbose = FALSE)
#P42 <- NormalizeData(P42, verbose = FALSE, assay = "Spatial")

P11 <- SCTransform(P11, assay = "Spatial", verbose = FALSE)

P12 <- SCTransform(P12, assay = "Spatial", verbose = FALSE)

P28 <- SCTransform(P28, assay = "Spatial", verbose = FALSE)

P43 <- SCTransform(P43, assay = "Spatial", verbose = FALSE)
#P43 <- NormalizeData(P43, verbose = FALSE, assay = "Spatial")

### due to the fact that there are stong batch effects betweeen different ST sections, so it may be a good idea to integrate the data across sections

st.list = list(D0 = P0, D1_1=P31, D1_2=P32, D3_1=P33, D3_2=P34, D5_1=P41,D5_2=P42,D7_1=P11, D7_2=P12, D28=P28, D56=P43 )

# run SCT on both datasets
st.list = lapply(st.list, SCTransform, assay = "Spatial", method = "poisson")

# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000000 * 1024^2)  # set allowed size to 2K MiB


# now we perform the actual integration
st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
                              verbose = FALSE)

int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = st.features)
Pyelonephritis.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                                  verbose = FALSE)

#rm(int.anchors, st.list)
#gc()

saveRDS(Pyelonephritis.integrated, file =sprintf("%s/PreliminaryIntegratedSptatil_Pyelonephritis.RDS",outdir))



### due to the integration take a long time, we ignore the previous and rea accordingly
Pyelonephritis.integrated<-readRDS(file = "PreliminaryIntegratedSptatil_Pyelonephritis.RDS")
## save the integrated 
## change the order:

levels(factor(Pyelonephritis.integrated$orig.ident))

Pyelonephritis.integrated$orig.ident<-ordered(factor(Pyelonephritis.integrated$orig.ident), levels=c("PyeloD0_1","PyeloD1_1","PyeloD1_2","PyeloD3_1","PyeloD3_2","PyeloD5_1","PyeloD5_2","PyeloD7_1","PyeloD7_2","PyeloD28_1","PyeloD56_1"))
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
Pyelonephritis.integrated <- FindClusters(Pyelonephritis.integrated, verbose = FALSE, resolution = 0.2)
Pyelonephritis.integrated <- RunUMAP(Pyelonephritis.integrated, dims = 1:15)


#### Run harmony to remove the 
# Pyelonephritis.integrated_harmony<- RunHarmony(Pyelonephritis.integrated, group.by.vars = "orig.ident")
# Pyelonephritis.integrated_harmony <- RunPCA(Pyelonephritis.integrated_harmony, verbose = FALSE)
# Pyelonephritis.integrated_harmony <- FindNeighbors(Pyelonephritis.integrated_harmony, dims = 1:15)
# Pyelonephritis.integrated_harmony <- FindClusters(Pyelonephritis.integrated_harmony, verbose = FALSE, resolution = 0.4)
# Pyelonephritis.integrated_harmony <- RunUMAP(Pyelonephritis.integrated_harmony,reduction = "harmony", dims = 1:15)
# 

### hereis to have low number of clusters

saveRDS(Pyelonephritis.integrated, file =sprintf("%s/PreliminaryIntegratedSptatil_Pyelonephritis_afterUMAP.RDS",outdir))


pdf("SpatialIntegrated_Pyeloneprhitis_Dimplot.pdf", height = 4, width = 9)
DimPlot(Pyelonephritis.integrated, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()

pdf("SpatialIntegrated_Pyeloneprhitis_Dimplot_Seperate.pdf", height = 4, width = 30)
DimPlot(Pyelonephritis.integrated, reduction = "umap", split.by ="orig.ident" )
dev.off()

pdf("SpatialIntegrated_Pyeloneprhitis_Dimplot_Seperate_rotated.pdf", height = 30, width = 4)
DimPlot(Pyelonephritis.integrated, reduction = "umap", split.by ="orig.ident",ncol = 1 )
dev.off()

?DimPlot
pdf("SpatialIntegrated_Pyeloneprhitis_SpatialDimplot.pdf", height = 8, width = 40)
SpatialDimPlot(Pyelonephritis.integrated)
dev.off()



pdf("SpatialIntegrated_Pyeloneprhitis_SpatialDimplot_withlabels.pdf", height = 8, width = 40)
SpatialDimPlot(Pyelonephritis.integrated, label = TRUE,label.size=2.5)
dev.off()

?SpatialDimPlot

pdf("SpatialIntegrated_Pyeloneprhitis_SpatialDimplot_withlabels_rotated.pdf", height = 40, width = 8)
SpatialDimPlot(Pyelonephritis.integrated, label = TRUE,label.size=2.5, ncol=1)
dev.off()

# 
# 
# #################################
# # then merge the objects
# #################################
# # merge did perform worse
# Pyelonephritis<-merge(P0,P31)
# Pyelonephritis<-merge(Pyelonephritis,P32)
# Pyelonephritis<-merge(Pyelonephritis,P33)
# Pyelonephritis<-merge(Pyelonephritis,P34)
# Pyelonephritis<-merge(Pyelonephritis,P41)
# Pyelonephritis<-merge(Pyelonephritis,P42)
# Pyelonephritis<-merge(Pyelonephritis,P11)
# Pyelonephritis<-merge(Pyelonephritis,P12)
# Pyelonephritis<-merge(Pyelonephritis,P28)
# Pyelonephritis<-merge(Pyelonephritis,P43)
# 
# ## save the combined prelimnary dataset in rds
# 
# saveRDS(Pyelonephritis, file =sprintf("%s/PreliminaryMergedSptatil_Pyelonephritis.RDS",outdir))
# 
# 
# ### loading the Pyelonephritis
# 
# Pyelonephritis<-readRDS(file = "PreliminaryMergedSptatil_Pyelonephritis.RDS")
# 

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
Pyelonephritis.integrated[['slice1']] # equivalent to P31@assays$slice1

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

## chagne the order the rowname

TotalSpot$name<-ordered(factor(TotalSpot$name), levels=c("PyeloD0_1","PyeloD1_1","PyeloD1_2","PyeloD3_1","PyeloD3_2","PyeloD5_1","PyeloD5_2","PyeloD7_1","PyeloD7_2","PyeloD28_1","PyeloD56_1"))

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
# Fig S. the cell-type cluster proportion across different time points
################################################################
Pyelonephritis.integrated@meta.data$seurat_clusters
### add the cluster proportion within all the datasets
table(Pyelonephritis.integrated@meta.data)
## 4.3 add the proportion of each cluster and show on the annotated clusters from all the datasets
length(Pyelonephritis.integrated@meta.data$seurat_clusters)
df<- data.frame(clu=names(table(Pyelonephritis.integrated@meta.data$seurat_clusters)), totalnumber=sprintf("%1.0f",table(Pyelonephritis.integrated@meta.data$seurat_clusters)), percentage= sprintf("%1.2f", 100*table(Pyelonephritis.integrated@meta.data$seurat_clusters)/length(Pyelonephritis.integrated@meta.data$seurat_clusters)))
names(table(Pyelonephritis.integrated@meta.data$seurat_clusters))
table(Pyelonephritis.integrated@meta.data$seurat_clusters)
df
## then add the proportion into the meta data
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
# 

#####################################################################################################################
### 4. Normalization the variance of molecular counts expresses, dimension reduction
#####################################################################################################################
### without the reference
## normalized bu SCTranform and check the structure
## SCTransform returns the Seurat object with a new assay called SCT, where its counts slot stores the corrected UMI counts, the data slot stores the log-normalized version of the corrected UMI counts, and scale.data slot stores the pearson residuals (normalized values) and is used as PCA input.

dim(Pyelonephritis.integrated@assays$SCT@counts)
dim(Pyelonephritis.integrated@assays$SCT@data) 
dim(Pyelonephritis.integrated@assays$SCT@scale.data) 
#

## Computes the correlation of the log normalized data and sctransform residuals with the
# # number of UMIs
# Pyelonephritis_norm <- GroupCorrelation(Pyelonephritis_norm, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
# Pyelonephritis_norm <- GroupCorrelation(Pyelonephritis_norm, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)
# 
# p1 <- GroupCorrelationPlot(Pyelonephritis_norm, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2 <- GroupCorrelationPlot(Pyelonephritis_norm, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
#   theme(plot.title = element_text(hjust = 0.5))
# p1 + p2

# ## Run PCA 
# 
# DefaultAssay(Pyelonephritis) <-"SCT"
# VariableFeatures(Pyelonephritis) <-c (VariableFeatures(P31),VariableFeatures(P32),VariableFeatures(P33),VariableFeatures(P34),VariableFeatures(P41),VariableFeatures(P42), VariableFeatures(P43))
# Pyelonephritis <- RunPCA(Pyelonephritis, assay = "SCT", verbose = FALSE)
# # compute K nearest neighbors (KNN)
# Pyelonephritis <- FindNeighbors(Pyelonephritis, reduction = "pca", dims = 1:20)
# # Leiden algorithm for community detection
# Pyelonephritis <- FindClusters(Pyelonephritis, verbose = FALSE)
# # PCA result is the default UMAP input, use dimensions 1:30 as input features
# Pyelonephritis <- RunUMAP(Pyelonephritis, reduction = "pca", dims = 1:20)
# 
# pdf(file = "Pyeonephritis_umap_clustered_Dimplot.pdf",  height = 12, width = 40)
# plot3 <- DimPlot(Pyelonephritis_norm, reduction = "umap", label = TRUE, split.by = "group") + NoLegend()
# plot3
# dev.off()
# pdf(file = "Pyeonephritis_Spactial_clustered_SpatialDimplot.pdf",  height = 12, width = 40)
# plot4 <- SpatialDimPlot(Pyelonephritis_norm, label = TRUE, label.size = 3) + NoLegend()
# plot4
# dev.off()

Idents(Pyelonephritis.integrated)

# save the final integrated into rds file
saveRDS(Pyelonephritis.integrated, file = "Pyelonephritis.integrated_Beforedeconvolution_V0325.rds")

# read the rds file

Pyelonephritis.integrated<- readRDS(file = "Pyelonephritis.integrated_Beforedeconvolution_V0325.rds")

Idents(Pyelonephritis.integrated)
## write the rds file into hd5

BiocManager::install(c("LoomExperiment", "SingleCellExperiment"))

library(SingleCellExperiment)
library(HDF5Array)



saveHDF5SummarizedExperiment(Pyelonephritis.integrated, dir = "./",
                             prefix = "", replace = FALSE,
                             chunkdim = NULL, level = NULL, as.sparse = NA,
                             verbose = NA)


# we have three options to perform deconvolution 
##############################################################################################
### First option: Using Seurat itself instead of cluster by themselves, we clusster based on the reference scRNA-seq using Seurat
##############################################################################################
# with reference
allen_reference <- readRDS("../../../../PIPseq/Results/Evaluation_V022723/CellEvaluation/Final_CombinedInformation_Pyelonephritis_V0321.RDS")

# here we set up 10000 normlizes the full datstet but learns noise models on 10k cells that speeds up SCTransform dramtially wiht no loss in performanace
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 10000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20)

levels(factor(allen_reference[[]]$old.ident))

dev.off()
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "old.ident", label = TRUE)

# After subsetting, we renormalize cortex
Pyelonephritis.integrated <- SCTransform(Pyelonephritis.integrated, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)


anchors <- FindTransferAnchors(reference = allen_reference, query = Pyelonephritis.integrated, normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$old.ident, prediction.assay = TRUE,
                                  weight.reduction = Pyelonephritis.integrated[["pca"]], dims = 1:20)
Pyelonephritis.integrated[["predictions"]] <- predictions.assay
DefaultAssay(Pyelonephritis.integrated) <- "predictions"



## check the row name and column name for the 
row.names(Pyelonephritis.integrated)
## the row name is the cell types.
## based on these preidction scores, we can also predict cell types whose location is patially restricted, we use the same methods based on marked point process to deinfe spatially varibale feautures, 
## but used the cell type preicton scores as the marks rather than gene expression
Pyelonephritis.integrated<- FindSpatiallyVariableFeatures(Pyelonephritis.integrated, assay = "predictions", selection.method = "markvariogram", features = rownames(Pyelonephritis.integrated), r.metric = 5, slot = "data")

top.clusters <- head(SpatiallyVariableFeatures(Pyelonephritis.integrated), 4)
SpatiallyVariableFeatures(Pyelonephritis.integrated)
top.clusters
SpatialPlot(object = Pyelonephritis.integrated, features = top.clusters, ncol = 2)

head(Pyelonephritis.integrated[["predictions"]]@meta.features)

Pyelonephritis.integrated[["predictions"]]@data
saveRDS(Pyelonephritis.integrated, file = "Pyelonephritis.integrated_withRefereceAnnotion_seurat_v0325.rds")

Pyelonephritis.integrated<-readRDS("Pyelonephritis.integrated_withRefereceAnnotion_seurat_v0325.rds")


### plot the 
Pyelonephritis.integrated[["predictions"]]@meta.features
#Finally, we show that our integrative procedure is capable of recovering the known spatial localization patterns of both kidney and non-neuronal subsets, 

# including laminar excitatory, layer-1 astrocytes, and the cortical grey matter.

#SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
DefaultAssay(Pyelonephritis.integrated)<-"SCT"
## check the cortex PT biomarkers
pdf("SpatialFeaturePlots_Pyeloneprhitis_CortexMarkers_Test.pdf", height = 40, width = 20)
CortexGenes<- c("Slc22a19","Serpina1f","Akr1c18","Mep1b","Napsa","Aadat","Kap","Havcr1")
SpatialFeaturePlot(Pyelonephritis.integrated, features = CortexGenes)
dev.off()
## generate the spatial dimplot

pdf("SpatialFeaturePlots_Pyeloneprhitis_Myeloblast_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Myeloblast"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_Endothelium_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Endothelium"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_Tcells_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("T cells"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()

pdf("SpatialFeaturePlots_Pyeloneprhitis_ProximaltubuleI_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Proximal tubule I"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_ProximaltubuleII_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Proximal tubule II"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()



pdf("SpatialFeaturePlots_Pyeloneprhitis_ProximaltubuleIII_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Proximal tubule III"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_ ProximaltubuleIV_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Proximal tubule IV"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()
# 

pdf("SpatialFeaturePlots_Pyeloneprhitis_ProximaltubuleV_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Proximal tubule V"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()

pdf("SpatialFeaturePlots_Pyeloneprhitis_ ProximaltubuleVI_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Proximal tubule VI"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()

# 

pdf("SpatialFeaturePlots_Pyeloneprhitis_Distaltubule_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Distal tubule"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_Intercalatedcells_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Intercalated cells"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()

pdf("SpatialFeaturePlots_Pyeloneprhitis_Principalcells_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Principal cells"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()



pdf("SpatialFeaturePlots_Pyeloneprhitis_ThinlimbofLOH_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Thin limb of LOH"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_ThickascendinglimbofLOH_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Thick ascending limb of LOH"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_Fibroblast_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Fibroblast"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_Bcells_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("B cells"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()



pdf("SpatialFeaturePlots_Pyeloneprhitis_NKcells_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("NK cells"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()



pdf("SpatialFeaturePlots_Pyeloneprhitis_LOH_Test.pdf", height = 8, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Loop of Henle"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


########################
### using the marker genes
### change the default assay to SCT
DefaultAssay(Pyelonephritis.integrated) <- "SCT"

pdf("SpatialFeaturePlots_Pyeloneprhitis_Proximal_tubule_WithGene_Test.pdf", height = 12, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Lrp2","Slc34a1","Acsm1"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()
pdf("SpatialFeaturePlots_Pyeloneprhitis_DitalTubule_WithGene_Test.pdf", height = 12, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Slc12a3","Pvalb","Wnk1"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


pdf("SpatialFeaturePlots_Pyeloneprhitis_IC_WithGene_Test.pdf", height = 12, width = 40)
SpatialFeaturePlot(Pyelonephritis.integrated, features = c("Atp6v1g3","Aqp6","Slc4a1"), pt.size.factor = 1.6, ncol = 11, crop = TRUE)
dev.off()


rownames(Pyelonephritis.integrated[["predictions"]]@meta.features)

str(Pyelonephritis.integrated)

######################################################################
# ### measure the proportion of each spots of cell types
### this stored in the matrix of Pyelonephritis.integrated[["predictions"]]@data
#################################################################

# the matrix for the predicted cell types and spot id
row.names(Pyelonephritis.integrated[["predictions"]]@data)

colnames(Pyelonephritis.integrated[["predictions"]]@data)

#Pyelonephritis.integrated[["predictions"]]@data[, ]
################
# the matrix for the cluster of the spot ID
Pyelonephritis.integrated@meta.data

# cluster ID Pyelonephritis.integrated$seurat_clusters
###

## merge the Pyelonephritis.integrated@meta.data and Pyelonephritis.integrated[["predictions"]]@data

t(Pyelonephritis.integrated[["predictions"]]@data)

Pyelonephritis.integrated@meta.data[c("AAACCGTTCGTCCAGG-1_1"),]
row.names(t(Pyelonephritis.integrated[["predictions"]]@data))

Table1 <- Pyelonephritis.integrated@meta.data[,c("orig.ident","seurat_clusters","TimepointCluster")]
## merge is way too troublesome, therefore
FinalIntergratedTable<- merge(Pyelonephritis.integrated@meta.data[,c("orig.ident","seurat_clusters","TimepointCluster")], t(Pyelonephritis.integrated[["predictions"]]@data), all=TRUE, by="row.names")


colnames(FinalIntergratedTable)

##########################################################################################
# The heat map figure to show the different compostion within timepoints in the cluster spots
###########################################################################################
### We combined all the dataset without group by orig.ident and measure the different composition

ClusterAllGroup_Median<-FinalIntergratedTable %>% group_by(seurat_clusters) %>% summarise(ProximaltubuleI=median(`Proximal tubule I`),  ProximaltubuleII=median(`Proximal tubule II`), ProximaltubuleIII=median(`Proximal tubule III`), ProximaltubuleIV= median(`Proximal tubule IV`),
                                                                                                           ProximaltubuleV=median(`Proximal tubule V`), ProximaltubuleVI=median(`Proximal tubule VI`), Distaltubule=median(`Distal tubule`),  ThinlimbofLOH=median(`Thin limb of LOH`), ThickascendinglimbofLOH=median(`Thick ascending limb of LOH`), 
                                                                                                           LoopofHenle=median(`Loop of Henle`), Intercalatedcells=median(`Intercalated cells`), Principalcells=median(`Principal cells`), 
                                                                                                           Fibroblast=median(`Fibroblast`), Endothelium=median(`Endothelium`), Myeloblast=median(`Myeloblast`), Tcells=median(`T cells`, na.rm=TRUE),  Bcells=median(`B cells`),  
                                                                                                           NKCells=median(`NK cells`),  .groups= 'drop') %>% as.data.frame()



ClusterAllGroup_Median

### we calculate the median value for each cluster and also based on the timepotins
### we used the median value to  scale the cell-type compositions within each niche. 
ClusterTimepointGroup_Median<-FinalIntergratedTable %>% group_by(orig.ident,seurat_clusters) %>% summarise(ProximaltubuleI=median(`Proximal tubule I`),  ProximaltubuleII=median(`Proximal tubule II`), ProximaltubuleIII=median(`Proximal tubule III`), ProximaltubuleIV= median(`Proximal tubule IV`),
                                                                                                           ProximaltubuleV=median(`Proximal tubule V`), ProximaltubuleVI=median(`Proximal tubule VI`), Distaltubule=median(`Distal tubule`),  ThinlimbofLOH=median(`Thin limb of LOH`), ThickascendinglimbofLOH=median(`Thick ascending limb of LOH`), 
                                                                                                           LoopofHenle=median(`Loop of Henle`), Intercalatedcells=median(`Intercalated cells`), Principalcells=median(`Principal cells`), 
                                                                                                           Fibroblast=median(`Fibroblast`), Endothelium=median(`Endothelium`), Myeloblast=median(`Myeloblast`), Tcells=median(`T cells`, na.rm=TRUE),  Bcells=median(`B cells`),  
                                                                                                           NKCells=median(`NK cells`),  .groups= 'drop') %>% as.data.frame()
#ClusterTimepointGroup_Median
### draw the heatmap for these combined clusters

Combined_Median_tidy <- ClusterAllGroup_Median %>% mutate_at( vars(-seurat_clusters),scale) %>% pivot_longer(cols = -c(seurat_clusters), names_to = "Property", values_to = "Value")

### reorder the levels
Combined_Median_tidy$Property<- ordered(factor(Combined_Median_tidy$Property), levels=c("ProximaltubuleI", "ProximaltubuleII", 
                                                                                  "ProximaltubuleIII"  ,  "ProximaltubuleIV",  "ProximaltubuleV", "ProximaltubuleVI",
                                                                                  "Principalcells", "Intercalatedcells",  "ThickascendinglimbofLOH", "ThinlimbofLOH",
                                                                                  "LoopofHenle","Distaltubule",  "Fibroblast", "Endothelium","Myeloblast" ,
                                                                                  "Tcells" ,"Bcells" , "NKCells"))
pdf("Combined_Median_tidy_Heatmap.pdf",height = 6, width = 8)
Combined_Median_tidyHeatmap<- na.omit(Combined_Median_tidy) %>% heatmap(Property, seurat_clusters, Value, scale="row",cluster_rows = FALSE, cluster_columns = TRUE ,palette_value=circlize::colorRamp2(seq(-4, 4, length.out = 11), RColorBrewer::brewer.pal(11, "RdBu") ) )
Combined_Median_tidyHeatmap
dev.off()


## setting up the conditions 
Condition<- levels(PyeloD0_1_Median_tidy$orig.ident)
## run the for liips to generate the heatmap figures

### generate multiple heatmap plots



#pdf("Combined_heatmaps_Seurat_0320.pdf", height = 4, width = 40)

## attention par() and layout() not work for ggplot and 
# par(mfrow=c(11,1))
## instead we should use grid.arrange()
## grid.arrange(plot1, plot2, ncol=2)

## create a new empty list and new number
i=0
my_vect <- c()
for (Pdate in Condition){
  ### generate the transfered objects for each condition
  ### 
 i =i+1
  Pdate_Median_tidy<- ClusterTimepointGroup_Median %>% filter(orig.ident==Pdate)  %>% mutate_at(vars (-orig.ident, -seurat_clusters),scale) %>% pivot_longer(cols = -c(orig.ident, seurat_clusters), names_to = "Property", values_to ="Value") 
  
  ### change the levels of Property at each condition
  ### tranfer the ClusterTimepointGroup_Median into a tide "element-feature-independent variables"  data frame, where the independent variable
  ### similar function to melt, but this will transfer the objects
  levels(factor(Pdate_Median_tidy$Property))
  Pdate_Median_tidy$Property<- ordered(factor(Pdate_Median_tidy$Property), levels=c("ProximaltubuleI", "ProximaltubuleII", 
                                                                                            "ProximaltubuleIII"  ,  "ProximaltubuleIV",  "ProximaltubuleV", "ProximaltubuleVI",
                                                                                            "Principalcells", "Intercalatedcells",  "ThickascendinglimbofLOH", "ThinlimbofLOH",
                                                                                            "LoopofHenle","Distaltubule",  "Fibroblast", "Endothelium","Myeloblast" ,
                                                                                            "Tcells" ,"Bcells" , "NKCells"))
  
  ## change the NA value to Zero
  #PyeloD0_1_Median_tidy$Value[is.nan(PyeloD0_1_Median_tidy$Value)]<-0
  # ignore the na.omit
  na.omit(PyeloD0_1_Median_tidy)
  
 pdf(paste0(Pdate,"_Median_tidyHeatmap.pdf"),height = 6, width = 8)
  #PyeloD3_1_Median_tidyHeatmap<- na.omit(PyeloD3_1_Median_tidy) %>% heatmap(Property, seurat_clusters, Value, scale="row",cluster_rows = FALSE, cluster_columns = TRUE ,palette_value = c( "#E6E6FA", "#FFFFFF","#F1948A")) 
  
  Pdate_Median_tidyHeatmap<- na.omit(Pdate_Median_tidy) %>% heatmap(Property, seurat_clusters, Value, scale="row",cluster_rows = FALSE, cluster_columns = FALSE ,palette_value=circlize::colorRamp2(seq(-4, 4, length.out = 11), RColorBrewer::brewer.pal(11, "RdBu") ) )
  
  #PyeloD3_1_Median_tidyHeatmap<- na.omit(PyeloD3_1_Median_tidy) %>% heatmap(Property, seurat_clusters, Value, scale="row",cluster_rows = FALSE, cluster_columns = FALSE , palette_value = c("blue","white","red"), color_legend_min=-4, color_legend_max=4)
  
  #my_vect<- append(my_vect,Pdate_Median_tidyHeatmap)
  #plots[[i]]<- Pdate_Median_tidyHeatmap
  print(Pdate_Median_tidyHeatmap)
  dev.off()
}

Condition
###############################################################################################
# spearman correlation of the celltype matrix 
###############################################################################################


for (Pdate in Condition){
  ### generate the transfered objects for each condition
  ### 
  #i =i+1
  #Pdate=c("PyeloD0_1")
  Pdate_Median_tidy_matrix <-ClusterTimepointGroup_Median %>% filter(orig.ident==Pdate) %>% select(-c('orig.ident','seurat_clusters'))
  
  #Pdate_Median_tidy_matrix
  #Pdate_Median_tidy_matrix<- ClusterTimepointGroup_Median %>% filter(orig.ident==Pdate) 
  
  #%>% mutate_at(vars (-orig.ident, -seurat_clusters),scale) %>% pivot_longer(cols = -c(orig.ident, seurat_clusters), names_to = "Property", values_to ="Value") 
  
 
  Pdate_Median_tidy_matrix_cor<-cor(as.matrix(Pdate_Median_tidy_matrix),method = "spearman")
  
  ## replace the NA into zerro
  Pdate_Median_tidy_matrix_cor <- replace(Pdate_Median_tidy_matrix_cor, is.na(Pdate_Median_tidy_matrix_cor),0)
  
  #testRes = cor.mtest(BiomarkerExpCor, conf.level = 0.95)
  Pdate_Median_tidy_matrix_cor_melt = melt(Pdate_Median_tidy_matrix_cor)
  
  ## create a dataframe
  
  #Pdate_Median_tidy_matrix_cor_melt$name= rep(Pdate, length(Pdate_Median_tidy_matrix_cor_melt$Var1))
 
  ## merge different dataframe
  
  pdf(paste0(Pdate,"_Correlation_tidyHeatmap.pdf"),height = 6, width = 8)
  #pdf("BioMarkerExpCorrelation.pdf", width = 10, height =10)
  ## generate the coeffecient figure
  corrplot(Pdate_Median_tidy_matrix_cor, order = 'hclust', type='upper',tl.col = 'black', insig='blank',
           addCoef.col ='black', number.cex = 0.5, tl.cex= 0.7, tl.srt = 45, diag=FALSE,col = colorRampPalette(c("LightBlue","white","FireBrick"))(100))
  dev.off()
  
  ## other way to plot
  # Pdate_Median_tidy_matrix_cor = rcorr(as.matrix(Pdate_Median_tidy_matrix), type='spearman')
  # 
  # Pdate_Median_tidy_matrix_cor_melt = melt(Pdate_Median_tidy_matrix_cor$r)
  # 
  # 
  # pdf(paste0(Pdate,"_Correlation_tidyHeatmap.pdf"),height = 6, width = 8)
  # 
  # #txtsize <- par('din')[2] / 2
  # ggplot(Pdate_Median_tidy_matrix_cor_melt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + 
  #   theme(axis.text.x = element_text(angle=90, hjust=TRUE))
  # dev.off()
   #+
   # xlab("") + ylab("") + 
   # geom_text(label=Pdate_Median_tidy_matrix_cor_melt$label, size=txtsize) + 
   # geom_text(label=Pdate_Median_tidy_matrix_cor_melt$strike, size=txtsize * 4, color="red", alpha=0.4) 

 
}



##############################################################################################
##############################################################################################
# we have three options to perform deconvolution 
##############################################################################################
### Second option: Using SCDC for deconvolution of cell types 
##############################################################################################
##############################################################################################
##############################################################################################


inst = installed.packages()

if (!("xbioc" %in% rownames(inst))) {
  remotes::install_github("renozao/xbioc", dependencies = FALSE)
}
if (!("SCDC" %in% rownames(inst))) {
  remotes::install_github("meichendong/SCDC", dependencies = FALSE)
}

suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))

#

# select gene for ddecovolution: 
# 1. use varible genes in the SC data
# 2. Use varible genes in both SC and ST data
# 3. DE genes between clusters in the SC data

# with reference
allen_reference <- readRDS("../../../../PIPseq/Results/Evaluation_V022723/CellEvaluation/Final_CombinedInformation_Pyelonephritis_V0321.RDS")

allen_reference[[]]
Idents(allen_reference)<- allen_reference$Annot

#Idents(allen_reference)
#strsplit(allen_reference$Annot,split=":")[2]
rownames(Pyelonephritis.integrated)
allen_reference@active.assay ="RNA"
markers_sc <- FindAllMarkers(allen_reference, only.pos = TRUE, logfc.threshold = 0.1,
                             test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
                             return.thresh = 0.05, assay = "RNA")

length(rownames(Pyelonephritis.integrated))

# Filter for genes that are also present in the ST data
markers_sc <- markers_sc[markers_sc$gene %in% rownames(Pyelonephritis.integrated), ]

table(markers_sc$cluster)
markers_sc
# Select top 20 genes per cluster, select top by first p-value, then absolute
# diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
                                                               1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(-100, p_val) %>%
  top_n(50, pct.diff) %>%
  top_n(20, log.pct.diff) -> top20

m_feats <- unique(as.character(top20$gene))
m_feats
### create expression sets
eset_SC <- ExpressionSet(assayData = as.matrix(allen_reference@assays$RNA@counts[m_feats,
]), phenoData = AnnotatedDataFrame(allen_reference@meta.data))

as.character(unique(eset_SC$Annot))

Pyelonephritis.integrated@assays$Spatial@counts
eset_ST <- ExpressionSet(assayData = as.matrix(Pyelonephritis.integrated@assays$Spatial@counts[m_feats,
]), phenoData = AnnotatedDataFrame(Pyelonephritis.integrated@meta.data))



deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "Annot",
                                     ct.sub = as.character(unique(eset_SC$Annot)))



deconvolution_crc

##############################################################################################
# we have three options to perform deconvolution 
##############################################################################################
### Second option: Using CARD for deconvolution of cell types 
##############################################################################################
##############################################################################################
##############################################################################################

devtools::install_github('YingMa0107/CARD')
# load package
library(SingleCellExperiment)
library(SummarizedExperiment)
library(concaveman)
library(sp)
library(Matrix)
library(methods)

library(ggcorrplot)
library(ggplot2)
install.packages("MuSiC")
library(MuSiC)
install.packages("CARD")
library(CARD)



#### the following is for timepoint samples




## tranfer the ClusterTimepointGroup_Median into a tide "element-feature-independent variables"  data frame, where the independent variable
## similar function to melt, but this will transfer the objects

PyeloD0_1_Median_tidy<- ClusterTimepointGroup_Median %>% filter(orig.ident=="PyeloD0_1")  %>% mutate_at(vars (-orig.ident, -seurat_clusters),scale) %>% pivot_longer(cols = -c(orig.ident, seurat_clusters), names_to = "Property", values_to ="Value") 

#factor(PyeloD0_1_Median_tidy)
levels(factor(PyeloD0_1_Median_tidy$Property))
PyeloD0_1_Median_tidy$Property<- ordered(factor(PyeloD0_1_Median_tidy$Property), levels=c("ProximaltubuleI", "ProximaltubuleII", 
                                                                                          "ProximaltubuleIII"  ,  "ProximaltubuleIV",  "ProximaltubuleV", "ProximaltubuleVI",
                                                                                          "Principalcells", "Intercalatedcells",  "ThickascendinglimbofLOH", "ThinlimbofLOH",
                                                                                          "LoopofHenle","Distaltubule",  "Fibroblast", "Endothelium","Myeloblast" ,
                                                                                          "Tcells" ,"Bcells" , "NKCells"))





                                                                                                                      
## change the NA value to Zero
#PyeloD0_1_Median_tidy$Value[is.nan(PyeloD0_1_Median_tidy$Value)]<-0
# ignore the na.omit
na.omit(PyeloD0_1_Median_tidy)

pdf("PyeloD0_1_Median_tidyHeatmap.pdf",height = 6, width = 8)
PyeloD0_1_Median_tidyHeatmap<- na.omit(PyeloD0_1_Median_tidy) %>% heatmap(Property, seurat_clusters, Value, scale="row",cluster_rows = FALSE, cluster_columns = TRUE ,palette_value=circlize::colorRamp2(seq(-4, 4, length.out = 11), RColorBrewer::brewer.pal(11, "RdBu") ) )
PyeloD0_1_Median_tidyHeatmap
dev.off()


### 


PyeloD3_1_Median_tidy<- ClusterTimepointGroup_Median %>% filter(orig.ident=="PyeloD3_1")  %>% mutate_at(vars (-orig.ident, -seurat_clusters),scale) %>% pivot_longer(cols = -c(orig.ident, seurat_clusters), names_to = "Property", values_to ="Value") 

levels(factor(PyeloD3_1_Median_tidy$Property))
PyeloD3_1_Median_tidy$Property<- ordered(factor(PyeloD3_1_Median_tidy$Property), levels=c("ProximaltubuleI", "ProximaltubuleII", 
                                                                                          "ProximaltubuleIII"  ,  "ProximaltubuleIV",  "ProximaltubuleV", "ProximaltubuleVI",
                                                                                          "Principalcells", "Intercalatedcells",  "ThickascendinglimbofLOH", "ThinlimbofLOH",
                                                                                          "LoopofHenle","Distaltubule",  "Fibroblast", "Endothelium","Myeloblast" ,
                                                                                          "Tcells" ,"Bcells" , "NKCells"))


## change the NA value to Zero
#PyeloD0_1_Median_tidy$Value[is.nan(PyeloD0_1_Median_tidy$Value)]<-0
# ignore the na.omit
na.omit(PyeloD0_1_Median_tidy)

pdf("PyeloD3_1_Median_tidyHeatmap.pdf",height = 6, width = 8)
#PyeloD3_1_Median_tidyHeatmap<- na.omit(PyeloD3_1_Median_tidy) %>% heatmap(Property, seurat_clusters, Value, scale="row",cluster_rows = FALSE, cluster_columns = TRUE ,palette_value = c( "#E6E6FA", "#FFFFFF","#F1948A")) 

PyeloD3_1_Median_tidyHeatmap<- na.omit(PyeloD3_1_Median_tidy) %>% heatmap(Property, seurat_clusters, Value, scale="row",cluster_rows = FALSE, cluster_columns = TRUE ,palette_value=circlize::colorRamp2(seq(-4, 4, length.out = 11), RColorBrewer::brewer.pal(11, "RdBu") ) )

#PyeloD3_1_Median_tidyHeatmap<- na.omit(PyeloD3_1_Median_tidy) %>% heatmap(Property, seurat_clusters, Value, scale="row",cluster_rows = FALSE, cluster_columns = FALSE , palette_value = c("blue","white","red"), color_legend_min=-4, color_legend_max=4)

PyeloD3_1_Median_tidyHeatmap
dev.off()

?heatmap

ClusterTimepointGroup_Median_tidy %>% filter(orig.ident=="PyeloD0_1") %>% drop_na() 
#print (ClusterTimepointGroup_Median_tidy %>% filter(orig.ident== "PyeloD0_1"),n=100)
  # Scale 




library(ggplot2)
hub (development)

devtools::install_github("stemangiola/tidySummarizedExperiment")
library(tidySummarizedExperiment)
## then reshape the cluster matrix into a varaible and value
ClusterTimepointGroup_MedianMelt<-reshape2::melt(ClusterTimepointGroup_Median,id.vars= c("orig.ident","seurat_clusters"))

### then use the tidyheatmap with mutlipe groups


TestDataframe<- ClusterTimepointGroup_MedianMelt |as_


ClusterTimepointGroup_MedianMeltHeatmap<-ClusterTimepointGroup_MedianMelt %>% group_by(orig.ident) %>% heatmap(seurat_clusters, variable, value, scale="row",palette_value = c("red", "white", "blue")) 

## filter the different and generate the heatmaps
PyeloD0_1<- ClusterTimepointGroup_MedianMelt %>% filter(orig.ident == "PyeloD0_1") 


## transfer the datasetinto a tide "element-feature-independent variable' data frame


PyeloD0_1|>heatmap(seurat_clusters, variable, value ,scale= "row")

heatmap(seurat_clusters, variable, value, scale="row",palette_value = c("red", "white", "blue")) 

#|>   add_tile(orig.ident, show_legend = FALSE)

devtools::install_github("stemangiola/tidyHeatmap")
library(tidyHeatmap)

## Here we generate multiple heatmaps grouped by the orig.ident

ClusterTimepointGroup_Median %>% filter(orig.ident == "PyeloD0_1") %>% mutate_at(vars(- `orig.ident`, -`seurat_clusters`), scale) %>% heatmap()




### we used the median value to scale the cellptype compositions without the timepoints, this is mainly for the main figure 

### we also calculate the P value using the wilcoxon test for each cell types within the niches


ClusterTimepointGroup_Median
df<- data.frame(clu=names(table(Pyelonephritis.integrated@meta.data$seurat_clusters)), totalnumber=sprintf("%1.0f",table(Pyelonephritis.integrated@meta.data$seurat_clusters)), percentage= sprintf("%1.2f", 100*table(Pyelonephritis.integrated@meta.data$seurat_clusters)/length(Pyelonephritis.integrated@meta.data$seurat_clusters)))
names(table(Pyelonephritis.integrated@meta.data$seurat_clusters))
table(Pyelonephritis.integrated@meta.data$seurat_clusters)
df
## then add the proportion into the meta data
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





##########################################################################################
### second option: deconvolve using the top DE genes per cluster, 
##########################################################################################


allen_reference@assays$RNA@counts[m_feats,]
## for SCDC both the SC and ST data neee to be in the format of an expression set with the countmatircs
eset_SC <- ExpressionSet(assayData = as.matrix(allen_reference@assays$RNA@counts[m_feats,]), phenoData = AnnotatedDataFrame(allen_reference@meta.data))



deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = Pyelonephritis.integrated, sc.eset = Pyelonephritis.integrated, ct.varname = "subclass",
                                     ct.sub = as.character(unique(Pyelonephritis.integrated$subclass)))

# now we have a mattrix with predict proportion of each cell types for each visium spot in pre.est.mvw

head(deconvolution_crc)



Pyelonephritis[['slice1']]@key
# 
levels(factor(allen_reference[[]]$old.ident))

Pyelonephritis.integrated@meta.data
Idents(Pyelonephritis.integrated)

pdf("SpatialIntegrated_Pyeloneprhitis_Dimplot_Seperate.pdf", height = 5, width = 30)
DimPlot(Pyelonephritis.integrated, reduction = "umap", split.by ="orig.ident" )
dev.off()

pdf("SpatialIntegrated_Pyeloneprhitis_SpatialDimplot.pdf", height = 12, width = 40)
SpatialDimPlot(Pyelonephritis.integrated)
dev.off()



# 
# # the tissue coordinate for single cell RNAs
# 
# coldata <- GetTissueCoordinates(controlsample,
#                                 cols = c("row", "col", "tissue"),
#                                 scale = NULL)
# coldata
# 
# # get the information whether on the tissue
# controlsample$tissue <- coldata[rownames(controlsample@meta.data), "tissue"]
# 
# controlsample$tissue <- ifelse(controlsample$tissue == 1, 
#                                "on_tissue", 
#                                "not_on_tissue")
# ## add the condition information (Four condtions)
# 
# controlsample$Condition<-replicate(n=length(controlsample$tissue), 'Control')
# 
# controlsample@meta.data
# # We do a comparison of profiles between spots in tissue and not in tissue
# 
# rm(coldata)
# 

### here also add the percentage of mt genes
controlsample[["percent.mt"]] <- PercentageFeatureSet(controlsample, 
                                                      pattern = "^MT-")

### only select the ones that located in the region
controlsample_seurat <- subset(controlsample, subset = tissue == "on_tissue")

SpatialFeaturePlot(controlsample_seurat,
                   features = c("nCount_Spatial"))
SpatialFeaturePlot
slide_qc_p_control <- SpatialFeaturePlot(controlsample_seurat,
                                 features = c("nCount_Spatial", 
                                              "nFeature_Spatial", 
                                              "percent.mt"),
                                 ncol = 3)

slide_qc_p_control

ControlQC<-controlsample@meta.data %>% select(-orig.ident)

#controlsample@meta.data %>% select(-orig.ident)
### Add a new column for the condition information

ControlQC
colnames(ControlQC)

#####################################################################################################################
### Generated a matrix for the acutepyelo, here we added the condition and tissue for the meta datatset
#####################################################################################################################
### 
# the tissue coordinate for single cell RNAs in acutepyelo
acutepyelocoldata <- GetTissueCoordinates(acutepyelo,
                                cols = c("row", "col", "tissue"),
                                scale = NULL)
acutepyelocoldata

# get the information whether on the tissue
acutepyelo$tissue <- acutepyelocoldata[rownames(acutepyelo@meta.data), "tissue"]

acutepyelo$tissue <- ifelse(acutepyelocoldata$tissue == 1, 
                               "on_tissue", 
                               "not_on_tissue")
acutepyelo$Condition<-replicate(n=length(acutepyelo$tissue), 'acutepyelo')

acutepyelo@meta.data

acutepyeloQC<-acutepyelo@meta.data %>% select(-orig.ident)

colnames(acutepyeloQC)

## add the condition information (Four condtions)

acutepyelo$Condition<-replicate(n=length(acutepyelo$tissue), 'acutepyelo')

acutepyelo@meta.data
# We do a comparison of profiles between spots in tissue and not in tissue
acutepyelo@meta.data %% select(tissue)
rm(coldata)

acutepyelo[["percent.mt"]] <- PercentageFeatureSet(acutepyelo, 
                                                      pattern = "^MT-")

### only select the ones that located in the region
acutepyelo_seurat <- subset(acutepyelo, subset = tissue == "on_tissue")

slide_qc_acutepyelo_control <- SpatialFeaturePlot(acutepyelo_seurat,
                                         features = c("nCount_Spatial", 
                                                      "nFeature_Spatial", 
                                                      "percent.mt"),
                                         ncol = 3)
?SpatialFeaturePlot
slide_qc_acutepyelo_control


#####################################################################################################################
### Generated a matrix for the acuteonchronicpyelo, here we added the condition and tissue for the meta datatset
#####################################################################################################################
### 
# the tissue coordinate for single cell RNAs in acutepyelo
acuteonchronicpyelocoldata <- GetTissueCoordinates(acuteonchronicpyelo,
                                          cols = c("row", "col", "tissue"),
                                          scale = NULL)
acuteonchronicpyelocoldata

# get the information whether on the tissue
acuteonchronicpyelo$tissue <- acuteonchronicpyelocoldata[rownames(acuteonchronicpyelo@meta.data), "tissue"]

acuteonchronicpyelo$tissue <- ifelse(acuteonchronicpyelocoldata$tissue == 1, 
                            "on_tissue", 
                            "not_on_tissue")
acuteonchronicpyelo$Condition<-replicate(n=length(acuteonchronicpyelo$tissue), 'acuteonchronicpyelo')

acuteonchronicpyelo@meta.data

## add the mitochondrion info

acuteonchronicpyelo[["percent.mt"]] <- PercentageFeatureSet(acuteonchronicpyelo, 
                                                   pattern = "^MT-")

### only select the ones that located in the region
acuteonchronicpyelo_seurat <- subset(acuteonchronicpyelo, subset = tissue == "on_tissue")

slide_qc_acuteonchronicpyelo <- SpatialFeaturePlot(acuteonchronicpyelo_seurat,
                                                  features = c("nCount_Spatial", 
                                                               "nFeature_Spatial", 
                                                               "percent.mt"),
                                                  ncol = 3)
slide_qc_acuteonchronicpyelo 

acuteonchronicpyeloQC<-acuteonchronicpyelo@meta.data %>% select(-orig.ident)

acuteonchronicpyeloQC
colnames(acuteonchronicpyeloQC)
#acutepyeloQC$Conditon<-replicate(n=nrow(acutepyeloQC), 'acutepyelo')

#####################################################################################################################
### Generated a matrix for the chronic sample, here we added the condition and tissue for the meta datatset
#####################################################################################################################
### 
# the tissue coordinate for single cell RNAs in acutepyelo
chronicsamplecoldata <- GetTissueCoordinates(chronicsample,
                                                   cols = c("row", "col", "tissue"),
                                                   scale = NULL)
chronicsamplecoldata

# get the information whether on the tissue
chronicsample$tissue <- chronicsamplecoldata[rownames(chronicsample@meta.data), "tissue"]

chronicsample$tissue <- ifelse(chronicsample$tissue == 1, 
                                     "on_tissue", 
                                     "not_on_tissue")
chronicsample$Condition<-replicate(n=length(chronicsample$tissue), 'chronic')

chronicsample@meta.data

chronicsample[["percent.mt"]] <- PercentageFeatureSet(chronicsample, 
                                                      pattern = "^MT-")
slide_qc_p <- SpatialFeaturePlot(chronicsample,
                                 features = c("nCount_Spatial", 
                                              "nFeature_Spatial", 
                                              "percent.mt"),
                                 ncol = 3)

slide_qc_p
chronicsampleQC<-chronicsample@meta.data %>% select(-orig.ident)

chronicsampleQC


##### Genreate the basic statistics of the spatial transcriptomesx
### combine multiple samples

samcol<-intersect(colnames(ControlQC),colnames(chronicsampleQC))
CombinedQC<-merge(ControlQC,chronicsampleQC, by =samcol, all=TRUE)[samcol]

samcol2<-intersect(colnames(CombinedQC),colnames(acutepyeloQC))
CombinedQC<-merge(CombinedQC,acutepyeloQC, by =samcol2, all=TRUE)[samcol2]

samcol3<-intersect(colnames(CombinedQC),colnames(acuteonchronicpyeloQC))

CombinedQC<-merge(CombinedQC,acuteonchronicpyeloQC, by =samcol2, all=TRUE)[samcol2]

## measure the number of 
nrow(CombinedQC[CombinedQC$Condition=="Control",])
nrow(CombinedQC[CombinedQC$Condition=="chronic",])

nrow(CombinedQC[CombinedQC$Condition=="acuteonchronicpyelo",])
nrow(CombinedQC[CombinedQC$Condition=="acutepyelo",])


### generated the figure for nFeature_Spatial

ggplot(CombinedQC,aes(x = Condition, fill = Condition, y = nCount_Spatial/nFeature_Spatial)) +
  geom_violin() +theme_classic()

  #facet_grid(. ~ nFeature_Spatial, scales = "free")

ggplot(CombinedQC,aes(x = Condition, y = nCount_Spatial,fill = Condition)) +
  geom_violin()+ geom_boxplot(width=0.1)+ theme_classic()

ggplot(CombinedQC,aes(x = Condition,  y =  nFeature_Spatial,fill = Condition)) +
  geom_violin() + geom_boxplot(width=0.1)+ theme_classic()

tail(CombinedQC)
head(ControlQC)
head(chronicsampleQC)
head(acutepyeloQC)
head(acuteonchronicpyeloQC)
### Generate a matrix for the acuteonchronicpyelo condition

#






### got the number of qc_features including number of Spatal number of Feature_Spatial
tissue_qc <- controlsample@meta.data %>%
  select(-orig.ident) %>%
  pivot_longer(-tissue, 
               names_to = "qc_features",
               values_to = "counts") %>%
  ggplot(aes(x = qc_features, fill = tissue, y = counts)) +
  geom_violin() +
  facet_grid(. ~ qc_features, scales = "free")

tissue_qc

# Filter useful spots 

sample_seurat <- subset(sample_seurat, subset = tissue == "on_tissue")

# Continue with analysis

sample_seurat[["orig.ident"]] <- sample_name

# Get mitochondrial genes percentage ------------------------------------------------
sample_seurat[["percent.mt"]] <- PercentageFeatureSet(sample_seurat, 
                                                      pattern = "^MT-")
# QC relationships -------------------------------------------------------------------
qc_p1 <- sample_seurat@meta.data %>%
  ggplot(aes(x = nCount_Spatial, y = nFeature_Spatial)) +
  geom_point() +
  theme_classic() +
  ggtitle(paste0("nspots ", ncol(sample_seurat)))

qc_p2 <- sample_seurat@meta.data %>%
  ggplot(aes(x = nCount_Spatial, y = percent.mt)) +
  geom_point() +
  theme_classic() +
  ggtitle(paste0("nspots ", ncol(sample_seurat)))

qc_panel <- cowplot::plot_grid(qc_p1, qc_p2, ncol = 2, align = "hv")

slide_qc_p <- SpatialFeaturePlot(sample_seurat,
                                 features = c("nCount_Spatial", 
                                              "nFeature_Spatial", 
                                              "percent.mt"),
                                 ncol = 3)

qc_panel_a <- cowplot::plot_grid(qc_panel, slide_qc_p, 
                                 nrow = 2, ncol = 1, 
                                 rel_heights = c(0.5, 0.5))





path<-

  
  
  

qc_feats <- c("sample_names",
              "Estimated Number of Cells",
              "Mean Reads per Cell",
              "Median Genes per Cell",
              "Fraction Reads in Cells",
              "Total Genes Detected",
              "Median UMI Counts per Cell")

sample_names <- list.files(path)

slide_files <- paste0(path,
                      sample_names,
                      "/outs/metrics_summary.csv")

qc_stats <- tibble(sample_names,
                   qc_stats = map(slide_files, read_csv)) %>%
  unnest() %>%
  dplyr::select(all_of(qc_feats))

qc_stats$`Fraction Reads in Cells` <- gsub("%", "", qc_stats$`Fraction Reads in Cells`) %>%
  as.numeric()

qc_stats_plts <- qc_stats %>%
  pivot_longer(-sample_names, names_to = "qc_feature") %>%
  group_by(qc_feature) %>%
  nest() %>%
  mutate(qc_plt = map2(qc_feature, data, function(dat_label, dat) {
    
    ggplot(dat, aes(x = sample_names,
                    y = value)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5),
            axis.text = element_text(size = 12),
            axis.title.y = element_text(size = 13))  +
      xlab("") + ylab(dat_label)
    
    
  }))

all_panels <- cowplot::plot_grid(plotlist = qc_stats_plts$qc_plt, align = "vh", ncol = 1)

pdf(height = 20, width = 17, file = "./processed_snrnaseq/initial_qc/all_qcs.pdf")

plot(all_panels)

dev.off()

qc_stats[, c("sample_names",
             "Estimated Number of Cells",
             "Mean Reads per Cell",
             "Median Genes per Cell")] %>% 
  write.table(row.names = F, col.names = T, quote = F, sep = ",",
              file = "./processed_snrnaseq/initial_qc/all_qcs.csv")
