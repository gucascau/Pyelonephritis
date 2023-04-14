##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@nationwidechildrens.org                         
# Copyright (c) 2023 Kidney and Urology Tract Center
# Nationwide Children's Hospital
# Description: 
#   The script is to perform the integration of PIP-seq single cells 
#     1. integrating four time points of pyelonephritis
#     2. Identifying the cell types that are present in four datasets
#     3. Obtaining the cell type markers taht conserved in four datasets, 
#     4. illustrating the feature figure using Dimpplot, feature plot, dotplot, vlnplot
#  
################################################################

library(dplyr)
library(Seurat)
library(patchwork)

#install.packages('devtools') #assuming it is not already installed

library(devtools)

#install_github('andreacirilloac/updateR')
# install.packages("sctransform")
# install.packages("ggthemes")
# loading the libarary
library(updateR)


library("sctransform")
library("ggthemes")
library(Seurat)
library(Matrix)
library(nichenetr)
library(tidyverse)
library(dplyr)

library(biomaRt)
library("data.table")


library(umap)
library(patchwork)
library(cowplot)

library(ggplot2)


### we used discre color palettes from ggsci

library("ggsci")
library("ggplot2")
library("gridExtra")

##########################################################################
# step 1: preparation for the working environment and the sample details 
##########################################################################

## reading the loacal file and creat the seurat object:

## set pathway:
setwd("/Users/XXW004/Documents/Projects/RuizRosado/PIPseq/Results/Evaluation_V022723/CellEvaluationChangedColor/")

### create a vector of convenient sample names, such as "D0","D3","D7" and "D28"

samples=c("D0","D3","D7","D28")

outdir ="/Users/XXW004/Documents/Projects/RuizRosado/PIPseq/Results/Evaluation_V022723/CellEvaluationChangedColor/"

##########################################################################
### step 2: read in the feature-barcode matrics generated from PIPseeker
##########################################################################

## generate a new list of elements of different types
DataPipseq = list() 

DataPipseq[[1]]<- ReadMtx(mtx = "0DPI/matrix.mtx.gz",features = "0DPI/features.tsv.gz",cells = "0DPI/barcodes.tsv.gz")
DataPipseq[[2]]<- ReadMtx(mtx = "3DPI/matrix.mtx.gz",features = "3DPI/features.tsv.gz",cells = "3DPI/barcodes.tsv.gz")
DataPipseq[[3]]<- ReadMtx(mtx = "7DPI/matrix.mtx.gz",features = "7DPI/features.tsv.gz",cells = "7DPI/barcodes.tsv.gz")
DataPipseq[[4]]<- ReadMtx(mtx = "28DPI/matrix.mtx.gz",features = "28DPI/features.tsv.gz",cells = "28DPI/barcodes.tsv.gz")

dim(DataPipseq[[4]])

##########################################################################
### step 3:Convert each feature-barcode matrix to a Seurat object
##########################################################################
# First create an empty list to hold the Seurat objects
scrna.list = list(); 
#scrna.list[[1]]=CreateSeuratObject(counts = Day0_Expression_Matrix, Project="Day0",min.cells = 3, min.features = 200)
#scrna.list[[1]][["DataSet"]] = samples[1];

#scrna.list[[2]]=CreateSeuratObject(counts = Day0_Expression_Matrix, Project="Day0",min.cells = 3, min.features = 200)
#scrna.list[[2]][["DataSet"]] = samples[2];

#scrna.list[[3]]=CreateSeuratObject(counts = Day0_Expression_Matrix, Project="Day0",min.cells = 3, min.features = 200)
#scrna.list[[3]][["DataSet"]] = samples[3];

#scrna.list[[4]]=CreateSeuratObject(counts = Day0_Expression_Matrix, Project="Day0",min.cells = 3, min.features = 200)
#scrna.list[[4]][["DataSet"]] = samples[4];
# a more efficient way
# here we created the list to build the seurat objects
for (i in 1:length(DataPipseq)){
  scrna.list[[i]]= CreateSeuratObject(counts = DataPipseq[[i]], Project= samples[i],min.cells = 3, min.features = 200);
  
  ### here we measure the mitochondrial proportion for each datasets
  scrna.list[[i]][["percent.mt"]] <- PercentageFeatureSet(scrna.list[[i]], pattern = "^mt-")
  
  ### filtered the low quality cells that have less than 200 RNA or larger than 4000 and percentage mt high than 15%
  scrna.list[[i]]<- subset(scrna.list[[i]],subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
  ## measure the rDNA proportion for each datasets
  scrna.list[[i]][["rDNA"]] <- PercentageFeatureSet(scrna.list[[i]], pattern = "^Rp[sl][[:digit:]]")
  
  scrna.list[[i]][["DataSet"]]=samples[i];
  ### define the default with RNA
  DefaultAssay(scrna.list[[i]]) <- "RNA"
}

##########################################################################
### step 4:Merge the Seurate objects into a single object
##########################################################################
# check the different objects:
scrna.list[[1]]@meta.data

## merge might not be a good way for the intergration, instead we used the integration by find the integration features
#scrna <- merge(x=scrna.list[[1]], y=c(scrna.list[[2]],scrna.list[[3]],scrna.list[[4]]), add.cell.ids = c("A","B","C","D"), project="Pyelonephritis");
## we used Intergration funciton

# normalize and identify variable features for each dataset independently, be attendtion we need to using the normalizatiopn for the datasets
ifnb.list <- lapply(X = scrna.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)
# We then identify anchors using the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
scrna<-  FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
scrna <- IntegrateData(anchorset = scrna)
DefaultAssay(scrna) <- "integrated"

# Run the standard workflow for visualization and clustering
scrna <- ScaleData(scrna, verbose = FALSE)
scrna <- RunPCA(scrna, npcs = 20, verbose = FALSE)
scrna <- RunUMAP(scrna, reduction = "pca", dims = 1:20)
scrna <- FindNeighbors(scrna, reduction = "pca", dims = 1:20)
scrna <- FindClusters(scrna, resolution = 0.6)


## change the order of days 
levels(factor(scrna@meta.data$DataSet))

scrna@meta.data$DataSet <- ordered(factor(scrna@meta.data$DataSet),levels = c("D0", "D3", "D7","D28"))

levels(scrna@meta.data$DataSet)


# Visualization
pdf("UmapDimPlot.pdf", width = 13, height = 6)
p1 <- DimPlot(scrna, reduction = "umap", group.by = "DataSet")
p2 <- DimPlot(scrna, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

dev.off()

pdf("UmapDimPlot_Withlabel.pdf", width = 15, height = 5)
p4 <- DimPlot(scrna, reduction = "umap", split.by ="DataSet",label = TRUE)

p4

dev.off()


pdf(sprintf("%s/UMAP.%d.pdf", outdir, nPC), width = 10, height = 5);
#p1 <- DimPlot(object = scrna, reduction = "tsne", group.by = "DataSet", split.by ="DataSet",pt.size=0.1)
p2 <- DimPlot(object = scrna, reduction = "umap", group.by = "DataSet", split.by ="DataSet",  pt.size=0.1)

?DimPlot
p2
dev.off()

##########################################################################
### step 4:Annotate the cluster and measure the different proportion
##########################################################################

## 4.1 check the structure of seurat integrated datassets

str(scrna@meta.data)
# Because it took such a long time, I will save all the intergerated object as an RDS file, next time I will not intergrate and read the RDS file diretly
saveRDS(scrna, file =sprintf("%s/MergedSeuratObject_Pyelonephritis.RDS",outdir))
scrna@meta.data
# For new datasets, we only read the RDS file 
#################
#################
## we can generate mergeseurate object from here
rm(allen_reference)
rm(Cluster0subset)
scrna<-readRDS("../CellEvaluation/MergedSeuratObject_Pyelonephritis.RDS")
# For new datasets, we generated more clusters
#scrna<-FindClusters(scrna, resolution = 0.6)

nrow(scrna@meta.data)

D0<-scrna@meta.data %>% filter(`DataSet` == "D0") %>% pull(DataSet)
length(D0)

D3<-scrna@meta.data %>% filter(`DataSet` == "D3") %>% pull(DataSet)
length(D3)


D7<-scrna@meta.data %>% filter(`DataSet` == "D7") %>% pull(DataSet)
length(D7)

D28<-scrna@meta.data %>% filter(`DataSet` == "D28") %>% pull(DataSet)
length(D28)


## 4.2 measure the number of cluster 3 in D0 and D7
D0Cluster3<- scrna@meta.data %>% filter(`DataSet` == "D0" & `seurat_clusters` == 3) %>% pull(DataSet)
length(D0Cluster3)
D3Cluster3<- scrna@meta.data %>% filter(`DataSet` == "D3" & `seurat_clusters` == 3) %>% pull(DataSet)
D7Cluster3<- scrna@meta.data %>% filter(`DataSet` == "D7" & `seurat_clusters` == 3) %>% pull(DataSet)
D28Cluster3<- scrna@meta.data %>% filter(`DataSet` == "D28" & `seurat_clusters` == 3) %>% pull(DataSet)

# different number of cells in each cluster
table(scrna@meta.data$seurat_clusters)
# total number of cells in all clusters
length(scrna@meta.data$seurat_clusters)
# check the name of different clusters

names(table(scrna@meta.data$seurat_clusters))
## 4.3 add the proportion of each cluster and show on the annotated clusters 
df<- data.frame(clu=names(table(scrna@meta.data$seurat_clusters)), percentage= sprintf("%1.2f", 100*table(scrna@meta.data$seurat_clusters)/length(scrna@meta.data$seurat_clusters)))
## then add the proportion into the meta data
scrna@meta.data$ClusterPercent <- df[match(scrna@meta.data$seurat_clusters,df$clu),2]


## we also added the proportion of each cluster based on the different cluster at different time points

# get the number of each cluster at different time points
table(scrna@meta.data$TimepointCluster)
# get the number of all cluster at different time points
table(scrna@meta.data$DataSet)
# Firstly, create a new column to combine Dataset and cluster
scrna@meta.data$TimepointCluster<-paste0(scrna@meta.data$DataSet,":",scrna@meta.data$seurat_clusters)
tail(scrna@meta.data)
## add the proportion of each cluster at different time points
TimeDf <- data.frame(Timecluster=names(table(scrna@meta.data$TimepointCluster)), percentage= sprintf("%1.2f", 100*table(scrna@meta.data$TimepointCluster)/length(scrna@meta.data %>% filter (`DataSet` == strsplit(scrna@meta.data$TimepointCluster,":") [[1]][1]) %>% pull(DataSet))))
TimeDf

## we then added the TimeClusterPercentage
scrna@meta.data$TimeClusterPercentage <- TimeDf[match(scrna@meta.data$TimepointCluster, TimeDf$Timecluster),2]
scrna@meta.data

## add the annotation with the proportion

#nrow(scrna@meta.data[scrna@meta.data$seurat_clusters==3& scrna@meta.data$DataSet == "D0",])

#nrow(scrna@meta.data[scrna@meta.data$seurat_clusters==3& scrna@meta.data$DataSet == "D7",])
### before adding the annotation, we added the number of cluster and generate the figure

######################################################################################################################################################################
### For Figure 1A: The Uniform Manifold Approximation and Projection (UMAP) of sequenced cells from different time points of bacterial infections on kidney
######################################################################################################################################################################
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=12)

color_list

### using weanderson
### only for three - five color combination
#install.packages("wesanderson")
library(wesanderson)
## names
names(wes_palettes)



### use Rcolor Breweer paletters 
library(RColorBrewer) 
# scale_color_brewer(palette = "Dark2")

###### the best way to get more than 20 color by interpolate existing ones with constructor fucntion colorRampPallet

## get the number of colors
colorCounts =length(unique(scrna@meta.data$seurat_clusters))
## get the color palette
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
clustercolor<-getPalette(colorCounts)

nPC = 20
pdf(sprintf("%s/UMAP_V0310_SCIcolor.%d.pdf", outdir, nPC), width =6, height = 5);
#p1 <- DimPlot(object = scrna, reduction = "tsne", group.by = "DataSet", split.by ="DataSet",pt.size=0.1)
p2 <- DimPlot(object = scrna, reduction = "umap",  pt.size=0.1, label = T,repel = TRUE ,cols=clustercolor) 
  #scale_fill_manual(values = getPalette(colorCounts))
#p3 <- DimPlot(object = scrna, reduction = "umap", group.by = "DataSet", split.by ="DataSet",  pt.size=0.1)

#DimPlot(scrna, label = FALSE,reduction = "umap",  pt.size=0.1)
p2
dev.off()


######################################################################################################################################################################
### For Figure 1: The Uniform Manifold Approximation and Projection (UMAP) of sequenced cells within the time bacterial infections on kidney
######################################################################################################################################################################
pdf("UmapDimPlot_split_V0310_SCIcolor.pdf", width = 15, height = 4)
DimPlot(scrna, label = TRUE, reduction = "umap", pt.size=0.1, split.by = "DataSet",repel = TRUE,cols =clustercolor ) 
dev.off()

### the seperated c

cluster_colors <- c(
 "cornflowerblue", "darkblue", "cyan",
 "springgreen", "olivedrab", "chartreuse", "darkgreen",
 "gray", "black", "saddlebrown",
 "slateblue", "purple", "pink", "lightpink3", "hotpink", "darkred",
 "lightsalmon"
, "orange", "tan", "firebrick1", "indianred"
)

## 4.4 based on the marker genes, we rename the identity to cluster name

scrna<-RenameIdents(scrna, `0`= "Proximal tubule I", `1`= "Proximal tubule II", `2`= "T cells", 
                    `3`= "Myeloblast", `4`= "Proximal tubule III", `5`= "Proximal tubule IV", 
                    `6`= "Thick ascending limb of LOH", `7`= "Principal cells", `8`= "Intercalated cells", `9`= "NK cells",
                    `10`= "Thin limb of LOH", `11`= "Endothelium", `12`= "Plasmacytoid dendritic cells", 
                    `13`= "Proximal tubule V", `14`= "B cells", `15`= "Distal tubule", `16`= "Urothelium",`17`= "Fibroblast")

## get the current id of the cells
Idents(object = scrna)
#DimPlot(object = scrna, reduction = "umap",  pt.size=0.1, label = T,repel = TRUE )

# reorder the levels
levels(x=scrna) <- c("Proximal tubule I", "Proximal tubule II", "Proximal tubule III" ,  "Proximal tubule IV", "Proximal tubule V", 
                     "Plasmacytoid dendritic cells", "Distal tubule",
                    "Intercalated cells",  "Principal cells", "Thin limb of LOH", "Thick ascending limb of LOH", 
                     "Urothelium" ,"Endothelium","Fibroblast",  "Myeloblast" , "T cells" ,"NK cells" , "B cells" )

scrna@meta.data$Annot <- paste0(scrna@meta.data$seurat_clusters,":",Idents(scrna)) 


## paste the identity cell type, proportion of cells
scrna$newAnnot<-paste0(scrna@meta.data$seurat_clusters,":",Idents(scrna)," (",scrna@meta.data$ClusterPercent,")")

scrna@meta.data
## change the order the annotation
levels(as.factor(scrna@meta.data$newAnnot))

# reorder the levels
scrna@meta.data$newAnnot<- ordered(factor(scrna@meta.data$newAnnot), levels=c("0:Proximal tubule I (30.60)","1:Proximal tubule II (29.48)",
                                                                              "2:T cells (8.11)","3:Myeloblast (7.86)","4:Proximal tubule III (4.55)",
                                                                              "5:Proximal tubule IV (4.54)","6:Thick ascending limb of LOH (4.47)","7:Principal cells (2.81)",
                                                                              "8:Intercalated cells (1.65)","9:NK cells (1.35)",
                                                                              "10:Thin limb of LOH (1.11)","11:Endothelium (0.75)",
                                                                              "12:Plasmacytoid dendritic cells (0.72)",
                                                                              "13:Proximal tubule V (0.61)","14:B cells (0.59)","15:Distal tubule (0.42)","16:Urothelium (0.31)","17:Fibroblast (0.08)" ))
levels(as.factor(scrna@meta.data$newAnnot))
pdf("UmapDimPlot_Withlabels_V0310_scicolor.pdf", width = 9, height = 6)
DimPlot(scrna, label = FALSE,reduction = "umap",  pt.size=0.1, group.by = "newAnnot",cols =clustercolor)

dev.off()
pdf("UmapDimPlot_V0310_SCIcolor.pdf", width = 8, height = 6)
DimPlot(scrna, label = FALSE, reduction = "umap", pt.size=0.1, group.by = "ident",repel = TRUE,cols =clustercolor)
dev.off()



######################################################################################################################################################################
### For Figure 1B: Frequency of different cluster proportion within the time points
######################################################################################################################################################################
scrna@meta.data
scrna@meta.data$TimeClusterPercentage


Idents(scrna)
as.numeric(scrna@meta.data$TimeClusterPercentage)

## using the sci color
pdf("Frequency_split_V0310_SCIcolor.pdf", width = 15, height = 4)
ggplot(scrna@meta.data, aes(fill=DataSet, y=as.numeric(scrna@meta.data$TimeClusterPercentage), x=Idents(scrna))) + 
  geom_bar(position="dodge", stat="identity")+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits=c(0,40)) 
dev.off()

### we only consider the PT cells
## only 

PTCells<-scrna@meta.data[scrna@meta.data$seurat_clusters %in% c(0,1,4,5,13),]
## change the levels of clusters
PTCells$Annot <- ordered(factor(PTCells$Annot), levels= c("0:Proximal tubule I","1:Proximal tubule II","4:Proximal tubule III","5:Proximal tubule IV","13:Proximal tubule V"))

#PTCells$TimeClusterPercentage
pdf("Frequency_split_V0310_PT_SCIcolor.pdf", width = 6, height = 4)
ggplot(PTCells, aes(fill=DataSet, y=as.numeric(PTCells$TimeClusterPercentage), x=PTCells$Annot)) + 
  geom_bar(position="dodge", stat="identity")+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ scale_fill_manual(values= c("#F1948A", "#F4D03F","#BB8FCE","#5DADE2","#76D7C4"))
dev.off()



scrna@meta.data
# then plot the numeber of counts and feature in each clusters
pdf("UmapDimPlot_Vcounts_Feature_V0310.pdf", width = 8, height = 12)

Ncount_RNA <- VlnPlot(scrna, features = "nCount_RNA",cols = cluster_colors, pt.size = 0, group.by = "newAnnot") +
  theme(axis.title.x = element_blank())+
  NoLegend()
NFeature_RNA <- VlnPlot(scrna, features = "nFeature_RNA",cols = cluster_colors, pt.size = 0, group.by = "newAnnot") +
  xlab("cluster_id") +
  NoLegend()
plot_grid(Ncount_RNA,NFeature_RNA, ncol = 1)

dev.off()
scrna@meta.data
### only draw the figures for D0 and D7
nrow(scrna@meta.data)

scrna@meta.data


######################################################################################################################################################################
### For Figure 1B: Dotplot for the marker genes
######################################################################################################################################################################
 
DefaultAssay(scrna) <- "RNA"
Clustermakrers<- c("Lrp2","Slc34a1","Acsm1","Ireb2","Alas2","C78334","Slc12a3","Pvalb","Wnk1","Atp6v1g3","Aqp6","Slc4a1",
                   "Aqp2","Hsd11b2","Scnn1g","Fst","Bst1","Slc14a2","Umod","Slc12a1","Cldn10","Krt5","Krt14","Upk1b",
                   "Pecam1","Cdh5","Nrp1","Col3a1","Lum","Fbn1","C1qa","C1qb","Aif1","Lst1","Lyz2","Cebpb","Cxcr6","Cd247","Nkg7","Cd7","Ccl5","Cd79a","Cd79b","Ms4a1")

pdf("DotPlot_Markers_V0310_Colorscaled.pdf", width =15, height = 6)
DotPlot (object = scrna, features = Clustermakrers,cols = c("lightgray","steelblue")) + RotatedAxis()  
#+ scale_fill_gradient(low = "white", high = "steelblue")
dev.off()


library(viridis)
pdf("DotPlot_Markers_split_V0310_Coloscale.pdf", width =15, height = 6)
scrna@meta.data
DotPlot (scrna, features = Clustermakrers,group.by  = "DataSet", cols = c("lightgray","darkblue"))+RotatedAxis() 
dev.off()

# library(ComplexHeatmap)
# Heatmap(scrna,
#         heatmap_legend_param=list(title="expression"),
#         column_title = "clustered dotplot", 
#         #col=col_fun,
#         rect_gp = gpar(type = "none"),
#         #cell_fun = cell_fun,
#         row_names_gp = gpar(fontsize = 5),
#         row_km = 4,
#         border = "black",
#         top_annotation = column_ha)
# 
# #ordered(levels(identity(scrna)), levels=c("Proximal tubule I", "Proximal tubule II", 
#  #                                                                             "Proximal tubule III"  ,  "Proximal tubule IV", "Proximal tubule VI", 
#                                                                               "Principal cells", "Intercalated cells",  "Thick ascending limb of LOH", "Thin limb of LOH",
#                                                                               "Urothelium" , "Myeloblast" , 
#                                                                               "T cells" ,"NK cells" , "B cells" , "Fibroblast", "Endothelium" ))
# 

## reorder the identity


###### infer cell types

###############################################################################################################################
### For several of interesting genes, such as Umod urotmodulin, sopdium blance and blood pressure regulation.
### The excretion of uromodulin in urine may provide definse against urinary tract infections cause by uropathogenic bacteria
###############################################################################################################################

## generate the Umod feature across conditions
pdf("Umod_VlnPlot.pdf", height = 4, width = 15)
VlnPlot(scrna, features = "Umod",pt.size = 0, split.by = "DataSet") 
dev.off()

pdf("Umod_VlnPlot_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Umod",pt.size = 0, split.by = "DataSet",cols = c("gray","Darkblue") ) 
#+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()



## generate the Umod feature across conditions
pdf("Lyz2_VlnPlot.pdf", height = 4, width = 15)
VlnPlot(scrna, features = "Lyz2",pt.size = 0, split.by = "DataSet") 
dev.off()

pdf("Lyz2_VlnPlot_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Lyz2",pt.size = 0, split.by = "DataSet") 
dev.off()



### check all the neutrophil makres
pdf("S100a8_neutrophil_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "S100a8",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("S100a8_neutrophil_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "S100a8", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()

pdf("S100a9_neutrophil_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "S100a9",pt.size = 0, split.by = "DataSet") 
dev.off()
scrna@meta.data
pdf("S100a9_neutrophil_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "S100a9", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()

pdf("Plac8_neutrophil_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Plac8",pt.size = 0, split.by = "DataSet") 
dev.off()

pdf("Plac8_neutrophil_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "Plac8", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()

pdf("Ly6g_neutrophil_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Ly6g",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("Ly6g_neutrophil_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "Ly6g", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()

pdf("Itgam_neutrophil_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Itgam",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("Itgam_neutrophil_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "Itgam", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()


pdf("Cd244a_neutrophil_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Cd244a",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("Cd244a_neutrophil_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "Cd244a", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()


pdf("Krt5_Urothelium_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Krt5",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("Krt5_Urothelium_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "Krt5", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()

pdf("Krt14_Urothelium_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Krt14",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("Krt14_Urothelium_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "Krt14", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()


pdf("Krt14_Urothelium_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Krt14",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("Krt14_Urothelium_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "Krt14", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()


pdf("Upk1b_Urothelium_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Upk1b",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("Upk1b_Urothelium_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = "Upk1b", pt.size=1, ncol = 2, group.by="seurat_clusters", split.by = "DataSet");
dev.off()


## using the published markers

## set diffult array to RNA
DefaultAssay(scrna) <- "RNA"
#pdf("Marker_Genes.pdf", height = 16, width = 8)

## for endothelium and Fibroblasts

pdf("EndotheliumFibroblast_Marker_Genes.pdf", height = 6, width = 9)
EndotheliumFibroblast<- FeaturePlot(object = scrna, features = c("Pecam1", "Cd34","Vwf", "Col3a1", "Lum","Fbn1"),   min.cutoff = 'q0', max.cutoff = 'q80', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = EndotheliumFibroblast, ncol = 1)
dev.off()

## for the Proximal Tubule cells and distal tubule cells
pdf("ProximalTubuleDistalTubule_Genes.pdf", height = 6, width = 9)
ProximalTubuleDistalTubule<- FeaturePlot(object = scrna, features = c("Slc34a1", "Lrp2","Acsm1","Slc12a3", "Pvalb", "Wnk1"),   min.cutoff = 'q0', max.cutoff = 'q80', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = ProximalTubuleDistalTubule, ncol = 1)
dev.off()

## for the Collecting duct cells IC and PC
pdf("CollectingDuct_IC_PC_Genes.pdf", height = 6, width = 9)
CollectingDuct<- FeaturePlot(object = scrna, features = c("Slc4a1","Aqp6","Atp6v1g3","Aqp2","Scnn1g","Hsd11b2"),   min.cutoff = 'q0', max.cutoff = 'q80', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = CollectingDuct, ncol = 1)
dev.off()

## for the Urothelium
pdf("LoopHenle_Genes.pdf", height = 6, width = 9)
LoopHenle<- FeaturePlot(object = scrna, features = c("Fst","Bst1","Slc14a2","Slc12a1","Umod","Cldn10"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = LoopHenle, ncol = 1)
dev.off()

### Ascending thin lib of LOH (cannot find)

pdf("AscendingthinlibofLOH_Genes.pdf", height = 6, width = 9)
LoopHenle<- FeaturePlot(object = scrna, features = c("Epha7","Mx2","Clcnka","Enox1","Thsd4","Nos1"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = LoopHenle, ncol = 1)
dev.off()


### Pericytes (No pericytes were found)

pdf("Pericytes_Genes.pdf", height = 6, width = 9)
LoopHenle<- FeaturePlot(object = scrna, features = c("Vim","Tagln","Myh11","Pdgfrb","Ren1","Gata3"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = LoopHenle, ncol = 1)
dev.off()

### Mesangia (No Mesangia were found)

pdf("Mesangial_Genes.pdf", height = 6, width = 9)
Mesangial<- FeaturePlot(object = scrna, features = c("Serpine2","Fhl2","Des","Prkca","Art3","Nt5e"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = Mesangial, ncol = 1)
dev.off()

### ProlifeartionTubule (No ProlifeartionTubule were found)

pdf("ProliferationTuble_Genes.pdf", height = 6, width = 9)
ProliferationTuble<- FeaturePlot(object = scrna, features = c("Agt","Rnf24","Slc22a7","Slc22a13","Art3","Nt5e"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = ProliferationTuble, ncol = 1)
dev.off()


### Podocyte (No ProlifeartionTubule were found)
pdf("Podocyte_Genes.pdf", height = 6, width = 9)
ProliferationTuble<- FeaturePlot(object = scrna, features = c("Nphs1","Nphs2","Mafb","Slc22a13","Art3","Nt5e"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = ProliferationTuble, ncol = 1)
dev.off()



## for the immune cells

pdf("ImmuneTotal_Genes.pdf", height = 6, width = 9)
Immunecells<- FeaturePlot(object = scrna, features = c("Fst","Bst1","Slc14a2","Slc12a1","Umod","Cldn10"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = Immunecells, ncol = 1)
dev.off()

pdf("Mcrophage_Neutrophil_Genes.pdf", height = 6, width = 9)
MicrophageNeutrophil<- FeaturePlot(object = scrna, features = c("C1qa","C1qb","Aif1","Lst1","Lyz2","Cebpb"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = MicrophageNeutrophil, ncol = 1)
dev.off()

pdf("Dendritic_neutrophil_Genes.pdf", height = 6, width = 9)
Dendritic_neutrophil<- FeaturePlot(object = scrna, features = c("Itgax","Cxcr4","Ly6d","Lst1","Lyz2","Cebpb"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = Dendritic_neutrophil, ncol = 1)
dev.off()

pdf("T_B_Genes.pdf", height = 6, width = 9)
T_B<- FeaturePlot(object = scrna, features = c("Cxcr6","Cd247","Nkg7","Cd79a","Cd79b","Ms4a1"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = T_B, ncol = 1)
dev.off()

FeaturePlot(object = scrna, features = c("Cd19","Ptprc","Slc14a2","Slc12a1","Umod","Cldn10"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());

pdf("CD45T_Cd8effectorT_Genes.pdf", height = 6, width = 9)
CD45T_Cd8effectorT<- FeaturePlot(object = scrna, features = c("Ccr7","Ms4a4b","Klf2","Ccl5","Cd8b1","Tnfrsf4"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = CD45T_Cd8effectorT, ncol = 1)
dev.off()

pdf("NkCell_Genes.pdf", height = 6, width = 9)
NK<- FeaturePlot(object = scrna, features = c("Nkg7","Cd7","Ccl5","Tmsb10","Ly6c2","Cxcr6"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = NK, ncol = 1)
dev.off()

pdf("Tregulatory_Genes.pdf", height = 6, width = 9)
Tregulatory<- FeaturePlot(object = scrna, features = c("S100a4","Rgs1","Cd3g","Ltb","Capg","Izumo1r"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = Tregulatory, ncol = 1)
dev.off()

pdf("Vasculature_Genes.pdf", height = 6, width = 6)
FeaturePlot(object = scrna, features = c("Plvap","Kdr","Pecam1"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
dev.off()

#### Basophgil (No Basophgil were found)
pdf("Basophgil_Genes.pdf", height = 6, width = 9)
Basophgil<- FeaturePlot(object = scrna, features = c("Ifitm1","Hdc","Mcmpt8","Fcer1a","Csrp3","Ms4a2"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = Basophgil, ncol = 1)
dev.off()

### check netrophil 


FeaturePlot(object = scrna, features = c("S100a8", "S100a9","Lyz2"),  ncol=3, reduction = "umap", ) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());

#### pacemaker cells
FeaturePlot(object = scrna, features = c("Shh", "Cd34","Vwf"),  ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());

### here we have already put the identity into the G1 G2M and S 

Idents(scrna)<- scrna@meta.data$seurat_clusters
table(Idents(scrna))
### identify the marker genes of cluster 12 and 15 cannot find

#table(scrna@meta.data$old.ident)
cluster12.markers <- FindMarkers(scrna, ident.1 = "12", min.pct = 0.25)
head(cluster12.markers, n = 18)

cluster16.markers <- FindMarkers(scrna, ident.1 = "16", min.pct = 0.25)
head(cluster16.markers, n = 18)
#cluster12.markers

## Plasmacytoid dendritic cells
## -- Erythroid-like and erythroid precursor cells
# https://www.science.org/doi/10.1126/science.aar2131
pdf("cluster12.markers.pdf", height = 6, width = 9) 
dev.off()
cluster12.markersFeature<-FeaturePlot(object = scrna, features = c("Ireb2","Alas2","Kit"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = cluster12.markersFeature, ncol = 1)
FeaturePlot(object = scrna, features = c("Alas2","Ireb2","C78334"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());

pdf("cluster16.markers.pdf", height = 6, width = 9) 
cluster16.markersFeature<-FeaturePlot(object = scrna, features = c("Duoxa2","S100a14","Sprr2a1","Sfn"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = cluster16.markersFeature, ncol = 1)
dev.off()

FeaturePlot(object = scrna, features = c("Alas2","Ireb2","C78334"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());


FeaturePlot(object = scrna, features = c("Nphs2","Vim","Myh11"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());

table((scrna@meta.data)$seurat_clusters)
table((scrna@meta.data)$new)
cluster15.markers <- FindMarkers(scrna, ident.1 = "15", min.pct = 0.25)

head(cluster15.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(scrna, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 10)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones

pbmc.markers <- FindAllMarkers(scrna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers[pbmc.markers$cluster=="1",]

write.csv(pbmc.markers, file="AllBiomarkers.csv")

pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

pdf(file = "Cluter1FeaturePlot.pdf",width =10, height = 5)
FeaturePlot(object = scrna, features = c("Cyp4b1", "Slc34a1"),  ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
dev.off()
### save the biomarkers in a files

pdf("cluster15.markers.pdf", height = 6, width = 9)
cluster12.markers<-FeaturePlot(object = scrna, features = c("Ltf","Sprr2a2","Sprr2a1"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = cluster12.markers, ncol = 1)
dev.off()


fp <- FeaturePlot(object = scrna, features = c("Pecam1", "Cd34","Vwf","Atp1b1","Slc34a1","Lrp2","Cldn8","Aqp2","Cldn8","Slc4a1", "Scnn1g","Slc12a1", "Aqp1","Col3a1", "Ptprc", "Cd74", "Lyz2"),  ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = fp, ncol = 1)

dev.off()
immuneplot<-VlnPlot(scrna, features = c("Ptprc", "Cd74", "Lyz2"), group.by = "seurat_clusters",
                    pt.size = 0, combine = FALSE)
wrap_plots(plots = immuneplot, ncol = 1)

## save the scrna into 


saveRDS(scrna, file =sprintf("%s/Final_CombinedInformation_Pyelonephritis_V0330.RDS",outdir))



## ddd the biomarkers
# idnetify conserved cell type markers


# Step 6. Calculate a cell cycle score for each cell

cell.cycle.tirosh <- read.csv("http://genomedata.org/rnaseq-tutorial/scrna/CellCycleTiroshSymbol2ID.csv", header=TRUE); # read in the list of genes

### convert the gene id from human into mouse
head(cell.cycle.tirosh)
dim(cell.cycle.tirosh)
cell.cycle.tirosh$Gene.Symbol<- cell.cycle.tirosh$Gene.Symbol %>% convert_human_to_mouse_symbols()
# ignore the NxA in the cell cycles
cell.cycle.tirosh = cell.cycle.tirosh %>% na.omit()

cell.cycle.tirosh
dim(cell.cycle.tirosh)

### create the S-phage and G2/M phage genes
s.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G1/S")]; # create a vector of S-phase genes
length(s.genes)
g2m.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G2/M")]; # create a vector of G2/M-phase genes
length(g2m.genes)

scrna <- CellCycleScoring(object=scrna, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
head(scrna[[]])
scrna@meta.data

## visualize the grouping by cell cycle phages
pdf("CellCylcle_phages_cluster.pdf", height = 5, width = 15)
scrna@meta.data
DimPlot(scrna, reduction = "umap", split.by = "Phase", group.by="seurat_clusters")
dev.off()

pdf("Upk1b_Urothelium_Feature.pdf", height = 4, width = 15)
FeaturePlot(object = scrna, features = "Upk1b",pt.size = 0, split.by = "DataSet") 
dev.off()
pdf("CellCylcle_VlnFeature.pdf", height = 4, width = 15)
VlnPlot(object = scrna, features = c("Nkg7", "Cd7", "Cd5", "Cxcr6"), pt.size=1, ncol = 2, group.by="Phase", split.by = "DataSet");
dev.off()

# Visualize the distribution of cell cycle markers across
pdf("CellCylcle_phages_RidgePlot.pdf", height = 12, width = 8)
RidgePlot(scrna, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
dev.off()

## generate the phages in different cluster
scrna[[]]
### Generate the phage frquency at different cluster

## add the proportion of each cluster at different time points for cell phage
## 4.3 add the proportion of each cluster and show on the annotated clusters 
PhageProportion<- data.frame(clustername=names(table(scrna@meta.data$old.ident)), Gpercentage= sprintf("%1.2f", 100*table(scrna@meta.data$old.ident)/length(scrna@meta.data$seurat_clusters)))
## then add the proportion into the meta data

## add the proportion of each cluster at different time points
#TimeDf <- data.frame(Timecluster=names(table(scrna@meta.data$old.ident)), Gpercentage= sprintf("%1.2f", 100*table(scrna@meta.data->filter(``)/length(scrna@meta.data %>% filter (`DataSet` == scrna@meta.data$old.ident[[1]][1]) %>% pull(DataSet)))

                     
              #       scrna@meta.data$old.ident[[]]

#scrna@meta.data %>% filter(`Phase`== "G1")-> pull(DataSet)

## we also added the proportion of each cluster based on the different cluster at different time points

# get the number of each cluster at different time points
table(scrna@meta.data$TimepointCluster)
# get the number of all cluster at different time points
table(scrna@meta.data$DataSet)
# Firstly, create a new column to combine Dataset and cluster
scrna@meta.data$TimepointCluster<-paste0(scrna@meta.data$DataSet,":",scrna@meta.data$seurat_clusters)
tail(scrna@meta.data)
## add the proportion of each cluster at different time points
TimeDf <- data.frame(Timecluster=names(table(scrna@meta.data$TimepointCluster)), percentage= sprintf("%1.2f", 100*table(scrna@meta.data$TimepointCluster)/length(scrna@meta.data %>% filter (`DataSet` == strsplit(scrna@meta.data$TimepointCluster,":") [[1]][1]) %>% pull(DataSet))))
TimeDf



### then I double check the clusters defined by any of the technical difference, Mitochondrial genes show certain dependency on clusters???
FeaturePlot(srna)


## save the object 

#saveRDS(scrna, file =sprintf("%s/Final_CombinedInformation_Pyelonephritis_V0321.RDS",outdir))








































### fully understand the seurate object meta data which is stored in scrna@meta.data
## summary stattiscs, sample name, cluster membership fo each cell, cell cylce phase for each cells, batch or sample for each cell, other custom annotaion for each cell

scrna[[]]
scrna@meta.data;

## examine structure and contens of meta data
tail(scrna@meta.data)

### access gene/ number of UMIs for each cell
head(scrna@meta.data$nFeature_RNA)
head(scrna@meta.data$nCount_RNA)
levels(x=scrna)

#### step 5: Quality control plots

# calcualte the mitochondrial transcript percentage for each cell:

# find the 
#mit.genes<- grep (pattern = "^mt-", x= rownames(x= scrna), value = TRUE)
#mit.genes

#rownames(x= scrna)

#grep (pattern = "^Clec", x= rownames(x= scrna), value = TRUE)

#percent.mito<- Matrix::colSums(x=GetAssayData(object = scrna, slot = 'counts')[mit.genes,])/Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts'))

#scrna[[]]
#GetAssayData(object = scrna, slot = 'counts')
# Calculate the ribosomal transcript percentage for each cell:
#ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = scrna), value = TRUE);
#ribo.genes
#percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna, slot = 'counts'));
#scrna[['percent.ribo']] <- percent.ribo;

#scrna$percent.ribo

#QC and selecting cells for further analysis
#scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^mt-")

#scrna@meta.data$DataSet<- ordered(factor(scrna@meta.data$DataSet), levels=c(1,3,4,2))
scrna@meta.data$DataSet
## Plot as violin plots
scrna@meta.data

pdf("MitoRibsomeVlnPlot.pdf", width = 13, height = 6)
vln <- VlnPlot(object = scrna, features = c("percent.mt", "rDNA"), pt.size=0, ncol = 2, group.by="DataSet");
print(vln);
dev.off();

pdf("GeneCountsPlot.pdf", width = 10, height = 10)
vln <- VlnPlot(object = scrna, features = "nCount_RNA", pt.size=0, group.by="DataSet", y.max=25000)
print(vln)
dev.off();

pdf("FeaturePlot.pdf", width = 10, height = 10)
vln <- VlnPlot(object = scrna, features = "nFeature_RNA", pt.size=0, group.by="DataSet")
print(vln)
dev.off()

pdf(sprintf("%s/Scatter1.pdf", outdir), width = 8, height = 6);
scatter <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size=0.1)
print(scatter);
dev.off();

pdf(sprintf("%s/Scatter2.pdf", outdir), width = 8, height = 6);
scatter <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1)
print(scatter);
dev.off();

pdf(sprintf("%s/Scatter3.pdf", outdir), width = 8, height = 6);
scatter <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1)
print(scatter);
dev.off();


#####
#Run the standard workflow for visualization and clustering


### Step 9. Normalize the data, detect variable genes, and scale the data

# normalize the data

scrna <- NormalizeData(object = scrna, normalization.method = "LogNormalize", scale.factor = 1e6);

## now identify and plot the most varible genes, which would be used for downstream analyses. This is the critical step to reduce the contirtbution noice,. Considering adjust the cutoffs if you think that important genes are being excluded.

scrna <- FindVariableFeatures(object = scrna, selection.method = 'vst', mean.cutoff = c(0.1,8), dispersion.cutoff = c(1, Inf))
print(paste("Number of Variable Features: ",length(x = VariableFeatures(object = scrna))));
pdf(sprintf("%s/VG.pdf", outdir), useDingbats=FALSE)
vg <- VariableFeaturePlot(scrna)
print(vg);
dev.off()

# Scale and center the data:

scrna <- ScaleData(object = scrna, features = rownames(x = scrna), verbose=FALSE);

# save the normalized,scaled seurat object
saveRDS(scrna, file = sprintf("%s/VST.rds", outdir))

# Step 10. Reduce the dimensionality of the data using Principal Component Analysis

scrna <- RunPCA(object = scrna, npcs = 100, verbose = FALSE);

# There are several easy ways to investigate these questions. First, visualize the PCA “loadings.” Each “component” identified by PCA is a linear combination, or weighted sum, of the genes in the data set. Here, the “loadings” represent the weights of the genes in any given component. These plots tell you which genes contribute most to each component:

pdf(sprintf("%s/VizDimLoadings.pdf", outdir), width = 8, height = 30);
vdl <- VizDimLoadings(object = scrna, dims = 1:3)
print(vdl);
dev.off();

# Second, use the DimHeatmap function to generate heatmaps that summarize the expression of the most highly weighted genes in each principal component. As noted in the Seurat documentation, “both cells and genes are ordered according to their PCA scores. Setting cells.use to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated gene sets.
pdf(sprintf("%s/PCA.heatmap.multi.pdf", outdir), width = 8.5, height = 24);
hm.multi <- DimHeatmap(object = scrna, dims = 1:12, cells = 500, balanced = TRUE);
print(hm.multi);
dev.off();


### Step 11. Generate 2-dimensional layouts of the data using two related algorithms, t-SNE and UMAP

scrna <- RunUMAP(object = scrna, reduction = "pca", dims = 1:nPC);
scrna <- RunTSNE(object = scrna, reduction = "pca", dims = 1:nPC);

pdf(sprintf("%s/UMAP.%d.pdf", outdir, nPC), width = 10, height = 10);
p1 <- DimPlot(object = scrna, reduction = "tsne", group.by = "DataSet", split.by ="DataSet",pt.size=0.1)
p2 <- DimPlot(object = scrna, reduction = "umap", group.by = "DataSet", split.by ="DataSet",  pt.size=0.1)

p1+p2
dev.off()

##
save(scrna,file = "scrna.Rdata")



# Step 6. Calculate a cell cycle score for each cell

cell.cycle.tirosh <- read.csv("http://genomedata.org/rnaseq-tutorial/scrna/CellCycleTiroshSymbol2ID.csv", header=TRUE); # read in the list of genes

### convert the gene id from human into mouse
head(cell.cycle.tirosh)
dim(cell.cycle.tirosh)
cell.cycle.tirosh$Gene.Symbol<- cell.cycle.tirosh$Gene.Symbol %>% convert_human_to_mouse_symbols()
# ignore the NxA in the cell cycles
cell.cycle.tirosh = cell.cycle.tirosh %>% na.omit()

cell.cycle.tirosh
dim(cell.cycle.tirosh)


### create the S-phage and G2/M phage genes
s.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G1/S")]; # create a vector of S-phase genes
length(s.genes)
g2m.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G2/M")]; # create a vector of G2/M-phase genes
length(g2m.genes)

scrna <- CellCycleScoring(object=scrna, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
head(scrna[[]])
scrna@meta.data

## visualize the grouping by cell cycle phages

DimPlot(scrna, reduction = "umap", split.by = "DataSet", group.by="Phase")

# Visualize the distribution of cell cycle markers across
RidgePlot(scrna, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)

### Calcuate a cell cycle socore for each cell

table(scrna[[]]$Phase)

### then I double check the clusters defined by any of the technical difference, Mitochondrial genes show certain dependency on clusters???
FeaturePlot(srna)



# step 7: Filter the cells to remove debris, dead cells and probable oublets
min <- min(scrna@meta.data$nFeature_RNA);
m <- median(scrna@meta.data$nFeature_RNA)
max <- max(scrna@meta.data$nFeature_RNA)    
s <- sd(scrna@meta.data$nFeature_RNA)
min1 <- min(scrna@meta.data$nCount_RNA)
max1 <- max(scrna@meta.data$nCount_RNA)
m1 <- mean(scrna@meta.data$nCount_RNA)
s1 <- sd(scrna@meta.data$nCount_RNA)
Count93 <- quantile(scrna@meta.data$nCount_RNA, 0.93) # calculate value in the 93rd percentile
print(paste("Feature stats:",min,m,max,s));
print(paste("UMI stats:",min1,m1,max1,s1,Count93));

# Now, filter the data using the subset function and your chosen thresholds. Note that for large data sets with diverse samples, it may be beneficial to use sample-specific thresholds for some parameters. If you are not sure what thresholds to use, the following will work well for the purposes of this course:

scrna <- subset(x = scrna, subset = nFeature_RNA > 200  & nCount_RNA < Count93 & percent.mito < 0.2)
scrna@meta.data
table(scrna@meta.data$DataSet)
length(scrna@meta.data$nCount_RNA)


### read the srna datasets.

load("scrna.Rdata")

#check the Rdatasets
scrna.list
scrna[[]]

scrna@meta.data$DataSet <- ordered(factor(scrna@meta.data$DataSet),levels = c("D0", "D3", "D7","D28"))

# Visualization

pdf("combined_intergration.pdf", height = 6, width = 6)
p1 <- DimPlot(scrna, reduction = "umap", group.by = "DataSet")
p2 <- DimPlot(scrna, reduction = "umap", label = TRUE, repel = TRUE)
p2
p1 + p2
dev.off()
pdf("combined_intergration_Split.pdf", height = 4, width = 12)
DimPlot(scrna, reduction = "umap", split.by = "DataSet")
dev.off()

###### infer cell types

## using the published markers

pdf("Marker_Genes.pdf", height = 16, width = 8)

## for endothelium and Fibroblasts

pdf("EndotheliumFibroblast_Marker_Genes.pdf", height = 6, width = 9)
EndotheliumFibroblast<- FeaturePlot(object = scrna, features = c("Pecam1", "Cd34","Vwf", "Col3a1", "Lum","Fbn1"),   min.cutoff = 'q0', max.cutoff = 'q80', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = EndotheliumFibroblast, ncol = 1)
dev.off()

## for the Proximal Tubule cells and distal tubule cells
pdf("ProximalTubuleDistalTubule_Genes.pdf", height = 6, width = 9)
ProximalTubuleDistalTubule<- FeaturePlot(object = scrna, features = c("Slc34a1", "Lrp2","Acsm1", "Slc12a1","Slc12a3", "Pvalb"),   min.cutoff = 'q0', max.cutoff = 'q80', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = ProximalTubuleDistalTubule, ncol = 1)
dev.off()

## for the Collecting duct cells IC and PC
pdf("CollectingDuct_IC_PC_Genes.pdf", height = 6, width = 9)
CollectingDuct<- FeaturePlot(object = scrna, features = c("Slc4a1","Aqp6","Atp6v1g3","Aqp2","Scnn1g","Hsd11b2"),   min.cutoff = 'q0', max.cutoff = 'q80', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = CollectingDuct, ncol = 1)
dev.off()

## for the Urothelium



############
pdf("LoopHenle_Genes.pdf", height = 6, width = 9)
LoopHenle<- FeaturePlot(object = scrna, features = c("Fst","Bst1","Slc14a2","Slc12a1","Umod","Cldn10"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = LoopHenle, ncol = 1)
dev.off()


## for the immune cells

pdf("ImmuneTotal_Genes.pdf", height = 6, width = 9)
Immunecells<- FeaturePlot(object = scrna, features = c("Fst","Bst1","Slc14a2","Slc12a1","Umod","Cldn10"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = Immunecells, ncol = 1)
dev.off()

pdf("Mcrophage_Neutrophil_Genes.pdf", height = 6, width = 9)
MicrophageNeutrophil<- FeaturePlot(object = scrna, features = c("C1qa","C1qb","Aif1","Lst1","Lyz2","Cebpb"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = MicrophageNeutrophil, ncol = 1)
dev.off()

pdf("Dendritic_neutrophil_Genes.pdf", height = 6, width = 9)
MicrophageNeutrophil<- FeaturePlot(object = scrna, features = c("Itgax","Cxcr4","Ly6d","Lst1","Lyz2","Cebpb"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = MicrophageNeutrophil, ncol = 1)
dev.off()

pdf("T_B_Genes.pdf", height = 6, width = 9)
T_B<- FeaturePlot(object = scrna, features = c("Cxcr6","Cd247","Nkg7","Cd79a","Cd79b","Fcrl1"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = T_B, ncol = 1)
dev.off()


pdf("CD45T_Cd8effectorT_Genes.pdf", height = 6, width = 9)
CD45T_Cd8effectorT<- FeaturePlot(object = scrna, features = c("Ccr7","Ms4a4b","Klf2","Ccl5","Cd8b1","Tnfrsf4"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = CD45T_Cd8effectorT, ncol = 1)
dev.off()

pdf("NkCell_Genes.pdf", height = 6, width = 9)
NK<- FeaturePlot(object = scrna, features = c("Nkg7","Cd7","Ncr1","Ccl5","Ly6c2","Cxcr6"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = NK, ncol = 1)
dev.off()

pdf("Tregulatory_Genes.pdf", height = 6, width = 9)
Tregulatory<- FeaturePlot(object = scrna, features = c("S100a4","Rgs1","Cd3g","Ltb","Capg","Izumo1r"),   min.cutoff = 'q0', max.cutoff = 'q90', ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = Tregulatory, ncol = 1)
dev.off()


fp <- FeaturePlot(object = scrna, features = c("Pecam1", "Cd34","Vwf","Atp1b1","Slc34a1","Lrp2","Cldn8","Aqp2","Cldn8","Slc4a1", "Scnn1g","Slc12a1", "Aqp1","Col3a1", "Ptprc", "Cd74", "Lyz2"),  ncol=3, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
wrap_plots(plots = fp, ncol = 1)

dev.off()
immuneplot<-VlnPlot(scrna, features = c("Ptprc", "Cd74", "Lyz2"), group.by = "seurat_clusters",
                    pt.size = 0, combine = FALSE)
wrap_plots(plots = immuneplot, ncol = 1)
## ddd the biomarkers
# idnetify conserved cell type markers

nk.markers <- FindConservedMarkers(scrna,  grouping.var = "DataSet", verbose = FALSE)
head(nk.markers)



# perform intergration

FeaturePlot(object = scrna, 
            features = c("S100a4"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

### add the 'Pdgfrb 

DefaultAssay(scrna) <- "RNA"

FeaturePlot(object = scrna, 
            features = c("Il17f"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

pdf(sprintf("%s/Fut4.df", outdir), width = 10, height = 10);
FeaturePlot(object = scrna, 
            features = c("S100a"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();


pdf(sprintf("%s/Ptprc_immuneCells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Ptprc"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();
pdf(sprintf("%s/Lyz2_immunecells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Lyz2"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();
pdf(sprintf("%s/Cd74_immunecells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd74"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();


pdf(sprintf("%s/NCr1_NKcells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Ncr1"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Cd8a_CD8Tcell.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd8a"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Cd79a_CD4Tcell.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd79a"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Cd79b_CD4Tcell.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd79b"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/F13a1_Neutrophil.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("F13a1"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Hp_Neutrophil.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Hp"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Cd209a_DCcells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd209a"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Cd79a_Microphagecells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd79a"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Slc12a1_DistaltubularCells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Slc12a1"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();
pdf(sprintf("%s/Aqp2_CollectingDuctCells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Aqp2"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();
pdf(sprintf("%s/Cldn8_CollectingDuctCells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cldn8"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Col3a1_Fibroblasts.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Col3a1"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Vwf_EndotheliumCells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Vwf"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Cd68_Macrophages.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd68"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();


pdf(sprintf("%s/Hdc_Granulocytes.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Hdc"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();


pdf(sprintf("%s/Lrp2_ProximalTubuleCells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Lrp2"),
            min.cutoff = 'q10', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Scnn1g_PC.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Scnn1g"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Slc4a1_IC.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Slc4a9"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Upk1a_Urothelium.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Shh"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Mafb_Podocytes.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Mafb"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Gzmb_DendriticCells.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Gzmb"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();
pdf(sprintf("%s/Cd3e_activatedTCell.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd3e"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off()

pdf(sprintf("%s/Aqp1_thinlimbofloopofHenle.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Aqp1"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)
dev.off()

pdf(sprintf("%s/Slc12a1_thickascendinglimbofloopofHenle.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Slc12a1"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)
dev.off()

pdf(sprintf("%s/Fst_descendinglimbofloopofHenle.pdf", outdir), width = 12, height = 4);

FeaturePlot(object = scrna, 
            features = c("Fst"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off()

FeaturePlot(object = scrna, 
            features = c("Cxcr6"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

pdf(sprintf("%s/Slc12a3_TubuleCells", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Slc12a3"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();
pdf(sprintf("%s/Atp1b1_TubuleCells", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Atp1b1"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

pdf(sprintf("%s/Atp1b1_TubuleCells", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Acta2"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

### cell death/ Atf3, Btg2, Cd14, Cdkn1a
#pdf(sprintf("%s/Atp1b1_TubuleCells", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Cd14"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

#dev.off();

#### Interfereon:
FeaturePlot(object = scrna, 
            features = c("Tmcc1"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)
### Tnfa signlaing via NFkb

FeaturePlot(object = scrna, 
            features = c("Cxcl2"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)
FeaturePlot(object = scrna, 
            features = c("Egr1"),
            min.cutoff = 'q5', 
            #max.cutoff = 'q90', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

## measure the vlnplot
VlnPlot(scrna, features = c("Ptprc"), group.by  ="DataSet")

scrna[[]]
pdf("Ptprc_vlnplot_immune.pdf",width = 15, height = 3 )

immuneplot<-VlnPlot(scrna, features = c("Ptprc", "Cd74", "Lyz2"), split.by = "DataSet", group.by = "seurat_clusters",
        pt.size = 0, combine = FALSE)
wrap_plots(plots = immuneplot, ncol = 1)
dev.off()






print(plot_grid(p1, p2,ncol = 1));
dev.off();

### add the 'Pdgfrb 

FeaturePlot(object = scrna, 
            features = c("Scnn1g"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

pdf(sprintf("%s/Fut4.df", outdir), width = 10, height = 10);
FeaturePlot(object = scrna, 
            features = c("Fut4"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();


pdf(sprintf("%s/Ptprc.pdf", outdir), width = 12, height = 4);
FeaturePlot(object = scrna, 
            features = c("Ptprc"),
            min.cutoff = 'q0', 
            reduction = "umap",
            split.by ="DataSet",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

dev.off();

### customize featureplot in mutiple condtions

#add some dummy metadata
scrna[[]]

# measure how many cells each conditon
table(scrna@meta.data$DataSet)




?DimPlot
# Color the t-SNE and UMAP plots by some potential confounding variables. Here’s an example in which we color each cell according to the number of UMIs it contains:

feature.pal = rev(colorRampPalette(brewer.pal(11,"Spectral"))(50)); # a useful color palette
pdf(sprintf("%s/umap.%d.colorby.UMI.pdf", outdir, nPC), width = 10, height = 8);
fp <- FeaturePlot(object = scrna, features = c("nCount_RNA"), cols = feature.pal, pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank()); # the text after the ‘+’ simply removes the axis using ggplot syntax
print(fp);
dev.off();

# To investigate this, plot several principal components on the t-SNE/UMAP, for example the following code plots the first principal component and prints the plot to a file:

pdf(sprintf("%s/UMAP.%d.colorby.PCs.pdf", outdir, nPC), width = 12, height = 6);
redblue=c("blue","gray","red"); # another useful color scheme
fp1 <- FeaturePlot(object = scrna, features = 'PC_1', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(fp1);
dev.off();







### Normalizing the data
mydata <- NormalizeData(Seurat_object_D0, normalization.method = "LogNormalize", scale.factor = 10000)


### identifcaiton of highly variable features
mydata <- FindVariableFeatures(mydata,selection.method = "vst", nfeatures = 2000)

#### identify the 10 most highly varible genes
top10<-head(VariableFeatures(mydata),10)
top10
### plot varible features with and without labels
plot1<-VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge=0, ynudge=0)
?LabelPoints
plot1+plot2

### SCale the data we here used the scaling
all.genes<-rownames(mydata)

mydata<-ScaleData(mydata, features = all.genes)
mydata

### Two method are used to save the expression information of each gene
#ScaledData<-as.matrix(GetAssayData(object = mydata[["RNA"]], slot = "data"))
#ScaledData

NormalData<-as.matrix(mydata@assays[["RNA"]]@counts)
### crea

### perform linear dimensial reduction

mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))
# Examine and visualize PCA results a few different ways
print(mydata[["pca"]], dims = 1:5, nfeatures = 5)

# plot these genes
VizDimLoadings(mydata, dims = 1:2, reduction = "pca")

### generate the visualizing both cells and features that define the PCA,
DimPlot(mydata, reduction = "pca")
### generate the  easy exploration of the primary sources of heterogeneity in a dataset
DimHeatmap(mydata, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mydata, dims = 1:15, cells = 500, balanced = TRUE)

### Determine the ‘dimensionality’ of the dataset
## To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set. 
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
mydata <- JackStraw(mydata, num.replicate = 100)
mydata <- ScoreJackStraw(mydata, dims = 1:20)

###The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).
###In this case it appears that there is a sharp drop-off in significance after the first 10-12 PCs.

JackStrawPlot(mydata, dims = 1:15)

# An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.

ElbowPlot(mydata)



# t-SNE and Clustering
mydata <- RunUMAP(mydata, reduction = "pca", dims = 1:8)
mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:8)
mydata <- FindClusters(mydata, resolution = 0.2)
#install.packages("reticulate")


### mecom

FeaturePlot(object = mydata, 
            features = c("Fut4"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)



Ens<-row.names(mydata)

dim(mydata)

DimPlot(mydata,reduction = "umap")

mydata

### remove the subcluster of endothelium
?subset
sub_obj <-subset(mydata, idents = c(0,1,2,4,5,6))

sub_obj <- subset(mydata, idents = 3, invert=TRUE)

DimPlot(sub_obj,reduction = "umap")


### Then Run clustering 
sub_obj <- RunUMAP(sub_obj, reduction = "pca", dims = 1:12)
sub_obj <- FindNeighbors(sub_obj, reduction = "pca", dims = 1:12)
sub_obj <- FindClusters(sub_obj, resolution = 0.2)

DimPlot(sub_obj,reduction = "umap")


### mecom

FeaturePlot(object = sub_obj, 
            features = c("ENSG00000085276"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)











### 1.2 strandard pre-processing workflow (removed MT genes <20)

## QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^MT-")
mean(mydata[["percent.mt"]]$percent.mt)
mean(mydata$nCount_RNA)
mean(mydata$nFeature_RNA)
mydata$orig.ident




### using the gene ID 
#raw_counts_1<- read.table("GSE151186_RAW/GSM4568142_filtered_feature_bc_matrix_h5.witMT.GeneID.csv",sep = ",", header = T,row.names = 1)

#head(raw_counts_1)


mydata<-CreateSeuratObject(counts=raw_counts_1, min.cells=3, min.genes=200, project ="scRNAseq")

### Compute the length of cells for mydata
row.names(mydata$nFeature_RNA)

length(mydata$orig.ident)

write.table(mydata$nCount_RNA,"ncout_RNA.tsv",sep = "\t")
length(mydata$nCount_RNA)

length(mydata$nCount_RNA)
length(unique(mydata$nFeature_RNA))
### 1.2 strandard pre-processing workflow (removed MT genes <20)

## QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^MT-")
mean(mydata[["percent.mt"]]$percent.mt)
mean(mydata$nCount_RNA)
mean(mydata$nFeature_RNA)
mydata$orig.ident




# Visualize QC metrics as a violin plot
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

?VlnPlot
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#### We filtered the We filter cells that have unique feature counts over 2,500 or less than 200
#### We filter cells that have >5% mitochondrial counts

mydata <- subset(mydata, subset = nFeature_RNA > 200  & percent.mt < 20)

### Compute the length of cells for mydata
length(unique(mydata@meta.data$nCount_RNA))

length(mydata$nCount_RNA)
mydata
mean(mydata$nCount_RNA)
mean(mydata$nFeature_RNA)
mean(mydata$percent.mt)
### Normalizing the data
mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000)

### identifcaiton of highly variable features
mydata <- FindVariableFeatures(mydata,selection.method = "vst", nfeatures = 2000)

#### identify the 10 most highly varible genes
top10<-head(VariableFeatures(mydata),10)
top10
### plot varible features with and without labels
plot1<-VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge=0, ynudge=0)
?LabelPoints
plot1+plot2

### SCale the data we here used the scaling
all.genes<-rownames(mydata)

mydata<-ScaleData(mydata, features = all.genes)
mydata

### Two method are used to save the expression information of each gene
#ScaledData<-as.matrix(GetAssayData(object = mydata[["RNA"]], slot = "data"))
#ScaledData

NormalData<-as.matrix(mydata@assays[["RNA"]]@counts)
### crea

### perform linear dimensial reduction

mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))
# Examine and visualize PCA results a few different ways
print(mydata[["pca"]], dims = 1:5, nfeatures = 5)

# plot these genes
VizDimLoadings(mydata, dims = 1:2, reduction = "pca")

### generate the visualizing both cells and features that define the PCA,
DimPlot(mydata, reduction = "pca")
### generate the  easy exploration of the primary sources of heterogeneity in a dataset
DimHeatmap(mydata, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mydata, dims = 1:15, cells = 500, balanced = TRUE)

### Determine the ‘dimensionality’ of the dataset
## To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set. 
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
mydata <- JackStraw(mydata, num.replicate = 100)
mydata <- ScoreJackStraw(mydata, dims = 1:20)

###The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line).
###In this case it appears that there is a sharp drop-off in significance after the first 10-12 PCs.

JackStrawPlot(mydata, dims = 1:15)

# An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.

ElbowPlot(mydata)



# t-SNE and Clustering
mydata <- RunUMAP(mydata, reduction = "pca", dims = 1:12)
mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:12)
mydata <- FindClusters(mydata, resolution = 0.2)
#install.packages("reticulate")

Ens<-row.names(mydata)

dim(mydata)

DimPlot(mydata,reduction = "umap")

mydata

### remove the subcluster of endothelium
?subset
sub_obj <-subset(mydata, idents = c(0,1,2,4,5,6))

sub_obj <- subset(mydata, idents = 3, invert=TRUE)

DimPlot(sub_obj,reduction = "umap")


### Then Run clustering 
sub_obj <- RunUMAP(sub_obj, reduction = "pca", dims = 1:12)
sub_obj <- FindNeighbors(sub_obj, reduction = "pca", dims = 1:12)
sub_obj <- FindClusters(sub_obj, resolution = 0.2)

DimPlot(sub_obj,reduction = "umap")


### mecom

FeaturePlot(object = sub_obj, 
            features = c("ENSG00000085276"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

VlnPlot(sub_obj, features = c("ENSG00000085276"),pt.size = 0, assay = "RNA") 

## 

## check BMP10 ENSG00000163217 (cardiogenic mesoderm)


FeaturePlot(object = sub_obj, 
            features = c("ENSG00000163217"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
## check BMP4 ENSG00000125378 (cardiogenic mesoderm)

FeaturePlot(object = sub_obj, 
            features = c("ENSG00000125378"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)


### Building trajectories with Monocle 3

install.packages("Signac")
install.packages("SeuratWrappers")
install.packages("monocle3")
install.packages("patchwork")

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)
### only select the subcluster 
sub_objEndocadium<- subset(mydata, idents = c(0,1,2))

## then perform the tragectories analyses



## find all markers of cluster 2
cluster2.markers <- FindMarkers(mydata, ident.1 = 2, min.pct =  0.25)
head(cluster2.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

pbmc.markers
### merge marker with gene id
genes<-read.table("../../Database/ENSEMBL_GENEID.txt", row.names = 1)
pbmc.markers$gene

merge()

write.table(pbmc.markers, "ClusterMarker.txt", sep="\t")

### BMP-2 ENSG00000125845 

FeaturePlot(object = mydata, 
            features = c("ENSG00000125845"),
            min.cutoff = 'q0',
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### CDH2 ENSG00000170558
FeaturePlot(object = mydata, 
            features = c("ENSG00000170558"),
            min.cutoff = 'q0',
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### OCT4 ENSG00000204531

FeaturePlot(object = mydata, 
            features = c("ENSG00000204531"),
            min.cutoff = 'q0',
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### SOX2 ENSG00000181449 
FeaturePlot(object = mydata, 
            features = c("ENSG00000181449"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

## BMP10 ENSG00000163217 
FeaturePlot(object = mydata, 
            features = c("ENSG00000163217"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### Check the MEP 



### CD90 ENSG00000154096
FeaturePlot(object = mydata, 
            features = c("ENSG00000154096"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### 1.  ENSG00000113721

FeaturePlot(object = mydata, 
            features = c("ENSG00000113721"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

## 3. SIRPA ENSG00000198053
FeaturePlot(object = mydata, 
            features = c("ENSG00000198053"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

FeaturePlot(object = mydata, 
            features = c("ENSG00000174059"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### FOXA2 ENSG00000125798

FeaturePlot(object = mydata, 
            features = c("ENSG00000125798"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)


## check BMP10 ENSG00000163217 (cardiogenic mesoderm)


FeaturePlot(object = mydata, 
            features = c("ENSG00000163217"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

VlnPlot(mydata, features = c("ENSG00000163217"),pt.size = 0, assay = "RNA") 

## check BMP4 ENSG00000125378 (cardiogenic mesoderm)

FeaturePlot(object = sub_obj, 
            features = c("ENSG00000125378"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### endocardial progenitors (NKX2-5, ISL1 and ETV2)

VlnPlot(mydata, features = c("ENSG00000125378"),pt.size = 0, assay = "RNA") 

## 1 NKX2-5 ENSG00000183072

FeaturePlot(object = mydata, 
            features = c("ENSG00000183072"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### 2. ISL1 ENSG00000016082

FeaturePlot(object = mydata, 
            features = c("ENSG00000016082"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### 3. ETV2 ENSG00000105672 
FeaturePlot(object = mydata, 
            features = c("ENSG00000105672"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

## check endocadium NPR3 ENSG00000113389
DotPlot(object = mydata, features = c("ENSG00000113389"), assay = "RNA")  
DotPlot(object = mydata, features = c("ENSG00000113389"), assay = "RNA")

FeaturePlot(object = mydata, 
            features = c("ENSG00000113389"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

## check endocadium GATA5 ENSG00000130700 
FeaturePlot(object = mydata, 
            features = c("ENSG00000130700"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
## check endocadium GATA4 ENSG00000136574 
FeaturePlot(object = mydata, 
            features = c("ENSG00000136574"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

## check endocadium NFATC1 ENSG00000131196 
FeaturePlot(object = mydata, 
            features = c("ENSG00000131196"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
## check endocadium NRG1 ENSG00000157168
FeaturePlot(object = mydata, 
            features = c("ENSG00000157168"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check cor CXCR4  ENSG00000121966 are markers of coronary endothelium 
FeaturePlot(object = mydata, 
            features = c("ENSG00000121966"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### check cor APLN  ENSG00000171388 are markers of coronary endothelium 
FeaturePlot(object = mydata, 
            features = c("ENSG00000171388"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check PECAM CD31  ENSG00000261371
FeaturePlot(object = mydata, 
            features = c("ENSG00000261371"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check CHD5/ CD144 ENSG00000179776 
FeaturePlot(object = mydata, 
            features = c("ENSG00000179776"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check CDH11 ENSG00000140937 (Mesenchymal marker)
FeaturePlot(object = mydata, 
            features = c("ENSG00000140937"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### check SOX9 ENSG00000125398  (Mesenchymal marker)
FeaturePlot(object = mydata, 
            features = c("ENSG00000125398"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### check POSTN ENSG00000133110   (Mesenchymal marker),
FeaturePlot(object = mydata, 
            features = c("ENSG00000133110"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

####MSX2 ENSG00000120149

FeaturePlot(object = mydata, 
            features = c("ENSG00000120149"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### VIM ENSG00000026025
FeaturePlot(object = mydata, 
            features = c("ENSG00000026025"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### check COL1A1 ENSG00000108821 (Mesenchymal marker),
FeaturePlot(object = mydata, 
            features = c("ENSG00000108821"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

#### NRG1 induce Trabecular Myocardium
### trabecular markers BMP10, NPPA, NPPB, and IRX3
## 1. BMP10
FeaturePlot(object = mydata, 
            features = c("ENSG00000163217"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
###2. NPPB,ENSG00000120937
FeaturePlot(object = mydata, 
            features = c("ENSG00000120937"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### 3. IRX3 ENSG00000177508

FeaturePlot(object = mydata, 
            features = c("ENSG00000177508"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
###4. NPPA, ENSG00000175206
FeaturePlot(object = mydata, 
            features = c("ENSG00000175206"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

## 4. GJA5, ENSG00000265107 distinguishing marker of the trabecular lineagein the early embryo
FeaturePlot(object = mydata, 
            features = c("ENSG00000265107"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check VIM gene ENSG00000026025 
FeaturePlot(object = mydata, 
            features = c("ENSG00000026025"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)


## check NKX2-5 ENSG00000183072 
FeaturePlot(object = mydata, 
            features = c("ENSG00000183072"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

# check  NKAIN2 ENSG00000188580
FeaturePlot(object = mydata, 
            features = c("ENSG00000188580"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check Endothelium PECOM1 ENSG00000261371 
FeaturePlot(object = mydata, 
            features = c("ENSG00000261371"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
### check Endotherlium HEY2 ENSG00000135547
FeaturePlot(object = mydata, 
            features = c("ENSG00000135547"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

##### check the proginitor markers:CD133 ENSG00000007062  

FeaturePlot(object = mydata, 
            features = c("ENSG00000007062"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check Endothelium CXCR4 ENSG00000121966
FeaturePlot(object = mydata, 
            features = c("ENSG00000128052"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)



### check Endothelium APLN ENSG00000171388 
FeaturePlot(object = mydata, 
            features = c("ENSG00000171388"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### Check cadiomycyte : SMPX  ENSG00000091482 ENSG00000122367  TNNI3:ENSG00000129991

FeaturePlot(object = mydata, 
            features = c("ENSG00000129991"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

## Check cadiomycyte : SIRPA ENSG00000198053
FeaturePlot(object = mydata, 
            features = c("ENSG00000198053"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)


DefaultAssay(mydata) <- "RNA"
## check ISL1 :ENSG00000016082
FeaturePlot(object = mydata, 
            features = c("ENSG00000016082"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)



### Visualizing MECOM
pdf("MECOM_undeveloped_com.pdf")
FeaturePlot(object = mydata, 
            features = c("ENSG00000085276"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

VlnPlot(mydata, features = c("ENSG00000085276"),pt.size = 0, assay = "RNA") 






















#py_install('umap-learn', pip = T, pip_ignore_installed = T) 

?DimPlot
p1 <- DimPlot(EC.combined_intergrated, reduction = "umap", group.by = "stim",pt.size = 0.1)
p2 <- DimPlot(EC.combined_intergrated, reduction = "umap", pt.size = 0.1,label = TRUE)

DimPlot(EC.combined_intergrated, reduction = "umap", group.by = "stim",pt.size = 0.1)

DimPlot(EC.combined_intergrated, reduction = "umap", pt.size = 0.1,label = TRUE)
p1
p2+p1

###
#p1 <- DimPlot(EC.combined, reduction = "pca", group.by = "D12")
#p2 <- DimPlot(EC.combined, reduction = "pca", label = TRUE)
pdf("Undeveloped_Control_comp.pdf", width=12, height=5)
plot_grid(p1, p2)

dev.off()
DimPlot(EC.combined_intergrated, reduction = "umap",pt.size = 0.1)

### Feature for FBL ENSG00000105202


FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000105202"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",
            pt.size = 0.2,
            repel = TRUE)
DotPlot(object = EC.combined_intergrated, features = c("ENSG00000105202"), split.by = "stim", assay = "RNA") 
VlnPlot(EC.combined_intergrated, features = c("ENSG00000105202"), split.by = "stim",split.plot = TRUE,pt.size = 0, assay = "RNA") 













## read the second files

### read the read counts for each gene in each cells
raw_counts_2<- read.table("../Single_Cells/iPSC/GSM4125588_d83_LVp_filtered_gene_bc_matrices_h5.withMT.csv",sep = ",", header = T,row.names = 1)


head(raw_counts_2)

mydata_2<-CreateSeuratObject(counts=raw_counts_2, min.cells=3, min.genes=200, project ="scRNAseq")


### strandard pre-processing workflow

## QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
mydata_2[["percent.mt"]] <- PercentageFeatureSet(mydata_2, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(mydata_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(mydata_2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mydata_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#### We filtered the We filter cells that have unique feature counts over 2,500 or less than 200
#### We filter cells that have >5% mitochondrial counts

mydata_2 <- subset(mydata_2, subset = nFeature_RNA > 200  & percent.mt < 20)


### Normalizing the data
mydata_2 <- NormalizeData(mydata_2, normalization.method = "LogNormalize", scale.factor = 10000)

###For clarity, in this previous line of code (and in future commands), we provide the default values for certain parameters in the function call. However, this isn't required and the same behavior can be achieved with
#mydata_2 <- NormalizeData(mydata_2)

### SCale the data we here used the scaling
all.genes_2<-rownames(mydata_2)

all.genes_2

mydata_2<-ScaleData(mydata_2, features = all.genes_2)


### Two method are used to save the expression information of each gene
ScaledData_2<-as.matrix(GetAssayData(object = mydata_2[["RNA"]], slot = "data"))
ScaledData_2

NormalData_2<-as.matrix(mydata_2@assays[["RNA"]]@counts)



### Then intergrate the 
mydata$stim <- "hPSC-derived"
mydata_2$stim<- "Fetal"

EC.anchors_integrated <- FindIntegrationAnchors(object.list = list(mydata, mydata_2), dims = 1:30)
EC.combined_intergrated <- IntegrateData(anchorset = EC.anchors_integrated, dims = 1:30)

DefaultAssay(EC.combined_intergrated) <- "integrated"

### Run the standard work flow visualizaiton and cluster



## scale data
EC.combined_intergrated <- ScaleData(EC.combined_intergrated, verbose = FALSE)
EC.combined_intergrated<-RunPCA(EC.combined_intergrated,npcs = 30, verbose = FALSE)



#DimHeatmap(EC.combined, dims = c(1,2,3,4,5,6))

# t-SNE and Clustering
EC.combined_intergrated <- RunUMAP(EC.combined_intergrated, reduction = "pca", dims = 1:30)
EC.combined_intergrated <- FindNeighbors(EC.combined_intergrated, reduction = "pca", dims = 1:30)
EC.combined_intergrated <- FindClusters(EC.combined_intergrated, resolution = 0.05)
#install.packages("reticulate")


#py_install('umap-learn', pip = T, pip_ignore_installed = T) 

?DimPlot
p1 <- DimPlot(EC.combined_intergrated, reduction = "umap", group.by = "stim",pt.size = 0.1)
p2 <- DimPlot(EC.combined_intergrated, reduction = "umap", pt.size = 0.1,label = TRUE)

DimPlot(EC.combined_intergrated, reduction = "umap", group.by = "stim",pt.size = 0.1)

DimPlot(EC.combined_intergrated, reduction = "umap", pt.size = 0.1,label = TRUE)
p1
p2+p1

###
#p1 <- DimPlot(EC.combined, reduction = "pca", group.by = "D12")
#p2 <- DimPlot(EC.combined, reduction = "pca", label = TRUE)
pdf("Undeveloped_Control_comp.pdf", width=12, height=5)
plot_grid(p1, p2)

dev.off()
DimPlot(EC.combined_intergrated, reduction = "umap",pt.size = 0.1)

### Feature for FBL ENSG00000105202


FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000105202"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",
            pt.size = 0.2,
            repel = TRUE)
DotPlot(object = EC.combined_intergrated, features = c("ENSG00000105202"), split.by = "stim", assay = "RNA") 
VlnPlot(EC.combined_intergrated, features = c("ENSG00000105202"), split.by = "stim",split.plot = TRUE,pt.size = 0, assay = "RNA") 


## check endocadium NPR3 ENSG00000113389
DotPlot(object = EC.combined_intergrated, features = c("ENSG00000113389"), split.by = "stim", assay = "RNA")  
DotPlot(object = EC.combined_intergrated, features = c("ENSG00000113389"), group.by  = "stim",assay = "RNA")


FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000113389"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check PECAM CD31  ENSG00000261371
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000261371"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check CHD5/ CD144 ENSG00000179776 
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000179776"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check CDH11 ENSG00000140937  
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000140937"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)
## check NKX2-5 ENSG00000183072 
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000183072"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

# check  NKAIN2 ENSG00000188580
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000188580"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check Endothelium PECOM1 ENSG00000261371 
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000261371"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)


##### check the proginitor markers:CD133 ENSG00000007062  

FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000007062"),
            min.cutoff = 'q5', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### check Endothelium CXCR4 ENSG00000121966
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000128052"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)



### check Endothelium APLN ENSG00000171388 
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000171388"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

### Check cadiomycyte : SMPX  ENSG00000091482 ENSG00000122367  TNNI3:ENSG00000129991

FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000129991"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)

## Check cadiomycyte : SIRPA ENSG00000198053
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000198053"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)


DefaultAssay(EC.combined_intergrated) <- "RNA"
## check ISL1 :ENSG00000016082
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000016082"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            pt.size = 0.2,
            repel = TRUE)



### Visualizing MECOM
pdf("MECOM_undeveloped_com.pdf")
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000085276"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)
VlnPlot(EC.combined_intergrated, features = c("ENSG00000085276"), split.by = "stim",split.plot = TRUE,pt.size = 0, assay = "RNA") 

dev.off()

?FeaturePlot
EC.combined_intergrated@meta.data$nCount_RNA

?VlnPlot
pdf ("MECOM_Cluster_Expression_withdot.pdf",width=12, height=5)

install.packages("RColorBrewer")
VlnPlot(EC.combined_intergrated, features = c("ENSG00000085276"), split.by = "stim",split.plot = TRUE,pt.size = 0.01, assay = "RNA") 
??Vlnplot

VlnPlot(EC.combined_intergrated, features = c("ENSG00000085276"), split.by = "stim",split.plot = TRUE,pt.size=0, assay = "RNA") 


### Other Srterial Markers, such as ENSG00000265107 -- GJA5
VlnPlot(EC.combined_intergrated, features = c("ENSG00000265107"), split.by = "stim",split.plot = TRUE,pt.size=0, assay = "RNA") 

FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000265107"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

### Second Arterial Marker,ENSG00000128917  DLL4 
VlnPlot(EC.combined_intergrated, features = c("ENSG00000128917"), split.by = "stim",split.plot = TRUE,pt.size=0, assay = "RNA") 

FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000128917"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)

### The third Aterial Marker , ENSG00000164683 HEY1

VlnPlot(EC.combined_intergrated, features = c("ENSG00000164683"), split.by = "stim",split.plot = TRUE,pt.size=0, assay = "RNA") 
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000164683"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)


### THe ENSG00000136715 SAP130

VlnPlot(EC.combined_intergrated, features = c("ENSG00000136715"), split.by = "stim",split.plot = TRUE,pt.size=0, assay = "RNA") 
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000136715"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)


### the ENSG00000204961 Pcdha9 

VlnPlot(EC.combined_intergrated, features = c("ENSG00000204961"), split.by = "stim",split.plot = TRUE,pt.size=0, assay = "RNA") 

FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000204961"),
            min.cutoff = 'q0', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",
            pt.size = 0.2, 
            # col= rev(brewer.pal(n = 11, name = "RdBu")),
            repel = TRUE)
### 
DotPlot(object = EC.combined_intergrated, features = c("ENSG00000085276"), split.by = "stim", assay = "RNA")  
DotPlot(object = EC.combined_intergrated, features = c("ENSG00000085276"), group.by  = "stim",assay = "RNA")  
VlnPlot(EC.combined_intergrated, features = c("ENSG00000085276"), group.by = "stim",
        pt.size = 0.01, combine = FALSE, log=FALSE,assay = "RNA") 

FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000099822"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",pt.size = 0.2,
            repel = TRUE)

dev.off()


### Check rbfox2 
FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000100320"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",pt.size = 0.2,
            repel = TRUE)
DotPlot(object = EC.combined_intergrated, features = c("ENSG00000100320"), split.by = "stim", assay = "RNA")  
################# The marker information ############


###Identification of conserved markers in all conditions

# Find markers for every cluster compared to all remaining cells, report only the positive ones

pbmc.markers <- FindAllMarkers(EC.combined_intergrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Final_markers<-pbmc.markers %>% group_by(cluster) %>% top_n(n = 30)

write.table(Final_markers,"Undeveloped_Control_Com_Endocadium_top30_v0511.txt",sep = "\t")


##### here is to find specific cluster markers 

#new.cluster.ids <- c("Endocardium","Endothelium", "Endothelium - Pericyte", "Gamma delta T cells", 
#                     "Gamma delta T cells & Fibroblasts", "Fibroblast I","Lymphatic EC","Endothelium_EC-X", "Macrophage",  "Fibroblast", "Erythroid-like and erythroid precursor cells",
#                   "Cardiomyocyte & Fibroblast","Cardiomyocyte")






#### Annotation of the cluster
new.cluster.ids <- c("Endocardium","Choriocapillaris ECs", "Arterial ECs", "Venular ECs", 
                     "Angiogenic ECs", "Lymphatic ECs","Interferon-activated ECs", "Endocardium II",  "Microphage", "Erythroid precursor cells", "Cardiomyocyte")
names(new.cluster.ids) <- levels(EC.combined_intergrated)
EC.combined_intergrated <- RenameIdents(EC.combined_intergrated, new.cluster.ids)
DimPlot(EC.combined_intergrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



### Draw the marker genes and bubble plot:


##### Then use DotPlot to draw the cluster using marker genes

Idents(EC.combined_intergrated)<-factors(Idents(EC.combined_intergrated), levels=c("Endocardium","Choriocapillaris ECs", "Arterial ECs", "Venular ECs", 
                                                                                   "Angiogenic ECs", "Lymphatic EC","Interferon-activated ECs", "Endocardium II",  "Microphage", "Erythroid precursor cells", "Cardiomyocyte"))

MakersToplot<- c("ENSG00000133110","ENSG00000140937","ENSG00000113389","ENSG00000167434","ENSG00000102760","ENSG00000179403","ENSG00000140092","ENSG00000164683","ENSG00000184113","ENSG00000148773","ENSG00000168078","ENSG00000137804","ENSG00000143248","ENSG00000074181","ENSG00000185633","ENSG00000133800","ENSG00000163993","ENSG00000117707","ENSG00000171848","ENSG00000166803","ENSG00000167900","ENSG00000185774","ENSG00000180730","ENSG00000188580",
                 "ENSG00000011600","ENSG00000173369","ENSG00000159189","ENSG00000213934","ENSG00000206172","ENSG00000244734","ENSG00000205678","ENSG00000122367","ENSG00000091482","ENSG00000085276")

DotPlot(EC.combined_intergrated,features =MakersToplot, cols = c("gray", "blue") ) +RotatedAxis()
?DotPlot


#### identify the conserved markers:
FindConservedMarkers(EC.combined_intergrated, ident.1 = 1,
                     grouping.var = "stim",
                     only.pos = TRUE,
                     min.diff.pct = 0.25,
                     min.pct = 0.25,
                     logfc.threshold = 0.25)

MarkerGenes<-c()

## GAPDH ENSG00000111640 have almost no diffeerne  : Selection of reference genes for gene expression studies in heart failure for left and right ventricles
DotPlot(object = EC.combined_intergrated, features = c("ENSG00000085276","ENSG00000100320", "ENSG00000179403","ENSG00000159189","ENSG00000075624","ENSG00000111640"), cols = c("gray", "blue"),split.by = "stim", assay = "RNA",scale.max = 100, scale.min = 0, col.max = 3, col.min=0)  



DotPlot(object = EC.combined_intergrated, features = c("ENSG00000085276","ENSG00000100320", "ENSG00000175792","ENSG00000183207","ENSG00000111640"), cols = c("gray", "blue"),split.by = "stim", assay = "RNA",scale.max = 100, scale.min = 0, col.max = 3, col.min=0)  


FeaturePlot(object = EC.combined_intergrated, 
            features = c("ENSG00000100320"),
            min.cutoff = 'q10', 
            reduction = "umap",
            label = FALSE,
            split.by = "stim",pt.size = 0.2,
            repel = TRUE)

### Draw the feature of ENSG00000100320





  
