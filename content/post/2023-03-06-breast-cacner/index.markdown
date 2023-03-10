---

title: 'Breast cancer'
author: 'Aviv Leemann'
date: '2023-03-06'
slug: breast-cancer
categories: ['single cell','breast-cancer']
subtitle: 'single cell Breast cancer'
summary:  'Breast cancer data was obtained from single cell portal, single cell analysis executed with **R** program and Seurat package, Pallad expression was examined in Breast cancer data'
lastmod: '2023-03-06T17:52:47+02:00'
featured: no
image:
  caption: ''
  focal_point: ''                                                                            
  preview_only: no
projects: ['single cell','Breast cancer']

---




# single cell Breast cancer -analysis

Breast cancer data was obtained from single cell portal. single cell analysis executed with **R** program and Seurat package, Pallad expression was examined in Breast cancer data.
our lab found PALLD express in breast cancr, PALLD expression was examined between different cell type , different cluster groups, alsotop 10 mrker for each cluster was obtained

## libreries


```r
library(pacman)
library(Seurat)
labrvec <- c("tidyverse","ggplot2","Matrix","AnnotationDbi","readxl",
             "Seurat","future","patchwork","stringr")
pacman::p_load(labrvec,install = TRUE,character.only = TRUE)
```

## files


```r
# reading RDS
BCdata <- readRDS("BCdatageneral.rds")
```

## pre-processing

wont be used at the creation of this files because already created and RDS will be use


```r
# define nultiprocess for seurat
# 
# Breast.cancer.SCdata<-Read10X(
#   data.dir = indir ,
#   gene.column = 1,
#   unique.features = TRUE,
#   strip.suffix = FALSE
# ) 
# #class(Breast.cancer.SCdata)
# print(object.size(Breast.cancer.SCdata), units="Gb")
# # head(Breast.cancer.SCdata)
# mini.atlas <- read.csv('C:/Users/user/Desktop/BCsinglecelldata/Whole_miniatlas_meta.csv')
# gene_targets <- read_excel("C:/Users/user/Desktop/BCsinglecelldata/spikesearchtargetgenes.xlsx",)
# # summary(mini.atlas)
# gene_targets <- pull(gene_targets,"Gene name")
# gene_targets <- gene_targets[!is.na(gene_targets)]
# gene_targets <- gene_targets[1:128]
# gene_targets <- gene_targets[!duplicated(gene_targets)]
# #gene_targets
# colnames(mini.atlas)[1] <- "NAMES"
# #mini.atlas
# rownames(mini.atlas) <- pull(mini.atlas,"NAMES") 
# #mini.atlas.C <- mini.atlas %>% filter(normal_cell_call == "cancer")
# cancer.cell.types <- pull(mini.atlas,"normal_cell_call") %>%  unique()
# #cancer.cell.types
# mini.atlas <- mini.atlas[-1,-1]
# # head(mini.atlas)
# BCdata <- CreateSeuratObject(
#   Breast.cancer.SCdata,
#   project = "BC",
#   assay = "RNA",
#   min.cells = 0,
#   min.features = 0,
#   names.field = 1,
#   names.delim = "_",
#   meta.data = mini.atlas
# )

# looking at the idents in seurat object
# head(BCdata[[]])
# #BCdata
# BCdata@meta.data$project <- "BC"
# BCdata <- StashIdent(object = BCdata, save.name = "old.ident")
# Idents(object = BCdata) <- "project"
# Idents(object = BCdata)
#BCdata <- subset(x = BCdata, idents = "cancer")
# row.names(BCdata)
# # BCdata[[]]
# # extract the mitochondrial DNA
# BCdata[["Percent_mito"]] <-unlist(BCdata[["Percent_mito"]]) %>% as.numeric()
# BCdata[["nCount_RNA"]] <-unlist(BCdata[["nCount_RNA"]]) %>% as.numeric()
# BCdata[["nFeature_RNA"]] <-unlist(BCdata[["nFeature_RNA"]]) %>% as.numeric()
# BCdata[["percent.mt"]] <- PercentageFeatureSet(object = BCdata, pattern = "^MT-")
# #Idents(object = BCdata) <- "celltype_subset"
# BCdata[["normal_cell_call"]] %>%  unique()
# fixing to column Percent_mito as numeric
```

## Feature-scatter & violin


```r
Idents(BCdata) <- "Breast cancer"

plot1 <- FeatureScatter(BCdata, feature1 = "nCount_RNA", feature2 = "Percent_mito",raster=FALSE)
plot2 <- FeatureScatter(BCdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE)
plot2 + plot1
```

<img src="{{< blogdown/postref >}}index_files/figure-html/feature scatter-1.jpeg" width="2100" />

```r
VlnPlot(BCdata, features = c("nFeature_RNA", "nCount_RNA","Percent_mito"), ncol = 3,pt.size  = 0)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/feature scatter-2.jpeg" width="2100" />

```r
VlnPlot(BCdata, features = c("nFeature_RNA", "nCount_RNA","Percent_mito"), ncol = 3)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/feature scatter-3.jpeg" width="2100" />

## variable Feature-plot

BCdata \<- NormalizeData(BCdata)

BCdata \<- FindVariableFeatures(BCdata, selection.method = "vst", nfeatures = 2000)


```r
## what are our variable genes ?
top10 <- head(VariableFeatures(BCdata), 10)
top10
```

```
##  [1] "SCGB2A2"   "IGHV1-69"  "SCGB1B2P"  "IGFBP7"    "IGHG1"     "CCL4"     
##  [7] "IGLV2-28"  "IGKV2D-28" "IGKV2D-29" "CD74"
```

```r
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BCdata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

<img src="{{< blogdown/postref >}}index_files/figure-html/variable Feature plot-1.jpeg" width="2100" />

## PCA-plot


```r
# all_genes <- rownames(BCdata)
# # all_genes
# BCdata <- ScaleData(BCdata)
# BCdata <- RunPCA(BCdata, features = VariableFeatures(object = BCdata))
VizDimLoadings(BCdata, dims = 1:2, reduction = "pca")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/PCA plot-1.jpeg" width="2100" />

```r
DimPlot(BCdata, reduction = "pca")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/PCA plot-2.jpeg" width="2100" />

## Heat-map


```r
# Examine and visualize PCA results a few different ways
#(BCdata[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(BCdata, dims = 1, cells = 500, balanced = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/heatmap-1.jpeg" width="2100" />

```r
DimHeatmap(BCdata, dims = 1:9, cells = 500, balanced = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/heatmap-2.jpeg" width="2100" />

## JackStrawPlot & ElbowPlot


```r
# BCdata <- JackStraw(BCdata, num.replicate = 100)
# BCdata <- ScoreJackStraw(BCdata, dims = 1:20)
JackStrawPlot(BCdata, dims = 1:15)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/JackStrawPlot & ElbowPlot-1.jpeg" width="2100" />

```r
ElbowPlot(BCdata)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/JackStrawPlot & ElbowPlot-2.jpeg" width="2100" />

## Umap & tsne


```r
# BCdata <- FindNeighbors(BCdata, dims = 1:10)
# BCdata <- FindClusters(BCdata, resolution = 0.4)
# BCdata <- RunUMAP(BCdata, dims = 1:10)
# BCdata <- RunTSNE(BCdata, dims = 1:10)
# Look at cluster IDs of the first 5 cells
head(Idents(BCdata), 5)                     
```

```
## CID3586_AAGACCTCAGCATGAG CID3586_AAGGTTCGTAGTACCT CID3586_ACCAGTAGTTGTGGCC 
##            Breast cancer            Breast cancer            Breast cancer 
## CID3586_ACCCACTAGATGTCGG CID3586_ACTTGTTAGGGAAACA 
##            Breast cancer            Breast cancer 
## Levels: Breast cancer
```

```r
DimPlot(BCdata, reduction = "umap",group.by = "celltype_major" ,label = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-1-1.jpeg" width="2100" />

```r
DimPlot(BCdata, reduction = "tsne",group.by = "celltype_major" ,label = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-1-2.jpeg" width="2100" />

```r
#saveRDS(BCdata, file = "C:/Users/user/Desktop/BCsinglecelldata/BCdata.rds")
```

## heat map markers for each cell types major


```r
DimPlot(BCdata, reduction = "umap",group.by ="celltype_major" , label = TRUE, pt.size = 0.5)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/heat map markers cell type major-1.jpeg" width="2100" />

```r
BCdata.markers <- readRDS('BCdata.markers.rds')

# Idents(BCdata) <- "celltype_major"
# BCdata.markers <- FindAllMarkers(BCdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BCdata.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(BCdata, features = top10$gene) + NoLegend()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/heat map markers cell type major-2.jpeg" width="2100" />

```r
# saveRDS(BCdata.markers, file = "BCdata.markers.rds")
```

## violin plot PALLD
               
## Idnets normal_cell_call


```r
Idents(BCdata) <- 'normal_cell_call'
DimPlot(BCdata, reduction = "umap" ,label = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.jpeg" width="2100" />

```r
DimPlot(BCdata, reduction = "tsne" ,label = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-2.jpeg" width="2100" />

```r
VlnPlot(BCdata, features = c("PALLD"), sort = TRUE) + NoLegend()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-3.jpeg" width="2100" />

## cancer subset


```r
BCdata_cancer <-  readRDS('BCdata_subset_cancer.rds')
Idents(BCdata_cancer) <- "celltype_minor"
```

## umap & tsne cancer cells


```r
DimPlot(BCdata_cancer, reduction = "umap",label = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/umap & tsne cancer-1.jpeg" width="2100" />

```r
DimPlot(BCdata_cancer, reduction = "tsne" ,label = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/umap & tsne cancer-2.jpeg" width="2100" />


```r
VlnPlot(BCdata_cancer, features = c("PALLD"), sort = TRUE) + NoLegend()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/violin cancer-1.jpeg" width="2100" />

## find all marker by cancer idents


```r
Idents(BCdata) <- 'normal_cell_call'
# BCdata_markers_normal_cell_call <- FindAllMarkers(BCdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BCdata_markers_normal_cell_call <- readRDS('BCdata_markers_normal_cell_call.rds')
BCdata_markers_normal_cell_call %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(BCdata, features = top10$gene) + NoLegend()
```

<img src="{{< blogdown/postref >}}index_files/figure-html/markers by cancer-1.jpeg" width="2100" />

```r
#saveRDS(BCdata_markers_normal_cell_call, file = "BCdata_markers_normal_cell_call.rds")
```


