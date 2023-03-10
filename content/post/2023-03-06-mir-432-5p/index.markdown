---
title: 'miR-432-5p'
author: 'Aviv leemann'
date: '2023-03-06'
slug: 'mir-432-5p'
tags: 
  -'mir-432-5p'
subtitle: 'single cell analysis hsa-mir432-5p'
summary: 'enrichment analysis purpose is to check if there are pathways, molecular function, biological process & cellular components that is connected or related to this bio-marker target.'
lastmod: '2023-03-06T20:06:52+02:00'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: ['enrichment','mir-432-5p']
---




## Enrichmet analysis: hsa-miR-432-5p targets

enrichment analysis purpose is to check if there are:
pathways, molecular function, biological process & cellular components that is connected or related to this bio-marker target.

## libreris

pacman library purpose is to load multiple libraries from a vector


```r
library(pacman)
lab <-  c("GO.db","GOstats","cRegulome","readxl","ggplot2","tidyverse","R.utils","MeSHDbi","AnnotationHub","KEGGREST","AnnotationDbi","org.Hs.eg.db","clusterProfiler","meshr","pathview","ReactomePA","meshes")
pacman::p_load(lab,install = FALSE,character.only = TRUE)
```

## files

loading miR targets files acquired from mirWalk


```r
mir432_5p_targets <-read.csv("miRWalk_miRNA_Targets_432_5p_3UTR.csv")
# annotations_ahb <- read.csv("annotations_ahb.csv")
```

## head


```r
mir432_5p_targets %>% head()
```

```
##          mirnaid     refseqid genesymbol
## 1 hsa-miR-432-5p NM_001353214        DYM
## 2 hsa-miR-432-5p NM_001374442        DYM
## 3 hsa-miR-432-5p NM_001291932    SEC14L2
## 4 hsa-miR-432-5p NM_001301068       ARL1
## 5 hsa-miR-432-5p    NM_003387      WIPF1
## 6 hsa-miR-432-5p    NM_012215        OGA
##                                                                                    duplex
## 1 TCTTGGAGTAGGTCATTGGGTGG#CCACCCGGACCCCTCCAGC#..((((((..((((..(((((((#))))))))))).)))))).
## 2 TCTTGGAGTAGGTCATTGGGTGG#CCACCCGGACCCCTCCAGC#..((((((..((((..(((((((#))))))))))).)))))).
## 3 TCTTGGAGTAGGTCATTGGGTGG#ACCCAGGATTGTCTCCAAC#..((((((..((((.((((((..#))))))))))..)))))).
## 4 TCTTGGAGTAGGTCATTGGGTGG#CACCCATTCCTACTTTAGT#..((((((((((....((((((.#))))))..)))))))))).
## 5 TCTTGGAGTAGGTCATTGGGTGG#GCCTAATGACCCATTTCAC#...((((((.(((((((((((..#))))))))))).)))))).
## 6   TCTTGGAGTAGGTCATTGGGTGG#TACCTAGTGACTTCCAAC#..(((((...((((((((((((.#))))))))))))))))).
##   start  end  bindingp energy seed accessibility    au phylopstem phylopflank
## 1  2491 2510 0.9500000  -24.3    1   0.005542395 0.441  0.0000000   0.0000000
## 2  1756 1775 0.9500000  -24.3    1   0.005542395 0.441  0.0000000   0.0000000
## 3  3555 3574 0.9529915  -21.4    1   0.002860247 0.456  0.0000000   0.0000000
## 4  2116 2135 0.9538462  -22.2    1   0.232450300 0.603  1.2249763   0.8797262
## 5  3702 3721 0.9538462  -22.9    0   0.053886420 0.647  2.1715973   2.1847662
## 6  4641 4659 0.9538462  -21.3    0   0.048461950 0.588  0.9082432   1.0060216
##          me number_of_pairings binding_region_length
## 1 -8.224395                 17                    19
## 2 -8.224395                 17                    19
## 3 -7.991346                 16                    19
## 4 -8.071470                 16                    19
## 5 -7.757047                 17                    19
## 6 -9.569532                 17                    18
##   longest_consecutive_pairings position validated TargetScan miRDB
## 1                           11     3UTR                    0     0
## 2                           11     3UTR                    0     0
## 3                           10     3UTR                    0     0
## 4                           10     3UTR                    0     0
## 5                           11     3UTR                    0     0
## 6                           17     3UTR                    0     0
```

```r
mir432_5p_targets %>% colnames()
```

```
##  [1] "mirnaid"                      "refseqid"                    
##  [3] "genesymbol"                   "duplex"                      
##  [5] "start"                        "end"                         
##  [7] "bindingp"                     "energy"                      
##  [9] "seed"                         "accessibility"               
## [11] "au"                           "phylopstem"                  
## [13] "phylopflank"                  "me"                          
## [15] "number_of_pairings"           "binding_region_length"       
## [17] "longest_consecutive_pairings" "position"                    
## [19] "validated"                    "TargetScan"                  
## [21] "miRDB"
```

## preprocessing


```r
colnames(mir432_5p_targets)[3] <- "SYMBOL"
colnames(mir432_5p_targets)[1] <- "miR"
colnames(mir432_5p_targets)[7] <- "score"
```

maps the gene to ENTREZID


```r
mir432_5p_targets$ENTREZ = mapIds(org.Hs.eg.db,
                                  key = mir432_5p_targets$SYMBOL,  # data to use for retrieval
                                  column =  "ENTREZID", # information to retreive for given data
                                  keytype = "SYMBOL",multiVals =  "first" )
```


```r
# mir_532_5_mapIDs <- bitr( mir432_5p_targets$SYMBOL,  "SYMBOL","ENTREZID", OrgDb = org.Hs.eg.db)
```


## KEGG

run KEGG enrichment via cluster Profiler package


```r
clusterProfiler_KEGG_432_5p <- enrichKEGG(mir432_5p_targets$ENTREZ)
```





```r
head(summary(clusterProfiler_KEGG_432_5p)[, -8])
```

```
##                ID              Description GeneRatio  BgRatio       pvalue
## hsa05032 hsa05032       Morphine addiction   35/1487  91/8292 2.819029e-06
## hsa04727 hsa04727        GABAergic synapse   31/1487  89/8292 9.720160e-05
## hsa04728 hsa04728     Dopaminergic synapse   41/1487 132/8292 1.570806e-04
## hsa04713 hsa04713    Circadian entrainment   32/1487  97/8292 2.406306e-04
## hsa04014 hsa04014    Ras signaling pathway   64/1487 236/8292 2.570959e-04
## hsa04371 hsa04371 Apelin signaling pathway   42/1487 139/8292 2.572457e-04
##              p.adjust       qvalue Count
## hsa05032 0.0009387366 0.0007744911    35
## hsa04727 0.0142771350 0.0117791441    31
## hsa04728 0.0142771350 0.0117791441    41
## hsa04713 0.0142771350 0.0117791441    32
## hsa04014 0.0142771350 0.0117791441    64
## hsa04371 0.0142771350 0.0117791441    42
```

```r
dotplot(clusterProfiler_KEGG_432_5p, showCategory=10) + labs(title ="hsa-miR-432-5p targets KEGG pathways")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/KEGG enrichment 432_5-1.jpeg" width="2100" />

## GO analysis

### BP analysis

no results


```r
clusterProfiler_bp_432_5p <- enrichGO(mir432_5p_targets$ENTREZ, ont="BP",OrgDb =org.Hs.eg.db)
#summary(clusterProfiler_bp_432_5p) %>% head()
#dotplot(clusterProfiler_bp_432_5p, showCategory=20)
```

no results

### CC analysis


```r
clusterProfiler_cc_432_5p <- enrichGO(mir432_5p_targets$ENTREZ, ont="CC",OrgDb =org.Hs.eg.db)
head(summary(clusterProfiler_cc_432_5p)[,-8])
```

```
##                    ID                           Description GeneRatio   BgRatio
## GO:0005769 GO:0005769                        early endosome  104/3317 414/19869
## GO:0031965 GO:0031965                      nuclear membrane   79/3317 305/19869
## GO:0005635 GO:0005635                      nuclear envelope  116/3317 488/19869
## GO:0031228 GO:0031228 intrinsic component of Golgi membrane   25/3317  69/19869
## GO:0000151 GO:0000151              ubiquitin ligase complex   78/3317 309/19869
## GO:0043197 GO:0043197                       dendritic spine   49/3317 175/19869
##                  pvalue    p.adjust      qvalue Count
## GO:0005769 6.841258e-06 0.005199356 0.004500828   104
## GO:0031965 2.694279e-05 0.007917530 0.006853817    79
## GO:0005635 3.125341e-05 0.007917530 0.006853817   116
## GO:0031228 7.110974e-05 0.011447826 0.009909821    25
## GO:0000151 7.705279e-05 0.011447826 0.009909821    78
## GO:0043197 1.164215e-04 0.011447826 0.009909821    49
```

```r
dotplot(clusterProfiler_cc_432_5p, showCategory=10) + labs(title ="hsa-miR-432-5p targets CC GO")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/CC 432_5-1.jpeg" width="2100" />


```r
goplot(clusterProfiler_cc_432_5p )
```

<img src="{{< blogdown/postref >}}index_files/figure-html/goplot cc-1.jpeg" width="2100" />

### MF analysis


```r
clusterProfiler_mf_432_5p <- enrichGO(mir432_5p_targets$ENTREZ, ont="MF",OrgDb =org.Hs.eg.db)
head(summary(clusterProfiler_mf_432_5p)[,-8] )
```

```
##                    ID                   Description GeneRatio   BgRatio
## GO:0019001 GO:0019001     guanyl nucleotide binding  111/3240 401/18432
## GO:0032561 GO:0032561 guanyl ribonucleotide binding  111/3240 401/18432
## GO:0003924 GO:0003924               GTPase activity   94/3240 331/18432
## GO:0005525 GO:0005525                   GTP binding  104/3240 379/18432
## GO:0019003 GO:0019003                   GDP binding   30/3240  73/18432
## GO:0003925 GO:0003925            G protein activity   18/3240  41/18432
##                  pvalue     p.adjust       qvalue Count
## GO:0019001 2.696462e-07 0.0001620574 0.0001558271   111
## GO:0032561 2.696462e-07 0.0001620574 0.0001558271   111
## GO:0003924 6.380072e-07 0.0002556282 0.0002458007    94
## GO:0005525 9.916314e-07 0.0002979852 0.0002865293   104
## GO:0019003 2.039762e-06 0.0004903587 0.0004715070    30
## GO:0003925 7.950900e-05 0.0159283036 0.0153159448    18
```

```r
dotplot(clusterProfiler_mf_432_5p, showCategory=10) + labs(title ="hsa-miR-432-5p targets MF GO")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/MF 432_5-1.jpeg" width="2100" />


```r
goplot(clusterProfiler_mf_432_5p )
```

<img src="{{< blogdown/postref >}}index_files/figure-html/goplot mf-1.jpeg" width="2100" />

## Reactome Pathways


```r
reactome_432_5<- enrichPathway(gene=mir432_5p_targets$ENTREZ, pvalueCutoff = 0.05, readable=TRUE)
dotplot(reactome_432_5,showCategory=10) + labs(title ="hsa-miR-432-5p targets reactome pathways")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/reactome pathways-1.jpeg" width="2100" />

## MESH


```r
# already created 
# creating meshDBi
# ah <- AnnotationHub(localHub=TRUE)
# hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
# file_hsa <- hsa[[1]]
# db <- MeSHDbi::MeSHDb(file_hsa)
```


```r
# saved & read the mehses 
# meshes432 <- enrichMeSH(mir432_5p_targets$ENTREZ, MeSHDb = db, database='gendoo', category = 'C')
# saveRDS(meshes432,'meshes432.rds')
meshes432 <- readRDS('meshes432.rds')
dotplot(meshes432) + labs(title ="hsa-miR-432-5p targets mesh gendoo Diseases")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/mesh gendoo-1.jpeg" width="2100" />


```r
# meshes432_pubmed <- enrichMeSH(mir432_5p_targets$ENTREZ, MeSHDb = db, database='gene2pubmed', category = 'C')
# dotplot(meshes432_pubmed) + labs(title ="hsa-miR-432-5p targets mesh gene2pubmed Diseases")
```

# session info


```r
sessionInfo()[6]
```

```
## $basePkgs
## [1] "stats4"    "stats"     "graphics"  "grDevices" "utils"     "datasets" 
## [7] "methods"   "base"
```

