library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(Signac)
options(warn=0)
setwd("~/gibbs/DOGMAMORPH/Ranalysis")

results<-readRDS( "Objects/20230518Completeobj.rds")

#basically we're going to compare each cluster to each other cluster and then calculate the DE overlap between each pair. 
#Similar clusters should have similar DE overlaps allowing us to merge them based on their clustering

DefaultAssay(results)<-"RNA"

#first going to try this an "easy way" using find all markers, then I can consider expanding it if needed
#17-24 (less 22) have too few cells to test, so those are auto merging with the cluster they most closely overlap with
table(results$seurat_clusters)
results$newclust<-case_when(results$seurat_clusters==17~8, results$seurat_clusters==18~10,results$seurat_clusters==19~3,
                            results$seurat_clusters==20~4,results$seurat_clusters==21~6,results$seurat_clusters==23~2,
                            results$seurat_clusters==24~6, results$seurat_clusters==0~0,T~as.numeric(results$seurat_clusters))
Idents(results)<-results$newclust
DEresult<-FindAllMarkers(results)

#now we split it and calc overlaps 
DEresults<-split(DEresult, DEresult$cluster)

Overlapvector<-c()

for(i in names(DEresults)){
  for(j in names(DEresults)){
    Overlapvector<-c(Overlapvector, length(intersect(DEresults[[i]]$gene, DEresults[[j]]$gene)))
  }
}

Overlapvector<-matrix(Overlapvector,dimnames = list(names(DEresults), names(DEresults)), nrow = length(DEresults) )

heatmap(Overlapvector, scale = "row")

#it didnt end up being super useful since small clusters are very noisy in the overall dataset, instead trying it with all the different datapoints 
DefaultAssay(results)<-"RNA"

todo<-unique(Idents(results))
DEresult<-list()
for (i in todo){
  cur<-list()
  for (j in todo){
    if(i==j){
      cur[[j]]<-1
      next
    }
    cur[[j]]<-FindMarkers(results, i,j)
  }
  DEresult[[i]]<-cur
}

saveRDS(DEresult, "Objects/20230522pairwiseDE.rds")
