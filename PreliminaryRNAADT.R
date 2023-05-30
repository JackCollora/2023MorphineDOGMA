library(Seurat)
library(dplyr)
library(Matrix)
library(SeuratWrappers)
library(harmony)
setwd("~/gibbs/DOGMAMORPH/Ranalysis")

results<-readRDS("Objects/dsbSeuratObjs_binary.rds")

#first taking each object with some very modest cutoffs, getting the singlets and then moving forward 
#this is a strict cutoff to get the good "RNA" cells. There are likely some cells that exist only in the ATAC which are also good enough to look at separately. 
for (i in names(results)){
  results[[i]]$mt<-PercentageFeatureSet(results[[i]], "^MT-")
  results[[i]]<-subset(results[[i]], subset= HTO_classification.global =="Singlet" & nFeature_RNA > 200 & mt<20)
  DefaultAssay(results[[i]])<-"RNA"
  results[[i]]<-RenameCells(results[[i]], new.names=paste(i, colnames(results[[i]]), sep="_"))
}


#normalizing and integrating using fastmnn
results<-lapply(results, NormalizeData)
results<-lapply(results, FindVariableFeatures)
results<-RunFastMNN(results, features=SelectIntegrationFeatures(results))

#now the ADT 
DefaultAssay(results)<-"ADT.dsb"
results<-ScaleData(results)
results<-FindVariableFeatures(results)
results<-RunPCA(results)
results<-RunHarmony(results,group.by.vars = "orig.ident")

#now we can do WNN

results <- FindMultiModalNeighbors(
  results, reduction.list = list("mnn", "harmony"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

results <- RunUMAP(results, nn.name="weighted.nn")
results <- FindClusters(results, graph.name = "wsnn", algorithm= 3, resolution = 0.8, verbose = FALSE)

DimPlot(results, label=TRUE)+NoLegend()
DimPlot(results, label=TRUE, group.by = "orig.ident")+NoLegend()
FeaturePlot(results, c("CD4-TotalA","CD56-TotalA","CD19-TotalA","CD14-TotalA", "CD8-TotalA"), min.cutoff = 'q5', max.cutoff='q95')

#throwing the metadata on to then do some testing 
hashkeymeta<-read.csv("metadata/Morphine.Meta.csv")
colnames(hashkeymeta)<-gsub("X...","", colnames(hashkeymeta))

results$Hashkey<-gsub("hashtag","", results$HTO_classification)
results$Hashkey<-gsub("-TotalA","", results$Hashkey)
meta<-results@meta.data
meta$cellid<-rownames(meta)
metamerge<-merge(meta, hashkeymeta, by.x=c("orig.ident", "Hashkey"),  by.y=c("Orig.Ident", "Hashkey"), all.x=TRUE)
rownames(metamerge)<-metamerge$cellid

results<-AddMetaData(results, metamerge)

table(results$Participant, results$Timepoint, results$Cell.Source)
table(results$Cell.Source)

saveRDS(results, "Objects/20220516RNAADT.rds")
