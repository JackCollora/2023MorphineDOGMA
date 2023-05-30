library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(SeuratWrappers)
library(harmony)
library(dplyr)
setwd("~/gibbs/DOGMAMORPH/Ranalysis")

results<-readRDS( "~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230518Completeobj.rds")
results$merged_clusters<-case_when(results$seurat_clusters==16~8,results$seurat_clusters==17~8, results$seurat_clusters==18~10,results$seurat_clusters==19~3,
                                   results$seurat_clusters==20~4,results$seurat_clusters==21~6,results$seurat_clusters==22~2,results$seurat_clusters==23~2,
                                   results$seurat_clusters==24~6, results$seurat_clusters==0~0,T~(as.numeric(results$seurat_clusters)-1))

results$merged_clusters<-case_when(results$merged_clusters==0~"Naive_CD4_T_1",results$merged_clusters==1~"Memory_CD4_Polarized_1",results$merged_clusters==2~"Memory_CD4_Polarized_2",results$merged_clusters==3~"Memory_CD4_1",
                                   results$merged_clusters==4~"Cytotoxic_T",results$merged_clusters==5~"Naive_B",results$merged_clusters==6~"Naive_CD4_T_2",results$merged_clusters==7~"Treg_proliferating",
                                   results$merged_clusters==8~"Naive_CD4_T_3",results$merged_clusters==9~"NK",results$merged_clusters==10~"Memory_B",results$merged_clusters==11~"CD14+_Mono",
                                   results$merged_clusters==12~"cDC",results$merged_clusters==13~"pDC",results$merged_clusters==14~"Memory_CD4_2",results$merged_clusters==15~"Plasma")

results$merged_clusters<-factor(results$merged_clusters, levels = c("Memory_CD4_Polarized_1","Memory_CD4_Polarized_2","Treg_proliferating","Cytotoxic_T","NK","Memory_CD4_1","Memory_CD4_2",
                                                                    "Naive_CD4_T_1","Naive_CD4_T_2","Naive_CD4_T_3","cDC","pDC","CD14+_Mono","Naive_B","Memory_B","Plasma"))

Idents(results)<-results$merged_clusters
DefaultAssay(results)<-"RNA"
for (i in c(500, 1000, 1500, 2000, 3500, 5000)){
results<-SplitObject(results, "orig.ident")
results<-lapply(results, NormalizeData)
results<-lapply(results, FindVariableFeatures, nfeatures = i)
results<-RunFastMNN(results, features=SelectIntegrationFeatures(results, nfeatures = i))


DefaultAssay(results)<-"ADT.dsb"
results<-ScaleData(results)
results<-FindVariableFeatures(results)
results<-RunPCA(results)
results<-RunHarmony(results,group.by.vars = "orig.ident", reduction.save="adt_harmony")

####ATACint####
DefaultAssay(results)<-"ATAC"

# basic preprocessing
results <- FindTopFeatures(results, min.cutoff = 20)
results <- RunTFIDF(results)
results <- RunSVD(results)

results<-RunHarmony(results,group.by.vars = "orig.ident", reduction.save="atac_harmony")


#now we can do WNN

results <- FindMultiModalNeighbors(
  results, reduction.list = list("mnn", "adt_harmony", "atac_harmony"), 
  dims.list = list(1:30, 1:18, 2:30), modality.weight.name = "RNA.weight"
)

results <- RunUMAP(results, nn.name="weighted.nn")
pdf(paste0("Parameterscan/","nfeature_",i,".pdf"))
print(DimPlot(results, label=TRUE)+NoLegend())
print(DimPlot(results, label=TRUE, group.by = "orig.ident")+NoLegend())
print(DimPlot(results, label=TRUE, group.by = "merged_clusters")+NoLegend())
dev.off()
}

for (i in c(10,20,30,40,50)){
  results<-SplitObject(results, "orig.ident")
  results<-lapply(results, NormalizeData)
  results<-lapply(results, FindVariableFeatures)
  results<-RunFastMNN(results, features=SelectIntegrationFeatures(results))
  
  
  DefaultAssay(results)<-"ADT.dsb"
  results<-ScaleData(results)
  results<-FindVariableFeatures(results)
  results<-RunPCA(results, npcs = i)
  results<-RunHarmony(results,group.by.vars = "orig.ident", reduction.save="adt_harmony")
  
  ####ATACint####
  DefaultAssay(results)<-"ATAC"
  
  # basic preprocessing
  results <- FindTopFeatures(results, min.cutoff = 20)
  results <- RunTFIDF(results)
  results <- RunSVD(results,n = i)
  
  results<-RunHarmony(results,group.by.vars = "orig.ident", reduction.save="atac_harmony")
  
  
  #now we can do WNN
  
  results <- FindMultiModalNeighbors(
    results, reduction.list = list("mnn", "adt_harmony", "atac_harmony"), 
    dims.list = list(1:i, 1:i, 2:i), modality.weight.name = "RNA.weight"
  )
  
  results <- RunUMAP(results, nn.name="weighted.nn")
  pdf(paste0("Parameterscan/","nPCsint_",i,".pdf"))
  print(DimPlot(results, label=TRUE)+NoLegend())
  print(DimPlot(results, label=TRUE, group.by = "orig.ident")+NoLegend())
  print(DimPlot(results, label=TRUE, group.by = "merged_clusters")+NoLegend())
  dev.off()
}

for (i in c(5,10,20,30,50, 100)){
  results<-SplitObject(results, "orig.ident")
  results<-lapply(results, NormalizeData)
  results<-lapply(results, FindVariableFeatures)
  results<-RunFastMNN(results, features=SelectIntegrationFeatures(results))
  
  
  DefaultAssay(results)<-"ADT.dsb"
  results<-ScaleData(results)
  results<-FindVariableFeatures(results)
  results<-RunPCA(results)
  results<-RunHarmony(results,group.by.vars = "orig.ident", reduction.save="adt_harmony")
  
  ####ATACint####
  DefaultAssay(results)<-"ATAC"
  
  # basic preprocessing
  results <- FindTopFeatures(results, min.cutoff = 20)
  results <- RunTFIDF(results)
  results <- RunSVD(results)
  
  results<-RunHarmony(results,group.by.vars = "orig.ident", reduction.save="atac_harmony")
  
  
  #now we can do WNN
  
  results <- FindMultiModalNeighbors(
    results, reduction.list = list("mnn", "adt_harmony", "atac_harmony"), 
    dims.list = list(1:30, 1:18, 2:30), modality.weight.name = "RNA.weight", k.nn=i
  )
  
  results <- RunUMAP(results, nn.name="weighted.nn", n.neighbors = i)
  pdf(paste0("Parameterscan/","nNeighbors_",i,".pdf"))
  print(DimPlot(results, label=TRUE)+NoLegend())
  print(DimPlot(results, label=TRUE, group.by = "orig.ident")+NoLegend())
  print(DimPlot(results, label=TRUE, group.by = "merged_clusters")+NoLegend())
  dev.off()
}

for (i in c(5,10,20,30,50, 100)){
  results<-SplitObject(results, "orig.ident")
  results<-lapply(results, NormalizeData)
  results<-lapply(results, FindVariableFeatures)
  results<-RunFastMNN(results, features=SelectIntegrationFeatures(results))
  
  
  DefaultAssay(results)<-"ADT.dsb"
  results<-ScaleData(results)
  results<-FindVariableFeatures(results)
  results<-RunPCA(results)
  results<-RunHarmony(results,group.by.vars = "orig.ident", reduction.save="adt_harmony")
  
  ####ATACint####
  DefaultAssay(results)<-"ATAC"
  
  # basic preprocessing
  results <- FindTopFeatures(results, min.cutoff = 20)
  results <- RunTFIDF(results)
  results <- RunSVD(results)
  
  results<-RunHarmony(results,group.by.vars = "orig.ident", reduction.save="atac_harmony")
  
  
  #now we can do WNN
  
  results <- FindMultiModalNeighbors(
    results, reduction.list = list("mnn", "adt_harmony", "atac_harmony"), 
    dims.list = list(1:30, 1:18, 2:30), modality.weight.name = "RNA.weight"
  )
  
  results <- RunUMAP(results, nn.name="weighted.nn", local.connectivity = i)
  pdf(paste0("Parameterscan/","nConnectivity_",i,".pdf"))
  print(DimPlot(results, label=TRUE)+NoLegend())
  print(DimPlot(results, label=TRUE, group.by = "orig.ident")+NoLegend())
  print(DimPlot(results, label=TRUE, group.by = "merged_clusters")+NoLegend())
  dev.off()
}

