library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(SeuratWrappers)
library(harmony)
library(dplyr)
setwd("~/gibbs/DOGMAMORPH/Ranalysis")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
results<-readRDS( "~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230518Completeobj.rds")

#cluster annotations from ADT 

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

#grabbing just CD4 and spliting things up by original identity 
t_cells<-c("Memory_CD4_Polarized_1","Memory_CD4_Polarized_2","Treg_proliferating","Cytotoxic_T","Memory_CD4_1","Memory_CD4_2",
           "Naive_CD4_T_1","Naive_CD4_T_2","Naive_CD4_T_3")
DefaultAssay(results)<-"RNA"
results_sub<-subset(results, idents=t_cells)
results_sub<-SplitObject(results_sub, "orig.ident")
results_sub<-lapply(results_sub, NormalizeData)
results_sub<-lapply(results_sub, FindVariableFeatures)
results_sub<-RunFastMNN(results_sub, features=SelectIntegrationFeatures(results_sub))


DefaultAssay(results_sub)<-"ADT.dsb"
results_sub<-ScaleData(results_sub)
results_sub<-FindVariableFeatures(results_sub)
results_sub<-RunPCA(results_sub)
results_sub<-RunHarmony(results_sub,group.by.vars = "orig.ident", reduction.save="adt_harmony")

####ATACint####
DefaultAssay(results_sub)<-"ATAC"

# basic preprocessing
results_sub <- FindTopFeatures(results_sub, min.cutoff = 20)
results_sub <- RunTFIDF(results_sub)
results_sub <- RunSVD(results_sub)

results_sub<-RunHarmony(results_sub,group.by.vars = "orig.ident", reduction.save="atac_harmony")


#now we can do WNN

results_sub <- FindMultiModalNeighbors(
  results_sub, reduction.list = list("mnn", "adt_harmony", "atac_harmony"), 
  dims.list = list(1:30, 1:18, 2:30), modality.weight.name = "RNA.weight"
)

results_sub <- RunUMAP(results_sub, nn.name="weighted.nn")
results_sub <- FindClusters(results_sub, graph.name = "wsnn", algorithm= 3, resolution = 0.8, verbose = FALSE)


DimPlot(results_sub, label=TRUE)+NoLegend()
DimPlot(results_sub, label=TRUE, group.by = "orig.ident")+NoLegend()
DimPlot(results_sub, label=TRUE, group.by = "merged_clusters")+NoLegend()

saveRDS(results_sub, "Objects/20230526CD4Obj.rds")

CD4_marks_Collora_RNA<-c("GZMB", "CCL5", "TBX21", "GZMK","GATA3","RORC", "CTSH","FOXP3", "IL2RA","MKI67", "CXCR5", "CCR7", "SELL", "TCF7", "mt")

DefaultAssay(results_sub)<-"RNA"

FeaturePlot(results_sub, CD4_marks_Collora_RNA)
DefaultAssay(results_sub)<-"ADT.dsb"

FeaturePlot(results_sub, features= c("CD3-TotalA","CD4-TotalA", "CD8-TotalA","CD45RA-TotalA","CD45RO-TotalA", "CD25-TotalA","CD69-TotalA","HLA-DR-DP-DQ-TotalA","CD185-TotalA",
                                                               "CD279-TotalA","TIGIT-TotalA","CD62L-TotalA","CD152-TotalA","KLRG1-TotalA", "TCRVa7.2-TotalA","TCRVd2-TotalA", "TCRab-TotalA", 
                                                               "CX3CR1-TotalA","CD194-TotalA","CD196-TotalA" , "CD278-TotalA"), min.cutoff = 'q5', max.cutoff = 'q95', reduction="umap")
