library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)
library(SeuratWrappers)
library(harmony)

setwd("~/gibbs/DOGMAMORPH/Ranalysis")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

#importing the ARCHR peakset here, going to call a separate MACS2 one after on all to see if there's something to be gained from it

readPeaks<-function(peakdir){
  data1<-readMM(paste0(peakdir,"/counts.mtx"))
  cols<-read.csv(paste0(peakdir,"/cells.csv"), header = TRUE)
  rows<-read.csv(paste0(peakdir,"/peaks.csv"), header = TRUE)

  rows<-paste(rows$seqnames, rows$start, rows$end, sep = "-")
  dimnames(data1)<-list(rows,cols$x)
  return(data1)
}

data<-readPeaks("DOGMA_ATAC/export/peak_counts/")
results<-readRDS("Objects/dsbSeuratObjs_binary.rds")


for (i in names(results)){
  data1<-data[,grep(paste0("^",i), colnames(data))]
  colnames(data1)<-gsub(i,"", colnames(data1))
  colnames(data1)<-gsub("#","", colnames(data1))
  
  #keeping only those in both assays
  data1<-data1[,colnames(data1)%in%colnames(results[[i]])]
  results[[i]]$keep<-colnames(results[[i]])%in% colnames(data1)
  results[[i]]<-subset(results[[i]],keep==TRUE)
  
  results[[i]][["ATAC"]] <- CreateChromatinAssay(
  counts = data1,
  sep = c("-", "-"),
  fragments = paste0("../data/",i, "/atac_fragments.tsv.gz"),
  annotation = annotation
)
  results[[i]]$mt<-PercentageFeatureSet(results[[i]], "^MT-")
  results[[i]]<-subset(results[[i]], subset= HTO_classification.global =="Singlet" & nFeature_RNA > 200 & mt<20)
  DefaultAssay(results[[i]])<-"RNA"
  results[[i]]<-RenameCells(results[[i]], new.names=paste(i, colnames(results[[i]]), sep="_"))
}

saveRDS(results,"Objects/trimodalcheckpoint1.")

#normalizing and integrating using fastmnn
results<-lapply(results, NormalizeData)
results<-lapply(results, FindVariableFeatures)
results<-RunFastMNN(results, features=SelectIntegrationFeatures(results))

#now the ADT 
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
results <- FindClusters(results, graph.name = "wsnn", algorithm= 3, resolution = 0.8, verbose = FALSE)

DimPlot(results, label=TRUE)+NoLegend()
DimPlot(results, label=TRUE, group.by = "orig.ident")+NoLegend()
DefaultAssay(results)<-"ADT.dsb"
FeaturePlot(results, c("CD4-TotalA","CD56-TotalA","CD19-TotalA","CD14-TotalA", "CD8-TotalA"), min.cutoff = 'q5',
            max.cutoff='q95', reduction = "umap")



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

saveRDS(results, "Objects/20230518Completeobj.rds")
