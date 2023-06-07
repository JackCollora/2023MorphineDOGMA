suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

#this is the object post chromvar. Should be the last "checkpoint" object that analysis is plotted from. Contains annotations and all assays. 
results<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230606FinalObj.rds")
DefaultAssay(results)<-"RNA"

genes<-read.csv("~/gibbs/DOGMAMORPH/refs/Shalek_Antiviral_genes.csv")
genes<-list(genes$X...GENE[genes$CLUSTER=="Id"], genes$X...GENE[genes$CLUSTER=="IIIc"], genes$X...GENE[genes$CLUSTER=="IIId"])

names(genes)<-c("Core_antiviral","Peaked_Inflam","Sustained_Inflam")

results<-AddModuleScore(results, features = genes,name = names(genes))

#doesnt appear these are different across conditions. 
FeaturePlot(results,"Core_antiviral1", min.cutoff='q5', max.cutoff = 'q95')
FeaturePlot(results,"Core_antiviral1", min.cutoff=0, max.cutoff = 0.2, split.by = "Treatment")

FeaturePlot(results,"Peaked_Inflam2", min.cutoff='q5', max.cutoff = 'q95')
FeaturePlot(results,"Peaked_Inflam2", min.cutoff=0, max.cutoff = 0.15, split.by = "Treatment")

FeaturePlot(results,"Sustained_Inflam3", min.cutoff='q5', max.cutoff = 'q95')
FeaturePlot(results,"Sustained_Inflam3", min.cutoff=0, max.cutoff = 0.2, split.by = "Treatment")
