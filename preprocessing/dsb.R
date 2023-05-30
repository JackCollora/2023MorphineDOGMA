library(Seurat)
library(dplyr)
library(Matrix)
library(dsb)
setwd("~/gibbs/DOGMAMORPH/Ranalysis")

readCITE<-function(citedir){
  data1<-readMM(paste0(citedir,"/matrix.mtx"))
  cols<-read.csv(paste0(citedir,"/barcodes.tsv"), header = FALSE)
  rows<-read.csv(paste0(citedir,"/features.tsv"), header = FALSE)
  data1<-as.data.frame(t(data1))
  colnames(data1)<-paste(cols$V1, 1, sep="-") 
  rownames(data1)<-rows$V1
  data1<-data1[!grepl("hash", rownames(data1)),]
  rownames(data1)<-gsub("_","-", rownames(data1))
  return(data1)
}

results<-readRDS("Objects/UnfilteredSeuratObjs.rds")
todo<-list.files("../data/")
raws=list()
for (i in todo){
  raws[[i]]<-readCITE(paste0("../data/",i,"/CITEmatrix/"))
  raws[[i]]<-raws[[i]][,setdiff(colnames(raws[[i]]), colnames(results[[i]]))]
  raws[[i]]<-raws[[i]][,colSums(raws[[i]])>0]
  isotypes<-grep("Mouse|Rat|Hamster",rownames(raws[[i]]), value = TRUE)
  raws[[i]]<-DSBNormalizeProtein(
    cell_protein_matrix = results[[i]]@assays$ADT@counts, 
    empty_drop_matrix = raws[[i]], 
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    isotype.control.name.vec = isotypes
  )
  results[[i]][['ADT.dsb']]<-CreateAssayObject(counts=raws[[i]])
}
saveRDS(results, "Objects/dsbSeuratObjs.rds")

#now doing the Z-scoring for each antibody for each


#biggest deal is making this translation matrix, to map each to their corresponding isotype, had to be made manually. 

key<-read.csv("metadata/totalA_iso.csv", sep=",")
key<-key[1:163,1:2]
#fixing the rats which lost their extra letter due to encoding issues 
key[c(85,117),2]<-"RatIgG1K-TotalA"
key[c(38,86),2]<-"RatIgG1Y-TotalA"
#fixing rest which are all kappa 
colnames(key)<-c("ab","iso")
key$iso<-gsub("Mouse IgG1","Mouse-IgG1-isotype-TotalA",key$iso)
key$iso<-gsub("Mouse IgG2a","Mouse-IgG2a-isotype-TotalA",key$iso)
key$iso<-gsub("Mouse IgG2b","Mouse-IgG2b-isotype-TotalA",key$iso)
key$iso<-gsub("Rat IgG2b","Rat-IgG2b-isotype-TotalA",key$iso)
key$iso<-gsub("Rat IgG2a","RatIgG2a-TotalA",key$iso)
key$iso<-gsub("Rat IgG2c","RatIgG2c-TotalA",key$iso)
key$iso<-gsub("Armenian Hamster IgG","ArmenianHamsterIgG-TotalA",key$iso)
table(key$iso%in%isotypes)
#now for each of these, we'll calculate a cutoff by taking the mean and adding 2Sds to it, 
#then go through and map each into a true/false 
#assay will basically be good only for fisher testing between groups but that's fine I guess

results<-readRDS("Objects/dsbSeuratObjs.rds")

for (i in todo){
  cutoffs<-list()
  for (j in isotypes){
    cutoffs[[j]]<-mean(results[[i]]@assays$ADT.dsb[j,])+2*sd(results[[i]]@assays$ADT.dsb[j,])
  }
  cutoffs<-t(as.data.frame(cutoffs))
  colnames(cutoffs)<-"cuts"
  cutoffs<-as.data.frame(cutoffs)
  cutoffs$iso<-gsub("\\.","-",rownames(cutoffs))
  
  cutoffs<-merge(key, cutoffs, all=TRUE,by.X="iso", by.Y="iso" )
  results[[i]][['ADT.dsb.bin']]<-CreateAssayObject(results[[i]]@assays$ADT.dsb@data>cutoffs$cuts)
  
}
saveRDS(results, "Objects/dsbSeuratObjs_binary.rds")
