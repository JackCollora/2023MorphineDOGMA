library(Seurat)
library(dplyr)
library(Matrix)
setwd("~/gibbs/DOGMAMORPH/Ranalysis")
options(warn = 0)

#quick function to load the cite data and add it as ADT/HTO
addCITEAssays<-function(citedir, seuratobj, assaynames=c("ADT","HTO")){
  data1<-readMM(paste0(citedir,"matrix.mtx"))
  cols<-read.csv(paste0(citedir,"barcodes.tsv"), header = FALSE)
  rows<-read.csv(paste0(citedir,"features.tsv"), header = FALSE)
  data1<-as.data.frame(t(data1))
  colnames(data1)<-paste(cols$V1, 1, sep="-") 
  rownames(data1)<-rows$V1
  
  #filtering had to be added to remove those with no cite data and put objects in agreement 
  data1<-data1[,colnames(data1)%in%colnames(seuratobj)]
  seuratobj$keep<-colnames(seuratobj)%in% colnames(data1)
  seuratobj<-subset(seuratobj,keep==TRUE)
  seuratobj[[assaynames[1]]]<-CreateAssayObject(data1[!grepl("hash", rownames(data1)),])
  #filter for hashs in the sample
  seuratobj[[assaynames[2]]]<-CreateAssayObject(data1[grepl("hash", rownames(data1))&(rowMeans(data1)>1),])
  return(seuratobj)
}


todo<-list.files("../data/")

#loop is just using the cite to define two groups
#the things that are not called doublets (includes the negative)
#things that are celled singlets
results<-list()
for(i in todo){
  curdata<-paste0( "../data/", i,"/")
  data<-Read10X(paste0(curdata, "RNAmatrix/"))
  data<-CreateSeuratObject(data$`Gene Expression`, project = i)
  data<-addCITEAssays(paste0(curdata,"CITEmatrix/"), data)
  data<-HTODemux(data)
  table(data$HTO_classification.global)
  write.table(data$HTO_classification[data$HTO_classification.global!="Doublet"], paste0(curdata,"notdoublet_cells.txt" ),quote = FALSE, row.names = TRUE, col.names = FALSE)
  write.table(data$HTO_classification[data$HTO_classification.global=="Singlet"], paste0(curdata,"singlet_cells.txt" ),quote = FALSE, row.names = TRUE, col.names = FALSE)
  results[[i]]<-data
  }
saveRDS(results, "Objects/UnfilteredSeuratObjs.rds")

