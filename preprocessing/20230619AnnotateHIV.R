library(dplyr)
library(tidyr)
library(Seurat)
library(Signac)
library(Biostrings)
library(ShortRead)

#need the reference barcodes for ATAC<->RNA conversion
bcs<-read.table(gzfile("/vast/palmer/apps/avx2/software/CellRanger-ARC/2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz"), header = FALSE)
bcs2<-read.table(gzfile("/vast/palmer/apps/avx2/software/CellRanger-ARC/2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz"), header = FALSE)
bcs$B<-bcs2$V1
colnames(bcs)<-c("RNA", "ATAC")

FixBCs<-function(querry, ref){
  #cellranger tests all 4 and uses the best one, so first we need to generate those 
  forward16<-sapply(querry, substr, 1,16)
  rev16<-sapply(querry, substr, 1,16)
  rev16<-DNAStringSet(rev16)
  rev16<-reverseComplement(rev16)
  rev16<-as.character(rev16)
  forward24<-sapply(querry, substr, 9,24)
  rev24<-sapply(querry, substr, 9,24)
  rev24<-DNAStringSet(rev24)
  rev24<-reverseComplement(rev24)
  rev24<-as.character(rev24)
  final<-list(forward16,rev16, forward24, rev24)
  final<-lapply(final, unname)
  #now we'll test the individual categories to go
  res<-list()
  for (i in 1:4){res[[i]]<-final[[i]][final[[i]] %in% ref$ATAC]}
  best<-max(unlist(lapply(res, length)))
  #return the one with the best subset
  best<-final[unlist(lapply(res, length))==best]
  best<-data.frame(best)
  colnames(best)<-"best"
  best$querry<-querry
  #merge with the ref
  best<-merge(best, ref, by.x="best", by.y="ATAC", all=TRUE)
  #subset to just those we care for
  best<-best[!is.na(best$querry),]
  return(best)
}

#function for RNA annotation
#just reads in the file and then does some fangeling to make it "correct" 
#filtering is on by default and requires at least 2 reads to support 
HIVRNAAnnotation<-function(RNAbarcodes){
  #read data, subset into just the reads 
  data<-readFastq(RNAbarcodes)
  data<-as.character(data@sread)
  data<-lapply(data, substring, first=c(1, 17), last=c(16, 28))
  #first we just get the unique cellids and generate readcounts. 
  cellid<-unique(sapply(data, "[[",1))
  final<-data.frame(cellid)
  final$HIV_reads_RNA<-table(sapply(data, "[[",1))[final$cellid]
  #then we do UMIs, so for each we'll first grab the indicies that are in that 
  UMIs<-list()
  for(i in final$cellid){
    counted<-c()
    ontarget<-grep(i, sapply(data, "[[",1))
    UMIs[[i]]<-0
    for(j in ontarget){
      if(!sapply(data, "[[",2)[[j]] %in% counted){
        counted<-c(counted, sapply(data, "[[",2)[[j]])
        UMIs[[i]]<-UMIs[[i]]+1
      }
    }
  }
  final$HIV_UMI_RNA<-unlist(UMIs)
  final<-final[final$HIV_UMI_RNA>1 | final$HIV_reads_RNA>1,]
  
  return(final)
}

#function for ATAC annotation 

HIVATACAnnotation<-function(ATACbarcodes, bcs){
  #read data, subset into just the reads 
  data<-readFastq(ATACbarcodes)
  data<-as.character(data@sread)
  
  #easiest if we just setup the altered bcs first 
  cellid2<-FixBCs(unique(data), bcs)
  cellid<-unique(cellid2$RNA)

  final<-data.frame(cellid)
  final<-final[!is.na(final$cellid),]
  #noUMI to add since we arent encoding the location in this representation 
  counts<-c()
  for(i in 1:length(final)){
    counts[[i]]<-table(data)[cellid2$querry[cellid2$RNA==final[[i]]]]
  }
  counts<-unlist(lapply(counts, sum, na.rm=TRUE))
  final<-data.frame(cbind(final, counts))  

  #filter
  final<-final[final$counts>1,]
  colnames(final)<-c("cellid","HIV_reads_ATAC")
  return(final)
}

#function to wrap both, return a unified dataframe
ProcessHIVAnnotation<-function(RNAbarcodes, ATACbarcodes, prefix="", bcs){
  RNA<-HIVRNAAnnotation(RNAbarcodes)
  ATAC<-HIVATACAnnotation(ATACbarcodes, bcs)
  Results<-merge(RNA, ATAC, by="cellid", all=TRUE)
  rownames(Results)<-paste(prefix, Results$cellid, sep="_")
  Results<-subset(Results, select = -cellid)
  return(Results)
}


#now actually doing the annotation for HIV, also going to add in the fine-mapped CD4 annotation while we're here. 

results<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230606FinalObj.rds")

meta<-results@meta.data
meta<-split(meta, meta$orig.ident)

todo<-names(meta)
final<-rbind(ProcessHIVAnnotation(paste0("~/gibbs/DOGMAMORPH/data/HIV_RNA2/", todo[[1]],".results.txt" ), paste0("~/gibbs/DOGMAMORPH/data/HIV_ATAC/", todo[[1]]), prefix = todo[[1]], bcs)
             , ProcessHIVAnnotation(paste0("~/gibbs/DOGMAMORPH/data/HIV_RNA2/", todo[[2]],".results.txt" ), paste0("~/gibbs/DOGMAMORPH/data/HIV_ATAC/", todo[[2]]), prefix = todo[[2]], bcs))
for(i in 3:length(todo)){
  final<-rbind(final, ProcessHIVAnnotation(paste0("~/gibbs/DOGMAMORPH/data/HIV_RNA2/", todo[[i]],".results.txt" ), paste0("~/gibbs/DOGMAMORPH/data/HIV_ATAC/", todo[[i]]), prefix = todo[[i]], bcs))
}

rownames(final)<-paste0(rownames(final), "-1")

#little over 1250 cell ids that meet the very permissive cutoff
#notably, many of these start with a polyG suggestive of a null read
#literally 0 cells have both RNA and DNA 
#360 cells pass QC apparently
table(rownames(final)%in%colnames(results))

results<-AddMetaData(results, final)
#17 RNA positive cells
table(results$HIV_reads_RNA)
#348 DNA positive cells
table(results$HIV_reads_ATAC)
results$HIV_reads_ATAC[is.na(results$HIV_reads_ATAC)]<-0
results$HIV_reads_RNA[is.na(results$HIV_reads_RNA)]<-0

DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>0])+NoLegend()
DimPlot(results, cells.highlight  = colnames(results)[results$HIV_reads_RNA>0])+NoLegend()


p0<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>0])+NoLegend()
p1<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>1])+NoLegend()
p2<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>2])+NoLegend()
p3<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>3])+NoLegend()
p4<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>4])+NoLegend()
p5<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>5])+NoLegend()
p6<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>6])+NoLegend()
p7<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>7])+NoLegend()
p8<-DimPlot(results, cells.highlight = colnames(results)[results$HIV_reads_ATAC>8])+NoLegend()


(p0+p1+p2)/(p3+p4+p5)/(p6+p7+p8)

#adding the annotation from CD4
results2<-readRDS("gibbs/DOGMAMORPH/Ranalysis/Objects/20230529CD4ObjAnno.rds")
meta<-results2@meta.data
meta$trash<-""
meta<-meta[,c("CD4anno", "trash")]

results<-AddMetaData(results, meta)
table(results$CD4anno, results$merged_clusters)
results@meta.data<-select(results@meta.data, -trash)
results$clusteranno<-case_when(!is.na(results$CD4anno)~as.character(results$CD4anno), TRUE ~ as.character(results$merged_clusters))

results$clusteranno<-factor(results$clusteranno, levels = unlist(c(levels(results$CD4anno), levels(results$merged_clusters))))
Idents(results)<-results$clusteranno
DimPlot(results, label = TRUE)

saveRDS(results, "gibbs/DOGMAMORPH/Ranalysis/Objects/202306192023FinalClusternames.rds")
