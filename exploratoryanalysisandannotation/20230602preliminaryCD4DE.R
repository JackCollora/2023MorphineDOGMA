library(UpSetR)
library(Seurat)
library(Signac) 
library(ggplot2)
library(dplyr)

#loading data
results<-readRDS( "~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230529CD4ObjAnno.rds")
#setting correct default
DefaultAssay(results)<-"RNA"
#generating a new variable that will make the timepoint 3s separate but the timepoint 0s together where in theory they are biologically equivalent
#will test this hypothesis later
results$Timepoint_3 <-case_when(results$Timepoint==0 ~ "0", T~as.character(results$Treatment))

#this comparison has no significant DE genes 
Prolif_3_B.N <-FindMarkers(results, "Bup.Nalo", 0, group.by = "Timepoint_3",subset.ident = "Proliferating")
Prolif_3_Meth <-FindMarkers(results, "Methadone", 0, group.by = "Timepoint_3",subset.ident = "Proliferating")
Prolif_3_Nal <-FindMarkers(results, "Naltrexone", 0, group.by = "Timepoint_3",subset.ident = "Proliferating")

#trying it more directly 
results.sub<-subset(results, idents = "Proliferating")
Idents(results.sub)<-results.sub$Treatment

#single DE gene here, basically just RunX1 in Bup.Nalo 
B.N <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Bup.Nalo")
Meth <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Methadone")
Nal <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Naltrexone")

#trying a larger population

results.sub<-subset(results, idents = "Exhausted_Th1")
Idents(results.sub)<-results.sub$Treatment

#small handful of differential genes here 
B.N <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Bup.Nalo")
Meth <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Methadone")
Nal <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Naltrexone")


results.sub<-subset(results, idents = "Effector_Th1_CD4_and_CD8")
Idents(results.sub)<-results.sub$Treatment

#another small handful here, in each case Bup.Nalo is having the largest effect size.
B.N <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Bup.Nalo")
Meth <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Methadone")
Nal <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Naltrexone")


results.sub<-subset(results, idents = "Polarized_Effector_T")
Idents(results.sub)<-results.sub$Treatment

#another small handful here, in each case Bup.Nalo is having the largest effect size.
B.N <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Bup.Nalo")
Meth <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Methadone")
Nal <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Naltrexone")

#checked if we're just filtering lowly expressed but differential genes, no difference
B.N <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Bup.Nalo", min.pct = -Inf)
Meth <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Methadone", min.pct = -Inf)
Nal <-FindMarkers(results.sub, 3,0, group.by = "Timepoint", subset.ident = "Naltrexone", min.pct = -Inf)

#trying to see if doing it as pseudobulk will work better 
library(DESeq2)

prepDESeq2<-function(seuratobj, group_col, rep_col){
  expression<-seuratobj@assays$RNA@counts
  results<-list()
  #first aggregate by participant
  for (i in names(table(seuratobj[[rep_col]]))){
    sub<-seuratobj[[rep_col]]==i
    results[[i]]<-expression[,sub[,1]]
    results[[i]]<-rowSums(results[[i]])
  }
  names(results)<-paste0("Psuedobulk_", names(results))
  expression<-as.data.frame(results)
  
  #now create the condition matrix based on that group_col 
  coldata<-group_by(seuratobj@meta.data, !!sym(group_col), !!sym(rep_col))%>%summarise()
  coldata<-as.data.frame(coldata)
  rownames(coldata)<-paste0("Psuedobulk_", coldata[[rep_col]])  
  coldata[,1]<-factor(coldata[,1])
  coldata<-coldata[colnames(expression),]
  #finally, we make a DESeq object
  expression<-DESeqDataSetFromMatrix(countData = expression, colData = coldata,
                         design = as.formula(paste("~ ", group_col) ) )
  return(expression)
}

dds<-prepDESeq2(results.sub, "Treatment","Participant")
dds<-dds[rowSums(counts(dds)) >= 10,]
dds <- DESeq(dds)

resultsNames(dds)

res_Meth_v_Nal<-results(dds, name = "Intercept")
res_Meth_v_Bup.Nal<-results(dds, name = "Treatment_Methadone_vs_Bup.Nalo")
res_Nal_v_Bup.Nal<-results(dds, name = "Treatment_Naltrexone_vs_Bup.Nalo")
# 12k genes up here
summary(res_Meth_v_Nal)
#basically no fold changes for either of these, over half have low counts. In all comparisons 15 are outliers 
summary(res_Meth_v_Bup.Nal)
summary(res_Nal_v_Bup.Nal)

#p values do backup that above observation 
hist(res_Meth_v_Nal$pvalue)
hist(res_Meth_v_Bup.Nal$pvalue)
hist(res_Nal_v_Bup.Nal$pvalue)

#trying this again but now we'll compare across groups
results.sub<-subset(results, idents = "Proliferating")

Idents(results.sub)<-paste(results.sub$Treatment, results.sub$Timepoint, sep="_")
table(Idents(results.sub))

#single DE gene here, basically just RunX1 in Bup.Nalo 
#still only a handful of genes different between conditions 
B.N_v_Meth_3 <-FindMarkers(results.sub, "Bup.Nalo_3","Methadone_3")
Meth_v_Nal_3 <-FindMarkers(results.sub, "Methadone_3","Naltrexone_3")
Nal_v_B.N_3 <-FindMarkers(results.sub, "Naltrexone_3","Bup.Nalo_3")


#overall it seems that the data are too sparse, we aren't getting differences in gene expression that are useful for downstream in the subset of clusters I've taken 
#two new approaches to test 
#1. Going to start less biased in what is interesting, will be generating distribution plots for clusters to find oens that are "perturbed" to focus further 
#2. Switching to ATAC-first. The data are more robust at an ATAC level 

VlnPlot(results, c("nCount_ATAC", "nFeature_ATAC"), pt.size = 0)


