library(Seurat)

results<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230626FinalClusternames.rds")
Idents(results)<-results$merged_clusters
results$T_Tp<-paste(results$Treatment, results$Timepoint, sep = "_")
results$T_Tp<-factor(results$T_Tp, levels = c("Methadone_0","Methadone_3","Bup.Nalo_0", "Bup.Nalo_3", "Naltrexone_0", "Naltrexone_3"))


clusters<-unique(results$merged_clusters)
treatmentcomps<-list(c("Methadone_3", "Methadone_0"), c("Bup.Nalo_3","Bup.Nalo_0"), c("Naltrexone_3", "Naltrexone_0"),
                     c("Methadone_3", "Bup.Nalo_3"), c("Bup.Nalo_3","Naltrexone_3"), c("Naltrexone_3", "Methadone_3"),
                     c("Methadone_0", "Bup.Nalo_0"), c("Bup.Nalo_0","Naltrexone_0"), c("Naltrexone_0", "Methadone_0"))
names(treatmentcomps)<-paste(c("Methadone_3", "Bup.Nalo_3", "Naltrexone_3", "Methadone_3", "Bup.Nalo_3", "Naltrexone_3", "Methadone_0","Bup.Nalo_0", "Naltrexone_0"),
                             c("Methadone_0", "Bup.Nalo_0", "Naltrexone_0", "Bup.Nalo_3", "Naltrexone_3", "Methadone_3", "Bup.Nalo_0", "Naltrexone_0","Methadone_0"), 
                             sep = "_vs_")
DefaultAssay(results)<-"RNA"
cluster_marks<-list()
for(i in clusters){
  current<-list()
  for (j in names(treatmentcomps)){
    current[[j]]<-FindMarkers(results, ident.1 = treatmentcomps[[j]][[1]],ident.2 = treatmentcomps[[j]][[2]], 
                              group.by = "T_Tp", subset.ident = i)
  }
  cluster_marks[[i]]<-current
}

saveRDS(cluster_marks, "~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230702AllClustersTreatmentconditionsRNA.rds")

DefaultAssay(results)<-"ATAC"
cluster_marks<-list()
for(i in clusters){
  current<-list()
  for (j in names(treatmentcomps)){
    current[[j]]<-FindMarkers(results, ident.1 = treatmentcomps[[j]][[1]],ident.2 = treatmentcomps[[j]][[2]], 
                              group.by = "T_Tp", subset.ident = i,
                              test.use = 'LR',
                              min.pct = 0.05,
                              latent.vars = 'nCount_ATAC')
  }
  cluster_marks[[i]]<-current
}

saveRDS(cluster_marks, "~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230702AllClustersTreatmentconditionsATAC.rds")

DefaultAssay(results)<-"chromvar"
cluster_marks<-list()
for(i in clusters){
  current<-list()
  for (j in names(treatmentcomps)){
    current[[j]]<-FindMarkers(results, ident.1 = treatmentcomps[[j]][[1]],ident.2 = treatmentcomps[[j]][[2]], 
                              group.by = "T_Tp", subset.ident = i, 
                              mean.fxn = rowMeans,
                              fc.name = "avg_diff")
  }
  cluster_marks[[i]]<-current
}

saveRDS(cluster_marks, "~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230702AllClustersTreatmentconditionschromVAR.rds")
