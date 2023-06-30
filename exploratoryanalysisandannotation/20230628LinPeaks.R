library(Seurat)
library(Signac)

results<-readRDS("gibbs/DOGMAMORPH/Ranalysis/Objects/20230626FinalClusternames.rds")
results$T_Tp<-paste(results$Treatment, results$Timepoint, sep = "_")
results$T_Tp<-factor(results$T_Tp, levels = c("Methadone_0","Methadone_3","Bup.Nalo_0", "Bup.Nalo_3", "Naltrexone_0", "Naltrexone_3"))
results<-SplitObject(results, "T_Tp")

results<-lapply(results, LinkPeaks, peak.assay="ATAC", expression.assay="RNA")
saveRDS(results, "gibbs/DOGMAMORPH/Ranalysis/Objects/20230628Linkpeakobj.rds")
