library(Seurat)
library(Signac)
results<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230606FinalObj.rds")
#gene set sources are msigdb 
#just testing here if opioid response or signaling genes are different between Timepoints or treatments, these results suggest no. 

GOBP_RESPONSE_TO_MORPHINE<-c("ADCY8","CNR1","DRD2","DRD3","AIF1","FOSB","GHR","GPI","OPRK1",
                             "OPRM1","PCSK1","PENK","PITX3","PPP2R2A","PPP5C","PRKCE","PRKCG",
                             "RGS4","SRR","SLC1A1","TACR1","TACR3","CXCR4","PPP1R1B","MAP1LC3A",
                             "PPP1R9B","PEA15","FADD")
DefaultAssay(results)<-"RNA"
results<-AddModuleScore(results, features = list(GOBP_RESPONSE_TO_MORPHINE))
FeaturePlot(results, "Cluster1", min.cutoff = 'q5', max.cutoff = 'q95')



FeaturePlot(results, "Cluster1", min.cutoff = 0, max.cutoff = 0.16, split.by = "Treatment")

results$T_Tp<-paste(results$Treatment, results$Timepoint, sep = "_")
results$T_Tp<-factor(results$T_Tp, levels = c("Methadone_0","Methadone_3","Bup.Nalo_0", "Bup.Nalo_3", "Naltrexone_0", "Naltrexone_3"))

FeaturePlot(results, "Cluster1", min.cutoff = 0, max.cutoff = 0.16, split.by = "T_Tp")
VlnPlot(results, "Cluster1",idents = "Naive_B", group.by = "T_Tp")
REACTOME_OPIOID_SIGNALLING<-c("CAMKK1","PRKAR2B","CAMK2B","GNA15","GNAI3",
                              "PDE4A","GNB5","CAMK2A","PRKACA","ADCY2","GNB1",
                              "GNA11","ITPR3","MAPK1","PDYN","PLCB4","PPP2CB",
                              "PPP2R1A","PDE4C","PPP3CB","PRKAR1A","CAMKK2","GNB3",
                              "OPRM1","PPP2R5D","PDE4D","PPP2CA","PRKAR2A","GNAI2",
                              "GNB4","POMC","PDE1A","PLA2G4A","CREB1","PPP3CC",
                              "ADCY7","ITPR2","PDE1B","PRKCG","GNG13","GNG11","GNGT1",
                              "GNAI1","ADCY4","PPP1R1B","PPP2R1B","PLCB2","ADCY3",
                              "PPP3CA","GNAL","PRKACB","CAMK2D","CAMK2G","PLCB3",
                              "ITPR1","CAMK4","PRKCA","PDE1C","ADCY8","GNA14","GNAQ",
                              "ADCY9","GNG3","PRKCD","ADCY1","CDK5","PRKACG","GNGT2",
                              "GNG8","GNG4","AHCYL1","GNB2","GNG12","PPP1CA","NBEA",
                              "GRK2","ADCY5","GNG5","ADCY6","GNG7","KPNA2","PLCB1",
                              "PRKX","PDE4B","GNG2","PRKAR1B","CALM1","GNAT3","PPP3R1",
                              "GNG10","GNG4","PRKAR2B","ADCY4","PRKACA")

results<-AddModuleScore(results, features = list(REACTOME_OPIOID_SIGNALLING))
FeaturePlot(results, "Cluster1", min.cutoff = 'q5', max.cutoff = 'q95')



FeaturePlot(results, "Cluster1", min.cutoff = 0, max.cutoff = 0.16, split.by = "Treatment")

results$T_Tp<-paste(results$Treatment, results$Timepoint, sep = "_")
results$T_Tp<-factor(results$T_Tp, levels = c("Methadone_0","Methadone_3","Bup.Nalo_0", "Bup.Nalo_3", "Naltrexone_0", "Naltrexone_3"))
FeaturePlot(results, "Cluster1", min.cutoff = 0, max.cutoff = 0.16, split.by = "T_Tp")
VlnPlot(results, "Cluster1",idents = "Naive_B", group.by = "T_Tp")
VlnPlot(results, "Cluster1",idents = "Treg_proliferating", group.by = "T_Tp")



