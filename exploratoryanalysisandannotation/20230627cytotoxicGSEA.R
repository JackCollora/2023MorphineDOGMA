library(dplyr)

gsea<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230622GSEACtyoTTrad.rds")

colnames(gsea)

#not super informative 
tops<-gsea %>% group_by(ID)%>%slice_max(order_by = NES, n=10)
bots<-gsea %>% group_by(ID)%>%slice_min(order_by = NES, n=10)

tops<-subset(tops, padj<0.05)
bots<-subset(bots, padj<0.05)

TNF<-subset(gsea, grepl("TNF", gsea$pathway)& padj < 0.05)
INF<-subset(gsea, grepl("INF", gsea$pathway)& padj < 0.05)
INT<-subset(gsea, grepl("INTERFERON", gsea$pathway)& padj < 0.05)
IL<-subset(gsea, grepl("_IL", gsea$pathway)& padj < 0.05)
IL6<-subset(gsea, grepl("IL6", gsea$pathway)& padj < 0.05)
TH<-subset(gsea, grepl("TH1|TH2|TH17|TH9|TH12|TREG|TFH", gsea$pathway)& padj < 0.05)
CD8<-subset(gsea, grepl("CD8", gsea$pathway)& padj < 0.05)


View(subset(gsea, ID=="Meth_vs_Nal_3" & padj < 0.05))
View(subset(gsea, ID=="Meth_vs_Nal_0" & padj < 0.05))
View(subset(gsea, ID=="Meth_vs_Bup.Nalo_3" & padj < 0.05))
View(subset(gsea, ID=="Meth_vs_Bup.Nalo_0" & padj < 0.05))
View(subset(gsea, ID=="Bup.Nalo_vs_Nal_3" & padj < 0.05))
View(subset(gsea, ID=="Bup.Nalo_vs_Nal_0" & padj < 0.05))

#picking out a few that sound interesting to plot for visual inspection, found 29 that are interesting between conditions. 
to_test<-c("PHONG_TNF_TARGETS_UP", "GSE18791_CTRL_VS_NEWCASTLE_VIRUS_DC_6H_DN","GSE6269_FLU_VS_E_COLI_INF_PBMC_DN","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
           "BOSCO_INTERFERON_INDUCED_ANTIVIRAL_MODULE", "GSE34205_HEALTHY_VS_RSV_INF_INFANT_PBMC_UP", "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM", "GSE22886_NAIVE_CD8_TCELL_VS_MEMORY_TCELL_UP",
           "DANG_BOUND_BY_MYC","GSE34205_HEALTHY_VS_FLU_INF_INFANT_PBMC_UP","MARSON_BOUND_BY_FOXP3_UNSTIMULATED","ZHENG_BOUND_BY_FOXP3","GSE26030_TH1_VS_TH17_RESTIMULATED_DAY5_POST_POLARIZATION_UP",
           "ICHIBA_GRAFT_VERSUS_HOST_DISEASE_D7_UP","GSE21670_UNTREATED_VS_TGFB_IL6_TREATED_CD4_TCELL_DN","WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA",
           "GSE24574_BCL6_HIGH_TFH_VS_TFH_CD4_TCELL_DN","GSE2770_TGFB_AND_IL4_ACT_VS_ACT_CD4_TCELL_2H_UP","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "GSE21063_WT_VS_NFATC1_KO_16H_ANTI_IGM_STIM_BCELL_UP","GSE41978_KLRG1_HIGH_VS_LOW_EFFECTOR_CD8_TCELL_DN","ICHIBA_GRAFT_VERSUS_HOST_DISEASE_35D_UP", "BLANCO_MELO_RESPIRATORY_SYNCYTIAL_VIRUS_INFECTION_A594_CELLS_UP",
           "BROWNE_HCMV_INFECTION_14HR_UP","WIELAND_UP_BY_HBV_INFECTION","REACTOME_HIV_INFECTION", "GSE43955_TH0_VS_TGFB_IL6_IL23_TH17_ACT_CD4_TCELL_60H_DN", 
           "GSE5960_TH1_VS_ANERGIC_TH1_UP", "GSE12392_IFNAR_KO_VS_IFNB_KO_CD8_NEG_SPLEEN_DC_UP")
