library(dplyr)

gsea<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230711GSEApolar2res.rds")

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
to_test<-c("HOUNKPE_HOUSEKEEPING_GENES","NAKAYA_PBMC_FLUMIST_AGE_18_50YO_3DY_DN","SHEN_SMARCA2_TARGETS_UP",
           "GSE21063_WT_VS_NFATC1_KO_16H_ANTI_IGM_STIM_BCELL_UP","GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_UP","GSE22886_CD8_VS_CD4_NAIVE_TCELL_DN",
           "GSE41978_KLRG1_HIGH_VS_LOW_EFFECTOR_CD8_TCELL_DN","GSE41978_ID2_KO_VS_ID2_KO_AND_BIM_KO_KLRG1_LOW_EFFECTOR_CD8_TCELL_DN","GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_UP",
           "GSE36888_UNTREATED_VS_IL2_TREATED_TCELL_6H_UP","GSE17974_0.5H_VS_72H_IL4_AND_ANTI_IL12_ACT_CD4_TCELL_UP","GSE23505_UNTREATED_VS_4DAY_IL6_IL1_IL23_TREATED_CD4_TCELL_UP",
           "GSE17974_0.5H_VS_72H_IL4_AND_ANTI_IL12_ACT_CD4_TCELL_UP","GSE2770_IL12_AND_TGFB_ACT_VS_ACT_CD4_TCELL_48H_DN","BROCKE_APOPTOSIS_REVERSED_BY_IL6",
           "GSE21670_TGFB_VS_IL6_TREATED_STAT3_KO_CD4_TCELL_DN","GSE21670_UNTREATED_VS_IL6_TREATED_STAT3_KO_CD4_TCELL_DN","REACTOME_INFECTIOUS_DISEASE",
           "THAKAR_PBMC_INACTIVATED_INFLUENZA_AGE_21_30YO_RESPONDERS_7DY_UP","DEBIASI_APOPTOSIS_BY_REOVIRUS_INFECTION_UP","HALLMARK_INTERFERON_GAMMA_RESPONSE",
           "GSE26030_TH1_VS_TH17_RESTIMULATED_DAY5_POST_POLARIZATION_UP","GSE27434_WT_VS_DNMT1_KO_TREG_DN","GSE24574_BCL6_HIGH_VS_LOW_TFH_CD4_TCELL_UP",
           "TIAN_TNF_SIGNALING_VIA_NFKB","PHONG_TNF_TARGETS_UP","SANA_TNF_SIGNALING_UP")


library(dplyr)

gsea<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230711GSEAnaive1res.rds")

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
TGFB<-subset(gsea, grepl("TGFB", gsea$pathway)& padj < 0.05)

checkforgene<-function(gene, leadingedge){
  return(gene%in%leadingedge)
}

res<-c()
for(i in 1:nrow(gsea)){
  res<-c(res, checkforgene("ETS1", gsea$leadingEdge[[i]])|checkforgene("ETS2", gsea$leadingEdge[[i]]))
}
ETS<-subset(gsea, res & padj < 0.05)


View(subset(gsea, ID=="Meth_vs_Nal_3" & padj < 0.05))
View(subset(gsea, ID=="Meth_vs_Nal_0" & padj < 0.05))
View(subset(gsea, ID=="Meth_vs_Bup.Nalo_3" & padj < 0.05))
View(subset(gsea, ID=="Meth_vs_Bup.Nalo_0" & padj < 0.05))
View(subset(gsea, ID=="Bup.Nalo_vs_Nal_3" & padj < 0.05))
View(subset(gsea, ID=="Bup.Nalo_vs_Nal_0" & padj < 0.05))

#picking out a few that sound interesting to plot for visual inspection, found 29 that are interesting between conditions. 
to_test<-c("REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION","REACTOME_G_ALPHA_12_13_SIGNALLING_EVENTS","GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_BCELL_DN",
           "GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_DN","HALLMARK_INTERFERON_GAMMA_RESPONSE","GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_DN",
           "GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_DN","GSE41978_ID2_KO_VS_ID2_KO_AND_BIM_KO_KLRG1_LOW_EFFECTOR_CD8_TCELL_DN","GSE10239_NAIVE_VS_KLRG1HIGH_EFF_CD8_TCELL_UP",
           "GSE7460_CD8_TCELL_VS_CD4_TCELL_ACT_UP","GSE9650_NAIVE_VS_EFF_CD8_TCELL_DN","GSE3039_NKT_CELL_VS_ALPHAALPHA_CD8_TCELL_DN",
           "GSE33424_CD161_INT_VS_NEG_CD8_TCELL_UP","GSE21670_UNTREATED_VS_IL6_TREATED_STAT3_KO_CD4_TCELL_DN","GSE21670_UNTREATED_VS_TGFB_IL6_TREATED_CD4_TCELL_DN",
           "GSE2770_TGFB_AND_IL4_ACT_VS_ACT_CD4_TCELL_2H_UP","GSE39820_CTRL_VS_IL1B_IL6_IL23A_CD4_TCELL_UP","GSE2770_IL12_AND_TGFB_ACT_VS_ACT_CD4_TCELL_48H_DN",
           "REACTOME_INFLUENZA_INFECTION","GSE34205_HEALTHY_VS_FLU_INF_INFANT_PBMC_UP","HOWARD_MONOCYTE_INACT_MONOV_INFLUENZA_A_INDONESIA_05_2005_H5N1_AGE_18_49YO_1DY_UP",
           "REACTOME_HIV_INFECTION","HALLMARK_INTERFERON_ALPHA_RESPONSE","ZHANG_INTERFERON_RESPONSE",
           "GSE24574_BCL6_HIGH_TFH_VS_TFH_CD4_TCELL_DN","GSE26030_TH1_VS_TH17_RESTIMULATED_DAY5_POST_POLARIZATION_UP","GSE3982_MEMORY_CD4_TCELL_VS_TH1_DN",
           "PHONG_TNF_RESPONSE_NOT_VIA_P38","WANG_TNF_TARGETS","ZHOU_TNF_SIGNALING_4HR","PHONG_TNF_TARGETS_UP","GSE11057_CD4_EFF_MEM_VS_PBMC_UP")
