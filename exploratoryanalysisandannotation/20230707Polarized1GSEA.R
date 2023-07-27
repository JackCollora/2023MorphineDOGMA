library(dplyr)

gsea<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230707GSEApolar1res.rds")

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
to_test<-c("GSE11057_CD4_EFF_MEM_VS_PBMC_UP","GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_UP","NAKAYA_PBMC_FLUMIST_AGE_18_50YO_3DY_DN","GSE11057_CD4_CENT_MEM_VS_PBMC_DN","GSE21670_UNTREATED_VS_TGFB_IL6_TREATED_CD4_TCELL_DN",
           "GSE2770_TGFB_AND_IL4_ACT_VS_ACT_CD4_TCELL_2H_UP","GSE8685_IL2_STARVED_VS_IL2_ACT_IL2_STARVED_CD4_TCELL_DN","GSE21670_IL6_VS_TGFB_AND_IL6_TREATED_STAT3_KO_CD4_TCELL_UP","GSE43955_10H_VS_60H_ACT_CD4_TCELL_WITH_TGFB_IL6_UP",
           "GSE43955_TH0_VS_TGFB_IL6_TH17_ACT_CD4_TCELL_52H_UP","GSE6269_FLU_VS_E_COLI_INF_PBMC_DN","GSE45365_CTRL_VS_MCMV_INFECTION_NK_CELL_DN","WIELAND_UP_BY_HBV_INFECTION","GSE36826_WT_VS_IL1R_KO_SKIN_STAPH_AUREUS_INF_UP",
           "BROWNE_INTERFERON_RESPONSIVE_GENES","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","GSE24574_BCL6_HIGH_TFH_VS_TFH_CD4_TCELL_DN","GSE24574_BCL6_HIGH_VS_LOW_TFH_CD4_TCELL_UP",
           "GSE43863_TH1_VS_TFH_MEMORY_CD4_TCELL_UP","GSE22886_NAIVE_CD4_TCELL_VS_12H_ACT_TH2_DN","GSE5960_TH1_VS_ANERGIC_TH1_UP","GSE14415_INDUCED_VS_NATURAL_TREG_UP","GSE14350_TREG_VS_TEFF_UP","GSE14308_TH1_VS_TH17_DN","ZHOU_TNF_SIGNALING_4HR",
           "ZHOU_TNF_SIGNALING_4HR","GSE16385_IFNG_TNF_VS_IL4_STIM_MACROPHAGE_UP","PHONG_TNF_TARGETS_UP")
