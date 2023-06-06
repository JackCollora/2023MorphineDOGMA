#working based on the tutorial from Signac

library(Signac)
library(Seurat)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
set.seed(1234)

results<-readRDS("~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230518Completeobj.rds")

#naming taken from the ADT notebook
results$merged_clusters<-case_when(results$seurat_clusters==16~8,results$seurat_clusters==17~8, results$seurat_clusters==18~10,results$seurat_clusters==19~3,
                                   results$seurat_clusters==20~4,results$seurat_clusters==21~6,results$seurat_clusters==22~2,results$seurat_clusters==23~2,
                                   results$seurat_clusters==24~6, results$seurat_clusters==0~0,T~(as.numeric(results$seurat_clusters)-1))

results$merged_clusters<-case_when(results$merged_clusters==0~"Naive_CD4_T_1",results$merged_clusters==1~"Memory_CD4_Polarized_1",results$merged_clusters==2~"Memory_CD4_Polarized_2",results$merged_clusters==3~"Memory_CD4_1",
                                   results$merged_clusters==4~"Cytotoxic_T",results$merged_clusters==5~"Naive_B",results$merged_clusters==6~"Naive_CD4_T_2",results$merged_clusters==7~"Treg_proliferating",
                                   results$merged_clusters==8~"Naive_CD4_T_3",results$merged_clusters==9~"NK",results$merged_clusters==10~"Memory_B",results$merged_clusters==11~"CD14+_Mono",
                                   results$merged_clusters==12~"cDC",results$merged_clusters==13~"pDC",results$merged_clusters==14~"Memory_CD4_2",results$merged_clusters==15~"Plasma")

results$merged_clusters<-factor(results$merged_clusters, levels = c("Memory_CD4_Polarized_1","Memory_CD4_Polarized_2","Treg_proliferating","Cytotoxic_T","NK","Memory_CD4_1","Memory_CD4_2",
                                                                    "Naive_CD4_T_1","Naive_CD4_T_2","Naive_CD4_T_3","cDC","pDC","CD14+_Mono","Naive_B","Memory_B","Plasma"))

Idents(results)<-results$merged_clusters

DefaultAssay(results)<-"ATAC"

pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)


results <- AddMotifs(
  object = results,
 pfm = pfm, 
 genome = BSgenome.Hsapiens.UCSC.hg38
)
BiocParallel::register(BPPARAM = BiocParallel::SerialParam())

results <- RunChromVAR(
  object = results,
  genome = BSgenome.Hsapiens.UCSC.hg38
)


saveRDS(results, "~/gibbs/DOGMAMORPH/Ranalysis/Objects/20230606FinalObj.rds")
