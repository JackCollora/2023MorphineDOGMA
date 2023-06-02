library(ArchR)
set.seed(1)
addArchRThreads(threads = 10) 
addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(
  inputFiles = "../data/JC4/atac_fragments.tsv.gz",
  sampleNames = c("JC4"),
  minTSS = 1, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles<-read
# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
#   LSIMethod = 1
# )

Proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = proj_name,
  copyArrows = FALSE
)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", 
                        name = "IterativeLSI", force=TRUE)

chrs <- getChromSizes(proj)
var_features <- proj@reducedDims[["IterativeLSI"]]$LSIFeatures
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)



# Peaks using NFR fragments
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addGroupCoverages(proj, maxFragmentLength=147)
proj <- addReproduciblePeakSet(proj)
# Counts
proj <- addPeakMatrix(proj, maxFragmentLength=147, ceiling=10^9)

# Save 
proj <- saveArchRProject(ArchRProj = proj)


