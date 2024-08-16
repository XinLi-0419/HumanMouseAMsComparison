library("Seurat")
library("patchwork")
library("ggplot2")
library("dplyr")
library("plyr")
library("ggthemes")
library("paletteer")
library("easyGgplot2")
library("rstatix")
library("ComplexHeatmap")
library("ggpubr")
library("scales")
library("gridExtra")
library("SummarizedExperiment")
library("MetaNeighbor")
library("circlize")
library("plotly")
library("shiny")
library("fmsb")
library("RColorBrewer")
library("ggrepel")
library("topGO")
library("pcaExplorer")
library("data.table")
library("SeuratWrappers")
library("raster")
library("tools")
library("viridis")
library("ggeasy")
library("cowplot")
library("patchwork")
library("scater")
library("org.Hs.eg.db")
library("Matrix")
library("DESeq2")
library("rcartocolor")
library("reshape2")
library("gplots")
library("cowplot")
library("flipPlots")
library("openxlsx")
# Other functions are in Extra.Functions.R

##hAMs##########################################################################################################################################################################################################################
# This letter is a reanalysis of previously reported human AM dataset and mouse AM dataset
# From DOI: 10.26508/lsa.202201458
# Available raw data: GSE193782; 10x counts: https://drive.google.com/drive/folders/1-YZ0Nd6cuKPk7VgHAOLhu5A4g-M_1yFU?usp=share_link
# Available originally processed data: https://cells.ucsc.edu/?ds=ams-supercluster
# Available process data for this letter: https://cells.ucsc.edu/?ds=ams-human-mouse+human-ams
# This letter referred to the previous AM subset classification but re-did a more proper classification with mininal changes

# Raw data processing #
# hBAL.list
hBAL.HC1.data <- Read10X(data.dir = "CJ17/count/sample_feature_bc_matrix") # https://drive.google.com/drive/folders/1-YZ0Nd6cuKPk7VgHAOLhu5A4g-M_1yFU?usp=share_link
hBAL.HC2.data <- Read10X(data.dir = "CJ18/count/sample_feature_bc_matrix") # https://drive.google.com/drive/folders/1-YZ0Nd6cuKPk7VgHAOLhu5A4g-M_1yFU?usp=share_link
hBAL.HC3.data <- Read10X(data.dir = "CJ23/count/sample_feature_bc_matrix") # https://drive.google.com/drive/folders/1-YZ0Nd6cuKPk7VgHAOLhu5A4g-M_1yFU?usp=share_link
hBAL.HC4.data <- Read10X(data.dir = "CJ24/count/sample_feature_bc_matrix") # https://drive.google.com/drive/folders/1-YZ0Nd6cuKPk7VgHAOLhu5A4g-M_1yFU?usp=share_link
hBAL.CF1.data <- Read10X(data.dir = "CJ19-3/count/sample_feature_bc_matrix") # https://drive.google.com/drive/folders/1-YZ0Nd6cuKPk7VgHAOLhu5A4g-M_1yFU?usp=share_link
hBAL.CF2.data <- Read10X(data.dir = "CJ21/count/sample_feature_bc_matrix") # https://drive.google.com/drive/folders/1-YZ0Nd6cuKPk7VgHAOLhu5A4g-M_1yFU?usp=share_link
hBAL.CF3.data <- Read10X(data.dir = "CJ22/count/sample_feature_bc_matrix") # https://drive.google.com/drive/folders/1-YZ0Nd6cuKPk7VgHAOLhu5A4g-M_1yFU?usp=share_link
hBAL.data <- list(hBAL.HC1.data, hBAL.HC2.data, hBAL.HC3.data, hBAL.HC4.data, hBAL.CF1.data, hBAL.CF2.data, hBAL.CF3.data)
hBAL.list.vector <- c("hBAL.HC1", "hBAL.HC2", "hBAL.HC3", "hBAL.HC4", "hBAL.CF1", "hBAL.CF2", "hBAL.CF3")
nFeature_RNA.min <- c(2000, 2750, 2500, 1300, 2000, 2500, 2500)
nFeature_RNA.max <- c(7500, 7000, 7000, 5000, 7500, 7000, 7000)
percent.mt.max <- c(12, 10, 12, 11, 11, 11, 11)
hBAL.list <- list()
for (i in 1:length(hBAL.list.vector)){
  hBAL.list[[i]] <- CreateSeuratObject(counts = hBAL.data[[i]][["Gene Expression"]], project = hBAL.list.vector[i], min.cells = 3, min.features = 200)
  hBAL.list[[i]]$"orig.ident" <- hBAL.list.vector[i]
  hBAL.list[[i]]$"orig.treatment" <- c(rep("hBAL.HC", 4), rep("hBAL.CF", 3))[i]
  hBAL.list[[i]]$"orig.tissue" <- "hBAL"
  hBAL.list[[i]]$"orig.species" <- "human"
  hBAL.list[[i]]$"Note" <- rep(rep(c("ADT:No.HLA-DR.No.CD169", "ADT:Yes.HLA-DR.Yes.CD169"), each = 2), 2)[-6][i]
  hBAL.list[[i]]$"percent.mt" <- PercentageFeatureSet(hBAL.list[[i]], pattern = "^MT-")
  hBAL.list[[i]] <- subset(hBAL.list[[i]], subset = nFeature_RNA > nFeature_RNA.min[i] & nFeature_RNA < nFeature_RNA.max[i] & percent.mt < percent.mt.max[i]) ##Cell Filter Adjustment##
  all.cells <- colnames(hBAL.list[[i]][["RNA"]])
  rownames(hBAL.data[[i]][["Antibody Capture"]]) <- mapvalues(rownames(hBAL.data[[i]][["Antibody Capture"]]), c("CD14.1", "CD93.1", "CD197", "XCR1.1", "CD163.1", "CD36.1"), c("CD14", "CD93", "CCR7", "XCR1", "CD163", "CD36"))
  hBAL.data[[i]][["ADT"]] <- hBAL.data[[i]][["Antibody Capture"]][, all.cells]
  hBAL.list[[i]][["ADT"]] <- CreateAssayObject(counts = hBAL.data[[i]][["ADT"]])
  hBAL.list[[i]] <- RenameCells(hBAL.list[[i]], add.cell.id = hBAL.list.vector[i])
  hBAL.list[[i]] <- hBAL.list[[i]] %>% SCTransform(vars.to.regress = "percent.mt", method = "glmGamPoi") %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50) %>% FindClusters()
  hBAL.list[[i]] <- NormalizeData(hBAL.list[[i]], normalization.method = "CLR", margin = 2, assay = "ADT")
}
# hBAL.combined
hBAL.features <- SelectIntegrationFeatures(object.list = hBAL.list, nfeatures = 10000)
hBAL.features.list <- PrepSCTIntegration(object.list = hBAL.list, anchor.features = hBAL.features)
hBAL.anchors <- FindIntegrationAnchors(object.list = hBAL.features.list, normalization.method = "SCT", anchor.features = hBAL.features)
hBAL.combined <- IntegrateData(anchorset = hBAL.anchors, normalization.method = "SCT", new.assay.name = "ITG")
hBAL.combined <- hBAL.combined %>% RunPCA(assay = "ITG") %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50, return.model=TRUE) %>% FindClusters()
# Based on the cell markers and top DEGs, defined the major cell types, including AMs, Mono, FOLR2.IMs, SPP1.IMs, DC1, DC2, Mig.DCs, pDCs, Lym, Cyc.Lym, Cyc.Mye, and Epi
# Detailed methods: https://www.life-science-alliance.org/content/5/11/e202201458

# hAMs #
hAMs <- hBAL.combined[, hBAL.combined$orig.treatment == "hBAL.HC"]
# Remove IFI27 and APOC2 to exclude the supercluster effect
hAMs[["ITG"]]@data[c(1, 3), ] <- 0 # IFI27, APOC2
hAMs[["ITG"]]@scale.data[c(1, 3), ] <- 0
DefaultAssay(hAMs) <- "ITG"
VariableFeatures(hAMs) <- setdiff(VariableFeatures(hAMs), c("IFI27", "APOC2"))
hAMs <- hAMs[, hAMs$Cell.Types == "AMs"]
hAMs <- hAMs %>% RunPCA(assay = "ITG") %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50, return.model=TRUE)
# ITG_snn_res.0.8
hAMs <- FindClusters(hAMs, graph.name = "ITG_snn", resolution = 0.8)
hAMs.0.8.DEGs <- FindAllMarkers(hAMs, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, recorrect_umi = FALSE)
# Based on hAMs.0.8.DEGs and the DEG expression pattern to classify AM subsets
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 0][1:100] # no
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 1][1:100] # no
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 2][1:100] # no
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 3][1:100] # ✓ CCL18
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 4][1:100] # no
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 5][1:100] # ✓ IGF1, HELLPAR
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 6][1:100] # ✓ all cholesterol
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 7][1:100] # FN1? not really, overlap with cluster 14
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 8][1:100] # no
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 9][1:100] # ✓ C15orf48, MARCKS, CCL23, CXCL9/10/11, CXCL5
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 10][1:100] # MME, SCD, FDX1? not obvious
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 11][1:100] # CXCR4, CD48, like cluster 9
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 12][1:100] # ✓ MT1X
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 13][1:100] # ✓ AREG, NR4A1, EGR1, FOSB, NR4A3
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 14][1:100] # ✓ CXCL5
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 15][1:100] # ✓ CCL4, CCL4L2, CCL20, CXCL8, CCL3, ICAM1
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 16][1:100] # no
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 17][1:100] # ✓ GDF15, CLGN, PSAT1 # RPL # SNHG/Starvation
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 18][1:100] # ✓ RSAD2, IFI44L, not really CXCL9/10/11
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 19][1:100] # CAMP
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 20][1:100] # C1orf56, LRRC75A, very un-obvious
hAMs.0.8.DEGs$gene[hAMs.0.8.DEGs$cluster == 21][1:100] # ✓ Cycling
# ITG_snn_res.0.8 cluster 9 is heterogeneous
# ITG_snn_res.0.8.0.4
hAMs <- FindSubCluster(hAMs, cluster = 9, graph.name = "ITG_snn", subcluster.name = "ITG_snn_res.0.8.0.4", resolution = 0.4)
DimPlot(hAMs, cells = colnames(hAMs)[hAMs$ITG_snn_res.0.8 == 9], group.by = "ITG_snn_res.0.8.0.4", label = TRUE)
hAMs.9 <- hAMs[, hAMs$ITG_snn_res.0.8 == 9]
Idents(hAMs.9) <- "ITG_snn_res.0.8.0.4"
hAMs.0.8.0.4.DEGs <- FindAllMarkers(hAMs.9, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, recorrect_umi = FALSE)
# Based on hAMs.0.8.0.4.DEGs and the DEG expression pattern to classify AM subsets
hAMs.0.8.0.4.DEGs$gene[hAMs.0.8.0.4.DEGs$cluster == "9_0"][1:100] # no
hAMs.0.8.0.4.DEGs$gene[hAMs.0.8.0.4.DEGs$cluster == "9_1"][1:100] # CCL18
hAMs.0.8.0.4.DEGs$gene[hAMs.0.8.0.4.DEGs$cluster == "9_2"][1:100] # CXCL5
hAMs.0.8.0.4.DEGs$gene[hAMs.0.8.0.4.DEGs$cluster == "9_3"][1:100] # CXCL10/11, CXCL9?
hAMs.0.8.0.4.DEGs$gene[hAMs.0.8.0.4.DEGs$cluster == "9_4"][1:100] # no

# AMs.Cell.Subtypes
hAMs$AMs.Cell.Subtypes <- gsub("SNHG.AMs", "GDF15.AMs", as.vector(hAMs$AMs.Cell.Types))
hAMs$AMs.Cell.Subtypes[grep("AMs.s", hAMs$AMs.Cell.Subtypes)] <- "Main.AMs"
hAMs$AMs.Cell.Subtypes[grep("CK.AMs.c", hAMs$AMs.Cell.Subtypes)] <- "Main.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$AMs.Cell.Subtypes == "CCL18.AMs" & hAMs@reductions$umap@cell.embeddings[, "umap_2"] > 2] <- "Main.AMs" # original CCL18.AMs not all TRUE
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 3] <- "CCL18.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 5] <- "IGF1.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 6] <- "Cholesterol.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 12] <- "MT.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 13] <- "GF.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 14] <- "CXCL5.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 15] <- "CK1.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 17] <- "GDF15.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 18] <- "IFN.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8 == 21] <- "Cyc.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8.0.4 == "9_1"] <- "CCL18.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8.0.4 == "9_2"] <- "CXCL5.AMs"
hAMs$AMs.Cell.Subtypes[hAMs$ITG_snn_res.0.8.0.4 == "9_3"] <- "CK2.AMs"
hAMs$AMs.Cell.Subtypes <- factor(hAMs$AMs.Cell.Subtypes, levels = 
                                               c("Main.AMs", "MT.AMs", "IFN.AMs", "CK1.AMs", "CK2.AMs", "CXCL5.AMs", "CCL18.AMs", "GDF15.AMs", "GF.AMs", "IGF1.AMs", "Cholesterol.AMs", "Cyc.AMs"))
DimPlot(hAMs, group.by = "AMs.Cell.Subtypes", label = TRUE)
Idents(hAMs) <- "AMs.Cell.Subtypes"
hAMs.Cell.Subtypes.DEGs <- FindAllMarkers(hAMs, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, recorrect_umi = FALSE)
# Table E1. AMs Subtypes DEGs.xlsx

##mAMs##########################################################################################################################################################################################################################
# This letter is a reanalysis of previously reported human AM dataset and mouse AM dataset
# From DOI: 10.1038/s41590-024-01826-9
# Available raw data: GSE225663
# Available originally processed data: https://cells.ucsc.edu/?ds=macrophage-atlas+mBAL
# Available process data for this letter: https://cells.ucsc.edu/?ds=ams-human-mouse+mouse-ams
# There is no reference for the AM subset classification for this letter

# Raw data processing #
# mBAL.list
mBAL.1.data <- Read10X(data.dir = "CJ41/outs/filtered_feature_bc_matrix") # GSE225663
mBAL.2.data <- Read10X(data.dir = "CJ43/outs/filtered_feature_bc_matrix") # GSE225663
mBAL.3.data <- Read10X(data.dir = "CJ45/outs/filtered_feature_bc_matrix") # GSE225663
mBAL.data <- list(mBAL.1.data, mBAL.2.data, mBAL.3.data)
mBAL.list.vector <- c("mBAL.1", "mBAL.2", "mBAL.3")
nFeature_RNA.min <- c(2600, 2400, 2600)
nFeature_RNA.max <- c(6000, 5500, 6100)
percent.mt.max <- c(7.5, 7.5, 7.5)
mBAL.list <- list()
for (i in 1:length(mBAL.list.vector)){
  mBAL.list[[i]] <- CreateSeuratObject(counts = mBAL.data[[i]], project = mBAL.list.vector[i], min.cells = 3, min.features = 200)
  mBAL.list[[i]]$"orig.ident" <- mBAL.list.vector[i]
  mBAL.list[[i]]$"orig.treatment" <- "Naive"
  mBAL.list[[i]]$"orig.tissue" <- "mBAL"
  mBAL.list[[i]]$"orig.species" <- "mouse"
  mBAL.list[[i]]$"percent.mt" <- PercentageFeatureSet(mBAL.list[[i]], pattern = "^mt-")
  mBAL.list[[i]] <- subset(mBAL.list[[i]], subset = nFeature_RNA > nFeature_RNA.min[i] & nFeature_RNA < nFeature_RNA.max[i] & percent.mt < percent.mt.max[i])
  mBAL.list[[i]] <- RenameCells(mBAL.list[[i]], add.cell.id = mBAL.list.vector[i])
  mBAL.list[[i]] <- mBAL.list[[i]] %>% SCTransform(vars.to.regress = "percent.mt", method = "glmGamPoi") %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50) %>% FindClusters()
}
##mBAL.combined
mBAL.features <- SelectIntegrationFeatures(object.list = mBAL.list, nfeatures = 10000)
mBAL.features.list <- PrepSCTIntegration(object.list = mBAL.list, anchor.features = mBAL.features)
mBAL.anchors <- FindIntegrationAnchors(object.list = mBAL.features.list, normalization.method = "SCT", anchor.features = mBAL.features)
mBAL.combined <- IntegrateData(anchorset = mBAL.anchors, normalization.method = "SCT", new.assay.name = "ITG")
mBAL.combined <- mBAL.combined %>% RunPCA(assay = "ITG") %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50, return.model=TRUE) %>% FindClusters()
# Based on the cell markers and top DEGs, defined the major cell types, including AMs, Hb.AMs, Cyc.AMs, IMs, Res.DCs, Mig.DCs, Lym, Epi
# Detailed methods: https://www.nature.com/articles/s41590-024-01826-9 Verification: https://cells.ucsc.edu/?ds=macrophage-atlas+mBAL

# mAMs #
mAMs <- mBAL.combined[, mBAL.combined$Cell.Types %in% c("AMs", "Cyc.AMs") & mBAL.combined$seurat_clusters != 43] # 43 is Hb.AMs
mAMs <- mAMs[, mAMs@reductions$umap@cell.embeddings[, "UMAP_1"] > 0.2 | mAMs@reductions$umap@cell.embeddings[, "UMAP_2"] > -7.5]
mAMs <- mAMs %>% RunPCA(assay = "ITG") %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50, return.model=TRUE)
# ITG_snn_res.0.8
mAMs <- FindClusters(mAMs, graph.name = "ITG_snn", resolution = 0.8)
mAMs.0.8.DEGs <- FindAllMarkers(mAMs, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, recorrect_umi = FALSE)
# Based on mAMs.0.8.DEGs and the DEG expression pattern to classify AM subsets
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 0][1:100] # no # but polarization with Krt79, Car4 and others # polarization
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 1][1:100] # no
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 2][1:100] # no # but polarization with Cd300a, Cd63 and others # polarization
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 3][1:100] # ✓ Car4
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 4][1:100] # ✓ Mt2/1 # metallothionein, not have as many gene as in human
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 5][1:100] # polarization with all the genes, including Slc7a11 # share with 17 and 12? # polarization
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 6][1:100] # ✓ Mcm6 and other Mcm # Cycling
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 7][1:100] # ✓ Rrm2, Pclaf and others # Cycling, should do cell-cycle score and compare m vs h
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 8][1:100] # ✓ Clec4e, Icam1, Cfb, Ccl9, Ccrl2, Cxcl16, Clec4d # CK (adhension) and NFkB
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 9][1:100] # ✓ Mt2/1, similar to 4
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 10][1:100] # ✓ Tppp3, Crip1, Fn1, Cav1 # cytoskeleton, or migration?
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 11][1:100] # ✓ H2-Aa, and other MHC
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 12][1:100] # ✓ Cxcr1, Gdf15, Cd200 # Trib3 # not very specific # Starvation
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 13][1:100] # ✓ Ccnb2, Ube2c, Cenpf, Birc5, Cdca3, Gm49980, Emp1, Lockd, Mki67 # Cycling
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 14][1:100] # ✓ Irf7, Isg15 # IFN
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 15][1:100] # ✓ Ccna2, Ube2c, Kif11 # Cycling
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 16][1:100] # ✓ AY036118 (ETS-related), Gm42418, maybe others including Ear2 # unknown # GO
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 17][1:100] # ✓ Hspa1a, Hspa1b, and others like Efr3b, Hfm1 # unknown # GO
mAMs.0.8.DEGs$gene[mAMs.0.8.DEGs$cluster == 18][1:100] # ✓ AREG, Egr1, Il1a/b, NFkB, Ccl3/4, TLR2, Cxcl1/2/3, not so much 10 # map to hAMs GF and CK
# Missing one cluster in ITG_snn_res.0.8, one cluster specific for GF including Egr1 was seen with mBAL.combined (cluster 41)
mAMs$ITG_snn_res.0.8 <- factor(mAMs$ITG_snn_res.0.8, levels = 0:18)
mAMs$ITG_snn_res.0.8[colnames(mAMs) %in% colnames(mBAL.combined)[mBAL.combined$seurat_clusters == 41]] <- 18 # Add this cluster as 18
Idents(mAMs) <- "ITG_snn_res.0.8"
mAMs.0.8.DEGs <- FindAllMarkers(mAMs, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, recorrect_umi = FALSE)
# not good enough to just use the cluster before # redo with module scores to create new cluster using scores
mAMs <- AddModuleScore(mAMs, features = list("features" = c18.DEGs[1:35][-c(13, 19, 29, 31, 32)]), name = "c18")
mAMs$AMs.Cell.Subtypes <- "Main.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 3] <- "Car4.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 8] <- "Ccl9.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 12] <- "Gdf15.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$c181 > 0.4] <- "CK.AMs" # Both CK and GF
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 4] <- "Mt.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 6] <- "Cyc.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 7] <- "Cyc.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 9] <- "Mt.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 10] <- "Tppp3.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 11] <- "Ag-presenting.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 13] <- "Cyc.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 14] <- "IFN.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 15] <- "Cyc.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 16] <- "AY036118.AMs"
mAMs$AMs.Cell.Subtypes[mAMs$ITG_snn_res.0.8 == 17] <- "HS.AMs"
mAMs$AMs.Cell.Subtypes <- factor(mAMs$AMs.Cell.Subtypes, levels = 
                                   c("Main.AMs", "IFN.AMs", "AY036118.AMs", "HS.AMs", "CK.AMs", "Ccl9.AMs", "Mt.AMs", "Ag-presenting.AMs", "Gdf15.AMs", "Car4.AMs", "Tppp3.AMs", "Cyc.AMs"))
DimPlot(mAMs, group.by = "AMs.Cell.Subtypes", label = TRUE)
Idents(mAMs) <- "AMs.Cell.Subtypes"
mAMs.Cell.Subtypes.DEGs <- FindAllMarkers(mAMs, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1, recorrect_umi = FALSE)
# Table E1. AMs Subtypes DEGs.xlsx

############################################################################################################################################################################################################################

# Figure.1A.hAMs.DimPlot
Figure.1.DimPlot <- DimPlot(hAMs, group.by = "AMs.Cell.Subtypes", cols = names(meta.data.list$hAMs.Cell.Subtypes), pt.size = 1.5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range1.01(hAMs@reductions$umap@cell.embeddings[, "umap_1"])) + ylim(range1.01(hAMs@reductions$umap@cell.embeddings[, "umap_2"])) + coord_cartesian(expand = FALSE)
Figure.1.DimPlot <- Figure.1.DimPlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
Figure.1.DimPlot.directory <- paste0("Figure.1A.hAMs", ".DimPlot.png")
png(Figure.1.DimPlot.directory, width = 1100, height = 1000, bg = "transparent")
print(Figure.1.DimPlot)
dev.off() 
# Figure.1E.mAMs.DimPlot
Figure.1.DimPlot <- DimPlot(mAMs, group.by = "AMs.Cell.Subtypes", cols = names(meta.data.list$mAMs.Cell.Subtypes), pt.size = 1.5) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range1.01(mAMs@reductions$umap@cell.embeddings[, "umap_1"])) + ylim(range1.01(mAMs@reductions$umap@cell.embeddings[, "umap_2"])) + coord_cartesian(expand = FALSE)
Figure.1.DimPlot <- Figure.1.DimPlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
Figure.1.DimPlot.directory <- paste0("Figure.1E.mAMs", ".DimPlot.png")
png(Figure.1.DimPlot.directory, width = 1100, height = 1000, bg = "transparent")
print(Figure.1.DimPlot)
dev.off() 

# Figure.1B.hAMs.DonutChart
Figure.1.DonutChart.Data <- as.data.frame(table(hAMs$AMs.Cell.Subtypes))[-12, ]
colnames(Figure.1.DonutChart.Data) <- c("AMs.Cell.Subtypes", "Number")
Figure.1.DonutChart <- Figure.1.DonutChart.Data %>% plot_ly(labels = ~AMs.Cell.Subtypes, values = ~Number, marker = list(colors = names(meta.data.list$hAMs.Cell.Subtypes)[-12]), sort = FALSE, direction = "clockwise", automargin = FALSE, textposition = "inside", textfont = list(size = 100))
Figure.1.DonutChart <- Figure.1.DonutChart %>% add_pie(hole = 0.5) %>% layout(showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
Figure.1.DonutChart.directory <- paste0("Figure.1B.hAMs", ".DonutChart.pdf")
save_image(Figure.1.DonutChart, Figure.1.DonutChart.directory, scale = 5)
# Figure.1F.mAMs.DonutChart
Figure.1.DonutChart.Data <- as.data.frame(table(mAMs$AMs.Cell.Subtypes))[-12, ]
colnames(Figure.1.DonutChart.Data) <- c("AMs.Cell.Subtypes", "Number")
Figure.1.DonutChart <- Figure.1.DonutChart.Data %>% plot_ly(labels = ~AMs.Cell.Subtypes, values = ~Number, marker = list(colors = names(meta.data.list$mAMs.Cell.Subtypes)[-12]), sort = FALSE, direction = "clockwise", automargin = FALSE, textposition = "inside", textfont = list(size = 100))
Figure.1.DonutChart <- Figure.1.DonutChart %>% add_pie(hole = 0.5) %>% layout(showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), plot_bgcolor = "transparent")
Figure.1.DonutChart.directory <- paste0("Figure.1F.mAMs", ".DonutChart.pdf")
save_image(Figure.1.DonutChart, Figure.1.DonutChart.directory, scale = 5)

# Figure.1C.hAMs.BarPlot
Figure.1.BarPlot.Data <- table(hAMs$AMs.Cell.Subtypes, hAMs$orig.ident)
Figure.1.BarPlot.Data <- Figure.1.BarPlot.Data[1:11, 1:4]
Figure.1.BarPlot.Data <- as.data.frame(t(Figure.1.BarPlot.Data) / colSums(Figure.1.BarPlot.Data))
colnames(Figure.1.BarPlot.Data) <- mapvalues(colnames(Figure.1.BarPlot.Data), c("Var1", "Var2"), c("orig.ident", "AMs.Cell.Subtypes"))
Figure.1.BarPlot.Data$orig.ident <- factor(Figure.1.BarPlot.Data$orig.ident, levels = levels(hAMs$orig.ident))
Figure.1.BarPlot.Data$AMs.Cell.Subtypes <- factor(Figure.1.BarPlot.Data$AMs.Cell.Subtypes, levels = rev(as.vector(meta.data.list$hAMs.Cell.Subtypes)))
Figure.1.BarPlot <- ggplot2.barplot(Figure.1.BarPlot.Data, xName = "orig.ident", yName= "Freq", width = 0.8, xTickLabelFont = c(22, "plain", color = "grey60"), yTickLabelFont = c(22, "plain", c(rep("black", 4), rep("red", 3))), groupName = "AMs.Cell.Subtypes", backgroundColor = "white", removePanelGrid = TRUE, removePanelBorder = TRUE, xShowTitle = FALSE, yShowTitle = FALSE, groupColors = rev(names(meta.data.list$hAMs.Cell.Subtypes))[-1]) + 
  scale_y_continuous(labels = scales::percent) + coord_flip(expand = FALSE) + ylab("Subject") + theme(axis.line = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(.25, "cm"), axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position = "none")
Figure.1.BarPlot.directory <- paste0("Figure.1C.hAMs", ".BarPlot.png")
png(Figure.1.BarPlot.directory, width = 1200, height = 200)
print(Figure.1.BarPlot)
dev.off()
# Figure.1G.mAMs.BarPlot
Figure.1.BarPlot.Data <- table(mAMs$AMs.Cell.Subtypes, mAMs$orig.ident)
Figure.1.BarPlot.Data <- Figure.1.BarPlot.Data[1:11, 1:3]
Figure.1.BarPlot.Data <- as.data.frame(t(Figure.1.BarPlot.Data) / colSums(Figure.1.BarPlot.Data))
colnames(Figure.1.BarPlot.Data) <- mapvalues(colnames(Figure.1.BarPlot.Data), c("Var1", "Var2"), c("orig.ident", "AMs.Cell.Subtypes"))
Figure.1.BarPlot.Data$orig.ident <- factor(Figure.1.BarPlot.Data$orig.ident, levels = paste0("mBAL.", 1:3))
Figure.1.BarPlot.Data$AMs.Cell.Subtypes <- factor(Figure.1.BarPlot.Data$AMs.Cell.Subtypes, levels = rev(as.vector(meta.data.list$mAMs.Cell.Subtypes)))
Figure.1.BarPlot <- ggplot2.barplot(Figure.1.BarPlot.Data, xName = "orig.ident", yName= "Freq", width = 0.8, xTickLabelFont = c(22, "plain", color = "grey60"), yTickLabelFont = c(22, "plain", c(rep("black", 4), rep("red", 3))), groupName = "AMs.Cell.Subtypes", backgroundColor = "white", removePanelGrid = TRUE, removePanelBorder = TRUE, xShowTitle = FALSE, yShowTitle = FALSE, groupColors = rev(names(meta.data.list$mAMs.Cell.Subtypes))[-1]) + 
  scale_y_continuous(labels = scales::percent) + coord_flip(expand = FALSE) + ylab("Subject") + theme(axis.line = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(.25, "cm"), axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position = "none")
Figure.1.BarPlot.directory <- paste0("Figure.1G.mAMs", ".BarPlot.png")
png(Figure.1.BarPlot.directory, width = 1200, height = 150)
print(Figure.1.BarPlot)
dev.off()

# Figure.1D.hAMs.VlnPlot
Figure.1.VlnPlot.Features.list <- list(
  "Main.AMs" = c("FABP4", "CSTA"),
  "IFN.AMs" = c("RSAD2", "IFIT1"),
  "MT.AMs" = c("MT1X", "MT1F"),
  "CK1.AMs" = c("CCL20", "CXCL8"),
  "CK2.AMs" = c("CXCL10", "CXCL11"),
  "CXCL5.AMs" = c("C15orf48", "CXCL5"),
  "CCL18.AMs" = c("CCL18", "CD99"),
  "GDF15.AMs" = c("GDF15", "TRIB3"),
  "GF.AMs" = c("AREG", "OSM"),
  "IGF1.AMs" = c("IGF1", "HELLPAR"),
  "Cholesterol.AMs" = c("CYP51A1", "MSMO1"),
  "Cyc.AMs" = c("PCLAF", "MKI67")
)
Figure.1.VlnPlot.Features <- as.vector(unlist(Figure.1.VlnPlot.Features.list))
Figure.1.VlnPlot.Features <- intersect(Figure.1.VlnPlot.Features, rownames(hAMs[["SCT"]]$data))
Figure.1.VlnPlot <- VlnPlot(hAMs, group.by = "AMs.Cell.Subtypes", stack = TRUE, features = Figure.1.VlnPlot.Features, fill.by = "ident", cols = names(meta.data.list$hAMs.Cell.Subtypes), flip = TRUE) + 
  theme(legend.position = "none", axis.title = element_blank(), strip.text.y = element_text(face = "italic", size = 25), axis.text.x = element_text(angle = 30, hjust = 0, size = 25), axis.ticks.x = element_blank()) + scale_x_discrete(position = "top")
Figure.1.VlnPlot.directory <- paste0("Figure.1D.hAMs", ".VlnPlot.png")
png(Figure.1.VlnPlot.directory, width = 800, height = 900) # 800 900 for original # angle 30 hjust 0
print(Figure.1.VlnPlot)
dev.off()
# Figure.1H.mAMs.VlnPlot
Figure.1.VlnPlot.Features.list <- list(
  "Main.AMs" = c("Fabp4", "Alox5ap"),
  "IFN.AMs" = c("Rsad2", "Ifit1"),
  "Mt.AMs" = c("Mt1", "Mt2"),
  "CK.AMs" = c("Cxcl2", "Ccl3"),
  "Ccl9.AMs" = c("Ccl9", "Cxcl16"),
  "Gdf15.AMs" = c("Gdf15", "Trib3"),
  "AY036118.AMs" = c("AY036118", "Gm42418"),
  "HS.AMs" = c("Hspa1a", "Hspa1b"),
  "Ag-presenting.AMs" = c("H2-Aa", "H2-Ab1"),
  "Car4.AMs" = c("Car4", "Tlr2"),
  "Tppp3.AMs" = c("Tppp3", "Crip1"),
  "Cyc.AMs" = c("Pclaf", "Mki67")
)
Figure.1.VlnPlot.Features <- as.vector(unlist(Figure.1.VlnPlot.Features.list))
Figure.1.VlnPlot.Features <- intersect(Figure.1.VlnPlot.Features, rownames(mAMs[["SCT"]]$data))
Figure.1.VlnPlot <- VlnPlot(mAMs, group.by = "AMs.Cell.Subtypes", stack = TRUE, features = Figure.1.VlnPlot.Features, fill.by = "ident", cols = names(meta.data.list$mAMs.Cell.Subtypes), flip = TRUE) +
  theme(legend.position = "none", axis.title = element_blank(), strip.text.y = element_text(face = "italic", size = 25), axis.text.x = element_text(angle = 30, hjust = 0, size = 25), axis.ticks.x = element_blank()) + scale_x_discrete(position = "top")
Figure.1.VlnPlot.directory <- paste0("Figure.1H.mAMs", ".VlnPlot.png")
png(Figure.1.VlnPlot.directory, width = 800, height = 915)
print(Figure.1.VlnPlot)
dev.off()

# Figure.2A.hAMs.FeaturePlot
Figure.2A.FeaturePlot.Features <- c("RSAD2", "IFIT1", "MT1X", "MT2A", "CCL3", "CCL4", "CXCL1", "CXCL2", "CXCL3", "CXCL9", "CXCL10", "CXCL5", "CCL18", "GDF15", "TRIB3", "EGR1", "OSM", "IGF1", "CYP51A1", "MSMO1", "PCLAF", "MKI67")
Figure.2A.FeaturePlot.Subtypes <- c(rep("IFN.AMs", 2), rep("MT.AMs", 2), rep("CK1.AMs", 5), rep("CK2.AMs", 2), "CXCL5.AMs", "CCL18.AMs", rep("GDF15.AMs", 2), rep("GF.AMs", 2), "IGF1.AMs", rep("Cholesterol.AMs", 2), rep("Cyc.AMs", 2))
for (i in 1:length(Figure.2A.FeaturePlot.Features)) {
  umap <- hAMs@reductions$umap@cell.embeddings[hAMs$AMs.Cell.Subtypes %in% Figure.2A.FeaturePlot.Subtypes[i], ]
  Figure.2.FeaturePlot <- FeaturePlot(hAMs, features = Figure.2A.FeaturePlot.Features[i], raster = FALSE, order = TRUE, pt.size = 0.5) + scale_colour_gradient(low = "#E5E5E5", high = "#990F0F") + 
    geom_density_2d(data = data.frame(umap), aes(x = umap_1, y = umap_2), bins = 15, linewidth = 2.5, alpha = 2, colour = "gray20") + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + # bins = ifelse(i == 2, ifelse(i == 5, 25, 8), 15) to make Mono tight, DCs bigger
    xlim(range.expand(hAMs@reductions$umap@cell.embeddings[, "umap_1"], expand.ratio = 0.05)) + ylim(range.expand(hAMs@reductions$umap@cell.embeddings[, "umap_2"], expand.ratio = 0.05)) + coord_cartesian(expand = FALSE)
  Figure.2.FeaturePlot <- ggplot_build(Figure.2.FeaturePlot)
  # highlight specific cell types
  umap_1 <- hAMs@reductions$umap@cell.embeddings[, "umap_1"][hAMs$AMs.Cell.Subtypes %in% meta.data.list$AMs.Cell.Subtypes[i]]
  Figure.2.FeaturePlot$data[[1]]$size[which(Figure.2.FeaturePlot$data[[1]]$x %in% umap[, "umap_1"])] <- 2
  # remove smaller contour
  Contour <- if(Figure.2A.FeaturePlot.Subtypes[i] == "CXCL5.AMs") c(6, 19, 14) else 
    if(Figure.2A.FeaturePlot.Subtypes[i] == "CCL18.AMs") c(1, 2) else 
      if(Figure.2A.FeaturePlot.Subtypes[i] == "GF.AMs") c(4) else c(which.max(table(Figure.2.FeaturePlot$data[[2]]$piece)))
  Figure.2.FeaturePlot$data[[2]] <- Figure.2.FeaturePlot$data[[2]][Figure.2.FeaturePlot$data[[2]]$piece == Contour, ] # to make Neu tight
  Figure.2.FeaturePlot.directory <- paste0("Figure.2A.hAMs.", Figure.2A.FeaturePlot.Features[i], ".FeaturePlot.png")
  png(Figure.2.FeaturePlot.directory, width = 1100/2, height = 1000/2)
  print(grid::grid.draw(ggplot_gtable(Figure.2.FeaturePlot)))
  dev.off() 
}
# Figure.2A.mAMs.FeaturePlot
Figure.2A.FeaturePlot.Features <- c("Rsad2", "Ifit1", "Mt1", "Mt2", "Ccl3", "Ccl4", "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl9", "Cxcl10", "Cxcl5", "Ccl3", "Gdf15", "Trib3", "Egr1", "Osm", "Igf1", "Cyp51", "Msmo1", "Pclaf", "Mki67")
Figure.2A.FeaturePlot.Subtypes <- c(rep("IFN.AMs", 2), rep("Mt.AMs", 2), rep("CK.AMs", 5), rep("IFN.AMs", 2), "CK1.AMs", "CK.AMs", rep("Gdf15.AMs", 2), rep("CK.AMs", 2), "Gdf15.AMs", rep("Car4.AMs", 2), rep("Cyc.AMs", 2))
for (i in 1:length(Figure.2A.FeaturePlot.Features)) {
  umap <- mAMs@reductions$umap@cell.embeddings[mAMs$AMs.Cell.Subtypes %in% Figure.2A.FeaturePlot.Subtypes[i], ]
  Figure.2.FeaturePlot <- FeaturePlot(mAMs, features = Figure.2A.FeaturePlot.Features[i], raster = FALSE, order = TRUE, pt.size = 0.5) + scale_colour_gradient(low = "#E5E5E5", high = "#990F0F") + 
    geom_density_2d(data = data.frame(umap), aes(x = umap_1, y = umap_2), bins = 15, linewidth = 2.5, alpha = 2, colour = "gray20") + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + # bins = ifelse(i == 2, ifelse(i == 5, 25, 8), 15) to make Mono tight, DCs bigger
    xlim(range.expand(mAMs@reductions$umap@cell.embeddings[, "umap_1"], expand.ratio = 0.05)) + ylim(range.expand(mAMs@reductions$umap@cell.embeddings[, "umap_2"], expand.ratio = 0.05)) + coord_cartesian(expand = FALSE)
  Figure.2.FeaturePlot <- Figure.2.FeaturePlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
  Figure.2.FeaturePlot <- ggplot_build(Figure.2.FeaturePlot)
  # highlight specific cell types
  umap_1 <- mAMs@reductions$umap@cell.embeddings[, "umap_1"][mAMs$AMs.Cell.Subtypes %in% meta.data.list$mAMs.Cell.Subtypes[i]]
  Figure.2.FeaturePlot$data[[1]]$size[which(Figure.2.FeaturePlot$data[[1]]$x %in% umap[, "umap_1"])] <- 2
  # remove smaller contour
  Contour <- if(Figure.2A.FeaturePlot.Subtypes[i] == "Gdf15.AMs") c(3, 7, 8) else c(which.max(table(Figure.2.FeaturePlot$data[[2]]$piece)))
  Figure.2.FeaturePlot$data[[2]] <- Figure.2.FeaturePlot$data[[2]][Figure.2.FeaturePlot$data[[2]]$piece == Contour, ] # to make Neu tight
  Figure.2.FeaturePlot.directory <- paste0("Figure.2A.mAMs.", Figure.2A.FeaturePlot.Features[i], ".FeaturePlot.png")
  png(Figure.2.FeaturePlot.directory, width = 1100/2, height = 1000/2, bg = "transparent")
  print(grid::grid.draw(ggplot_gtable(Figure.2.FeaturePlot)))
  dev.off() 
}

# Figure.2B.SankeyDiagram
# hAMs.Cell.Subtypes.DEGs # Table E1. AMs Subtypes DEGs.xlsx
# mAMs.Cell.Subtypes.DEGs # Table E1. AMs Subtypes DEGs.xlsx
mAMs.Subtypes.DEGs <- c() # named vector for DEGs
for (i in meta.data.list$mAMs.Cell.Subtypes) {
  DEGs <- mAMs.Cell.Subtypes.DEGs$gene[mAMs.Cell.Subtypes.DEGs$cluster == i & mAMs.Cell.Subtypes.DEGs$p_val_adj < 0.05]
  names(DEGs) <- rep(i, sum(mAMs.Cell.Subtypes.DEGs$cluster == i & mAMs.Cell.Subtypes.DEGs$p_val_adj < 0.05))
  mAMs.Subtypes.DEGs <- c(mAMs.Subtypes.DEGs, DEGs)
}
Human <- c()
Mouse <- c()
for (i in meta.data.list$hAMs.Cell.Subtypes) { # top 50 of hAMs, map to mAMs
  hAMs.Subtypes.DEGs <- hAMs.Cell.Subtypes.DEGs$gene[hAMs.Cell.Subtypes.DEGs$cluster == i & hAMs.Cell.Subtypes.DEGs$p_val_adj < 0.05]
  Mouse.Subset <- names(mAMs.Subtypes.DEGs[mAMs.Subtypes.DEGs %in% HtoM.melt(hAMs.Subtypes.DEGs[1:50])])
  Mouse <- c(Mouse, Mouse.Subset)
  Human <- c(Human, rep(i, length(Mouse.Subset)))
}
Figure.2.SankeyDiagram.Data <- data.frame("Human" = factor(Human, levels = meta.data.list$hAMs.Cell.Subtypes), "Mouse" = factor(Mouse, levels = meta.data.list$mAMs.Cell.Subtypes))
Cell.Subypes <- c(meta.data.list$hAMs.Cell.Subtypes, meta.data.list$mAMs.Cell.Subtypes)
colors <- c(names(meta.data.list$hAMs.Cell.Subtypes), names(meta.data.list$mAMs.Cell.Subtypes))
SankeyDiagram(Figure.2.SankeyDiagram.Data, max.categories = 100, font.size = 0, node.width = 15, node.padding = 5, colors = colors, link.color = "Source", node.position.automatic = FALSE, label.show.varname = FALSE)

# Figure.2C.Heatmap
Mouse.data <- mAMs[["SCT"]]@data # Humanize mAMs Matrix
Mouse.gene.ref <- c()
Mouse.gene.add <- c()
for (i in 1:nrow(Mouse.data)) {
  Human.Gene <- MtoH(rownames(Mouse.data)[i])
  if (length(Human.Gene) == 1) {
    if (Human.Gene != rownames(Mouse.data)[i]) {
      rownames(Mouse.data)[i] <- Human.Gene
    } else {
      rownames(Mouse.data)[i] <- paste0("REMOVE.", rownames(Mouse.data)[i])
    }
  } else {
    rownames(Mouse.data)[i] <- Human.Gene[1]
    add.number <- length(Human.Gene) - 1
    Mouse.gene.ref <- c(Mouse.gene.ref, rep(Human.Gene[1], add.number))
    Mouse.gene.add <- c(Mouse.gene.add, Human.Gene[-1]) 
  }
}
Mouse.data.add <- Mouse.data[Mouse.gene.ref, ] # add the genes with one to many match
rownames(Mouse.data.add) <- Mouse.gene.add
Mouse.data <- rbind(Mouse.data, Mouse.data.add) # Humanized mAMs Matrix
# downsample Mouse.data to max 1000 cells per AMs.Cell.Subtypes
downsample <- c()
for (i in meta.data.list$mAMs.Cell.Subtypes) {
  if (table(mAMs$AMs.Cell.Subtypes)[i] > 1000) {
    downsample <- c(downsample, which(sample.boolean(mAMs$AMs.Cell.Subtypes == i, 1000)))
  } else {
    downsample <- c(downsample, which(mAMs$AMs.Cell.Subtypes == i))
  }
}
Mouse.data.ds <- Mouse.data[, downsample]
# downsample hAMs to max 1000 cells per AMs.Cell.Subtypes
Idents(hAMs) <- "AMs.Cell.Subtypes"
hAMs.ds <- subset(hAMs, downsample = 1000)

# MetaNeighbor
# mAMs.ds
Gene.Matrix_M <- Mouse.data.ds
colData_M <- mAMs@meta.data[downsample, c("orig.species", "AMs.Cell.Subtypes")]
colData_M <- tibble::rownames_to_column(colData_M, "Cell_ID")
# hAMs.ds
Gene.Matrix_H <- hAMs.ds[["SCT"]]@data
colData_H <- hAMs.ds@meta.data[, c("orig.species", "AMs.Cell.Subtypes")]
colnames(colData_H) <- gsub("AMs.Cell.Subtypes", "AMs.Cell.Subtypes", colnames(colData_H))
colData_H <- tibble::rownames_to_column(colData_H, "Cell_ID")

Gene.Matrix.Gene.Order_M_H <- intersect(rownames(Gene.Matrix_M), rownames(Gene.Matrix_H))
Gene.Matrix.Combined_M_H <- cbind(Gene.Matrix_M[Gene.Matrix.Gene.Order_M_H, ], Gene.Matrix_H[Gene.Matrix.Gene.Order_M_H, ])
colData.Combined_M_H <- rbind(colData_M, colData_H)
SE.Combined_M_H <- SummarizedExperiment(assays = SimpleList(gene_matrix = Gene.Matrix.Combined_M_H), colData = colData.Combined_M_H)
Var.Genes_M_H <- variableGenes(dat = SE.Combined_M_H, exp_labels = SE.Combined_M_H$orig.species)
AMs.Cell.Subtypes.NV_M_H <- MetaNeighborUS(var_genes = Var.Genes_M_H, dat = SE.Combined_M_H, study_id = SE.Combined_M_H$orig.species, cell_type = SE.Combined_M_H$AMs.Cell.Subtypes)
Figure.2.Heatmap.Data <- AMs.Cell.Subtypes.NV_M_H[grep("human", rownames(AMs.Cell.Subtypes.NV_M_H)), grep("mouse", colnames(AMs.Cell.Subtypes.NV_M_H))]
dimnames(Figure.2.Heatmap.Data) <- list(str_split_i(rownames(Figure.2.Heatmap.Data), "\\|", 2), str_split_i(colnames(Figure.2.Heatmap.Data), "\\|", 2))
Figure.2.Heatmap.Data <- Figure.2.Heatmap.Data[match(meta.data.list$hAMs.Cell.Subtypes, rownames(Figure.2.Heatmap.Data)), match(meta.data.list$mAMs.Cell.Subtypes, colnames(Figure.2.Heatmap.Data))]
col_fun <- colorRamp2(range(Figure.2.Heatmap.Data), c("white", "#990F0F"))
left_annotation <- rowAnnotation(hAMs.Cell.Subtypes = meta.data.list$hAMs.Cell.Subtypes, col = list(hAMs.Cell.Subtypes = setNames(names(meta.data.list$hAMs.Cell.Subtypes), meta.data.list$hAMs.Cell.Subtypes)), simple_anno_size = unit(0.25, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
bottom_annotation <- columnAnnotation(mAMs.Cell.Subtypes = meta.data.list$mAMs.Cell.Subtypes, col = list(mAMs.Cell.Subtypes = setNames(names(meta.data.list$mAMs.Cell.Subtypes), meta.data.list$mAMs.Cell.Subtypes)), simple_anno_size = unit(0.25, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.2.Heatmap <- Heatmap(Figure.2.Heatmap.Data, cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE, col = col_fun, row_names_gp = gpar(fontsize = 21), column_names_gp = gpar(fontsize = 21), column_names_rot = 45, row_names_side = "left", show_heatmap_legend = FALSE, 
                            width = unit(7, "in"), height = unit(7, "in"), left_annotation = left_annotation, bottom_annotation = bottom_annotation)
Figure.2.Heatmap.directory <- paste0("Figure.2C.", "Heatmap.pdf")
pdf(Figure.2.Heatmap.directory, width = 9.5, height = 9.5)
print(Figure.2.Heatmap)
dev.off() 

# Figure.2D.hAMs.DimPlot
Figure.2.DimPlot <- DimPlot(hAMs, group.by = "orig.tissue", cols = "grey80", pt.size = 2) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range1.01(hAMs@reductions$umap@cell.embeddings[, "umap_1"])) + ylim(range1.01(hAMs@reductions$umap@cell.embeddings[, "umap_2"])) + coord_cartesian(expand = FALSE)
Figure.2.DimPlot <- Figure.2.DimPlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
Figure.2.DimPlot.directory <- paste0("Figure.2D.hAMs", ".DimPlot.png")
png(Figure.2.DimPlot.directory, width = 1100, height = 1000, bg = "transparent")
print(Figure.2.DimPlot)
dev.off() 
# Figure.2D.mAMs.DimPlot
Figure.2.DimPlot <- DimPlot(mAMs, group.by = "orig.tissue", cols = "grey80", pt.size = 2) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range1.01(mAMs@reductions$umap@cell.embeddings[, "umap_1"])) + ylim(range1.01(mAMs@reductions$umap@cell.embeddings[, "umap_2"])) + coord_cartesian(expand = FALSE)
Figure.2.DimPlot <- Figure.2.DimPlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
Figure.2.DimPlot.directory <- paste0("Figure.2D.mAMs", ".DimPlot.png")
png(Figure.2.DimPlot.directory, width = 1100, height = 1000, bg = "transparent")
print(Figure.2.DimPlot)
dev.off() 
# Figure.2D.DimPlot
# Humanized mAMs SeuratObject
mAMs.h <- CreateSeuratObject(counts = Mouse.data, assay = "SCT", meta.data = mAMs@meta.data) 
mAMs.h[["SCT"]]$data <- mAMs.h[["SCT"]]$counts
mAMs.h <- mAMs.h %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()
# Mapping.AMs
AMs_anchors <- FindTransferAnchors(reference = hAMs, query = mAMs.h, dims = 1:50, reference.reduction = "pca", normalization.method = "SCT")
predictions.Subtypes <- TransferData(anchorset = AMs_anchors, refdata = hAMs$AMs.Cell.Subtypes, dims = 1:50)
mAMs.h@meta.data <- cbind(mAMs.h@meta.data, predictions.Subtypes)
mAMs.h$predicted.id <- factor(mAMs.h$predicted.id, levels = meta.data.list$hAMs.Cell.Subtypes)
mAMs.h <- MapQuery(anchorset = AMs_anchors, reference = hAMs, query = mAMs.h, refdata = list(orig.tissue = "orig.tissue"), reference.reduction = "pca", reduction.model = "umap", new.reduction.name = "umap_mapping")
mAMs.h$hAMs.Cell.Subtypes <- mAMs.h$predicted.id
Figure.2.DimPlot <- DimPlot(mAMs.h, group.by = "hAMs.Cell.Subtypes", reduction = "ref.umap", cols = names(meta.data.list$hAMs.Cell.Subtypes), pt.size = 2, shuffle = TRUE) + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + 
  xlim(range1.01(mAMs.h@reductions$ref.umap@cell.embeddings[, "refUMAP_1"])) + ylim(range1.01(mAMs.h@reductions$ref.umap@cell.embeddings[, "refUMAP_2"])) + coord_cartesian(expand = FALSE)
Figure.2.DimPlot <- Figure.2.DimPlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
Figure.2.DimPlot.directory <- paste0("Figure.2D", ".DimPlot.png")
png(Figure.2.DimPlot.directory, width = 1100, height = 1000, bg = "transparent")
print(Figure.2.DimPlot)
dev.off() 

# Figure.2E.DonutChart
Figure.2.DonutChart.Data <- as.data.frame(t(table(mAMs.h$hAMs.Cell.Subtypes, mAMs.h$AMs.Cell.Subtypes)))
Figure.2.DonutChart.Data$mAMs.Cell.Subtypes <- NaturalLevel(rownames(Figure.2.DonutChart.Data))
reticulate::py_run_string("import sys")
for (i in 1:(length(colnames(Figure.2.DonutChart.Data))-1)){
  Figure.2.DonutChart <- Figure.2.DonutChart.Data %>% plot_ly(labels = ~mAMs.Cell.Subtypes, values = ~Figure.2.DonutChart.Data[[colnames(Figure.2.DonutChart.Data)[i]]], marker = list(colors = names(meta.data.list$mAMs.Cell.Subtypes)), sort = FALSE, direction = "clockwise", automargin = FALSE, textposition = "inside", textfont = list(size = 100))
  Figure.2.DonutChart <- Figure.2.DonutChart %>% add_pie(hole = 0.5) %>% layout(showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  Figure.2.DonutChart.directory <- paste0("Figure.2E", ".DonutChart.", colnames(Figure.2.DonutChart.Data)[i], ".pdf")
  save_image(Figure.2.DonutChart, Figure.2.DonutChart.directory, height = 1000, width = 1000, scale = 5)
}

# Figure.E1A.hAMs.DotPlot.Data
Figure.E1.DotPlot.Features <- hAMs.Cell.Subtypes.DEGs
Figure.E1.DotPlot.Features <- Figure.E1.DotPlot.Features %>% group_by(cluster) %>% slice_head(n = 5)
Figure.E1.DotPlot.Features <- setNames(Figure.E1.DotPlot.Features$gene, Figure.E1.DotPlot.Features$cluster)
Figure.E1.DotPlot.Features <- Figure.E1.DotPlot.Features[!duplicated(Figure.E1.DotPlot.Features)]
Figure.E1.DotPlot.Features <- lapply(split(Figure.E1.DotPlot.Features, names(Figure.E1.DotPlot.Features)), unname)
Figure.E1.DotPlot.Features <- Figure.E1.DotPlot.Features[meta.data.list$hAMs.Cell.Subtypes]
Figure.E1.DotPlot <- DotPlot(hAMs, group.by = "AMs.Cell.Subtypes", assay = "SCT", features = Figure.E1.DotPlot.Features, cols = c("#E5E5E5", "#990F0F"), dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), strip.text = element_text(size = 15), legend.position = "none")
Figure.E1.DotPlot.directory <- paste0("Figure.E1A.", "hAMs", ".DotPlot.pdf")
pdf(Figure.E1.DotPlot.directory, width = 20, height = 4.5)
print(Figure.E1.DotPlot)
dev.off() 
# Figure.E2A.mAMs.DotPlot.Data
Figure.E2.DotPlot.Features <- mAMs.Cell.Subtypes.DEGs
Figure.E2.DotPlot.Features <- Figure.E2.DotPlot.Features %>% group_by(cluster) %>% slice_head(n = 5)
Figure.E2.DotPlot.Features <- setNames(Figure.E2.DotPlot.Features$gene, Figure.E2.DotPlot.Features$cluster)
Figure.E2.DotPlot.Features <- Figure.E2.DotPlot.Features[!duplicated(Figure.E2.DotPlot.Features)]
Figure.E2.DotPlot.Features <- lapply(split(Figure.E2.DotPlot.Features, names(Figure.E2.DotPlot.Features)), unname)
Figure.E2.DotPlot.Features <- Figure.E2.DotPlot.Features[meta.data.list$mAMs.Cell.Subtypes]
Figure.E2.DotPlot <- DotPlot(mAMs, group.by = "AMs.Cell.Subtypes", assay = "SCT", features = Figure.E2.DotPlot.Features, cols = c("#E5E5E5", "#990F0F"), dot.scale = 6) + RotatedAxis() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), strip.text = element_text(size = 15), legend.position = "none")
Figure.E2.DotPlot.directory <- paste0("Figure.E2A.", "mAMs", ".DotPlot.pdf")
pdf(Figure.E2.DotPlot.directory, width = 20, height = 4.5)
print(Figure.E2.DotPlot)
dev.off() 

# Figure.E1B.hAMs.UpSetPlot
Figure.E1.DotPlot.Features <- hAMs.Cell.Subtypes.DEGs
Figure.E1.UpSetPlot.Data <- list()
for (i in 1:length(meta.data.list$hAMs.Cell.Subtypes)){
  Figure.E1.UpSetPlot.Data[[i]] <- Figure.E1.DotPlot.Features$gene[Figure.E1.DotPlot.Features$cluster == meta.data.list$hAMs.Cell.Subtypes[i] & Figure.E1.DotPlot.Features$p_val_adj < 0.05]
}
Figure.E1.UpSetPlot.Data.list <- Figure.E1.UpSetPlot.Data # for Figure.E1C.hAMs.Heatmap
Figure.E1.UpSetPlot.Data <- make_comb_mat(Figure.E1.UpSetPlot.Data, mode = "intersect")
top40 <- order(comb_size(Figure.E1.UpSetPlot.Data), decreasing = TRUE)[1:40] # top 40 overlap
comb_name_order <- character(0) # all Subtype comb_name
for (i in 1:length(meta.data.list$hAMs.Cell.Subtypes)) {
  vec <- rep("0", length(meta.data.list$hAMs.Cell.Subtypes))
  vec[i] <- "1"
  comb_name_order <- c(comb_name_order, paste(vec, collapse = ""))
}
Figure.E1.UpSetPlot.Data <- Figure.E1.UpSetPlot.Data[union(which(comb_name(Figure.E1.UpSetPlot.Data) %in% comb_name_order), top40)]
comb_col <- ifelse(comb_name(Figure.E1.UpSetPlot.Data) %in% comb_name_order, comb_name(Figure.E1.UpSetPlot.Data), "grey60")
comb_col <- mapvalues(comb_col, from = comb_name_order, to = gsub("#C7C7C7", "grey60", names(meta.data.list$hAMs.Cell.Subtypes)))
Figure.E1.UpSetPlot.Data <- t(Figure.E1.UpSetPlot.Data) # transpose to horizontal
Figure.E1.UpSetPlot <- UpSet(Figure.E1.UpSetPlot.Data, comb_order = order(comb_size(Figure.E1.UpSetPlot.Data), decreasing = TRUE), pt_size = unit(4, "mm"), column_names_gp = gpar(fontsize = 10), heatmap_width = unit(3, "in"), heatmap_height = unit(8, "in"), comb_col = comb_col, 
                             top_annotation = upset_top_annotation(Figure.E1.UpSetPlot.Data, show_annotation_name = FALSE, axis_param = list(side = "left", gp = gpar(fontsize = 8)), width = unit(2, "cm"), gp = gpar(fill = names(meta.data.list$hAMs.Cell.Subtypes))),
                             right_annotation = upset_right_annotation(Figure.E1.UpSetPlot.Data, show_annotation_name = FALSE, axis_param = list(side = "top", gp = gpar(fontsize = 8)), height = unit(3.75, "cm"), anno_size = unit(2, "in"), gp = gpar(fill = gsub("grey60", "grey90", comb_col))))
Figure.E1.UpSetPlot.directory <- paste0("Figure.E1B.", "hAMs", ".UpSetPlot.pdf")
pdf(Figure.E1.UpSetPlot.directory, width = 3.5, height = 8.5)
print(Figure.E1.UpSetPlot)
dev.off()
# Figure.E2B.mAMs.UpSetPlot
Figure.E2.DotPlot.Features <- mAMs.Cell.Subtypes.DEGs
Figure.E2.UpSetPlot.Data <- list()
for (i in 1:length(meta.data.list$mAMs.Cell.Subtypes)){
  Figure.E2.UpSetPlot.Data[[i]] <- Figure.E2.DotPlot.Features$gene[Figure.E2.DotPlot.Features$cluster == meta.data.list$mAMs.Cell.Subtypes[i] & Figure.E2.DotPlot.Features$p_val_adj < 0.05]
}
Figure.E2.UpSetPlot.Data.list <- Figure.E2.UpSetPlot.Data # for Figure.E2C.hAMs.Heatmap
Figure.E2.UpSetPlot.Data <- make_comb_mat(Figure.E2.UpSetPlot.Data, mode = "intersect")
top40 <- order(comb_size(Figure.E2.UpSetPlot.Data), decreasing = TRUE)[1:40] # top 40 overlap
comb_name_order <- character(0) # all Subtype comb_name
for (i in 1:length(meta.data.list$mAMs.Cell.Subtypes)) {
  vec <- rep("0", length(meta.data.list$mAMs.Cell.Subtypes))
  vec[i] <- "1"
  comb_name_order <- c(comb_name_order, paste(vec, collapse = ""))
}
Figure.E2.UpSetPlot.Data <- Figure.E2.UpSetPlot.Data[union(which(comb_name(Figure.E2.UpSetPlot.Data) %in% comb_name_order), top40)]
comb_col <- ifelse(comb_name(Figure.E2.UpSetPlot.Data) %in% comb_name_order, comb_name(Figure.E2.UpSetPlot.Data), "grey60")
comb_col <- mapvalues(comb_col, from = comb_name_order, to = gsub("#C7C7C7", "grey60", names(meta.data.list$mAMs.Cell.Subtypes)))
Figure.E2.UpSetPlot.Data <- t(Figure.E2.UpSetPlot.Data) # transpose to horizontal
Figure.E2.UpSetPlot <- UpSet(Figure.E2.UpSetPlot.Data, comb_order = order(comb_size(Figure.E2.UpSetPlot.Data), decreasing = TRUE), pt_size = unit(4, "mm"), column_names_gp = gpar(fontsize = 10), heatmap_width = unit(3, "in"), heatmap_height = unit(8, "in"), comb_col = comb_col, 
                             top_annotation = upset_top_annotation(Figure.E2.UpSetPlot.Data, show_annotation_name = FALSE, axis_param = list(side = "left", gp = gpar(fontsize = 8)), width = unit(2, "cm"), gp = gpar(fill = names(meta.data.list$mAMs.Cell.Subtypes))),
                             right_annotation = upset_right_annotation(Figure.E2.UpSetPlot.Data, show_annotation_name = FALSE, axis_param = list(side = "top", gp = gpar(fontsize = 8)), height = unit(3.75, "cm"), anno_size = unit(2, "in"), gp = gpar(fill = gsub("grey60", "grey90", comb_col))))
Figure.E2.UpSetPlot.directory <- paste0("Figure.E2B.", "mAMs", ".UpSetPlot.pdf")
pdf(Figure.E2.UpSetPlot.directory, width = 3.5, height = 8.5)
print(Figure.E2.UpSetPlot)
dev.off()

# Figure.E1C.hAMs.Heatmap
AMs.bg <- rownames(hAMs)[rowSums(hAMs[["SCT"]]@data) > 0]
Figure.E1.UpSetPlot.Data <- Figure.E1.UpSetPlot.Data.list
# AMs.Cell.Subtypes.topGO.Terms
AMs.Cell.Subtypes.GOs <- list()
AMs.Cell.Subtypes.topGO.Terms <- c()
for (i in 1:length(meta.data.list$hAMs.Cell.Subtypes)){
  AMs.DEG <- Figure.E1.UpSetPlot.Data[[i]]
  if (length(AMs.DEG) != 0) {
    AMs.Cell.Subtypes.GOs[[i]] <- topGOtable(AMs.DEG, AMs.bg, mapping = "org.Hs.eg.db", ontology = "BP", writeOutput = FALSE, topTablerows = 500, outputFile = AMs.Cell.Subtypes.GOs.directory)
  }
  AMs.Cell.Subtypes.topGO.Terms <- c(AMs.Cell.Subtypes.topGO.Terms, AMs.Cell.Subtypes.GOs[[i]]$GO.ID[1:5])
}
names(AMs.Cell.Subtypes.topGO.Terms) <- rep(meta.data.list$hAMs.Cell.Subtypes, each = 5)
# Table E2. Human AMs Subsets GOs.xlsx
# Figure.E1.Heatmap.Data
AMs.Cell.Subtypes.topGO.Terms <- AMs.Cell.Subtypes.topGO.Terms[!duplicated(AMs.Cell.Subtypes.topGO.Terms)]
Figure.E1.Heatmap.Data <- data.frame(matrix(ncol = 0, nrow = length(AMs.Cell.Subtypes.topGO.Terms)))
AMs.Cell.Subtypes.targetGOs <- list()
for (i in 1:length(meta.data.list$hAMs.Cell.Subtypes)){
  AMs.DEG <- Figure.E1.UpSetPlot.Data[[i]]
  if (length(AMs.DEG) != 0) {
    AMs.Cell.Subtypes.targetGOs[[i]] <- topGOtable.targetGO(AMs.DEG, AMs.bg, mapping = "org.Hs.eg.db", ontology = "BP", TargetGO = AMs.Cell.Subtypes.topGO.Terms, writeOutput = FALSE, topTablerows = 500, outputFile = AMs.Cell.Subtypes.targetGOs.directory)
    setnafill(AMs.Cell.Subtypes.targetGOs[[i]], "nocb", cols = "p.value_elim")
    AMs.Cell.Subtypes.targetGOs.List <- AMs.Cell.Subtypes.targetGOs[[i]][match(AMs.Cell.Subtypes.topGO.Terms, AMs.Cell.Subtypes.targetGOs[[i]]$GO.ID), ]
    Figure.E1.Heatmap.Data <- cbind(Figure.E1.Heatmap.Data, AMs.Cell.Subtypes.targetGOs.List$p.value_elim)
  }
}
rownames(Figure.E1.Heatmap.Data) <- as.vector(AMs.Cell.Subtypes.topGO.Terms)
colnames(Figure.E1.Heatmap.Data) <- meta.data.list$hAMs.Cell.Subtypes
Figure.E1.Heatmap.Data <- data.matrix(-log10(Figure.E1.Heatmap.Data), rownames.force = TRUE)
rownames(Figure.E1.Heatmap.Data) <- AMs.Cell.Subtypes.targetGOs[[1]]$Term[match(AMs.Cell.Subtypes.topGO.Terms, AMs.Cell.Subtypes.targetGOs[[1]]$GO.ID)]
col_fun = colorRamp2(range(Figure.E1.Heatmap.Data), c("white", "#990F0F"))
AMs.Cell.Subtypes <- names(AMs.Cell.Subtypes.topGO.Terms)
colors <- names(meta.data.list$hAMs.Cell.Subtypes)
cells <- unique(AMs.Cell.Subtypes)
column_annotation <- columnAnnotation(Cell.Subtypes = cells, col = list(Cell.Subtypes = setNames(colors, cells)), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
right_annotation <- rowAnnotation(AMs.Cell.Subtypes = AMs.Cell.Subtypes, col = list(AMs.Cell.Subtypes = setNames(colors, cells)), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.E1.Heatmap <- Heatmap(Figure.E1.Heatmap.Data, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE, column_names_gp = gpar(fontsize = 20), column_names_side = "top", column_names_rot = 35, row_names_gp = gpar(fontsize = 20), show_heatmap_legend = FALSE, 
                             heatmap_width = unit(10, "in"), heatmap_height = unit(15.4, "in"), col = col_fun, right_annotation = right_annotation, top_annotation = column_annotation)
Figure.E1.Heatmap.directory <- paste0("Figure.E1C.", "hAMs", ".Heatmap.png")
png(Figure.E1.Heatmap.directory, res = 150, width = 2600, height = 2400, bg = "transparent")
draw(Figure.E1.Heatmap, background = "transparent", padding = unit(c(0, 0, 0, 17), "cm"))
dev.off()
# Figure.E2C.mAMs.Heatmap
AMs.bg <- rownames(mAMs)[rowSums(mAMs[["SCT"]]@data) > 0]
Figure.E2.UpSetPlot.Data <- Figure.E2.UpSetPlot.Data.list
# AMs.Cell.Subtypes.topGO.Terms
AMs.Cell.Subtypes.GOs <- list()
AMs.Cell.Subtypes.topGO.Terms <- c()
for (i in 1:length(meta.data.list$mAMs.Cell.Subtypes)){
  AMs.DEG <- Figure.E2.UpSetPlot.Data[[i]]
  if (length(AMs.DEG) != 0) {
    AMs.Cell.Subtypes.GOs[[i]] <- topGOtable(AMs.DEG, AMs.bg, mapping = "org.Mm.eg.db", ontology = "BP", writeOutput = FALSE, topTablerows = 500, outputFile = AMs.Cell.Subtypes.GOs.directory)
  }
  AMs.Cell.Subtypes.topGO.Terms <- c(AMs.Cell.Subtypes.topGO.Terms, AMs.Cell.Subtypes.GOs[[i]]$GO.ID[1:5])
}
names(AMs.Cell.Subtypes.topGO.Terms) <- rep(meta.data.list$mAMs.Cell.Subtypes, each = 5)
# Table E3. Mouse AMs Subsets GOs.xlsx
# Figure.E2.Heatmap.Data
AMs.Cell.Subtypes.topGO.Terms <- AMs.Cell.Subtypes.topGO.Terms[!duplicated(AMs.Cell.Subtypes.topGO.Terms)]
Figure.E2.Heatmap.Data <- data.frame(matrix(ncol = 0, nrow = length(AMs.Cell.Subtypes.topGO.Terms)))
AMs.Cell.Subtypes.targetGOs <- list()
for (i in 1:length(meta.data.list$mAMs.Cell.Subtypes)){
  AMs.DEG <- Figure.E2.UpSetPlot.Data[[i]]
  if (length(AMs.DEG) != 0) {
    AMs.Cell.Subtypes.targetGOs[[i]] <- topGOtable.targetGO(AMs.DEG, AMs.bg, mapping = "org.Mm.eg.db", ontology = "BP", TargetGO = AMs.Cell.Subtypes.topGO.Terms, writeOutput = FALSE, topTablerows = 500, outputFile = AMs.Cell.Subtypes.targetGOs.directory)
    setnafill(AMs.Cell.Subtypes.targetGOs[[i]], "nocb", cols = "p.value_elim")
    AMs.Cell.Subtypes.targetGOs.List <- AMs.Cell.Subtypes.targetGOs[[i]][match(AMs.Cell.Subtypes.topGO.Terms, AMs.Cell.Subtypes.targetGOs[[i]]$GO.ID), ]
    Figure.E2.Heatmap.Data <- cbind(Figure.E2.Heatmap.Data, AMs.Cell.Subtypes.targetGOs.List$p.value_elim)
  }
}
rownames(Figure.E2.Heatmap.Data) <- as.vector(AMs.Cell.Subtypes.topGO.Terms)
colnames(Figure.E2.Heatmap.Data) <- meta.data.list$mAMs.Cell.Subtypes
Figure.E2.Heatmap.Data <- data.matrix(-log10(Figure.E2.Heatmap.Data), rownames.force = TRUE)
rownames(Figure.E2.Heatmap.Data) <- AMs.Cell.Subtypes.targetGOs[[1]]$Term[match(AMs.Cell.Subtypes.topGO.Terms, AMs.Cell.Subtypes.targetGOs[[1]]$GO.ID)]
col_fun = colorRamp2(range(Figure.E2.Heatmap.Data), c("white", "#990F0F"))
AMs.Cell.Subtypes <- names(AMs.Cell.Subtypes.topGO.Terms)
colors <- names(meta.data.list$mAMs.Cell.Subtypes)
cells <- unique(AMs.Cell.Subtypes)
column_annotation <- columnAnnotation(Cell.Subtypes = cells, col = list(Cell.Subtypes = setNames(colors, cells)), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
right_annotation <- rowAnnotation(AMs.Cell.Subtypes = AMs.Cell.Subtypes, col = list(AMs.Cell.Subtypes = setNames(colors, cells)), simple_anno_size = unit(0.4, "cm"), show_legend = FALSE, show_annotation_name = FALSE)
Figure.E2.Heatmap <- Heatmap(Figure.E2.Heatmap.Data, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = TRUE, column_names_gp = gpar(fontsize = 20), column_names_side = "top", column_names_rot = 35, row_names_gp = gpar(fontsize = 20), show_heatmap_legend = FALSE, 
                             heatmap_width = unit(10, "in"), heatmap_height = unit(15.4, "in"), col = col_fun, right_annotation = right_annotation, top_annotation = column_annotation)
Figure.E2.Heatmap.directory <- paste0("Figure.E2C.", "mAMs", ".Heatmap.png")
png(Figure.E2.Heatmap.directory, res = 150, width = 2600, height = 2400, bg = "transparent")
draw(Figure.E2.Heatmap, background = "transparent", padding = unit(c(0, 0, 0, 17), "cm"))
dev.off()

# Figure.E3A.mAMs.FeaturePlot & Figure.E3B.mAMs.FeaturePlot
Figure.E3.FeaturePlot.Features <- c("Ccl9", "Ccl9", "C3", "Relb", "Mmp14", "Clec4e", "Cxcl16", "Tlr2", "AY036118", "Gm42418", "Hspa1a", "Hspa1b", "H2-Aa", "H2-Ab1", "Car4", "Hivep3", "Tppp3", "Crip1")
Figure.E3.FeaturePlot.Subtypes <- c(rep("Ccl9.AMs", 7), "CK.AMs", rep(as.vector(meta.data.list$mAMs.Cell.Subtypes[7:11]), each = 2))
for (i in 1:length(Figure.E3.FeaturePlot.Features)) {
  umap <- mAMs@reductions$umap@cell.embeddings[mAMs$AMs.Cell.Subtypes %in% Figure.E3.FeaturePlot.Subtypes[i], ]
  Figure.E3.FeaturePlot <- FeaturePlot(mAMs, features = Figure.E3.FeaturePlot.Features[i], raster = FALSE, order = TRUE, pt.size = 0.5) + scale_colour_gradient(low = "#E5E5E5", high = "#990F0F") + 
    geom_density_2d(data = data.frame(umap), aes(x = umap_1, y = umap_2), bins = 15, linewidth = 2.5, alpha = 2, colour = "gray20") + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + # bins = ifelse(i == 2, ifelse(i == 5, 25, 8), 15) to make Mono tight, DCs bigger
    xlim(range.expand(mAMs@reductions$umap@cell.embeddings[, "umap_1"], expand.ratio = 0.05)) + ylim(range.expand(mAMs@reductions$umap@cell.embeddings[, "umap_2"], expand.ratio = 0.05)) + coord_cartesian(expand = FALSE)
  Figure.E3.FeaturePlot <- Figure.E3.FeaturePlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
  Figure.E3.FeaturePlot <- ggplot_build(Figure.E3.FeaturePlot)
  # highlight specific cell types
  umap_1 <- mAMs@reductions$umap@cell.embeddings[, "umap_1"][mAMs$AMs.Cell.Subtypes %in% meta.data.list$AMs.Cell.Subtypes[i]]
  Figure.E3.FeaturePlot$data[[1]]$size[which(Figure.E3.FeaturePlot$data[[1]]$x %in% umap[, "umap_1"])] <- 2
  # remove smaller contour
  Contour <- if(Figure.E3.FeaturePlot.Subtypes[i] == "Gdf15.AMs") c(3, 7, 8) else c(which.max(table(Figure.E3.FeaturePlot$data[[2]]$piece)))
  Figure.E3.FeaturePlot$data[[2]] <- Figure.E3.FeaturePlot$data[[2]][Figure.E3.FeaturePlot$data[[2]]$piece == Contour, ] # to make Neu tight
  Figure.E3.FeaturePlot.directory <- paste0("Figure.E3.mAMs.", Figure.E3.FeaturePlot.Features[i], ".FeaturePlot.png")
  png(Figure.E3.FeaturePlot.directory, width = 1100/2, height = 1000/2, bg = "transparent")
  print(grid::grid.draw(ggplot_gtable(Figure.E3.FeaturePlot)))
  dev.off() 
}
# Figure.E3A.hAMs.FeaturePlot & Figure.E3B.hAMs.FeaturePlot
# Figure.E3.FeaturePlot.Features <- c("CCL15", "CCL23", "C3", "RELB", "MMP14", "CLEC4E", "CXCL16", "TLR2", "AY036118", "Gm42418", "HSPA1A", "HSPA1B", "HLA-DQA2", "HLA-DQB1", "CA4", "HIVEP3", "TPPP3", "CRIP1")
Figure.E3.FeaturePlot.Features <- c("CCL15", "CCL23", "C3", "RELB", "MMP14", "CLEC4E", "CXCL16", "TLR2", "HSPA1A", "HSPA1B", "HLA-DQA2", "HLA-DQB1", "HIVEP3", "TPPP3", "CRIP1") # "AY036118", "Gm42418", "CA4" do not have expression
Figure.E3.FeaturePlot.Subtypes <- list(c("CK1.AMs", "CK2.AMs", "CXCL5.AMs", "CCL18.AMs"), c("CK1.AMs", "CK2.AMs", "CXCL5.AMs", "CCL18.AMs"), c("CK1.AMs", "CK2.AMs", "CXCL5.AMs", "CCL18.AMs"), 
                                        c("CK1.AMs", "CK2.AMs"), c("CK1.AMs", "CK2.AMs"), c("CK1.AMs", "CK2.AMs"), c(), c("CK1.AMs", "CK2.AMs"), 
                                        c(), c(), c(), c(), c(), c(), c())
for (i in 1:length(Figure.E3.FeaturePlot.Features)) {
  umap <- hAMs@reductions$umap@cell.embeddings[hAMs$AMs.Cell.Subtypes %in% Figure.E3.FeaturePlot.Subtypes[[i]], ]
  Figure.E3.FeaturePlot <- FeaturePlot(hAMs, features = Figure.E3.FeaturePlot.Features[i], raster = FALSE, order = TRUE, pt.size = 0.5) + scale_colour_gradient(low = "#E5E5E5", high = "#990F0F") + 
    geom_density_2d(data = data.frame(umap), aes(x = umap_1, y = umap_2), bins = 15, linewidth = 2.5, alpha = 2, colour = "gray20") + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + # bins = ifelse(i == 2, ifelse(i == 5, 25, 8), 15) to make Mono tight, DCs bigger
    xlim(range.expand(hAMs@reductions$umap@cell.embeddings[, "umap_1"], expand.ratio = 0.05)) + ylim(range.expand(hAMs@reductions$umap@cell.embeddings[, "umap_2"], expand.ratio = 0.05)) + coord_cartesian(expand = FALSE)
  Figure.E3.FeaturePlot <- Figure.E3.FeaturePlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
  Figure.E3.FeaturePlot <- ggplot_build(Figure.E3.FeaturePlot)
  # highlight specific cell types
  umap_1 <- hAMs@reductions$umap@cell.embeddings[, "umap_1"][hAMs$AMs.Cell.Subtypes %in% meta.data.list$AMs.Cell.Subtypes[i]]
  Figure.E3.FeaturePlot$data[[1]]$size[which(Figure.E3.FeaturePlot$data[[1]]$x %in% umap[, "umap_1"])] <- 2
  # remove smaller contour
  Contour <- which.max(table(Figure.E3.FeaturePlot$data[[2]]$piece))
  Figure.E3.FeaturePlot$data[[2]] <- Figure.E3.FeaturePlot$data[[2]][Figure.E3.FeaturePlot$data[[2]]$piece == Contour, ]
  Figure.E3.FeaturePlot.directory <- paste0("Figure.E3.hAMs.", Figure.E3.FeaturePlot.Features[i], ".FeaturePlot.png")
  png(Figure.E3.FeaturePlot.directory, width = 1100/2, height = 1000/2, bg = "transparent")
  print(grid::grid.draw(ggplot_gtable(Figure.E3.FeaturePlot)))
  dev.off() 
}
# for "AY036118", "Gm42418", "CA4", which do not have expression
umap <- hAMs@reductions$umap@cell.embeddings[hAMs$AMs.Cell.Subtypes %in% Figure.E3.FeaturePlot.Subtypes[[i]], ]
Figure.E3.DimPlot <- DimPlot(hAMs, group.by = "orig.tissue", cols = "lightgrey", raster = FALSE, order = TRUE, pt.size = 0.5) + 
  geom_density_2d(data = data.frame(umap), aes(x = umap_1, y = umap_2), bins = 15, linewidth = 2.5, alpha = 2, colour = "gray20") + theme_void() + theme(plot.title = element_blank(), legend.position = "none") + # bins = ifelse(i == 2, ifelse(i == 5, 25, 8), 15) to make Mono tight, DCs bigger
  xlim(range.expand(hAMs@reductions$umap@cell.embeddings[, "umap_1"], expand.ratio = 0.05)) + ylim(range.expand(hAMs@reductions$umap@cell.embeddings[, "umap_2"], expand.ratio = 0.05)) + coord_cartesian(expand = FALSE)
Figure.E3.DimPlot <- Figure.E3.DimPlot & theme(panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) # transparent
Figure.E3.DimPlot <- ggplot_build(Figure.E3.DimPlot)
Figure.E3.DimPlot.directory <- paste0("Figure.E3.hAMs.DimPlot.png")
png(Figure.E3.DimPlot.directory, width = 1100/2, height = 1000/2, bg = "transparent")
print(grid::grid.draw(ggplot_gtable(Figure.E3.DimPlot)))
dev.off() 

# Figure.E3C.hAMs.DonutChart
hAMs.Cyc <- hAMs[, hAMs$AMs.Cell.Subtypes == "Cyc.AMs"]
hAMs.Cyc <- CellCycleScoring(hAMs.Cyc, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
Figure.E3.DonutChart.Data <- data.frame("Phase" = names(table(hAMs.Cyc$Phase)), "Number" = as.numeric(table(hAMs.Cyc$Phase)))
Figure.E3.DonutChart <- Figure.E3.DonutChart.Data %>% plot_ly(labels = ~Phase, values = ~Number, marker = list(colors = col2hex(c("grey30", "grey60", "grey95"))), sort = FALSE, direction = "clockwise", automargin = FALSE, textposition = "inside", textfont = list(size = 200)) # plotly default colors
Figure.E3.DonutChart <- Figure.E3.DonutChart %>% add_pie(hole = 0.5) %>% layout(showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), plot_bgcolor = "transparent")
Figure.E3.DonutChart.directory <- paste0("Figure.E3C.hAMs", ".DonutChart.pdf")
save_image(Figure.E3.DonutChart, Figure.E3.DonutChart.directory, height = 1000, width = 1000, scale = 5)
# Figure.E3C.mAMs.DonutChart
mAMs.Cyc <- mAMs[, mAMs$AMs.Cell.Subtypes == "Cyc.AMs"]
mAMs.Cyc <- CellCycleScoring(mAMs.Cyc, s.features = HtoM.melt(cc.genes$s.genes), g2m.features = HtoM.melt(cc.genes$g2m.genes), set.ident = TRUE)
Figure.E3.DonutChart.Data <- data.frame("Phase" = names(table(mAMs.Cyc$Phase)), "Number" = as.numeric(table(mAMs.Cyc$Phase)))
Figure.E3.DonutChart <- Figure.E3.DonutChart.Data %>% plot_ly(labels = ~Phase, values = ~Number, marker = list(colors = col2hex(c("grey30", "grey60", "grey95"))), sort = FALSE, direction = "clockwise", automargin = FALSE, textposition = "inside", textfont = list(size = 200)) # plotly default colors
Figure.E3.DonutChart <- Figure.E3.DonutChart %>% add_pie(hole = 0.5) %>% layout(showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), plot_bgcolor = "transparent")
Figure.E3.DonutChart.directory <- paste0("Figure.E3C.mAMs", ".DonutChart.pdf")
save_image(Figure.E3.DonutChart, Figure.E3.DonutChart.directory, height = 1000, width = 1000, scale = 5)

# Figure.E3D.hAMs.BarPlot
Figure.E3.BarPlot.Data <- table(hAMs.Cyc$Phase, hAMs.Cyc$orig.ident)
Figure.E3.BarPlot.Data <- Figure.E3.BarPlot.Data[, 1:4]
Figure.E3.BarPlot.Data <- as.data.frame(t(Figure.E3.BarPlot.Data) / colSums(Figure.E3.BarPlot.Data))
colnames(Figure.E3.BarPlot.Data) <- mapvalues(colnames(Figure.E3.BarPlot.Data), c("Var1", "Var2"), c("orig.ident", "Phase"))
Figure.E3.BarPlot.Data$orig.ident <- factor(Figure.E3.BarPlot.Data$orig.ident, levels = paste0("HC", 1:4))
Figure.E3.BarPlot.Data$Phase <- factor(Figure.E3.BarPlot.Data$Phase, levels = rev(c("G1", "G2M", "S")))
Figure.E3.BarPlot <- ggplot2.barplot(Figure.E3.BarPlot.Data, xName = "orig.ident", yName= "Freq", width = 0.8, xTickLabelFont = c(22, "plain", color = "grey60"), yTickLabelFont = c(22, "plain", c(rep("black", 4), rep("red", 3))), groupName = "Phase", backgroundColor = "white", removePanelGrid = TRUE, removePanelBorder = TRUE, xShowTitle = FALSE, yShowTitle = FALSE, groupColors = rev(c("grey30", "grey60", "grey95"))) + 
  scale_y_continuous(labels = scales::percent) + coord_flip(expand = FALSE) + ylab("Subject") + theme(axis.line = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(.25, "cm"), axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position = "none")
Figure.E3.BarPlot.directory <- paste0("Figure.E3D.hAMs", ".BarPlot.png")
png(Figure.E3.BarPlot.directory, width = 1200, height = 200)
print(Figure.E3.BarPlot)
dev.off()
# Figure.E3D.mAMs.BarPlot
Figure.E3.BarPlot.Data <- table(mAMs.Cyc$Phase, mAMs.Cyc$orig.ident)
Figure.E3.BarPlot.Data <- as.data.frame(t(Figure.E3.BarPlot.Data) / colSums(Figure.E3.BarPlot.Data))
colnames(Figure.E3.BarPlot.Data) <- mapvalues(colnames(Figure.E3.BarPlot.Data), c("Var1", "Var2"), c("orig.ident", "Phase"))
Figure.E3.BarPlot.Data$orig.ident <- factor(Figure.E3.BarPlot.Data$orig.ident, levels = paste0("mBAL.", 1:3))
Figure.E3.BarPlot.Data$Phase <- factor(Figure.E3.BarPlot.Data$Phase, levels = rev(c("G1", "G2M", "S")))
Figure.E3.BarPlot <- ggplot2.barplot(Figure.E3.BarPlot.Data, xName = "orig.ident", yName= "Freq", width = 0.8, xTickLabelFont = c(22, "plain", color = "grey60"), yTickLabelFont = c(22, "plain", c(rep("black", 4), rep("red", 3))), groupName = "Phase", backgroundColor = "white", removePanelGrid = TRUE, removePanelBorder = TRUE, xShowTitle = FALSE, yShowTitle = FALSE, groupColors = rev(c("grey30", "grey60", "grey95"))) + 
  scale_y_continuous(labels = scales::percent) + coord_flip(expand = FALSE) + ylab("Subject") + theme(axis.line = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(.25, "cm"), axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), legend.position = "none")
Figure.E3.BarPlot.directory <- paste0("Figure.E3D.mAMs", ".BarPlot.png")
png(Figure.E3.BarPlot.directory, width = 1200, height = 150)
print(Figure.E3.BarPlot)
dev.off()

# Figure.E3E.SankeyDiagram
Figure.2.DonutChart.Data <- round(Figure.2.DonutChart.Data/10)
mAMs.Cell.Subtypes <- c()
hAMs.Cell.Subtypes <- c()
for (j in 1:12) {
  mAMs.Cell.Subtypes <- c(mAMs.Cell.Subtypes, rep(as.vector(meta.data.list$mAMs.Cell.Subtypes)[j], colSums(Figure.2.DonutChart.Data)[j]))
  for (i in 1:12) {
    hAMs.Cell.Subtypes <- c(hAMs.Cell.Subtypes, rep(as.vector(meta.data.list$hAMs.Cell.Subtypes)[i], Figure.2.DonutChart.Data[i, j]))
  }
}
Figure.E3.SankeyDiagram.Data <- data.frame("mAMs.Cell.Subtypes" = mAMs.Cell.Subtypes, "hAMs.Cell.Subtypes" = hAMs.Cell.Subtypes)
Cell.Subypes <- c(meta.data.list$mAMs.Cell.Subtypes, meta.data.list$hAMs.Cell.Subtypes)
colors <- c(names(meta.data.list$mAMs.Cell.Subtypes), names(meta.data.list$hAMs.Cell.Subtypes))
SankeyDiagram(Figure.E3.SankeyDiagram.Data, max.categories = 100, font.size = 0, node.width = 15, node.padding = 5, colors = colors, link.color = "Source", node.position.automatic = FALSE, label.show.varname = FALSE)

# Figure.E3F.hAMs.DotPlot
hChemokine <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0008009", columns = c('SYMBOL'), keytype = "GOALL")$SYMBOL)
hAMs.CK.AE <- AverageExpression(hAMs, group.by = "AMs.Cell.Subtypes", features = hChemokine, layer = "data")
hAMs.CK <- rbind(hAMs.CK.AE$SCT, hAMs.CK.AE$RNA[setdiff(rownames(hAMs.CK.AE$RNA), rownames(hAMs.CK.AE$SCT)), ]) # add genes from RNA (absent in SCT)
hAMs.CK.Heatmap <- Heatmap(as.matrix(hAMs.CK), cluster_columns = FALSE, row_dend_reorder = TRUE, row_dend_side = "right")
CK.Genes <- rownames(hAMs.CK)[row_order(hAMs.CK.Heatmap)]
CK.RNA.Genes <- setdiff(rownames(hAMs.CK.AE$RNA), rownames(hAMs.CK.AE$SCT))
Substitute.Genes <- rownames(hAMs[["SCT"]]@data)[grep("\\.", rownames(hAMs[["SCT"]]@data))][1:length(CK.RNA.Genes)]
hAMs[["SCT"]]@data[Substitute.Genes, ] <- hAMs[["RNA"]]@data[CK.RNA.Genes, ] # replace Substitute.Genes with CK RNA genes
Figure.E3.DotPlot <- ggplot_build(DotPlot(hAMs, features = mapvalues(CK.Genes, from = CK.RNA.Genes, to = Substitute.Genes), group.by = "AMs.Cell.Subtypes", cols = c("lightgrey", "#990F0F")) + coord_flip() + theme_linedraw() + RotatedAxis() + 
                                    theme(axis.title = element_blank(), axis.text.x = element_text(size = 8, angle = 30, vjust = 1, hjust = 1), axis.text.y = element_text(size = 8), axis.ticks.x = element_line(size = 9, color = names(meta.data.list$hAMs.Cell.Subtypes)), panel.grid.major = element_line(color = "grey", size = 0.3, linetype = 1), legend.position = "none") + 
                                    scale_x_discrete(labels = CK.Genes) + scale_size_continuous(range = c(0.1, 5)) + expand_limits(x = c(0, length(CK.Genes) + 0.25)))
Figure.E3.DotPlot$data[[1]]$shape <- 22
Figure.E3.DotPlot$data[[1]]$fill <- Figure.E3.DotPlot$data[[1]]$colour
Figure.E3.DotPlot$data[[1]]$colour <- "black"
Figure.E3.DotPlot$data[[1]]$alpha <- 0.85
Figure.E3.DotPlot.directory <- paste0("Figure.E3F.hAMs.DotPlot.pdf")
pdf(Figure.E3.DotPlot.directory, width = 4, height = 5)
print(grid::grid.draw(ggplot_gtable(Figure.E3.DotPlot)))
dev.off() 
# Figure.E3G.mAMs.DotPlot
mChemokine <- unique(AnnotationDbi::select(org.Mm.eg.db, keys = "GO:0008009", columns = c('SYMBOL'), keytype = "GOALL")$SYMBOL)
mAMs.CK.AE <- AverageExpression(mAMs, group.by = "AMs.Cell.Subtypes", features = mChemokine, layer = "data")
mAMs.CK <- rbind(mAMs.CK.AE$SCT, mAMs.CK.AE$RNA[setdiff(rownames(mAMs.CK.AE$RNA), rownames(mAMs.CK.AE$SCT)), ]) # add genes from RNA (absent in SCT)
mAMs.CK.Heatmap <- Heatmap(as.matrix(mAMs.CK), cluster_columns = FALSE, row_dend_reorder = TRUE, row_dend_side = "right")
CK.Genes <- rownames(mAMs.CK)[row_order(mAMs.CK.Heatmap)]
CK.RNA.Genes <- setdiff(rownames(mAMs.CK.AE$RNA), rownames(mAMs.CK.AE$SCT))
Substitute.Genes <- rownames(mAMs[["SCT"]]@data)[grep("Rik", rownames(mAMs[["SCT"]]@data))][1:length(CK.RNA.Genes)] # not as many genes with dot, change to Rik
mAMs[["SCT"]]@data[Substitute.Genes, ] <- mAMs[["RNA"]]@data[CK.RNA.Genes, ] # replace Substitute.Genes with CK RNA genes
Figure.E3.DotPlot <- ggplot_build(DotPlot(mAMs, features = mapvalues(CK.Genes, from = CK.RNA.Genes, to = Substitute.Genes), group.by = "AMs.Cell.Subtypes", cols = c("lightgrey", "#990F0F")) + coord_flip() + theme_linedraw() + RotatedAxis() + 
                                    theme(axis.title = element_blank(), axis.text.x = element_text(size = 8, angle = 30, vjust = 1, hjust = 1), axis.text.y = element_text(size = 8), axis.ticks.x = element_line(size = 9, color = names(meta.data.list$mAMs.Cell.Subtypes)), panel.grid.major = element_line(color = "grey", size = 0.3, linetype = 1), legend.position = "none") + 
                                    scale_x_discrete(labels = CK.Genes) + scale_size_continuous(range = c(0.1, 5)) + expand_limits(x = c(0, length(CK.Genes) + 0.25)))
Figure.E3.DotPlot$data[[1]]$shape <- 22
Figure.E3.DotPlot$data[[1]]$fill <- Figure.E3.DotPlot$data[[1]]$colour
Figure.E3.DotPlot$data[[1]]$colour <- "black"
Figure.E3.DotPlot$data[[1]]$alpha <- 0.85
Figure.E3.DotPlot.directory <- paste0("Figure.E3G.mAMs.DotPlot.pdf")
pdf(Figure.E3.DotPlot.directory, width = 4, height = 4)
print(grid::grid.draw(ggplot_gtable(Figure.E3.DotPlot)))
dev.off() 
