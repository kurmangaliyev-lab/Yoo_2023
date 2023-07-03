library(tidyverse)
library(Seurat)

#rm(list = ls())
#setwd("~/projects2/projLP1/transcriptome0/")

data0 <- readRDS("data/data_V1.0a.rds")

# 1. Load & Prep ----------------------------------------------------------
pal_t45 <- c(T4a="#1F78B4", T5a="#A6CEE3", T4b="#E31A1C",T5b="#FB9A99", T4c="#FF7F00",T5c="#FDBF6F", T4d="#33A02C", T5d="#B2DF8A")
pal_llpc <- c(LPC1="#E31A1C", LPC2="#E6B648", LLPC1="#1F78B4", LLPC2="#FF7F00", LLPC3="#33A02C")

# 2. Subset & Re-cluster LPC/LLPC -----------------------------------------
data1 <- subset(data0, idents = c("LLPC1", "LPC1"))
data1@active.assay <- "integrated"

# 2.1 Re-cluster
data1 <- data1 %>%
  FindVariableFeatures(nfeatures = 1000) %>%
  ScaleData() %>%
  RunPCA()
ElbowPlot(data1)

dims1 <- 1:9
data1 <- data1 %>%
  RunTSNE(dims = dims1) %>%
  FindNeighbors(dims = dims1) %>%
  FindClusters(resolution = 0.1)

# 2.2 Rename
VlnPlot(data1, features = c("dmrt99B", "br", "pdm3", "beat-IIa", "beat-IIb","beat-VI", "DIP-beta", "DIP-delta", "tey"))
data1 <- RenameIdents(data1, `0` = "LLPC1", `1` = "LLPC2", `2` = "LLPC3", `3` = "LPC1", `4` = "LPC2")
data1$type2 <- Idents(data1)

# 2.3 Plot tSNE
th0 <- theme(axis.title = element_blank(), axis.text = element_blank(), plot.title = element_blank()) + NoLegend()
cowplot::plot_grid(
  DimPlot(data1, group.by = "type1", pt.size = 2, label = T, raster = T) + th0,
  DimPlot(data1, group.by = "type2", pt.size = 2, cols = pal_llpc, label = T, raster = T) + th0
)
ggsave(filename = "results/tsne_1.pdf", width = 6, height = 3, units = "in")

# 3. Atlas V1.1 -----------------------------------------------------------
# 3.1 Rename & reorder
data0$type2 <- as.character(data0$type1)
data0$subtype2 <- as.character(data0$subtype1)

data0@meta.data[colnames(data1), "type2"] <- as.character(data1@meta.data[,"type2"])
data0@meta.data[colnames(data1), "subtype2"] <- as.character(data1@meta.data[,"type2"])

data0@meta.data$type2 = recode(data0@meta.data$type2, "Lpi3.4"="LPi3-4")
data0@meta.data$subtype2 = recode(data0@meta.data$subtype2, "Lpi3.4"="LPi3-4", "Tm9a"="Tm9", "Tm9b"="Tm9")

data0@meta.data %>%
  distinct(class1, type2, subtype2)

order0 <- data0@meta.data %>%
  distinct(class1, type2, subtype2) %>%
  mutate(type2 = factor(type2, gtools::mixedsort(unique(type2))), subtype2 = factor(subtype2, gtools::mixedsort(unique(subtype2)))) %>%
  arrange(class1, type2, subtype2)

data0$type2 <- factor(data0$type2, unique(as.character(order0$type2)))
data0$subtype2 <- factor(data0$subtype2, unique(as.character(order0$subtype2)))

# 3.2 Save Seurat-object
Idents(data0) <- data0$subtype2
saveRDS(data0, file = "data/data_V1.1a.rds", compress = F)

# 3.3 Main tSNE
mapped0 <- c(
  "Mi1", "Mi4", "Mi9", "Tm3", "Tm1", "Tm2", "Tm4", "Tm9", 
  "LPC1", "LPC2", "LLPC1", "LLPC2", "LLPC3", "LPLC1", "LPLC2", "LPi3-4", "TmY5a"
)
data0$tmp_label <- as.character(Idents(data0))
data0$tmp_label[!(data0$tmp_label %in% c(names(pal_t45), mapped0))] <- "."

DimPlot(data0, group.by = "tmp_label", cols = c(pal_t45, setNames(rep("lightslateblue", 17), mapped0)), label = T, repel = T, label.size = 3) +
  theme(axis.title = element_blank(), axis.text = element_blank(), plot.title = element_blank()) + 
  NoLegend()
ggsave(filename = "results/tsne_2.pdf", width = 6, height = 5.5, units = "in")

# 4. Hierarchical clustering at 48h APF -----------------------------------
# 4.1 Subset & var.genes
data2 <- subset(data0, sample == "W1118.B_48h" & class1 == "N")
data2 <- DietSeurat(data2, counts = F, data = T, scale.data = F, assays = "RNA", dimreducs = "tsne")

features0 <- FindVariableFeatures(data2, nfeatures = 1000) %>% VariableFeatures()
mat0 <- AverageExpression(data2, assays = "RNA", features = features0)[["RNA"]]

# 4.2 Cluster & plot
dist0 <- as.dist(1 - cor(log1p(mat0), method = "pearson"))
hclust0 <- hclust(dist0, method = "ward.D2")
tree0 <- ape::as.phylo(hclust0)

pdf(file = "results/dendro_1.pdf", width = 7, height = 7)
ape::plot.phylo(tree0, type = "fan",cex = 0.7, font = 1, label.offset = 0.005, rotate.tree = -60, no.margin = T)
dev.off()

# 5. Average Expression ---------------------------------------------------
# 5.1 Calculate
data0@meta.data <- data0@meta.data %>% unite(col = "group0",set, rep, time, sample, class1, type2, subtype2, remove = F, sep = "-")
mat0 <- AverageExpression(data0, assays = "RNA", group.by = "group0")[["RNA"]]

# 5.2 Filter lowly expressed genes
features0 <- rownames(data0)[which(Matrix::rowSums(GetAssayData(data0, assay = "RNA", slot = "data") > 0) >= 50)]
mat0 <- mat0[features0,]

# 5.3 Transform to long format
expr0 <- mat0 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "group0", values_to = "expr") %>%
  mutate(gene = factor(gene, levels = gtools::mixedsort(rownames(mat0))))

groups0 <- data0@meta.data %>%
  group_by(group0, set, rep, time, sample, class1, type2, subtype2) %>%
  count() %>%
  ungroup() %>%
  mutate(
    set = factor(set, levels = gtools::mixedsort(set) %>% unique),
    rep = factor(rep, levels = gtools::mixedsort(rep) %>% unique),
    time = factor(time, levels = gtools::mixedsort(time) %>% unique),
    sample = factor(sample, levels = c(
                      "DGRP.A_24h", "DGRP.B_24h", "W1118.A_24h", "DGRP.A_36h", "DGRP.B_36h", "DGRP.A_48h", "DGRP.B_48h", "W1118.A_48h", "W1118.B_48h", 
                      "DGRP.A_60h", "DGRP.B_60h","DGRP.A_72h", "DGRP.B_72h", "W1118.A_72h","DGRP.A_84h", "DGRP.B_84h", "DGRP.A_96h", "DGRP.B_96h", "W1118.A_96h"
                    ))
  )

expr0 <- inner_join(groups0, expr0, by = "group0") %>%
  select(-group0) %>%
  arrange(set, rep, time, sample, class1, type2, subtype2, gene) %>%
  as.data.frame()

# 5.4 Save
saveRDS(expr0, file = "data/expr_V1.1a.rds", compress = F)

