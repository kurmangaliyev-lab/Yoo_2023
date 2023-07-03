library(tidyverse)
library(Seurat)

#rm(list = ls())
#setwd("~/projects2/projLP1/transcriptome0/")

# 1. Load & Prep ----------------------------------------------------------
data0 <- readRDS("data/data_V1.1a.rds")
expr0 <- readRDS("data/expr_V1.1a.rds")

# 2. LPC/LLPC markers ---------------------------------------------------
types1 <- c("LPC1", "LPC2", "LLPC1", "LLPC2", "LLPC3")

# 2.1 Common & cell-type-specific markers
markers1 <- FindMarkers(data0, types1, max.cells.per.ident = 5000, logfc.threshold = 1, min.pct = 0.35, pseudocount.use = 0.01, only.pos = T)

data1 <- subset(data0, type2 %in% types1)
markers2 <- lapply(types1, function(type0){
  FindMarkers(data1, type0, max.cells.per.ident = 1000, logfc.threshold = 1, min.pct = 0.35, pseudocount.use = 0.01, only.pos = T) %>% 
    rownames_to_column("gene") %>%
    mutate(type2 = type0) %>%
    relocate(type2)
}) %>% bind_rows()

# 2.2 Plot markers: LPC/LLPC
head(markers1, n = 20)
markers2 %>%
  filter(p_val_adj < 0.01, avg_log2FC > 2) %>%
  group_by(type2) %>%
  slice(1:10) %>%
  as.data.frame()

pgenes0 <- c(
  "DIP-theta", "CG14459", "CG43689", "disco", "disco-r", "dmrt99B", "ham", "fru",
  "Sox21b", "lov", "tey", "vn", "pdm3", "bab1", "bab2", "br", "zfh2", "Dop1R2", "dpr5", "sdk", "DIP-kappa", "DIP-delta", "DIP-beta"
)

expr0 %>%
  filter(type2 %in% types1, gene %in% pgenes0) %>%
  mutate(type2 = factor(type2, types1), gene = factor(gene, rev(pgenes0)), expr = ifelse(expr > 20, 20, expr)) %>%
  ggplot(aes(x = time, y = gene, fill = log1p(expr))) +
  geom_tile(col = "grey20") + 
  scale_fill_distiller(palette = "RdBu") +
  scale_x_discrete(labels = c("24h","","48h","","72h","","96h")) +
  geom_hline(yintercept = c(15.5), linewidth = 1) +
  facet_grid(. ~ type2) +
  NoLegend() +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.spacing = unit(0.1, "lines"))

ggsave(filename = "results/markers_1.pdf", width = 4.5, height = 3.5, units = "in")

# 2.3 Plot markers: other neuron types
expr0 %>%
  filter(class1 %in% "N", !type2 %in% types1, gene %in% pgenes0) %>%
  droplevels() %>%
  group_by(type2, gene) %>%
  mutate(gene = factor(gene, rev(pgenes0)), expr = ifelse(expr > 20, 20, expr)) %>%
  summarize(expr = mean(expr)) %>%
  ggplot(aes(x = type2, y = gene, fill = log1p(expr))) +
  geom_tile(col = "grey20") +
  scale_fill_distiller(palette = "RdBu") +
  geom_hline(yintercept = c(15.5), size = 1) +
  NoLegend() +
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 6, face = "bold", angle = 90, hjust = 1, vjust = 0.5), axis.ticks.x = element_blank())

ggsave(filename = "results/markers_2.pdf", width = 12, height = 3.5, units = "in")

# 2.3: Validation: line plots
expr0 %>%
  filter(gene %in% c("br", "pdm3", "DIP-beta"), type2 %in% c("LLPC1", "LLPC2", "LLPC3", "LPC2")) %>%
  mutate(gene = factor(gene, c("br", "pdm3", "DIP-beta"))) %>%
  ggplot(aes(x = time, y = expr)) +
  stat_summary(aes(group = type2), fun = "mean", geom = "line") +
  geom_point(col = "grey40") +
  scale_x_discrete(labels = c("24h","", "48h", "", "72h", "","96h")) +
  facet_grid(gene ~ type2, scales = "free_y") +
  theme_bw() + NoLegend() +
  theme(
    axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12)
  )
ggsave(filename = "results/markers_3.pdf", width = 5, height = 3.5, units = "in")

# 3. Pairwise DEGs at 48h APF -----------------------------------------------------------
# 3.1 Subset & Prep
data2 <- subset(data0, time == "48h")
expr2 <- expr0 %>%
  filter(time == "48h") %>%
  group_by(subtype2, gene) %>%
  summarize(expr = mean(expr)) %>%
  ungroup()

# 3.2 (T4c vs T4d) AND (T5c vs T5d)
degT <- inner_join(
  FindMarkers(data2, "T4c", "T4d", pseudocount.use = 0.01) %>% rownames_to_column("gene"),
  FindMarkers(data2, "T5c", "T5d", pseudocount.use = 0.01) %>% rownames_to_column("gene"),
  by = "gene", suffix = c("_T4", "_T5")
)

exprT <- expr2 %>%
  filter(subtype2 %in% c("T4c", "T5c", "T4d", "T5d"), gene %in% degT$gene) %>%
  select(subtype2, gene, expr) %>%
  pivot_wider(names_from = subtype2, values_from = expr)

sigT <- inner_join(exprT, degT, by = "gene") %>% 
  select(-starts_with("pct")) %>%
  filter(T4c > 1 | T4d > 1 | T5c > 1 | T5d > 1) %>%
  filter((abs(avg_log2FC_T4) > log2(3) & p_val_adj_T4 < 0.01) | (abs(avg_log2FC_T5) > log2(3) & p_val_adj_T5 < 0.01)) %>%
  arrange((avg_log2FC_T4 + avg_log2FC_T5)/2) %>%
  mutate(sig_both = ifelse(p_val_adj_T4 < 0.01 & abs(avg_log2FC_T4) > log2(3) & p_val_adj_T5 < 0.01 & abs(avg_log2FC_T5) > log2(3), "yes", "no"))

## scatter
ggplot(sigT, aes(x = avg_log2FC_T4, y = avg_log2FC_T5, col = sig_both)) +
  geom_point(size = 1.5) + scale_color_manual(values = c("black", "red")) +
  lemon::scale_x_symmetric() + lemon::scale_y_symmetric() +
  geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  xlab("log2FC(T4c/T4d)") + ylab("log2FC(T5c/T5d)") +
  theme_classic() + NoLegend()

ggsave(filename = "results/deg_1.pdf", width = 3.5, height = 3.5, units = "in")

## heatmap
ptypes0 <- c("T4c", "T5c", "T4d", "T5d")
pgenes0 <- sigT$gene[sigT$sig_both == "yes"]
pdata0 <- expr0 %>%
  filter(subtype2 %in% ptypes0, gene %in% pgenes0) %>%
  select(subtype2, gene, time, expr) %>%
  group_by(subtype2, gene, time) %>%
  summarize(expr = mean(expr)) %>%
  ungroup() %>%
  mutate(subtype2 = factor(subtype2, ptypes0), gene = factor(gene, pgenes0), expr2 = ifelse(expr >= 20, 20, expr))

ggplot(pdata0, aes(x = gene, y = time, fill = log1p(expr2))) +
  geom_tile(col = "grey60", size = 0.05) + coord_fixed() +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_y_discrete(labels = c("24h","", "48h", "", "72h", "","96h")) +
  facet_grid(subtype2 ~ .) +
  NoLegend() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12), axis.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1), panel.spacing = unit(0.05, "cm")
  )
ggsave(filename = "results/deg_2.pdf", width = 2, height = 5, units = "in")

# 3. LLPC2 vs LLPC3
degL <- FindMarkers(data2, "LLPC2", "LLPC3", pseudocount.use = 0.01) %>% rownames_to_column("gene")

exprL <- expr2 %>%
  filter(subtype2 %in% c("LLPC2", "LLPC3"), gene %in% degL$gene) %>%
  select(subtype2, gene, expr) %>%
  pivot_wider(names_from = subtype2, values_from = expr)

sigL <- inner_join(exprL, degL, by = "gene") %>% 
  select(-starts_with("pct")) %>% 
  filter(abs(avg_log2FC) > log2(3) &  p_val_adj < 0.01) %>%
  filter(LLPC2 > 1 | LLPC3 > 1) %>%
  arrange(avg_log2FC)

## heatmap
pdata0 <- expr0 %>%
  filter(subtype2 %in% c("LLPC2", "LLPC3"), gene %in% sigL$gene) %>%
  select(subtype2, gene, time, expr) %>%
  group_by(subtype2, gene, time) %>%
  summarize(expr = mean(expr)) %>%
  ungroup() %>%
  mutate(gene = factor(gene, sigL$gene), expr2 = ifelse(expr >= 20, 20, expr))

ggplot(pdata0, aes(x = gene, y = time, fill = log1p(expr2))) +
  geom_tile(col = "grey60", size = 0.05) + coord_fixed() +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_y_discrete(labels = c("24h","", "48h", "", "72h", "","96h")) +
  facet_grid(subtype2 ~ .) +
  NoLegend() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12), axis.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1), panel.spacing = unit(0.05, "cm")
  )
ggsave(filename = "results/deg_3.pdf", width = 4, height = 3, units = "in")