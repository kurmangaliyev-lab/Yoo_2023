library(tidyverse)

#rm(list = ls())
#setwd("~/projects2/projLP1/transcriptome0/")

# 1. Load & Prep ----------------------------------------------------------
expr0 <- readRDS("data/expr_V1.1a.rds")

# 2. IgSF CAMs ------------------------------------------------------------
ptypes0 <- c("T4c", "T5c", "T4d", "T5d", "LLPC2", "LLPC3")
pgenes0 <- read.table("extra/IgSF.txt")[,1]

pdata0 <- expr0 %>%
  filter(subtype2 %in% ptypes0, gene %in% pgenes0) %>%
  select(subtype2, gene, time, expr) %>%
  group_by(subtype2, gene, time) %>%
  summarize(expr = mean(expr)) %>%
  ungroup() %>%
  mutate(subtype2 = factor(subtype2, ptypes0), expr2 = ifelse(expr >= 20, 20, expr))

ggplot(pdata0, aes(x = gene, y = time, fill = log1p(expr2))) +
  geom_tile(col = "grey60", size = 0.05) + coord_fixed() +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_y_discrete(labels = c("24h","", "48h", "", "72h", "","96h")) +
  facet_grid(subtype2 ~ .) +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 10), axis.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.spacing = unit(0.05, "cm"),
    strip.text = element_text(size = 12), legend.position = "none"
  )
ggsave(filename = "results/cams_1.pdf", width = 15, height = 5.5, units = "in")

# 3. Beat/Sides: T4/T5 & LLPC2/3 ------------------------------------------
lineplot0 <- function(pdata0){
  ggplot(pdata0, aes(x = time, y = expr, col = layer0)) +
    stat_summary(aes(group = subtype2), fun = "mean", geom = "line") +
    geom_point() +
    scale_color_manual(values = c("#FF7F00","#33A02C")) +
    scale_x_discrete(labels = c("24h","", "48h", "", "72h", "","96h")) +
    facet_grid(gene ~ subtype2, scales = "free_y") +
    theme_bw() +
    theme(
      axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12), legend.position = "none"
    ) 
}

##
ptypes1 <- c("T4c", "T5c", "T4d", "T5d")
pdata1 <- expr0 %>%
  filter(gene %in% c("side-II", "side-IV"), subtype2 %in% ptypes1) %>%
  mutate(layer0 = str_extract(subtype2, "c|d")) %>%
  mutate(gene = factor(gene, c("side-IV", "side-II")), subtype2 = factor(subtype2, ptypes1))
lineplot0(pdata1)
ggsave(filename = "results/cams_2a.pdf", width = 6, height = 3, units = "in")

##
ptypes2 <- c("LPi3-4", "LPC2", "LLPC2", "LLPC3")
pdata2 <- expr0 %>%
  filter(gene %in% c("beat-IIa", "beat-IIb", "beat-VI"), subtype2 %in% ptypes2) %>%
  mutate(layer0 = ifelse(subtype2 != "LLPC3", "c", "d")) %>%
  mutate(subtype2 = factor(subtype2, ptypes2))
lineplot0(pdata2)
ggsave(filename = "results/cams_2b.pdf", width = 6, height = 4, units = "in")

# 3. Beat/Sides: neurons & glia -------------------------------------------
pgenes0 <- grep("side|beat", levels(expr0$gene), value = T)
ptypes0 <- expr0 %>%
  filter(class1 %in% c("N", "G")) %>%
  distinct(class1, subtype2) %>%
  filter(!subtype2 %in% c(paste0("N", 143:195), paste0("G", 91:193))) %>%
  pull(subtype2) %>%
  as.character()
ptypes0 <- c(setdiff(ptypes0, grep("N|G", ptypes0, value = T)), grep("N|G", ptypes0, value = T))

pdata0 <- expr0 %>% 
  filter(subtype2 %in% ptypes0, time %in% c("24h", "48h", "72h", "96h"), gene %in% pgenes0) %>%
  group_by(subtype2, gene, time) %>%
  summarize(expr = mean(expr)) %>%
  ungroup() %>%
  mutate(gene = factor(gene, rev(pgenes0)), expr2 = ifelse(expr >= 20, 20, expr))

ggplot(pdata0, aes(x = subtype2, y = time, fill = log1p(expr2))) +
  geom_tile(col = "grey60", size = 0.05) + coord_fixed() +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  scale_y_discrete(labels = c("24h","","","96h"), position = "right") +
  facet_grid(gene ~ ., switch = "y") +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8), axis.title = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.spacing = unit(0.05, "cm"),
    strip.text = element_text(size = 8), legend.position = "none"
  )
ggsave(filename = "results/cams_3.pdf", width = 15.5, height = 11, units = "in")
