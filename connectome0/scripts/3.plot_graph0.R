library(tidyverse)
library(igraph)

#rm(list = ls())
#setwd("~/projects2/projLP1/connectome0")

# 1. Load & Prep
seed0 <- read.delim("data/seed0.txt")
con1 <- readRDS("data/con1.rds")
con2 <- readRDS("data/con2.rds")

t45 <- c("T4a","T5a", "T4b","T5b", "T4c","T5c", "T4d","T5d")

# 2. Input/Outputs: partners from 2.plotHeatmap0.R
inputs2 <- c("Mi1", "Mi4", "Mi9", "Tm3", "Tm1", "Tm2", "Tm4", "Tm9", "TmY15", "CT1")
outputs2 <- c(
  "LPC1", "LPC2", "LLPC1", "LLPC2", "LLPC3", "LPLC1", "LPLC2",
  "LPi1-2", "LPi2-1", "LPi3-4", "LPi4-3","LPi34-12",
  "Y3", "Y11", "Y12", "Tlp11", "Tlp12", "Tlp13", "Tlp14", "TmY4", "TmY5a", "TmY14", "TmY16", "TmY20",
  "HS", "DCH", "LPT/LPi2a", "H2","Am1", "LPi2b", "LPi3a", "LPT3b", "VS"
)
mapped0 <- c(t45, "Mi1","Tm3","Mi4","Mi9","Tm1","Tm2","Tm4","Tm9","LLPC1","LPC1","LPC2","LLPC2","LLPC3","LPLC1","LPLC2","LPi3-4","TmY5a")

# 3. Cytoscape: graph layout
# 3.1 Export
pdata0 <- con2 %>%
  filter((is_post == 1 & pre_type %in% inputs2 & roi %in% c("ME", "LO")) | (is_post == 0 & post_type %in% outputs2 & roi == "LOP")) %>%
  filter(weight >= 7)
write.table(pdata0, file = "extra/cytoscape/input0.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 3.2 Modify groph in Cytoscape: LP1.cys

# 3.3 Import
cyto0 <- jsonlite::fromJSON(txt = "extra/cytoscape/output0.cyjs")
cyto0 <- bind_cols(cyto0[["elements"]][["nodes"]][["data"]], cyto0[["elements"]][["nodes"]][["position"]])

layout0 <- data.frame(type = cyto0$name, x = cyto0$x, y = -cyto0$y)
layout0 <- layout0 %>% column_to_rownames("type") %>% select(x,y) %>% as.matrix()

# 4. Graph
# 4.1 Prep data
pdata0 <- pdata0 %>%
  mutate(is_mapped = ifelse(pre_type %in% mapped0 & post_type %in% mapped0, 1, 0)) %>%
  arrange(is_mapped)
net0 <- graph_from_data_frame(pdata0)
net0 <- set_graph_attr(net0, name = "layout", value = layout0[V(net0)$name,])
write.table(pdata0, file = "results/graph_data.txt", col.names = T, row.names = F, sep = "\t", quote = F)

# 4.2 Nodes/Edges
pal3 <- setNames(rep("darkgrey", length(V(net0))), V(net0)$name)
pal3[mapped0] <- "lightslateblue"
pal3[t45] <- c("#1F78B4","#A6CEE3", "#E31A1C","#FB9A99", "#FF7F00","#FDBF6F", "#33A02C","#B2DF8A")

V(net0)$color <- pal3
V(net0)$label.cex <- ifelse(V(net0)$name %in% t45 , 1.35, 1)
V(net0)$label.cex[V(net0)$name %in% c("LPT/LPi2a", "LPi34-12")] <- 0.8

E(net0)$color <- ifelse(pdata0$pre_type %in% mapped0 & pdata0$post_type %in% mapped0, "black", "grey")
E(net0)$width <- pdata0$weight/10

# 4.3 Plot
pdf(file = "results/con_graph_1.pdf", width = 12.5, height = 7.5)
plot(
  net0, edge.arrow.size = 0.3, asp = 0.6,
  vertex.size = 12, vertex.label.font = 2, vertex.label.color = "black", vertex.frame.color = NA
)
dev.off()
