library(tidyverse)

#rm(list = ls())
#setwd("~/projects2/projLP1/connectome0")

# 1. Load & Prep
con1 <- readRDS("data/con1.rds")
con2 <- readRDS("data/con2.rds")

seed0 <- read.delim("data/seed0.txt")
t45 <- c("T4a","T5a", "T4b","T5b", "T4c","T5c", "T4d","T5d")
seed0 <- seed0 %>%
  mutate(type = factor(type, t45), bodyID = as.character(bodyID)) %>%
  select(bodyID, type) %>%
  arrange(type)

# 2. Select Partners: mean weight>=10 to any T4/T5 subtype
# 2.1 Inputs (manual order)
con2 %>%
  filter(post_type %in% t45, !pre_type %in% t45, roi %in% c("ME", "LO"), weight >= 10) %>%
  pull(pre_type) %>% unique() %>% sort()

inputs1 <- c("Mi1", "Mi4", "Mi9", "Tm3", "Tm1", "Tm2", "Tm4", "Tm9", "TmY15", "CT1")

# 2.2 Outputs (manual order)
con2 %>%
  filter(pre_type %in% t45, !post_type %in% t45, roi == "LOP", weight >= 10) %>%
  pull(post_type) %>% unique() %>% sort()

outputs1 <- c(
  "LPC1", "LPC2", "LLPC1", "LLPC2", "LLPC3", "LPLC1", "LPLC2",
  "LPi1-2", "LPi2-1", "LPi3-4", "LPi4-3","LPi34-12",
  "Y3", "Y11", "Y12", "Tlp11", "Tlp12", "Tlp13", "Tlp14", "TmY4", "TmY5a", "TmY14", "TmY16", "TmY20",
  "HS", "DCH", "LPT/LPi2a", "H2","Am1", "LPi2b", "LPi3a", "LPT3b", "VS"
)

# 3. Connectomes: individual T4/T5 ~ average partner
# 3.1 Inputs
pdata1 <- con1 %>%
  filter(post_bodyID %in% seed0$bodyID, pre_type %in% inputs1, roi %in% c("ME","LO")) %>%
  group_by(pre_type, post_bodyID) %>%
  summarize(weight = sum(weight)) %>%
  select(neuron_x = post_bodyID, neuron_y = pre_type, weight) %>%
  ungroup()
  
# 3.2 Outputs
pdata2 <- con1 %>%
  filter(pre_bodyID %in% seed0$bodyID, post_type %in% outputs1, roi %in% c("LOP")) %>%
  group_by(pre_bodyID, post_type) %>%
  summarize(weight = sum(weight)) %>%
  select(neuron_x = pre_bodyID, neuron_y = post_type, weight) %>%
  ungroup()

# 3.3 Plot
pdata0 <- bind_rows(pdata1, pdata2) %>%
  complete(neuron_x, neuron_y, fill = list(weight = 0)) %>%
  mutate(neuron_y = factor(neuron_y, c(outputs1, inputs1))) %>%
  mutate(dir0 = ifelse(neuron_y %in% inputs1, "input", "output")) %>%
  as.data.frame()

ggplot(pdata0, aes(x = neuron_x, y = neuron_y, fill = weight)) +
  geom_tile(col = "grey60", size = 0.05) + coord_fixed() +
  scale_x_discrete(limits = seed0$bodyID, labels = seed0$type) +
  scale_fill_gradient2(low = "white", mid = "deepskyblue3", high = "darkblue", midpoint = 50, limits = c(0,100)) +
  geom_hline(yintercept = 33.5, col = "black", size = 0.5) +
  geom_vline(xintercept = c(seq(5.5, 35.5, 5)), col = "grey20", size = 0.25) +
  theme(
    axis.text.x = element_text(size = 7, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 7, face = "bold"),
    axis.title = element_blank(), axis.ticks = element_blank(), 
    legend.title = element_blank(), legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )
ggsave(filename = "results/con_heatmap_1.pdf", width = 7, height = 8, units = "in")

