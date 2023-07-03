library(tidyverse)
library(neuprintr)

#rm(list = ls())
#setwd("~/projects2/projLP1/connectome0")

# 1. Connect & Prep
npr0 <- neuprint_login(server = "https://neuprint.janelia.org", dataset = "fib19:v1.0", Cache = FALSE)
seed0 <- read.delim("data/seed0.txt")

# 2. Get Synapses
# 2.1 Inputs
inputs1 <- neuprint_connection_table(conn = npr0, bodyids = seed0$bodyID, partners = "inputs", by.roi = TRUE, details = TRUE)
inputs2 <- neuprint_connection_table(conn = npr0, bodyids = seed0$bodyID, partners = "inputs", by.roi = FALSE, details = TRUE)
inputs2 <- right_join(inputs1, inputs2) %>%
  filter(is.na(ROIweight)) %>%
  mutate(roi = "NotPrimary", ROIweight = weight)

np0 <- inputs1 %>%
  group_by(bodyid, partner, prepost, weight, name, type) %>%
  summarize(ROItotal = sum(ROIweight)) %>%
  mutate(roi = "NotPrimary", ROIweight = weight - ROItotal) %>%
  filter(ROIweight > 0) %>%
  select(-ROItotal)
inputs0 <- bind_rows(inputs1, inputs2, np0) %>% arrange(desc(weight))
inputs0 <- inner_join(seed0, inputs0, by = c("bodyID" = "bodyid"))

# 2.2 Outputs
outputs1 <- neuprint_connection_table(conn = npr0, bodyids = seed0$bodyID, partners = "outputs", by.roi = TRUE, details = TRUE)
outputs2 <- neuprint_connection_table(conn = npr0, bodyids = seed0$bodyID, partners = "outputs", by.roi = FALSE, details = TRUE)
outputs2 <- right_join(outputs1, outputs2) %>%
  filter(is.na(ROIweight)) %>%
  mutate(roi = "NotPrimary", ROIweight = weight)

np0 <- outputs1 %>%
  group_by(bodyid, partner, prepost, weight, name, type) %>%
  summarize(ROItotal = sum(ROIweight)) %>%
  mutate(roi = "NotPrimary", ROIweight = weight - ROItotal) %>%
  filter(ROIweight > 0) %>%
  select(-ROItotal)

outputs0 <- bind_rows(outputs1, outputs2, np0) %>% arrange(desc(weight))
outputs0 <- inner_join(seed0, outputs0, by = c("bodyID" = "bodyid"))

# 3. Combine Inputs/Outputs
con1 <- bind_rows(
  inputs0 %>%
    select(
      pre_bodyID = partner, pre_type = type.y, pre_name = name,
      post_bodyID = bodyID, post_type = type.x, post_name = instance,
      roi, weight = ROIweight
    ),
  outputs0 %>%
    select(
      pre_bodyID = bodyID, pre_type = type.x, pre_name = instance,
      post_bodyID = partner, post_type = type.y, post_name = name,
      roi, weight = ROIweight
    )
) %>% distinct() %>% drop_na(weight)

# 4. Neurons
neu0 <- neuprint_get_meta(conn = npr0, bodyids = unique(c(con1$pre_bodyID, con1$post_bodyID)))
neu0 <- neu0 %>% select(bodyID = bodyid, status, statusLabel, name, type)

# 5. Rename
rename0 <- c("CT1_part" = "CT1", "LPT-LPi2a" = "LPT/LPi2a", "HSN" = "HS", "HSS" = "HS", "HSE" = "HS")
neu0 <- neu0 %>% 
  mutate(
    bodyID = as.character(bodyID),
    type = recode(type, !!!rename0)
  )
con1 <- con1 %>%
  mutate(
    pre_bodyID = as.character(pre_bodyID),
    post_bodyID = as.character(post_bodyID),
    roi = factor(roi, c("ME","LO","LOP","NotPrimary")),
    pre_type = recode(pre_type, !!!rename0),
    post_type = recode(post_type, !!!rename0)
  )

# 6. Summarize Connectomes: T4/T5-centric
inputs0 <- con1 %>%
  filter(post_bodyID %in% seed0$bodyID, pre_bodyID != post_bodyID, !is.na(pre_type)) %>%
  group_by(pre_type, post_bodyID, post_type, roi) %>%
  summarize(weight = sum(weight)) %>%
  ungroup() %>%
  group_by(pre_type, post_type, roi) %>%
  arrange(desc(weight)) %>%
  summarize(weight = round(mean(weight))) %>%
  ungroup() %>%
  arrange(post_type, roi, desc(weight)) %>%
  mutate(is_post = 1) %>%
  as.data.frame()

outputs0 <- con1 %>%
  filter(pre_bodyID %in% seed0$bodyID, pre_bodyID != post_bodyID, !is.na(post_type)) %>%
  group_by(pre_bodyID, pre_type, post_type, roi) %>%
  summarize(weight = sum(weight)) %>%
  ungroup() %>%
  group_by(pre_type, post_type, roi) %>%
  arrange(desc(weight)) %>%
  summarize(weight = round(mean(weight))) %>%
  ungroup() %>%
  arrange(pre_type, roi, desc(weight)) %>%
  mutate(is_post = 0) %>%
  as.data.frame()

con2 <- bind_rows(inputs0, outputs0)

# 8. Save
saveRDS(neu0, file = "data/neu0.rds")
saveRDS(con1, file = "data/con1.rds")
saveRDS(con2, file = "data/con2.rds")

