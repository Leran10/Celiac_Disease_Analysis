library(tidyverse)
library(dplyr)

metadata2 <- read.csv("~/Desktop/metadata_final_update_11262024.csv")

metadata1 <- read.csv("~/Desktop/celiacs_metadata_final_110624_leran.csv") %>%
             mutate(sampleName = paste0(celiacs_group,"_",PIDu)) %>%
             mutate(sampleName = gsub("6Y","72M",sampleName)) %>%
             mutate(sampleName = gsub("7Y","84M",sampleName)) %>%
             select(sampleName,Name)


final <- left_join(metadata1,metadata2,by = "sampleName")

write.csv(final,"~/Handley Lab Dropbox/16S/Celiac/metadata/metadata_final_update_forDhoha_03132025.csv")

View(final)
View(metadata2)
  