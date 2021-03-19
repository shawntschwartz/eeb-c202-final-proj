## setup
wd <- "~/Dropbox/School/UCLA MS Biology/Coursework/EEB C202/eeb_c202_final_project"
setwd(wd)

library(tidyverse)

## load in data
fish_human <- read.csv("datasets/chaets-human-calls.csv", header = T, sep = ",")
fish_charisma <- read.csv("datasets/chaets-charisma-calls.csv", header = T, sep = ",")

birds_human <- read.csv("datasets/tanager-human-calls.csv", header = T, sep = ",")
birds_charisma <- read.csv("datasets/tanager-charisma-calls.csv", header = T, sep = ",")

## preview data
glimpse(fish_human)
glimpse(fish_charisma)
glimpse(birds_human)
glimpse(birds_charisma)

## combine datasets
fish_human <- fish_human %>%
  mutate(method = "human")

fish_charisma <- fish_charisma %>%
  mutate(method = "charisma")

birds_human <- birds_human %>%
  mutate(method = "human")

birds_charisma <- birds_charisma %>%
  mutate(method = "charisma")

fish <- rbind(fish_human, fish_charisma)
birds <- rbind(birds_human,  birds_charisma)

glimpse(fish)
glimpse(birds)

write.csv(fish, "datasets/chaets_all_cleaned.csv", row.names = F)
write.csv(birds, "datasets/tanager_all_cleaned.csv", row.names = F)