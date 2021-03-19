## setup
wd <- "~/Dropbox/School/UCLA MS Biology/Coursework/EEB C202/eeb_c202_final_project"
setwd(wd)

library(tidyverse)

## load in data
fish <- read.csv("datasets/chaets_all_cleaned.csv", header = T, sep = ",")
birds <- read.csv("datasets/tanager_all_cleaned.csv", header = T, sep = ",")

glimpse(fish)
glimpse(birds)

## fish
fish_human <- subset(fish, method == "human")
fish_charisma <- subset(fish, method == "charisma")

birds_human <- subset(birds, method == "human")
birds_charisma <- subset(birds, method == "charisma")

## permutation test (fish: human vs. charisma)
fish_observed_human_k <- mean(fish_human$k)
fish_observed_human_k
fish_observed_charisma_k <- mean(fish_charisma$k)
fish_observed_charisma_k
fish_observed_k <- fish_observed_human_k - fish_observed_charisma_k
fish_observed_k

iter <- 10000
perm_means <- rep(NA, iter)
for(ii in 1:iter) {
  shuffled_data <- sample(fish$method)
  perm_diff <- mean(fish$k[shuffled_data == "human"]) - mean(fish$k[shuffled_data == "charisma"])
  perm_means[ii] <- perm_diff
}
hist(perm_means, xlim = c(-1,1))
abline(v = fish_observed_k, lwd = 3, col = "blue")

#two-tailed test
p_two <- (sum(perm_means>=abs(fish_observed_k))+sum(perm_means<(-abs(fish_observed_k)))) / iter
p_two


## birds
birds_observed_human_k <- mean(birds_human$k)
birds_observed_human_k
birds_observed_charisma_k <- mean(birds_charisma$k)
birds_observed_charisma_k
birds_observed_k <- birds_observed_human_k - birds_observed_charisma_k
birds_observed_k

iter <- 10000
perm_means <- rep(NA, iter)
for(ii in 1:iter) {
  shuffled_data <- sample(birds$method)
  perm_diff <- mean(birds$k[shuffled_data == "human"]) - mean(birds$k[shuffled_data == "charisma"])
  perm_means[ii] <- perm_diff
}

hist(perm_means, xlim = c(-1,1))
abline(v = birds_observed_k, lwd = 3, col = "blue")

#two-tailed test
p_two <- (sum(perm_means>=abs(birds_observed_k))+sum(perm_means<(-abs(birds_observed_k)))) / iter
p_two

run_perm_t.test <- function(human, charisma, alldata, iter = 10000) {
  obs_mean_human_k <- mean(human$k)
  obs_mean_charisma_k <- mean(charisma$k)
  obs_mean_diff_k <- obs_mean_human_k - obs_mean_charisma_k
  
  perm_means <- rep(NA, iter)
  for(ii in 1:iter) {
    shuffled_data <- sample(alldata$method)
    perm_diff <- mean(alldata$k[shuffled_data == "human"]) - mean(alldata$k[shuffled_data == "charisma"])
    perm_means[ii] <- perm_diff
  }
  #two-tailed test
  p_two <- (sum(perm_means >= abs(obs_mean_diff_k)) + sum(perm_means < (-abs(obs_mean_diff_k)))) / iter
  return(p_two)
}

fish_pvals <- rep(NA, 100)
for(ii in 1:length(fish_pvals)) {
  print(paste0(ii, " of 100"))
  fish_pvals[ii] <- run_perm_t.test(fish_human, fish_charisma, fish)  
}
hist(fish_pvals)
mean(fish_pvals) #p <.01


birds_pvals <- rep(NA, 100)
for(ii in 1:length(birds_pvals)) {
  print(paste0(ii, " of 100"))
  birds_pvals[ii] <- run_perm_t.test(birds_human, birds_charisma, birds)  
}
hist(birds_pvals)
mean(birds_pvals) # p < .05
quantile(birds_pvals, .025)
quantile(birds_pvals, .975)



### distribution of colors
hist(fish_human$black)
hist(fish_charisma$black)
t.test(fish_human$black, fish_charisma$black, paired = T)

hist(fish_human$blue)
hist(fish_charisma$blue)
t.test(fish_human$blue, fish_charisma$blue, paired = T)

hist(fish_human$brown)
hist(fish_charisma$brown)

hist(fish_human$green)
hist(fish_charisma$green)

hist(fish_human$grey)
hist(fish_charisma$grey)

hist(fish_human$orange)
hist(fish_charisma$orange)

hist(fish_human$purple)
hist(fish_charisma$purple)

hist(fish_human$red)
hist(fish_charisma$red)

hist(fish_human$white)
hist(fish_charisma$white)

hist(fish_human$yellow)
hist(fish_charisma$yellow)


hist(birds_human$black)
hist(birds_charisma$black)

hist(birds_human$blue)
hist(birds_charisma$blue)

hist(birds_human$brown)
hist(birds_charisma$brown)

hist(birds_human$green)
hist(birds_charisma$green)

hist(birds_human$grey)
hist(birds_charisma$grey)

hist(birds_human$orange)
hist(birds_charisma$orange)

hist(birds_human$purple)
hist(birds_charisma$purple)

hist(birds_human$red)
hist(birds_charisma$red)

hist(birds_human$white)
hist(birds_charisma$white)

hist(birds_human$yellow)
hist(birds_charisma$yellow)











