#################################################################
##           Color Classification Pipeline Functions           ##
##    Adapted from Alfaro, Karan, Schwartz, & Shultz (2019)    ##
#################################################################
wd <- "~/Dropbox/School/UCLA MS Biology/Coursework/EEB C202/eeb_c202_final_project"
setwd(wd)

classifyByUniqueK <- function(path, kdf)
{
  image_paths <- kdf$img
  images_list <- rep(NA, length(image_paths))
  for(image in 1:length(images_list))
  {
    if(Sys.info()['sysname'] != "Windows")
    {
      images_list[image] <- tail(strsplit(image_paths[image], "/")[[1]], 1)
    }
    else
    {
      images_list[image] <- tail(strsplit(as.character(image_paths[image]), "\\\\")[[1]], 1)
    }
  }
  
  classifications <- list()
  
  for(image in 1:length(images_list))
  {
    pic <- pavo::getimg(file.path(path, images_list[image]), max.size = 3)
    cat(paste0("Image (", image, "/", length(images_list),"): ", images_list[image], 
               " with (k = ", kdf$k[image], ")\n"))
    classifications[[image]] <- pavo::classify(pic, kcols = kdf$k[image])
    cat("\n")
  }
  
  names(classifications) <- images_list
  return(classifications)
}

rgbEucDist <- function(rgb_table_altered, c1, c2) 
{
  euc_dist <- sqrt((rgb_table_altered[c1,"col1"]-rgb_table_altered[c2,"col1"])^2+(rgb_table_altered[c1,"col2"]-rgb_table_altered[c2,"col2"])^2) %>%
    .[1,1]
  return(euc_dist)
}

rgbLumDist <- function(rgb_table_altered, c1, c2)
{
  lum_dist <- sqrt((rgb_table_altered[c1,"lum"]-rgb_table_altered[c2,"lum"])^2) %>%
    .[1,1]
  return(lum_dist)
}

#input is a single classified image
calcEucLumDists <- function(classified_image)
{
  #extract RGB values for n colors
  class_rgb <- attr(classified_image, 'classRGB')
  class_rgb_altered <- class_rgb %>%
    rownames_to_column(var = "col_num") %>%
    as_tibble %>%
    mutate(col1 = (R-G)/(R+G), col2 = (G-B)/(G+B), lum = R+G+B) %>%
    select(col1, col2, lum)
  
  #create a matrix to hold colors based on the number of possible color comparisons
  euc_dists <- matrix(nrow=choose(nrow(class_rgb),2),ncol=4)
  
  combos_simple <- t(combn(rownames(class_rgb),2)) %>%
    as_tibble %>%
    transmute(c1 = as.numeric(V1), c2 = as.numeric(V2)) %>%
    as.data.frame()
  
  combos <- matrix(nrow=nrow(combos_simple),ncol=4)
  for(i in 1:nrow(combos_simple))
  {
    combos[i,1] <- combos_simple[i,1]
    combos[i,2] <- combos_simple[i,2]
    combos[i,3] <- rgbEucDist(class_rgb_altered,combos_simple[i,1],combos_simple[i,2])
    combos[i,4] <- rgbLumDist(class_rgb_altered,combos_simple[i,1],combos_simple[i,2])
  }
  
  combos <- combos %>%
    as.data.frame %>%
    as_tibble %>%
    dplyr::rename(c1 = V1,
                  c2 = V2,
                  dS = V3,
                  dL = V4) %>%
    as.data.frame
  
  return(combos)
}

#get distance data frame for each picture
getImgClassKDists <- function(classifications, euclidean_lum_dists) 
{
  return(map(.x=classifications,.f=euclidean_lum_dists))
}

#calculate the adjacency stats for each image, using the calculated distances as proxies for dS and dL
getAdjStats <- function(classifications, img_class_k_dists, xpts=100, xscale=100) 
{
  adj_k_dists_list <- list()
  
  for(i in 1:length(classifications)) 
  {
    adj_k_dists_list[[i]] <- pavo::adjacent(classimg = classifications[[i]],coldists=img_class_k_dists[[i]],xpts=xpts,xscale=xscale)
    cat("\n")
  }
  
  return(adj_k_dists_list)
}

#clean up and select relevant stats
getCleanedupStats <- function(adj_k_dists_list) 
{
  img_adj_k_dists <- Reduce(plyr::rbind.fill,adj_k_dists_list) %>%
    rownames_to_column(var = "name") %>%
    as_tibble()
  
  img_adj_k_dists_select <- img_adj_k_dists %>%
    dplyr::select(name,m,m_r,m_c,A,Sc,St,Jc,Jt,m_dS,s_dS,cv_dS,m_dL,s_dL,cv_dL)
  
  return(img_adj_k_dists_select)
}

classifyColorPipeline <- function(path, kdf)
{
  classifications <- classifyByUniqueK(path, kdf)
  classified_k_dists <- suppressWarnings(getImgClassKDists(classifications, calcEucLumDists))
  adj_stats_raw <- suppressWarnings(getAdjStats(classifications, classified_k_dists, 100, 100))
  adj_stats <- suppressWarnings(getCleanedupStats(adj_stats_raw))
  
  image_paths <- kdf$img
  images_list <- rep(NA, length(image_paths))
  for(image in 1:length(images_list))
  {
    if(Sys.info()['sysname'] != "Windows")
    {
      images_list[image] <- tail(strsplit(image_paths[image], "/")[[1]], 1)
    }
    else
    {
      images_list[image] <- tail(strsplit(as.character(image_paths[image]), "\\\\")[[1]], 1)
    }
  }
  
  adj_stats <- adj_stats %>%
    mutate(name = images_list)
  
  return(adj_stats)
}

runColorPCA <- function(adj_stats)
{
  #pca_input <- select(adj_stats, name, m, m_r, m_c, Sc, St, A, m_dS,m_dL) %>% as.data.frame %>%
  pca_input <- select(adj_stats, name, m, A, Jc, Jt, m_dS, m_dL) %>% as.data.frame %>%
    column_to_rownames("name")
  pca_res <- prcomp(pca_input, center = T, scale = T)
  return(pca_res)
}

getColorPCASummary <- function(pca_res)
{
  return(summary(pca_res))
}

## load in data
fish <- read.csv("datasets/chaets_all_cleaned.csv", header = T, sep = ",")
birds <- read.csv("datasets/tanager_all_cleaned.csv", header = T, sep = ",")

fish_human <- subset(fish, method == "human")
fish_charisma <- subset(fish, method == "charisma")
fish_charisma <- subset(fish_charisma, k != 1)
  
bird_human <- subset(birds, method == "human")
bird_charisma <- subset(birds, method == "charisma")

constdf_fish <- data.frame(img = fish_human$img, k = rep(4, nrow(fish_human)))
constdf_birds <- data.frame(img = birds_human$img, k = rep(4, nrow(birds_human)))

fish_path <- "images/chaets/"
birds_path <- "images/tangara/"

## Run color pattern analysis
### fish
fish_classifications_human <- classifyColorPipeline(fish_path, fish_human)
fish_color_pca_human <- runColorPCA(fish_classifications_human)
fish_color_pca_human_summary <- getColorPCASummary(fish_color_pca_human)
write.csv(fish_classifications_human, "scripts/pca_results/fish/human/fish_classifications_human.csv", row.names = F)
write.csv(fish_color_pca_human$rotation, "scripts/pca_results/fish/human/fish_color_pca_human_rotation.csv", row.names = F)
sink("scripts/pca_results/fish/human/fish_color_pca_human_summary.txt")
fish_color_pca_human_summary
sink()
saveRDS(fish_classifications_human, "scripts/pca_results/fish/human/fish_classifications_human.RDS")
saveRDS(fish_color_pca_human, "scripts/pca_results/fish/human/fish_color_pca_human.RDS")
fish_color_pca_human <- readRDS("scripts/pca_results/fish/human/fish_color_pca_human.RDS")
saveRDS(fish_color_pca_human_summary, "scripts/pca_results/fish/human/fish_color_pca_human_summary.RDS")

fish_classifications_charisma <- classifyColorPipeline(fish_path, fish_charisma)
fish_color_pca_charisma <- runColorPCA(fish_classifications_charisma)
fish_color_pca_charisma_summary <- getColorPCASummary(fish_color_pca_charisma)
write.csv(fish_classifications_charisma, "scripts/pca_results/fish/charisma/fish_classifications_charisma.csv", row.names = F)
write.csv(fish_color_pca_charisma$rotation, "scripts/pca_results/fish/charisma/fish_color_pca_charisma_rotation.csv", row.names = F)
sink("scripts/pca_results/fish/charisma/fish_color_pca_charisma_summary.txt")
fish_color_pca_charisma_summary
sink()
saveRDS(fish_classifications_charisma, "scripts/pca_results/fish/charisma/fish_classifications_charisma.RDS")
saveRDS(fish_color_pca_charisma, "scripts/pca_results/fish/charisma/fish_color_pca_charisma.RDS")
fish_color_pca_charisma <- readRDS("scripts/pca_results/fish/charisma/fish_color_pca_charisma.RDS")
saveRDS(fish_color_pca_charisma_summary, "scripts/pca_results/fish/charisma/fish_color_pca_charisma_summary.RDS")

fish_classifications_const <- classifyColorPipeline(fish_path, constdf_fish)
fish_color_pca_const <- runColorPCA(fish_classifications_const)
fish_color_pca_const_summary <- getColorPCASummary(fish_color_pca_const)
write.csv(fish_classifications_const, "scripts/pca_results/fish/const/fish_classifications_const.csv", row.names = F)
write.csv(fish_color_pca_const$rotation, "scripts/pca_results/fish/const/fish_color_pca_const_rotation.csv", row.names = F)
sink("scripts/pca_results/fish/const/fish_color_pca_const_summary.txt")
fish_color_pca_const_summary
sink()
saveRDS(fish_classifications_const, "scripts/pca_results/fish/const/fish_classifications_const.RDS")
saveRDS(fish_color_pca_const, "scripts/pca_results/fish/const/fish_color_pca_const.RDS")
fish_color_pca_const <- readRDS("scripts/pca_results/fish/const/fish_color_pca_const.RDS")
saveRDS(fish_color_pca_const_summary, "scripts/pca_results/fish/const/fish_color_pca_const_summary.RDS")

### bird
birds_classifications_human <- classifyColorPipeline(birds_path, birds_human)
birds_color_pca_human <- runColorPCA(birds_classifications_human)
birds_color_pca_human_summary <- getColorPCASummary(birds_color_pca_human)
write.csv(birds_classifications_human, "scripts/pca_results/birds/human/birds_classifications_human.csv", row.names = F)
write.csv(birds_color_pca_human$rotation, "scripts/pca_results/birds/human/birds_color_pca_human_rotation.csv", row.names = F)
sink("scripts/pca_results/birds/human/birds_color_pca_human_summary.txt")
birds_color_pca_human_summary
sink()
saveRDS(birds_classifications_human, "scripts/pca_results/birds/human/birds_classifications_human.RDS")
saveRDS(birds_color_pca_human, "scripts/pca_results/birds/human/birds_color_pca_human.RDS")
birds_color_pca_human <- readRDS("scripts/pca_results/birds/human/birds_color_pca_human.RDS")
saveRDS(birds_color_pca_human_summary, "scripts/pca_results/birds/human/birds_color_pca_human_summary.RDS")

birds_classifications_charisma <- classifyColorPipeline(birds_path, birds_charisma)
birds_color_pca_charisma <- runColorPCA(birds_classifications_charisma)
birds_color_pca_charisma_summary <- getColorPCASummary(birds_color_pca_charisma)
write.csv(birds_classifications_charisma, "scripts/pca_results/birds/charisma/birds_classifications_charisma.csv", row.names = F)
write.csv(birds_color_pca_charisma$rotation, "scripts/pca_results/birds/charisma/birds_color_pca_charisma_rotation.csv", row.names = F)
sink("scripts/pca_results/birds/charisma/birds_color_pca_charisma_summary.txt")
birds_color_pca_charisma_summary
sink()
saveRDS(birds_classifications_charisma, "scripts/pca_results/birds/charisma/birds_classifications_charisma.RDS")
saveRDS(birds_color_pca_charisma, "scripts/pca_results/birds/charisma/birds_color_pca_charisma.RDS")
birds_color_pca_charisma <- readRDS("scripts/pca_results/birds/charisma/birds_color_pca_charisma.RDS")
saveRDS(birds_color_pca_charisma_summary, "scripts/pca_results/birds/charisma/birds_color_pca_charisma_summary.RDS")

birds_classifications_const <- classifyColorPipeline(birds_path, constdf_birds)
birds_color_pca_const <- runColorPCA(birds_classifications_const)
birds_color_pca_const_summary <- getColorPCASummary(birds_color_pca_const)
write.csv(birds_classifications_const, "scripts/pca_results/birds/const/birds_classifications_const.csv", row.names = F)
write.csv(birds_color_pca_const$rotation, "scripts/pca_results/birds/const/birds_color_pca_const_rotation.csv", row.names = F)
sink("scripts/pca_results/birds/const/birds_color_pca_const_summary.txt")
birds_color_pca_const_summary
sink()
saveRDS(birds_classifications_const, "scripts/pca_results/birds/const/birds_classifications_const.RDS")
saveRDS(birds_color_pca_const, "scripts/pca_results/birds/const/birds_color_pca_const.RDS")
birds_color_pca_const <- readRDS("scripts/pca_results/birds/const/birds_color_pca_const.RDS")
saveRDS(birds_color_pca_const_summary, "scripts/pca_results/birds/const/birds_color_pca_const_summary.RDS")

### analyses
## run permutation tests for (birds & fish)
run_perm_t.test <- function(human, charisma, alldata, compvar, comp1 = "human", comp2 = "charisma", iter = 10000) {
  obs_mean_human_k <- mean(human[[compvar]])
  obs_mean_charisma_k <- mean(charisma[[compvar]])
  obs_mean_diff_k <- obs_mean_human_k - obs_mean_charisma_k
  
  #print(paste0("Mean ", comp1, " (", compvar, ") = ", obs_mean_human_k))
  #print(paste0("Mean ", comp2, " (", compvar, ") = ", obs_mean_charisma_k))
  
  perm_means <- rep(NA, iter)
  for(ii in 1:iter) {
    shuffled_data <- sample(alldata$method)
    perm_diff <- mean(alldata[[compvar]][shuffled_data == comp1]) - mean(alldata[[compvar]][shuffled_data == comp2])
    perm_means[ii] <- perm_diff
  }
  #two-tailed test
  p_two <- (sum(perm_means >= abs(obs_mean_diff_k)) + sum(perm_means < (-abs(obs_mean_diff_k)))) / iter
  return(p_two)
}
#### FISH ####
## a) human vs. charsima
fish_classifications_human_lab <- fish_classifications_human %>%
  select(name, m, A, Jc, Jt, m_dS, m_dL) %>%
  mutate(method = "human")
fish_classifications_charisma_lab <- fish_classifications_charisma %>%
  select(name, m, A, Jc, Jt, m_dS, m_dL) %>%
  mutate(method = "charisma")
perm_human_charisma_dat <- rbind(fish_classifications_human_lab, fish_classifications_charisma_lab)
perm_human_charisma_dat

##compute difference of means
head(fish_classifications_charisma_lab)
head(fish_classifications_human_lab)

round(mean(fish_classifications_charisma_lab$m) - mean(fish_classifications_human_lab$m), 2)
round(mean(fish_classifications_charisma_lab$A) - mean(fish_classifications_human_lab$A), 2)
round(mean(fish_classifications_charisma_lab$Jc) - mean(fish_classifications_human_lab$Jc), 2)
round(mean(fish_classifications_charisma_lab$Jt) - mean(fish_classifications_human_lab$Jt), 2)
round(mean(fish_classifications_charisma_lab$m_dS) - mean(fish_classifications_human_lab$m_dS), 2)
round(mean(fish_classifications_charisma_lab$m_dL) - mean(fish_classifications_human_lab$m_dL), 2)

#constants
iter_runs <- 100
iter <- 10000

sem <- function(df) {
  return(sd(df)/sqrt(length(df)))
}
## m
perm_human_charisma_pruns_m <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_m[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                     "m", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_m, "scripts/perm_results/perm_human_charisma_pruns_m.RDS")
perm_human_charisma_pruns_m <- readRDS("scripts/perm_results/fish/perm_human_charisma_pruns_m.RDS")
mean(perm_human_charisma_pruns_m)
upper_CI_emp <- quantile(perm_human_charisma_pruns_m, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_m, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_m.pdf")
hist(perm_human_charisma_pruns_m)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_m), lwd = 3, col = "blue")
dev.off()

## A
perm_human_charisma_pruns_A <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_A[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                     "A", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_A, "scripts/perm_results/perm_human_charisma_pruns_A.RDS")
perm_human_charisma_pruns_A <- readRDS("scripts/perm_results/perm_human_charisma_pruns_A.RDS")
mean(perm_human_charisma_pruns_A)
upper_CI_emp <- quantile(perm_human_charisma_pruns_A, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_A, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_A.pdf")
hist(perm_human_charisma_pruns_A)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_A), lwd = 3, col = "blue")
dev.off()

## Jc
perm_human_charisma_pruns_Jc <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_Jc[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                     "Jc", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_Jc, "scripts/perm_results/perm_human_charisma_pruns_Jc.RDS")
perm_human_charisma_pruns_Jc <- readRDS("scripts/perm_results/perm_human_charisma_pruns_Jc.RDS")
mean(perm_human_charisma_pruns_Jc)
upper_CI_emp <- quantile(perm_human_charisma_pruns_Jc, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_Jc, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_Jc.pdf")
hist(perm_human_charisma_pruns_Jc)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_Jc), lwd = 3, col = "blue")
dev.off()

## Jt
perm_human_charisma_pruns_Jt <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_Jt[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                     "Jt", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_Jt, "scripts/perm_results/perm_human_charisma_pruns_Jt.RDS")
perm_human_charisma_pruns_Jt <- readRDS("scripts/perm_results/perm_human_charisma_pruns_Jt.RDS")
mean(perm_human_charisma_pruns_Jt)
upper_CI_emp <- quantile(perm_human_charisma_pruns_Jt, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_Jt, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_Jt.pdf")
hist(perm_human_charisma_pruns_Jt)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_Jt), lwd = 3, col = "blue")
dev.off()


## m_dS
perm_human_charisma_pruns_m_dS <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_m_dS[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                     "m_dS", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_m_dS, "scripts/perm_results/perm_human_charisma_pruns_m_dS.RDS")
perm_human_charisma_pruns_m_dS <- readRDS("scripts/perm_results/perm_human_charisma_pruns_m_dS.RDS")
mean(perm_human_charisma_pruns_m_dS)
upper_CI_emp <- quantile(perm_human_charisma_pruns_m_dS, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_m_dS, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_m_dS.pdf")
hist(perm_human_charisma_pruns_m_dS)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_m_dS), lwd = 3, col = "blue")
dev.off()

## m_dL
perm_human_charisma_pruns_m_dL <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_m_dL[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                     "m_dL", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_m_dL, "scripts/perm_results/perm_human_charisma_pruns_m_dL.RDS")
perm_human_charisma_pruns_m_dL <- readRDS("scripts/perm_results/perm_human_charisma_pruns_m_dL.RDS")
mean(perm_human_charisma_pruns_m_dL)
upper_CI_emp <- quantile(perm_human_charisma_pruns_m_dL, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_m_dL, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_m_dL.pdf")
hist(perm_human_charisma_pruns_m_dL)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_m_dL), lwd = 3, col = "blue")
dev.off()


## b) human vs. const
fish_classifications_const_lab <- fish_classifications_const %>%
  select(name, m, A, Jc, Jt, m_dS, m_dL) %>%
  mutate(method = "const")
perm_human_const_dat <- rbind(fish_classifications_human_lab, fish_classifications_const_lab)
perm_human_const_dat

## m
perm_human_const_pruns_m <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_m[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                     "m", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_m, "scripts/perm_results/perm_human_const_pruns_m.RDS")
perm_human_const_pruns_m <- readRDS("scripts/perm_results/perm_human_const_pruns_m.RDS")
mean(perm_human_const_pruns_m)
upper_CI_emp <- quantile(perm_human_const_pruns_m, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_m, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_m.pdf")
hist(perm_human_const_pruns_m)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_m), lwd = 3, col = "blue")
dev.off()

## A
perm_human_const_pruns_A <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_A[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                     "A", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_A, "scripts/perm_results/perm_human_const_pruns_A.RDS")
perm_human_const_pruns_A <- readRDS("scripts/perm_results/perm_human_const_pruns_A.RDS")
mean(perm_human_const_pruns_A)
upper_CI_emp <- quantile(perm_human_const_pruns_A, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_A, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_A.pdf")
hist(perm_human_const_pruns_A)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_A), lwd = 3, col = "blue")
dev.off()

## Jc
perm_human_const_pruns_Jc <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_Jc[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                      "Jc", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_Jc, "scripts/perm_results/perm_human_const_pruns_Jc.RDS")
perm_human_const_pruns_Jc <- readRDS("scripts/perm_results/perm_human_const_pruns_Jc.RDS")
mean(perm_human_const_pruns_Jc)
upper_CI_emp <- quantile(perm_human_const_pruns_Jc, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_Jc, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_Jc.pdf")
hist(perm_human_const_pruns_Jc)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_Jc), lwd = 3, col = "blue")
dev.off()

## Jt
perm_human_const_pruns_Jt <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_Jt[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                      "Jt", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_Jt, "scripts/perm_results/perm_human_const_pruns_Jt.RDS")
perm_human_const_pruns_Jt <- readRDS("scripts/perm_results/perm_human_const_pruns_Jt.RDS")
mean(perm_human_const_pruns_Jt)
upper_CI_emp <- quantile(perm_human_const_pruns_Jt, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_Jt, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_Jt.pdf")
hist(perm_human_const_pruns_Jt)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_Jt), lwd = 3, col = "blue")
dev.off()


## m_dS
perm_human_const_pruns_m_dS <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_m_dS[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                        "m_dS", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_m_dS, "scripts/perm_results/perm_human_const_pruns_m_dS.RDS")
perm_human_const_pruns_m_dS <- readRDS("scripts/perm_results/perm_human_const_pruns_m_dS.RDS")
mean(perm_human_const_pruns_m_dS)
upper_CI_emp <- quantile(perm_human_const_pruns_m_dS, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_m_dS, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_m_dS.pdf")
hist(perm_human_const_pruns_m_dS)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_m_dS), lwd = 3, col = "blue")
dev.off()

## m_dL
perm_human_const_pruns_m_dL <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_m_dL[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                        "m_dL", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_m_dL, "scripts/perm_results/perm_human_const_pruns_m_dL.RDS")
perm_human_const_pruns_m_dL <- readRDS("scripts/perm_results/perm_human_const_pruns_m_dL.RDS")
mean(perm_human_const_pruns_m_dL)
upper_CI_emp <- quantile(perm_human_const_pruns_m_dL, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_m_dL, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_m_dL.pdf")
hist(perm_human_const_pruns_m_dL)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_m_dL), lwd = 3, col = "blue")
dev.off()

## c) charisma vs. const
perm_charisma_const_dat <- rbind(fish_classifications_charisma_lab, fish_classifications_const_lab)
perm_charisma_const_dat

round(mean(fish_classifications_charisma_lab$m) - mean(fish_classifications_human_lab$m), 2)
round(mean(fish_classifications_charisma_lab$A) - mean(fish_classifications_human_lab$A), 2)
round(mean(fish_classifications_charisma_lab$Jc) - mean(fish_classifications_human_lab$Jc), 2)
round(mean(fish_classifications_charisma_lab$Jt) - mean(fish_classifications_human_lab$Jt), 2)
round(mean(fish_classifications_charisma_lab$m_dS) - mean(fish_classifications_human_lab$m_dS), 2)
round(mean(fish_classifications_charisma_lab$m_dL) - mean(fish_classifications_human_lab$m_dL), 2) ##check perm test

round(mean(fish_classifications_charisma_lab$m) - mean(fish_classifications_const_lab$m), 2)
round(mean(fish_classifications_charisma_lab$A) - mean(fish_classifications_const_lab$A), 2)
round(mean(fish_classifications_charisma_lab$Jc) - mean(fish_classifications_const_lab$Jc), 2)
round(mean(fish_classifications_charisma_lab$Jt) - mean(fish_classifications_const_lab$Jt), 2)
round(mean(fish_classifications_charisma_lab$m_dS) - mean(fish_classifications_const_lab$m_dS), 2)
round(mean(fish_classifications_charisma_lab$m_dL) - mean(fish_classifications_const_lab$m_dL), 2)

round(mean(fish_classifications_human_lab$m) - mean(fish_classifications_const_lab$m), 2)
round(mean(fish_classifications_human_lab$A) - mean(fish_classifications_const_lab$A), 2)
round(mean(fish_classifications_human_lab$Jc) - mean(fish_classifications_const_lab$Jc), 2)
round(mean(fish_classifications_human_lab$Jt) - mean(fish_classifications_const_lab$Jt), 2)
round(mean(fish_classifications_human_lab$m_dS) - mean(fish_classifications_const_lab$m_dS), 2)
round(mean(fish_classifications_human_lab$m_dL) - mean(fish_classifications_const_lab$m_dL), 2)

## m
perm_charisma_const_pruns_m <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_m[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                  "m", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_m, "scripts/perm_results/perm_charisma_const_pruns_m.RDS")
perm_charisma_const_pruns_m <- readRDS("scripts/perm_results/perm_charisma_const_pruns_m.RDS")
mean(perm_charisma_const_pruns_m)
upper_CI_emp <- quantile(perm_charisma_const_pruns_m, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_m, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_m.pdf")
hist(perm_charisma_const_pruns_m)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_m), lwd = 3, col = "blue")
dev.off()

## A
perm_charisma_const_pruns_A <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_A[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                  "A", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_A, "scripts/perm_results/perm_charisma_const_pruns_A.RDS")
perm_charisma_const_pruns_A <- readRDS("scripts/perm_results/perm_charisma_const_pruns_A.RDS")
mean(perm_charisma_const_pruns_A)
upper_CI_emp <- quantile(perm_charisma_const_pruns_A, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_A, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_A.pdf")
hist(perm_charisma_const_pruns_A)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_A), lwd = 3, col = "blue")
dev.off()

## Jc
perm_charisma_const_pruns_Jc <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_Jc[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                   "Jc", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_Jc, "scripts/perm_results/perm_charisma_const_pruns_Jc.RDS")
perm_charisma_const_pruns_Jc <- readRDS("scripts/perm_results/perm_charisma_const_pruns_Jc.RDS")
mean(perm_charisma_const_pruns_Jc)
upper_CI_emp <- quantile(perm_charisma_const_pruns_Jc, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_Jc, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_Jc.pdf")
hist(perm_charisma_const_pruns_Jc)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_Jc), lwd = 3, col = "blue")
dev.off()

## Jt
perm_charisma_const_pruns_Jt <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_Jt[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                   "Jt", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_Jt, "scripts/perm_results/perm_charisma_const_pruns_Jt.RDS")
perm_charisma_const_pruns_Jt <- readRDS("scripts/perm_results/perm_charisma_const_pruns_Jt.RDS")
mean(perm_charisma_const_pruns_Jt)
upper_CI_emp <- quantile(perm_charisma_const_pruns_Jt, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_Jt, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_Jt.pdf")
hist(perm_charisma_const_pruns_Jt)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_Jt), lwd = 3, col = "blue")
dev.off()


## m_dS
perm_charisma_const_pruns_m_dS <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_m_dS[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                        "m_dS", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_m_dS, "scripts/perm_results/perm_charisma_const_pruns_m_dS.RDS")
perm_charisma_const_pruns_m_dS <- readRDS("scripts/perm_results/perm_charisma_const_pruns_m_dS.RDS")
mean(perm_charisma_const_pruns_m_dS)
upper_CI_emp <- quantile(perm_charisma_const_pruns_m_dS, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_m_dS, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_m_dS.pdf")
hist(perm_charisma_const_pruns_m_dS)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_m_dS), lwd = 3, col = "blue")
dev.off()

## m_dL
perm_charisma_const_pruns_m_dL <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_m_dL[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                        "m_dL", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_m_dL, "scripts/perm_results/perm_charisma_const_pruns_m_dL.RDS")
perm_charisma_const_pruns_m_dL <- readRDS("scripts/perm_results/perm_charisma_const_pruns_m_dL.RDS")
mean(perm_charisma_const_pruns_m_dL)
upper_CI_emp <- quantile(perm_charisma_const_pruns_m_dL, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_m_dL, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_m_dL.pdf")
hist(perm_charisma_const_pruns_m_dL)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_m_dL), lwd = 3, col = "blue")
dev.off()



round(mean(fish_classifications_charisma_lab$m) - mean(fish_classifications_const_lab$m), 2)
round(mean(fish_classifications_charisma_lab$A) - mean(fish_classifications_const_lab$A), 2)
round(mean(fish_classifications_charisma_lab$Jc) - mean(fish_classifications_const_lab$Jc), 2)
round(mean(fish_classifications_charisma_lab$Jt) - mean(fish_classifications_const_lab$Jt), 2)
round(mean(fish_classifications_charisma_lab$m_dS) - mean(fish_classifications_const_lab$m_dS), 2)
round(mean(fish_classifications_charisma_lab$m_dL) - mean(fish_classifications_const_lab$m_dL), 2)

round(mean(fish_classifications_human_lab$m) - mean(fish_classifications_const_lab$m), 2)
round(mean(fish_classifications_human_lab$A) - mean(fish_classifications_const_lab$A), 2)
round(mean(fish_classifications_human_lab$Jc) - mean(fish_classifications_const_lab$Jc), 2)
round(mean(fish_classifications_human_lab$Jt) - mean(fish_classifications_const_lab$Jt), 2)
round(mean(fish_classifications_human_lab$m_dS) - mean(fish_classifications_const_lab$m_dS), 2)
round(mean(fish_classifications_human_lab$m_dL) - mean(fish_classifications_const_lab$m_dL), 2)


#### BIRDS ####
## a) human vs. charsima
fish_classifications_human_lab <- birds_classifications_human %>%
  select(name, m, A, Jc, Jt, m_dS, m_dL) %>%
  mutate(method = "human")
fish_classifications_charisma_lab <- birds_classifications_charisma %>%
  select(name, m, A, Jc, Jt, m_dS, m_dL) %>%
  mutate(method = "charisma")
perm_human_charisma_dat <- rbind(fish_classifications_human_lab, fish_classifications_charisma_lab)
perm_human_charisma_dat

#constants
iter_runs <- 100
iter <- 10000

## m
perm_human_charisma_pruns_m <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_m[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                     "m", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_m, "scripts/perm_results/perm_human_charisma_pruns_m.RDS")
perm_human_charisma_pruns_m <- readRDS("scripts/perm_results/perm_human_charisma_pruns_m.RDS")
mean(perm_human_charisma_pruns_m)
upper_CI_emp <- quantile(perm_human_charisma_pruns_m, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_m, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_m.pdf")
hist(perm_human_charisma_pruns_m)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_m), lwd = 3, col = "blue")
dev.off()

## A
perm_human_charisma_pruns_A <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_A[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                     "A", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_A, "scripts/perm_results/perm_human_charisma_pruns_A.RDS")
perm_human_charisma_pruns_A <- readRDS("scripts/perm_results/perm_human_charisma_pruns_A.RDS")
mean(perm_human_charisma_pruns_A)
upper_CI_emp <- quantile(perm_human_charisma_pruns_A, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_A, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_A.pdf")
hist(perm_human_charisma_pruns_A)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_A), lwd = 3, col = "blue")
dev.off()

## Jc
perm_human_charisma_pruns_Jc <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_Jc[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                      "Jc", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_Jc, "scripts/perm_results/perm_human_charisma_pruns_Jc.RDS")
perm_human_charisma_pruns_Jc <- readRDS("scripts/perm_results/perm_human_charisma_pruns_Jc.RDS")
mean(perm_human_charisma_pruns_Jc)
upper_CI_emp <- quantile(perm_human_charisma_pruns_Jc, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_Jc, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_Jc.pdf")
hist(perm_human_charisma_pruns_Jc)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_Jc), lwd = 3, col = "blue")
dev.off()

## Jt
perm_human_charisma_pruns_Jt <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_Jt[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                      "Jt", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_Jt, "scripts/perm_results/perm_human_charisma_pruns_Jt.RDS")
perm_human_charisma_pruns_Jt <- readRDS("scripts/perm_results/perm_human_charisma_pruns_Jt.RDS")
mean(perm_human_charisma_pruns_Jt)
upper_CI_emp <- quantile(perm_human_charisma_pruns_Jt, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_Jt, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_Jt.pdf")
hist(perm_human_charisma_pruns_Jt)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_Jt), lwd = 3, col = "blue")
dev.off()


## m_dS
perm_human_charisma_pruns_m_dS <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_m_dS[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                        "m_dS", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_m_dS, "scripts/perm_results/perm_human_charisma_pruns_m_dS.RDS")
perm_human_charisma_pruns_m_dS <- readRDS("scripts/perm_results/perm_human_charisma_pruns_m_dS.RDS")
mean(perm_human_charisma_pruns_m_dS)
upper_CI_emp <- quantile(perm_human_charisma_pruns_m_dS, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_m_dS, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_m_dS.pdf")
hist(perm_human_charisma_pruns_m_dS)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_m_dS), lwd = 3, col = "blue")
dev.off()

## m_dL
perm_human_charisma_pruns_m_dL <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_charisma_pruns_m_dL[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_charisma_lab, perm_human_charisma_dat,
                                                        "m_dL", "human", "charisma", iter)
}
saveRDS(perm_human_charisma_pruns_m_dL, "scripts/perm_results/perm_human_charisma_pruns_m_dL.RDS")
perm_human_charisma_pruns_m_dL <- readRDS("scripts/perm_results/perm_human_charisma_pruns_m_dL.RDS")
mean(perm_human_charisma_pruns_m_dL)
upper_CI_emp <- quantile(perm_human_charisma_pruns_m_dL, 0.975)
lower_CI_emp <- quantile(perm_human_charisma_pruns_m_dL, 0.025)
pdf("scripts/perm_results/perm_human_charisma_pruns_m_dL.pdf")
hist(perm_human_charisma_pruns_m_dL)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_charisma_pruns_m_dL), lwd = 3, col = "blue")
dev.off()


## b) human vs. const
fish_classifications_const_lab <- birds_classifications_const %>%
  select(name, m, A, Jc, Jt, m_dS, m_dL) %>%
  mutate(method = "const")
perm_human_const_dat <- rbind(fish_classifications_human_lab, fish_classifications_const_lab)
perm_human_const_dat




## m
perm_human_const_pruns_m <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_m[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                  "m", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_m, "scripts/perm_results/perm_human_const_pruns_m.RDS")
perm_human_const_pruns_m <- readRDS("scripts/perm_results/perm_human_const_pruns_m.RDS")
mean(perm_human_const_pruns_m)
upper_CI_emp <- quantile(perm_human_const_pruns_m, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_m, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_m.pdf")
hist(perm_human_const_pruns_m)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_m), lwd = 3, col = "blue")
dev.off()

## A
perm_human_const_pruns_A <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_A[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                  "A", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_A, "scripts/perm_results/perm_human_const_pruns_A.RDS")
perm_human_const_pruns_A <- readRDS("scripts/perm_results/perm_human_const_pruns_A.RDS")
mean(perm_human_const_pruns_A)
upper_CI_emp <- quantile(perm_human_const_pruns_A, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_A, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_A.pdf")
hist(perm_human_const_pruns_A)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_A), lwd = 3, col = "blue")
dev.off()

## Jc
perm_human_const_pruns_Jc <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_Jc[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                   "Jc", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_Jc, "scripts/perm_results/perm_human_const_pruns_Jc.RDS")
perm_human_const_pruns_Jc <- readRDS("scripts/perm_results/perm_human_const_pruns_Jc.RDS")
mean(perm_human_const_pruns_Jc)
upper_CI_emp <- quantile(perm_human_const_pruns_Jc, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_Jc, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_Jc.pdf")
hist(perm_human_const_pruns_Jc)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_Jc), lwd = 3, col = "blue")
dev.off()

## Jt
perm_human_const_pruns_Jt <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_Jt[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                   "Jt", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_Jt, "scripts/perm_results/perm_human_const_pruns_Jt.RDS")
perm_human_const_pruns_Jt <- readRDS("scripts/perm_results/perm_human_const_pruns_Jt.RDS")
mean(perm_human_const_pruns_Jt)
upper_CI_emp <- quantile(perm_human_const_pruns_Jt, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_Jt, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_Jt.pdf")
hist(perm_human_const_pruns_Jt)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_Jt), lwd = 3, col = "blue")
dev.off()


## m_dS
perm_human_const_pruns_m_dS <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_m_dS[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                        "m_dS", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_m_dS, "scripts/perm_results/perm_human_const_pruns_m_dS.RDS")
perm_human_const_pruns_m_dS <- readRDS("scripts/perm_results/perm_human_const_pruns_m_dS.RDS")
mean(perm_human_const_pruns_m_dS)
upper_CI_emp <- quantile(perm_human_const_pruns_m_dS, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_m_dS, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_m_dS.pdf")
hist(perm_human_const_pruns_m_dS)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_m_dS), lwd = 3, col = "blue")
dev.off()

## m_dL
perm_human_const_pruns_m_dL <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_human_const_pruns_m_dL[ii] <- run_perm_t.test(fish_classifications_human_lab, fish_classifications_const_lab, perm_human_const_dat,
                                                        "m_dL", "human", "const", iter)
}
saveRDS(perm_human_const_pruns_m_dL, "scripts/perm_results/perm_human_const_pruns_m_dL.RDS")
perm_human_const_pruns_m_dL <- readRDS("scripts/perm_results/perm_human_const_pruns_m_dL.RDS")
mean(perm_human_const_pruns_m_dL)
upper_CI_emp <- quantile(perm_human_const_pruns_m_dL, 0.975)
lower_CI_emp <- quantile(perm_human_const_pruns_m_dL, 0.025)
pdf("scripts/perm_results/perm_human_const_pruns_m_dL.pdf")
hist(perm_human_const_pruns_m_dL)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_human_const_pruns_m_dL), lwd = 3, col = "blue")
dev.off()

## c) charisma vs. const
perm_charisma_const_dat <- rbind(fish_classifications_charisma_lab, fish_classifications_const_lab)
perm_charisma_const_dat

## m
perm_charisma_const_pruns_m <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_m[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                     "m", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_m, "scripts/perm_results/perm_charisma_const_pruns_m.RDS")
perm_charisma_const_pruns_m <- readRDS("scripts/perm_results/perm_charisma_const_pruns_m.RDS")
mean(perm_charisma_const_pruns_m)
upper_CI_emp <- quantile(perm_charisma_const_pruns_m, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_m, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_m.pdf")
hist(perm_charisma_const_pruns_m)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_m), lwd = 3, col = "blue")
dev.off()

## A
perm_charisma_const_pruns_A <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_A[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                     "A", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_A, "scripts/perm_results/perm_charisma_const_pruns_A.RDS")
perm_charisma_const_pruns_A <- readRDS("scripts/perm_results/perm_charisma_const_pruns_A.RDS")
mean(perm_charisma_const_pruns_A)
upper_CI_emp <- quantile(perm_charisma_const_pruns_A, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_A, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_A.pdf")
hist(perm_charisma_const_pruns_A)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_A), lwd = 3, col = "blue")
dev.off()

## Jc
perm_charisma_const_pruns_Jc <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_Jc[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                      "Jc", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_Jc, "scripts/perm_results/perm_charisma_const_pruns_Jc.RDS")
perm_charisma_const_pruns_Jc <- readRDS("scripts/perm_results/perm_charisma_const_pruns_Jc.RDS")
mean(perm_charisma_const_pruns_Jc)
upper_CI_emp <- quantile(perm_charisma_const_pruns_Jc, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_Jc, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_Jc.pdf")
hist(perm_charisma_const_pruns_Jc)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_Jc), lwd = 3, col = "blue")
dev.off()

## Jt
perm_charisma_const_pruns_Jt <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_Jt[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                      "Jt", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_Jt, "scripts/perm_results/perm_charisma_const_pruns_Jt.RDS")
perm_charisma_const_pruns_Jt <- readRDS("scripts/perm_results/perm_charisma_const_pruns_Jt.RDS")
mean(perm_charisma_const_pruns_Jt)
upper_CI_emp <- quantile(perm_charisma_const_pruns_Jt, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_Jt, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_Jt.pdf")
hist(perm_charisma_const_pruns_Jt)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_Jt), lwd = 3, col = "blue")
dev.off()


## m_dS
perm_charisma_const_pruns_m_dS <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_m_dS[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                        "m_dS", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_m_dS, "scripts/perm_results/perm_charisma_const_pruns_m_dS.RDS")
perm_charisma_const_pruns_m_dS <- readRDS("scripts/perm_results/perm_charisma_const_pruns_m_dS.RDS")
mean(perm_charisma_const_pruns_m_dS)
upper_CI_emp <- quantile(perm_charisma_const_pruns_m_dS, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_m_dS, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_m_dS.pdf")
hist(perm_charisma_const_pruns_m_dS)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_m_dS), lwd = 3, col = "blue")
dev.off()

## m_dL
perm_charisma_const_pruns_m_dL <- rep(NA, iter_runs)
for(ii in 1:iter_runs) {
  print(paste0("Run #", ii, "/", iter_runs))
  perm_charisma_const_pruns_m_dL[ii] <- run_perm_t.test(fish_classifications_charisma_lab, fish_classifications_const_lab, perm_charisma_const_dat,
                                                        "m_dL", "charisma", "const", iter)
}
saveRDS(perm_charisma_const_pruns_m_dL, "scripts/perm_results/perm_charisma_const_pruns_m_dL.RDS")
perm_charisma_const_pruns_m_dL <- readRDS("scripts/perm_results/perm_charisma_const_pruns_m_dL.RDS")
mean(perm_charisma_const_pruns_m_dL)
upper_CI_emp <- quantile(perm_charisma_const_pruns_m_dL, 0.975)
lower_CI_emp <- quantile(perm_charisma_const_pruns_m_dL, 0.025)
pdf("scripts/perm_results/perm_charisma_const_pruns_m_dL.pdf")
hist(perm_charisma_const_pruns_m_dL)
abline(v = lower_CI_emp, lwd = 3, col = "red")
abline(v = upper_CI_emp, lwd = 3, col = "red")
abline(v = mean(perm_charisma_const_pruns_m_dL), lwd = 3, col = "blue")
dev.off()



## make PCA plot using PCA summaries above
## check to make sure that charisma isn't different from human, then use charisma
##  in all downstream analyses
library(factoextra)
library(ggimage)
library(ggpubr)

## FISH ##
pdf("scripts/pca_results/fish/plots/fish_charisma_pca.pdf")
fviz_pca_var(fish_color_pca_charisma,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

pdf("scripts/pca_results/fish/plots/fish_charisma_pca_screeplot.pdf")
fviz_screeplot(fish_color_pca_charisma, addlabels = TRUE, ylim = c(0, 50))
dev.off()

# Contributions of variables to PC1
pc1_contrib <- fviz_contrib(fish_color_pca_charisma, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
pc2_contrib <- fviz_contrib(fish_color_pca_charisma, choice = "var", axes = 2, top = 10)
contrib <- ggarrange(pc1_contrib, pc2_contrib, 
          labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("scripts/pca_results/fish/plots/fish_charisma_pca_contrib.pdf")

img_paths_fish_charisma <- paste0("images/chaets/", rownames(as.data.frame(fish_color_pca_charisma$x)))
scores <- as.data.frame(fish_color_pca_charisma$x)
loadings_fish_charisma <- as.data.frame(fish_color_pca_charisma$rotation) 
p_fish_pca_charisma <- ggplot(as.data.frame(fish_color_pca_charisma$x), aes(PC1, PC2)) +
  geom_image(aes(image = img_paths_fish_charisma), size = .05, alpha = .5) +
  xlab("PC1 (43.4%)") +
  ylab("PC2 (20.2%)") +
  ggtitle("charisma-k") +
  theme_minimal() +
  geom_segment(data = loadings_fish_charisma, aes(x=0,y=0,xend=PC1,yend=PC2),
    arrow=arrow(length=unit(0.1,"cm")), color = "#ff0000") +
  geom_text(data = loadings_fish_charisma, aes(x=PC1, y=PC2, label=rownames(loadings_fish_charisma)),color="#ff0000", vjust = -0.25)
ggsave("scripts/pca_results/fish/plots/fish_charisma_pca_pc1pc2_2.pdf")






pdf("scripts/pca_results/fish/plots/fish_human_pca.pdf")
fviz_pca_var(fish_color_pca_human,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

pdf("scripts/pca_results/fish/plots/fish_human_pca_screeplot.pdf")
fviz_screeplot(fish_color_pca_human, addlabels = TRUE, ylim = c(0, 50))
dev.off()

# Contributions of variables to PC1
pc1_contrib <- fviz_contrib(fish_color_pca_human, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
pc2_contrib <- fviz_contrib(fish_color_pca_human, choice = "var", axes = 2, top = 10)
contrib <- ggarrange(pc1_contrib, pc2_contrib, 
                     labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("scripts/pca_results/fish/plots/fish_human_pca_contrib.pdf")

img_paths_fish_human <- paste0("images/chaets/", rownames(as.data.frame(fish_color_pca_human$x)))
scores <- as.data.frame(fish_color_pca_human$x)
loadings_fish_human <- as.data.frame(fish_color_pca_human$rotation) 
p_fish_pca_human <- ggplot(as.data.frame(fish_color_pca_human$x), aes(PC1, PC2)) +
  geom_image(aes(image = img_paths_fish_human), size = .05, alpha = .5) +
  xlab("PC1 (38.6%)") +
  ylab("PC2 (22.5%)") +
  ggtitle("human-k") +
  theme_minimal() +
  geom_segment(data = loadings_fish_human, aes(x=0,y=0,xend=PC1,yend=PC2),
               arrow=arrow(length=unit(0.1,"cm")), color = "#ff0000") +
  geom_text(data = loadings_fish_human, aes(x=PC1, y=PC2, label=rownames(loadings_fish_human)),color="#ff0000", vjust = -0.25)
ggsave("scripts/pca_results/fish/plots/fish_human_pca_pc1pc2_2.pdf")





pdf("scripts/pca_results/fish/plots/fish_const_pca.pdf")
fviz_pca_var(fish_color_pca_const,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

pdf("scripts/pca_results/fish/plots/fish_const_pca_screeplot.pdf")
fviz_screeplot(fish_color_pca_const, addlabels = TRUE, ylim = c(0, 50))
dev.off()

# Contributions of variables to PC1
pc1_contrib <- fviz_contrib(fish_color_pca_const, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
pc2_contrib <- fviz_contrib(fish_color_pca_const, choice = "var", axes = 2, top = 10)
contrib <- ggarrange(pc1_contrib, pc2_contrib, 
                     labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("scripts/pca_results/fish/plots/fish_const_pca_contrib.pdf")

img_paths_fish_const <- paste0("images/chaets/", rownames(as.data.frame(fish_color_pca_const$x)))
scores <- as.data.frame(fish_color_pca_const$x)
loadings_fish_const <- as.data.frame(fish_color_pca_const$rotation) 
p_fish_pca_const <- ggplot(as.data.frame(fish_color_pca_const$x), aes(PC1, PC2)) +
  geom_image(aes(image = img_paths_fish_const), size = .05, alpha = .5) +
  xlab("PC1 (31.0%)") +
  ylab("PC2 (26.9%)") +
  ggtitle("constant-k") +
  theme_minimal() +
  geom_segment(data = loadings_fish_const, aes(x=0,y=0,xend=PC1,yend=PC2),
               arrow=arrow(length=unit(0.1,"cm")), color = "#ff0000") +
  geom_text(data = loadings_fish_const, aes(x=PC1, y=PC2, label=rownames(loadings_fish_const)),color="#ff0000", vjust = -0.25)
ggsave("scripts/pca_results/fish/plots/fish_const_pca_pc1pc2.pdf")





## BIRDS ##
pdf("scripts/pca_results/birds/plots/birds_charisma_pca.pdf")
fviz_pca_var(birds_color_pca_charisma,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

pdf("scripts/pca_results/birds/plots/birds_charisma_pca_screeplot.pdf")
fviz_screeplot(birds_color_pca_charisma, addlabels = TRUE, ylim = c(0, 50))
dev.off()

# Contributions of variables to PC1
pc1_contrib <- fviz_contrib(birds_color_pca_charisma, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
pc2_contrib <- fviz_contrib(birds_color_pca_charisma, choice = "var", axes = 2, top = 10)
contrib <- ggarrange(pc1_contrib, pc2_contrib, 
                     labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("scripts/pca_results/birds/plots/birds_charisma_pca_contrib.pdf")

img_paths_birds_charisma <- paste0("images/tangara/", rownames(as.data.frame(birds_color_pca_charisma$x)))
scores <- as.data.frame(birds_color_pca_charisma$x)
loadings_birds_charisma <- as.data.frame(birds_color_pca_charisma$rotation) 
p_birds_pca_charisma <- ggplot(as.data.frame(birds_color_pca_charisma$x), aes(PC1, PC2)) +
  geom_image(aes(image = img_paths_birds_charisma), size = .1, alpha = .5) +
  xlab("PC1 (68.3%)") +
  ylab("PC2 (19.7%)") +
  ggtitle("charisma-k") +
  theme_minimal() +
  geom_segment(data = loadings_birds_charisma, aes(x=0,y=0,xend=PC1,yend=PC2),
               arrow=arrow(length=unit(0.1,"cm")), color = "#ff0000") +
  geom_text(data = loadings_birds_charisma, aes(x=PC1, y=PC2, label=rownames(loadings_birds_charisma)),color="#ff0000", vjust = -0.25)
ggsave("scripts/pca_results/birds/plots/birds_charisma_pca_pc1pc2.pdf")







pdf("scripts/pca_results/birds/plots/birds_human_pca.pdf")
fviz_pca_var(birds_color_pca_human,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

pdf("scripts/pca_results/birds/plots/birds_human_pca_screeplot.pdf")
fviz_screeplot(birds_color_pca_human, addlabels = TRUE, ylim = c(0, 50))
dev.off()

# Contributions of variables to PC1
pc1_contrib <- fviz_contrib(birds_color_pca_human, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
pc2_contrib <- fviz_contrib(birds_color_pca_human, choice = "var", axes = 2, top = 10)
contrib <- ggarrange(pc1_contrib, pc2_contrib, 
                     labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("scripts/pca_results/birds/plots/birds_human_pca_contrib.pdf")

img_paths_birds_human <- paste0("images/tangara/", rownames(as.data.frame(birds_color_pca_human$x)))
scores <- as.data.frame(birds_color_pca_human$x)
loadings_birds_human <- as.data.frame(birds_color_pca_human$rotation) 
p_birds_pca_human <- ggplot(as.data.frame(birds_color_pca_human$x), aes(PC1, PC2)) +
  geom_image(aes(image = img_paths_birds_human), size = .1, alpha = .5) +
  xlab("PC1 (59.0%)") +
  ylab("PC2 (23.9%)") +
  ggtitle("human-k") +
  theme_minimal() +
  geom_segment(data = loadings_birds_human, aes(x=0,y=0,xend=PC1,yend=PC2),
               arrow=arrow(length=unit(0.1,"cm")), color = "#ff0000") +
  geom_text(data = loadings_birds_human, aes(x=PC1, y=PC2, label=rownames(loadings_birds_human)),color="#ff0000", vjust = -0.25)
ggsave("scripts/pca_results/birds/plots/birds_human_pca_pc1pc2.pdf")








pdf("scripts/pca_results/birds/plots/birds_const_pca.pdf")
fviz_pca_var(birds_color_pca_const,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
dev.off()

pdf("scripts/pca_results/birds/plots/birds_const_pca_screeplot.pdf")
fviz_screeplot(birds_color_pca_const, addlabels = TRUE, ylim = c(0, 50))
dev.off()

# Contributions of variables to PC1
pc1_contrib <- fviz_contrib(birds_color_pca_const, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
pc2_contrib <- fviz_contrib(birds_color_pca_const, choice = "var", axes = 2, top = 10)
contrib <- ggarrange(pc1_contrib, pc2_contrib, 
                     labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("scripts/pca_results/birds/plots/birds_const_pca_contrib.pdf")

img_paths_birds_const <- paste0("images/tangara/", rownames(as.data.frame(birds_color_pca_const$x)))
scores <- as.data.frame(birds_color_pca_const$x)
loadings_birds_const <- as.data.frame(birds_color_pca_const$rotation) 
p_birds_pca_const <- ggplot(as.data.frame(birds_color_pca_const$x), aes(PC1, PC2)) +
  geom_image(aes(image = img_paths_birds_const), size = .1, alpha = .5) +
  xlab("PC1 (44.2%)") +
  ylab("PC2 (29.6%)") +
  ggtitle("constant-k") +
  theme_minimal() +
  geom_segment(data = loadings_birds_const, aes(x=0,y=0,xend=PC1,yend=PC2),
               arrow=arrow(length=unit(0.1,"cm")), color = "#ff0000") +
  geom_text(data = loadings_birds_const, aes(x=PC1, y=PC2, label=rownames(loadings_birds_const)),color="#ff0000", vjust = -0.25)
ggsave("scripts/pca_results/birds/plots/birds_const_pca_pc1pc2.pdf")


combo_pca_plots <- ggarrange(p_fish_pca_charisma, p_fish_pca_human, p_fish_pca_const, 
                             p_birds_pca_charisma, p_birds_pca_human, p_birds_pca_const,
                              labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)
combo_pca_plots
ggsave("figures/all_pca_plots_combined.pdf", width = 11, height = 8.5)
ggsave("figures/all_pca_plots_combined.jpeg", width = 11.5, height = 7.5)


## color pattern distributions
all_fish_classifications <- rbind(fish_classifications_human, fish_classifications_charisma, fish_classifications_const)
all_fish_classifications

hist(all_fish_classifications$m)
hist(all_fish_classifications$A)
hist(all_fish_classifications$Jc)
hist(all_fish_classifications$Jt)
hist(all_fish_classifications$m_dS)
hist(all_fish_classifications$m_dL)

library(tidyverse)
library(ggpubr)
p1 <-  ggplot(all_fish_classifications, aes(x=m)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_fish_classifications$m), sd = sd(all_fish_classifications$m))) +
  theme_classic()
p2 <-  ggplot(all_fish_classifications, aes(x=A)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_fish_classifications$A), sd = sd(all_fish_classifications$A))) +
  theme_classic()
p3 <-  ggplot(all_fish_classifications, aes(x=Jc)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_fish_classifications$Jc), sd = sd(all_fish_classifications$Jc))) +
  theme_classic()
p4 <-  ggplot(all_fish_classifications, aes(x=Jt)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_fish_classifications$Jt), sd = sd(all_fish_classifications$Jt))) +
  theme_classic()
p5 <-  ggplot(all_fish_classifications, aes(x=m_dS)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_fish_classifications$m_dS), sd = sd(all_fish_classifications$m_dS))) +
  theme_classic()
p6 <-  ggplot(all_fish_classifications, aes(x=m_dL)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_fish_classifications$m_dL), sd = sd(all_fish_classifications$m_dL))) +
  theme_classic()

combo_hist_fish <- ggarrange(p1, p2, p3, p4, p5, p6,
                        labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)
combo_hist_fish
ggsave("figures/fish_density_hist_combined.pdf", width = 10, height = 5)

all_birds_classifications <- rbind(birds_classifications_human, birds_classifications_charisma, birds_classifications_const)
all_birds_classifications

hist(all_birds_classifications$m)
hist(all_birds_classifications$A)
hist(all_birds_classifications$Jc)
hist(all_birds_classifications$Jt)
hist(all_birds_classifications$m_dS)
hist(all_birds_classifications$m_dL)

p1 <-  ggplot(all_birds_classifications, aes(x=m)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_birds_classifications$m), sd = sd(all_birds_classifications$m))) +
  theme_classic()
p2 <-  ggplot(all_birds_classifications, aes(x=A)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_birds_classifications$A), sd = sd(all_birds_classifications$A))) +
  theme_classic()
p3 <-  ggplot(all_birds_classifications, aes(x=Jc)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_birds_classifications$Jc), sd = sd(all_birds_classifications$Jc))) +
  theme_classic()
p4 <-  ggplot(all_birds_classifications, aes(x=Jt)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_birds_classifications$Jt), sd = sd(all_birds_classifications$Jt))) +
  theme_classic()
p5 <-  ggplot(all_birds_classifications, aes(x=m_dS)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_birds_classifications$m_dS), sd = sd(all_birds_classifications$m_dS))) +
  theme_classic()
p6 <-  ggplot(all_birds_classifications, aes(x=m_dL)) +
  geom_histogram(aes(y = ..density..), colour="black", fill="white") +
  stat_function(fun = dnorm, lwd = 2, col = 'red', 
                args = list(mean = mean(all_birds_classifications$m_dL), sd = sd(all_birds_classifications$m_dL))) +
  theme_classic()

combo_hist_birds <- ggarrange(p1, p2, p3, p4, p5, p6,
                             labels = c("A", "B", "C", "D", "E", "F"), ncol = 3, nrow = 2)
combo_hist_birds
ggsave("figures/birds_density_hist_combined.pdf", width = 10, height = 5)
