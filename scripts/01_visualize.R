## setup
wd <- "~/Dropbox/School/UCLA MS Biology/Coursework/EEB C202/eeb_c202_final_project"
setwd(wd)

library(tidyverse)
library(ggpubr)

## load in data
fish <- read.csv("datasets/chaets_all_cleaned.csv", header = T, sep = ",")
birds <- read.csv("datasets/tanager_all_cleaned.csv", header = T, sep = ",")

glimpse(fish)
glimpse(birds)

## visualize
### https://stackoverflow.com/questions/6957549/overlaying-histograms-with-ggplot2-in-r
plot_multi_histogram <- function(df, feature, label_column, x_lab, color_pal) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="dodge", aes(y = ..density..), color="black", binwidth = .5) +
    #geom_density(alpha=0.7) +
    #geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    #labs(x=feature, y = "Frequency") +
    labs(x=x_lab, y = "Frequency") +
    scale_color_manual(values=color_pal) +
    scale_fill_manual(values=color_pal) +
    scale_x_continuous(breaks = seq(1, 8, by = 1), limits = c(0,9)) +
    scale_y_continuous(breaks = seq(0, 0.9, by = 0.1), limits = c(0,0.9)) +
    theme_classic()
  plt + guides(fill=guide_legend(title=label_column))
}

## make plots and save
fish_hist <- plot_multi_histogram(fish, 'k', 'method', "Color Classes (k)", c("#999999", "#E69F00"))
fish_hist
ggsave("figures/fish_k_histogram.pdf", fish_hist)

birds_hist <- plot_multi_histogram(birds, 'k', 'method', "Color Classes (k)", c("#d795ed", "#34785f"))
birds_hist
ggsave("figures/birds_k_histogram.pdf", birds_hist)


combo_hist <- ggarrange(fish_hist, birds_hist, 
                     labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("figures/fish_birds_hist_combined.pdf", width = 10, height = 5)
