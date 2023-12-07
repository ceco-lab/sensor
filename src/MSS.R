# Load necessary libraries
library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)

# Set working directory
setwd("./data")

# Set seed for reproducibility
set.seed(6)

# Specify the variety
Variety <- "Aventicum"  # Change to "Aventicum" or "Delprim" if needed

# Read data based on variety
if (Variety == "Delprim") {
  MSS <- read.csv("MSS_Delprim.csv", sep = ",") # Delprim data
} else if (Variety == "Aventicum") {
  MSS <- read.csv("MSS_Aventicum.csv", sep = ",") # Aventicum data
}

#Keep only features
MSS_x <- MSS[, 3:ncol(MSS)]
str(MSS)

MSS$treatment
cols <-c("#D8B70A","darkgreen","#972D15","royalblue4")

#NMDS
nmds <- metaMDS(MSS_x,distance="gower")

nmds_plot<- vegan::scores(nmds) %>%
  cbind(MSS) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.8,s_shape=0.8) +
  geom_point(aes(color = treatment)) +
  xlim(-0.5, 0.5) +
  ylim(-0.2, 0.3) +
  annotate("text", x = -0.5, y = 0.3, label = paste0("stress: ", format(nmds$stress, digits = 3)), hjust = 0) +
    theme_bw() + scale_fill_manual(values=cols) + scale_color_manual(values=cols)  +
  theme(legend.position="none",aspect.ratio=1)

nmds_plot
ggsave(paste0("NMDS_MSS_", Variety, ".pdf"), dpi = 600)

# PERMANOVA 
sink(paste0("permanova_MSS_",Variety,".txt"))  # Redirect output to a file
adonis2(MSS_x ~ MSS$treatment, permutations = 999, method = "gower")
sink()  

#Random forest 
MSS_x$treatment <- as.factor(MSS$treatment)
sink(paste0("rf_MSS_",Variety,".txt"))  # Redirect output to a file
rfPermute(treatment ~ ., data = MSS_x, na.action = na.omit, ntree = 1000, num.rep = 150)
sink()  
