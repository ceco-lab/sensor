library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)


setwd("./data")

set.seed(6)

#Select plant variety and load corresponding data 

Variety="Delprim"
#Variety="Aventicum"

if (Variety == "Delprim") {
  GCMS <- read.csv("GCMS_Delprim.csv", sep = ",") # Delprim data
} else if (Variety == "Aventicum") {
  GCMS <- read.csv("GCMS_Aventicum.csv", sep = ",") # Aventicum data
}

#Keep only features
GCMS_x <- GCMS[, 2:ncol(GCMS)]

#NMDS
cols <-c("#D8B70A","darkgreen","#972D15","royalblue4")

nmds <- metaMDS(GCMS_x, distance = "gower", autotransform = FALSE)
nmds_plot <- vegan::scores(nmds)[1] %>%  
  cbind(GCMS) %>%
  ggplot(aes(x = sites.NMDS1, y = sites.NMDS2)) +
  xlab(label = "NMDS1") +
  ylab(label = "NMDS2") +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.6, s_shape = 0.8,expand=0.03) +
  geom_point(aes(color = treatment)) +
  annotate("text", x = -0.35, y = 0.3, label = paste0("stress: ", round(nmds$stress, digits = 4)), hjust = 0) +
  xlim(-0.4, 0.5) + ylim(-0.3, 0.3) +
  theme_bw() + scale_fill_manual(values = cols) + scale_color_manual(values = cols)+
  theme(legend.position="none")

nmds_plot
ggsave(paste0("NMDS_GC_", Variety, ".pdf"), width=8, height=8,dpi = 600)

# PERMANOVA 
sink(paste0("permanova_GC_",Variety,".txt"))  
adonis2(GCMS_x ~ GCMS$treatment, permutations = 999, method = "gower")
sink()  

#Random forest 
GCMS_x$treatment <- as.factor(GCMS$treatment)
sink(paste0("rf_GC_",Variety,".txt"))  # Redirect output to a file
rfPermute(treatment ~ ., data = GCMS_x, na.action = na.omit, ntree = 1000, num.rep = 150)
sink()  

