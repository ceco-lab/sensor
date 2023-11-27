library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)


set.seed(6)

setwd("./data")

Variety="Delprim"
#Variety="Aventicum"


if (Variety == "Delprim") {
  MSS <- read.csv("MSS_Delprim.csv", sep = ",") # Delprim data
} else if (Variety == "Aventicum") {
  MSS <- read.csv("MSS_Aventicum.csv", sep = ",") # Aventicum data
}

#Keep only features
MSS_x <- MSS[, 3:ncol(MSS)]
row.names(MSS_x) <- MSS$PrimaryID
#MSS_x  <-MSS_x[,seq(1,ncol(MSS_x),10)]


cols <-c("#D8B70A","darkgreen","#972D15","royalblue4")

nmds <- metaMDS(MSS_x,distance="gower")

nmds_plot<- vegan::scores(nmds) %>%
  cbind(MSS) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.8,s_shape=0.8) +
  geom_point(aes(color = treatment)) +
  annotate("text", x = -1, y = 0.5, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
    theme_bw() + scale_fill_manual(values=cols) + scale_color_manual(values=cols)  +
  theme(legend.position="none")

nmds_plot
ggsave(paste0("NMDS_MSS", Variety, ".pdf"), dpi = 600)

# PERMANOVA 
sink(paste0("permanova_MSS_",Variety,".txt"))  # Redirect output to a file
adonis2(MSS_x ~ MSS$treatment, permutations = 999, method = "gower")
sink()  

#Random forest 
MSS_x$treatment <- as.factor(MSS$treatment)
sink(paste0("rf_MSS_",Variety,".txt"))  # Redirect output to a file
rfPermute(treatment ~ ., data = MSS_x, na.action = na.omit, ntree = 1000, num.rep = 150)
sink()  
