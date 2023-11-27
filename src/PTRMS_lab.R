library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)

setwd("./data")

set.seed(6)

#Variety="Delprim"
Variety="Aventicum"

if (Variety == "Delprim") {
  PTRMS <- read.csv("PTR_Lab_Delprim.csv", sep = ",") # Delprim data
} else if (Variety == "Aventicum") {
  PTRMS <- read.csv("PTR_Lab_Aventicum.csv", sep = ",") # Aventicum data
}


Sample_ID <- unique(PTRMS$Sample_ID)
matt_90 <-   data.frame(aggregate(PTRMS[,3:38],
                                  by = list(PTRMS$Sample_ID,PTRMS$treatment), 
                                  FUN = function(x) quantile(x, probs = 0.9))) # function(x) quantile(x, probs = 0.95)


matt_90_stat <- matt_90
matt_90_stat$sample <- matt_90_stat$Group.1
colnames(matt_90_stat)[2] <- "treatment"
matt_90_x  <-matt_90_stat[,3:(ncol(matt_90_stat)-1)]
row.names(matt_90_x) <- matt_90_stat$sample


cols <-c("#D8B70A","darkgreen","#972D15","royalblue4")

nmds <- metaMDS(matt_90_x,distance="gower")

nmds_plot <- vegan::scores(nmds)[1] %>%  
  cbind(matt_90_stat) %>%
  ggplot(aes(x = sites.NMDS1, y = sites.NMDS2)) +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.6, s_shape = 0.8,expand=0.03) +
  geom_point(aes(color = treatment)) +
  xlab(label = "NMDS1") +
  ylab(label = "NMDS2") +
  annotate("text", x = -0.65, y = 0.21, label = paste0("stress: ", round(nmds$stress, digits = 4)), hjust = 0) +
  xlim(-0.7,0.8) + ylim(-0.26,0.26) +
  theme_bw() + scale_fill_manual(values=cols) + scale_color_manual(values=cols)  +
  theme(legend.position="none")

nmds_plot
ggsave(paste0("NMDS_PTRMS_", Variety, ".pdf"), dpi = 600)

# PERMANOVA 
sink(paste0("permanova_GC_",Variety,".txt"))  # Redirect output to a file
adonis2(matt_90_x ~ matt_90_stat$treatment, permutations = 999, method = "gower")
sink()  

#Random forest 
matt_90_x$treatment <- as.factor(matt_90_stat$treatment)
sink(paste0("rf_PTR_",Variety,".txt"))  # Redirect output to a file
rfPermute(treatment ~ ., data = matt_90_x, na.action = na.omit, ntree = 1000, num.rep = 150)
sink()  


