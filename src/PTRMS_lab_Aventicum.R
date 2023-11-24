library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)

setwd("./data")


Variety="Delprim"
#Variety="Aventicum"


if (Variety == "Delprim") {
  PTRMS <- read.csv("Delprim_PTRlab_new.csv", sep = ";") # Delprim data
} else if (Variety == "Aventicum") {
  PTRMS <- read.csv("Aventicum_PTRlab_new.csv", sep = ";") # Aventicum data
}


Sample_ID <- unique(PTRMS$Sample_ID)

View(PTRMS)
matt_90 <-   data.frame(aggregate(PTRMS[,4:39],
                                  by = list(PTRMS$Sample_ID,PTRMS$Inducer), 
                                  FUN = function(x) quantile(x, probs = 0.9))) # function(x) quantile(x, probs = 0.95)



matt_90_stat <- matt_90
matt_90_stat$sample <- matt_90_stat$Group.1
colnames(matt_90_stat)[2] <- "treatment"
matt_90_x  <-matt_90_stat[,3:(ncol(matt_90_stat)-1)]
row.names(matt_90_x) <- matt_90_stat$sample


cols2<-c("#D8B70A","darkgreen","#972D15","royalblue4")

nmds <- metaMDS(matt_90_x,distance="gower")

matt_plot_nmds<- vegan::scores(nmds)[1] %>%
  cbind(matt_90_stat) 

nmds_plot<- ggplot(data=matt_plot_nmds,aes(x = sites.NMDS1, y = sites.NMDS2)) +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.8,s_shape=0.7) +
  geom_point(aes(color = treatment)) +
  xlab(label = "NMDS1") +
  ylab(label = "NMDS2") +
  annotate("text", x = -0.65, y = 0.5, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  xlim(-0.7,0.8) + ylim(-0.5,0.5) +
  theme_bw() + scale_fill_manual(values=cols2) + scale_color_manual(values=cols2)  

nmds_plot
ggsave(paste0("NMDS_PTRMS", Variety, ".pdf"), dpi = 600)

# PERMANOVA 
sink(paste0("permanova_GC_",Variety,".txt"))  # Redirect output to a file
adonis2(matt_90_x ~ matt_90_stat$treatment, permutations = 999, method = "gower")
sink()  

str(matt_90_x)
#Random forest 
PTR_rf <- matt_90_x 
PTR_rf$treatment <- as.factor(matt_90_stat$treatment)
sink(paste0("rf_PTR_",Variety,".txt"))  # Redirect output to a file
rfPermute(treatment ~ ., data = PTR_rf, na.action = na.omit, ntree = 500, num.rep = 150)
sink()  


