library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)
library(wesanderson)

setwd("./data")

#Variety="Delprim"
Variety="Aventicum"
Variety

if (Variety == "Delprim") {
  GCMS <- read.csv("GCMS_all_data_delprim.csv", sep = ",") # Delprim data
} else if (Variety == "Aventicum") {
  GCMS <- read.csv("GCMS_all_data_aventicum.csv", sep = ",") # Aventicum data
}

GCMS_x <- GCMS[, 3:ncol(GCMS)]
row.names(GCMS_x) <- GCMS$sample

cols <- wes_palette("Cavalcanti1")
cols <- cols[c(4, 1, 5, 2)]
cols[4] <- "royalblue4"
cols[1] <- "darkgreen"
cols2 <- cols 


nmds <- metaMDS(GCMS_x, distance = "gower", autotransform = FALSE)

nmds_plot <- vegan::scores(nmds)[1] %>% # pcoa %>%# 
  cbind(GCMS) %>%
  ggplot(aes(x = sites.NMDS1, y = sites.NMDS2)) +
  xlab(label = "NMDS1") +
  ylab(label = "NMDS2") +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.8, s_shape = 0.8) +
  geom_point(aes(color = treatment)) +
  annotate("text", x = -0.45, y = 0.4, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  xlim(-0.5, 0.6) + ylim(-0.4, 0.4) +
  theme_bw() + scale_fill_manual(values = cols2) + scale_color_manual(values = cols2)

nmds_plot
ggsave(paste0("NMDS_GC", Variety, ".pdf"), dpi = 600)

# PERMANOVA 
sink(paste0("permanova_GC_",Variety,".txt"))  # Redirect output to a file
adonis2(GCMS_x ~ GCMS$treatment, permutations = 999, method = "gower")
sink()  


#Random forest 
GCMS_rf <- GCMS[, 2:ncol(GCMS)]
GCMS_rf$treatment <- as.factor(GCMS_rf$treatment)
sink(paste0("rf_GC_",Variety,".txt"))  # Redirect output to a file
rfPermute(treatment ~ ., data = GCMS_rf, na.action = na.omit, ntree = 1000, num.rep = 150)
sink()  

