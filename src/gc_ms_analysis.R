
library(xgboost)
library(caret)
library(rpart)
library(sp)
library(rayshader)
library(rgl)
library(readobj)
library(vegan)
library(rfPermute)
library(randomForest)
library(ggalt)
library(ggplot2)
library(dplyr)
library(ggtree)
library(ape)
library(vegan)
library(dendextend)
library(wesanderson)
######################################################################
#######################GC_ms_delprim
setwd("./data")


GCMS_delprim <- read.csv("GCMS_all_data_delprim.csv",sep=",") ## ad path var


GCMS_delprim_x  <-GCMS_delprim[,3:ncol(GCMS_delprim)]
row.names(GCMS_delprim_x) <- GCMS_delprim$sample





cols <- wes_palette("Cavalcanti1")
cols <-cols[c(4,1,5,2)]
cols[4] <- "royalblue4"
cols[1] <- "darkgreen"

####3 meta mds

cols2 <- cols#sample(hcl.colors(10, "Vik"))


nmds <- metaMDS(GCMS_delprim_x,distance="gower", autotransform = FALSE)

nmds_plot<-vegan::scores(nmds)[1] %>% # pcoa %>%# 
  cbind(GCMS_delprim) %>%
  ggplot(aes(x = sites.NMDS1, y = sites.NMDS2)) +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.8,s_shape=0.8) +
  geom_point(aes(color = treatment)) +
  annotate("text", x = -0.45, y = 0.4, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  xlim(-0.5,0.6) + ylim(-0.4,0.4) +
  theme_bw() + scale_fill_manual(values=cols2) + scale_color_manual(values=cols2)

nmds_plot
### data supervizer

GCMS_delprim_rf<-GCMS_delprim[,2:ncol(GCMS_delprim)]
#GCMS_delprim_rf <- GCMS_delprim_rf[,c(1,10,8)]
GCMS_delprim_rf$treatment <- as.factor(GCMS_delprim_rf$treatment)



#### randomforest 

model_1 = randomForest(treatment~., data = GCMS_delprim_rf, importance = TRUE,ntree=10000)

ozone.rp <- rfPermute(treatment ~ ., data = GCMS_delprim_rf, na.action = na.omit, ntree = 500, num.rep = 100)
print(ozone.rp)
