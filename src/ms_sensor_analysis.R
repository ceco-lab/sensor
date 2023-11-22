
library(ggtree)
library(vegan)
library(ggplot2)
library(wesanderson)
library(ggalt)
library(caret)
library(randomForest)
library(rfPermute)
library(ape)

######################################################################
#######################MSS sensors data

setwd("./data")


module1 <- read.csv("module1_Adventicum.csv",sep=",")

module1x  <-module1[,3:ncol(module1)]
row.names(module1x) <- module1$PrimaryID
module1x  <-module1x[,seq(1,ncol(module1x),10)]


cols <- wes_palette("Cavalcanti1")
cols <-cols[c(4,1,5,2)]
cols[4] <- "royalblue4"
cols[1] <- "darkgreen"
cols <- cols[c(2,1,3,4)]


####3 meta mds

cols2 <- cols#sample(hcl.colors(10, "Vik"))


nmds <- metaMDS(module1x,distance="gower")

nmds_plot<- vegan::scores(nmds) %>%
  cbind(module1) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_encircle(aes(group = ClassID, color = ClassID, fill = ClassID), alpha = 0.8,s_shape=0.8) +
  geom_point(aes(color = ClassID)) +
  annotate("text", x = -1, y = 0.5, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw() + scale_fill_manual(values=cols2) + scale_color_manual(values=cols2)

nmds_plot

### data supervizer

module1rf  <-module1x[,seq(1,ncol(module1x),10)]

module1rf_rf<-module1rf[,2:ncol(module1rf)]
module1rf_rf$treatment <- as.factor(substr(row.names(module1rf_rf), start = 1, stop = 2))



#### randomforest 

model_1 = randomForest(treatment~., data = module1rf_rf, importance = TRUE,ntree=10000)

ozone.rp <- rfPermute(treatment ~ ., data = module1rf_rf, na.action = na.omit, ntree = 100, num.rep = 50)
ozone.rp

