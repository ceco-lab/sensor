library(xgboost)
library(caret)
library(rpart)
library(sp)
library(rayshader)
library(rgl)
library(readobj)
library(geoR)
library(ape)
library(ggalt)
library(vegan)
library(ggtree)
library(wesanderson)
library(rfPermute)

myColorRamp <- function(colors, values) { 
  v <- (values - min(values))/diff(range(values)) 
  x <- colorRamp(colors)(v) 
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255) 
} 

# Set working directory
setwd("./data")

# Load data
Aventicum_PTRlab <- read.csv("Delprim_PTRlab_new.csv", sep = ";")


Sample_ID <- unique(Aventicum_PTRlab$Sample_ID)
matt_max_bck = NULL


matt_max <-   data.frame(aggregate(Aventicum_PTRlab[,4:39],
                                  by = list(Aventicum_PTRlab$Sample_ID,Aventicum_PTRlab$Inducer), # Data_ted1$Plant
                                  FUN = function(x) quantile(x, probs = 0.9))) # function(x) quantile(x, probs = 0.95)


matt_max_bck <- matt_max

matt_max_bck_stat <- matt_max_bck
matt_max_bck_stat$sample <- matt_max_bck_stat$Group.1
colnames(matt_max_bck_stat)[2] <- "treatment"
matt_max_x  <-matt_max_bck_stat[,3:(ncol(matt_max_bck_stat)-1)]
row.names(matt_max_x) <- matt_max_bck_stat$sample


cols <- wes_palette("Cavalcanti1")
cols <-cols[c(4,1,5,2)]
cols[4] <- "royalblue4"
cols[1] <- "darkgreen"
cols2 <- cols#sample(hcl.colors(10, "Vik"))


nmds <- metaMDS(matt_max_x,distance="gower")


matt_plot_nmds<- vegan::scores(nmds)[1] %>%
  cbind(matt_max_bck_stat) 

nmds_plot<- ggplot(data=matt_plot_nmds,aes(x = sites.NMDS1, y = sites.NMDS2)) +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.8,s_shape=0.7) +
  geom_point(aes(color = treatment)) +
  annotate("text", x = -0.65, y = 0.5, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  xlim(-0.7,0.8) + ylim(-0.5,0.5) +
  theme_bw() + scale_fill_manual(values=cols2) + scale_color_manual(values=cols2)  

nmds_plot


### data supervizer

matt_max_bck_stat_rf<-matt_max_bck_stat[,4:(ncol(matt_max_bck_stat)-1)]
#GCMS_delprim_rf <- GCMS_delprim_rf[,c(1,10,8)]
matt_max_bck_stat_rf$treatment <- as.factor(matt_max_bck_stat$treatment)


##################rf 

model_1 = randomForest(treatment~., data = matt_max_bck_stat_rf, importance = TRUE,ntree=10000)

ozone.rp <- rfPermute(treatment ~ ., data = matt_max_bck_stat_rf, na.action = na.omit, ntree = 1000, num.rep = 200)

print(summary(ozone.rp))



###########################################################################
##########################################################################