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
library(ggplot)
library(dplyr)
library(ape)
library(tidytree)
library(ggtree)


myColorRamp <- function(colors, values) { 
  v <- (values - min(values))/diff(range(values)) 
  x <- colorRamp(colors)(v) 
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255) 
} 

setwd("./data")

Data_PTRMS_outdoor <- read.csv("Data_field_PTRMS.csv",sep=",")


 vec_date <- as.Date(Data_PTRMS_outdoor$Date, format = "%m/%d/%y")
 vec_date_unique <- unique(vec_date)
 


 
 #### loop sec 
### ad chrono in sec 
Sample_ID <- unique(Data_PTRMS_outdoor$Sample_ID)
Data_PTRMS_outdoorx=NULL
 
 for (i in c(1:length(Sample_ID))) {
   
   
   Data_PTRMS_outdoorx_inter <- Data_PTRMS_outdoor[Data_PTRMS_outdoor$Sample_ID == Sample_ID[i],]
   Data_PTRMS_outdoorx_inter$chrono <- c(1:nrow(Data_PTRMS_outdoorx_inter))
   Data_PTRMS_outdoorx_inter <- Data_PTRMS_outdoorx_inter[,c(1:5,ncol(Data_PTRMS_outdoorx_inter),(6):(ncol(Data_PTRMS_outdoorx_inter)-1))]
   
   Data_PTRMS_outdoorx <- rbind(Data_PTRMS_outdoorx,Data_PTRMS_outdoorx_inter)
   
 }
 
 
 
 ##### substract background
 matt_quant_bck = NULL
 
 for (j in c(1:length(vec_date_unique))) {
   
   
   Data_PTRMS_outdoorx1 <- Data_PTRMS_outdoorx[vec_date == vec_date_unique[j],]
   Data_PTRMS_outdoorx1 <- Data_PTRMS_outdoorx1[Data_PTRMS_outdoorx1$chrono <26,]
   
   matt_quant <-   data.frame(aggregate(Data_PTRMS_outdoorx1[,7:34],
                                     by = list(Data_PTRMS_outdoorx1$Treatment,Data_PTRMS_outdoorx1$Plant), # 
                                     FUN = function(x) quantile(x, probs = 0.9))) # 
   
   background <- matt_quant[matt_quant$Group.1 == "background",] 
   
   matt_quant <- matt_quant[!(matt_quant$Group.1 == "background"),] 
   
   matt_quant2 <- matt_quant
   
   for ( i in c(3:ncol(matt_quant2))) { 
     
     matt_quant2[,i] <- matt_quant2[,i] - background[nrow(background),i]
     
   }
   
   date <- rep(j,nrow(matt_quant2))
   
   matt_quant2 <- data.frame(date,matt_quant2)
   
   matt_quant_bck <- rbind(matt_quant_bck,matt_quant2)
   
 }
 
 ## stop here 



 ###
 matt_quant_bck <- matt_quant_bck[!(matt_quant_bck$Group.2 == 29),] # !!!!! to check
 ########## heat map 
 
 
 matt_max_quant_stat <- matt_quant_bck
 matt_max_quant_stat$sample <- paste(matt_max_quant_stat$Group.1,matt_max_quant_stat$Group.2)
 colnames(matt_max_quant_stat)[2] <- "treatment"
 matt_quant_x  <-matt_max_quant_stat[,4:(ncol(matt_max_quant_stat)-1)]
 row.names(matt_quant_x) <- matt_max_quant_stat$sample
 #GCMS_delprim <- GCMS_delprim[,c(1,2,16,21)]
 matt_quant_x[matt_quant_x < 0] <- 0

d_matt_quant_x<- vegdist(matt_quant_x,method="bray") # method="man" # is a bit better
hc_matt_quant_x<- hclust(d_matt_quant_x, method = "complete")
matt_quant_treatment <- rev(levels(as.factor(matt_max_quant_stat$treatment)))

dend <- as.phylo(hc_matt_quant_x)
# order it the closest we can to the order of the observations:


g2 <- split(as.factor(matt_max_quant_stat$sample), as.factor(matt_max_quant_stat$treatment))
tree_plot2 <- groupOTU(dend, g2)

cols <- c("khaki","darkred")


circ <- ggtree(tree_plot2, aes(fill=group),size=1)+ #,layout='circular'
  geom_tiplab(size=2, offset=0.08)  +
  scale_color_manual(values=cols) + theme(legend.position="none")

df <- data.frame(matt_max_quant_stat$treatment)

rownames(df) <- tree_plot2$tip.label


ggheat_plot<-gheatmap(circ, df[, "matt_max_quant_stat.treatment", drop=F], offset=0, width=0.1,colnames = FALSE,
                      colnames_angle=90, colnames_offset_y = 0) + scale_fill_manual(values=cols)

ggheat_plot



#### nmds


cols2 <- cols#sample(hcl.colors(10, "Vik"))


nmds <- metaMDS(matt_quant_x,distance="gower")


matt_plot_nmds<- vegan::scores(nmds) %>%
  cbind(matt_max_quant_stat) 
matt_plot_nmds <- matt_plot_nmds[-40,]


pcoa <- cmdscale(d_matt_quant_x,k=2)
colnames(pcoa) <- c("NMDS1","NMDS2")


  matt_plot_nmds<- pcoa %>%
  cbind(matt_max_quant_stat) 
matt_plot_nmds <- matt_plot_nmds[-40,]





  
  nmds_plot<- ggplot(data=matt_plot_nmds,aes(x = NMDS1, y = NMDS2)) +
  geom_encircle(data = subset(matt_plot_nmds, treatment == "healthy"),color = cols2[2], fill = cols2[2], alpha = 0.5,s_shape=0.8) +
  geom_encircle(data = subset(matt_plot_nmds, treatment == "induced"),color = cols2[1], fill = cols2[1], alpha = 0.5,s_shape=0.8) +
  geom_point(data = subset(matt_plot_nmds, treatment == "healthy"),color = cols2[2], fill = cols2[2]) +
    geom_point(data = subset(matt_plot_nmds, treatment == "induced"),color = cols2[1], fill = cols2[1]) +
    xlim(-0.6,0.7) + ylim(-0.5,0.5) +
 # annotate("text", x = -1, y = 0.5, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw() + scale_fill_manual(values=cols2) + scale_color_manual(values=cols2)  




  
  
  ### data supervizer
  
  matt_max_quant_stat_rf<-matt_max_quant_stat[,4:(ncol(matt_max_quant_stat)-1)]
  #GCMS_delprim_rf <- GCMS_delprim_rf[,c(1,10,8)]
  matt_max_quant_stat_rf$treatment <- as.factor(matt_max_quant_stat$treatment)
  
  
  
  ##################### plsda 
  
  liney <- MASS::lda(treatment~., data = matt_max_quant_stat_rf)
  result <- predict(liney, matt_max_quant_stat_rf)
  lhat <- 1 - sum(result$class == matt_max_quant_stat_rf$treatment)/length(matt_max_quant_stat_rf$treatment)
  
  data <- data.frame(x1=result$x[,1], y=matt_max_quant_stat_rf$treatment)
  data$y <- factor(data$y)
  LDA <- ggplot(data, aes(x=x1, fill=y)) +
    geom_density(adjust=5, alpha=0.6) +
    xlab("x1") +
    ylab("Density") +
    ggtitle(sprintf("PLS-LDA, L = %.2f", lhat)) + scale_fill_manual(values=cols2)+
    theme_bw()
  
  
  
  
##################rf 
  
  
  model_1 = randomForest(treatment~., data = matt_max_quant_stat_rf, importance = TRUE,ntree=10000)
  
  ozone.rp <- rfPermute(treatment ~ ., data = matt_max_quant_stat_rf, na.action = na.omit, ntree = 1000, num.rep = 150)
  var_imp <-  plotImportance(ozone.rp, scale = TRUE,size = 3)
  
  

  
  #importance =  randomForest::importance(model_1)
  #varImportance = data.frame(Variables = row.names(importance),
   #                          Importance =round(importance[, "MeanDecreaseAccuracy"],2))
  
  #rankImportance= varImportance %>% mutate(Rank=paste("#",dense_rank(dplyr::desc(Importance))))
  
  #var_imp<- ggplot(rankImportance,aes(x=reorder(Variables,Importance),y=Importance,fill=Importance))+ 
  #  geom_bar(stat="identity") + 
  #  geom_text(aes(x = Variables, y = 0.5, label = Rank),hjust=0, vjust=0.55, size = 4, colour = "white") +
  #  labs(x = "Variables") +
   # coord_flip() + 
  #  theme_classic()
  
  
  
###################################################################################
  ################################################################################
  
  
  
  
  
  setwd("G:/My Drive/taf/postdoc neuchatel/ted_turling/ERC/PTR-MS field data")
  sink("table.txt")
  print(summary(ozone.rp))
  sink() 
  
  
  
  
  Rr_perm <- readLines("G:/My Drive/taf/postdoc neuchatel/ted_turling/ERC/PTR-MS field data/table.txt")
  Rr_perm <- Rr_perm[-length(Rr_perm)]
  
  
  pdf("PTRMS_field_result_new_pcoa.pdf")
  title <- "PTRMS field result"
  grid::grid.text(title,x = (0.5), y = (0.6))
  
  ggheat_plot
  nmds_plot
  LDA
  plotImportance(ozone.rp, scale = TRUE,size = 3)
  par(mar = c(0.1,1,1,0.1))
  plot(x=c(1,21),y=c(1,29.7),type="n", axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
  liner <- rev(seq(1,29,0.4)) ## interline
  for (i in c(1:length(Rr_perm))) {
    text(x=0, y=liner[i], labels=Rr_perm[i],cex = 0.5,pos = 4) 
  }
  dev.off()
  
  
  
  
  
  