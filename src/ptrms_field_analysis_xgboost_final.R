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


myColorRamp <- function(colors, values) { 
  v <- (values - min(values))/diff(range(values)) 
  x <- colorRamp(colors)(v) 
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255) 
} 

setwd("G:/My Drive/taf/postdoc neuchatel/ted_turling/ERC/PTR-MS field data")

Data_ted <- read.csv("Data_25_sec_field_corrected.csv",sep=",")



### SELCTED ORGAN 
#Data_ted_bck <- Data_ted[grep("#N/A",Data_ted$Plant_part) ,]
#Data_ted_org<- Data_ted[grep("systemic",Data_ted$Plant_part),]

#Data_ted <- rbind(Data_ted_bck,Data_ted_org)


p <- ggplot(Data_ted, aes(Treatment,C3H5O.))
p + geom_boxplot()

#### xgboost calssification 
 
#vec_date_unique <- unique(Data_ted$Tech_replicate)
Data_ted2 <-  Data_ted[-grep("leaf_wounded",Data_ted$Plant_part),]
Data_ted2 <-  Data_ted2[-grep("whorle",Data_ted2$Plant_part),]
#Data_ted2 <-  Data_ted2[-grep("T2S2",Data_ted2$Sample_ID),]
#Data_ted2 <-  Data_ted2[-grep("T3S2",Data_ted2$Sample_ID),]
#Data_ted2 <-  Data_ted2[-grep("T2S3",Data_ted2$Sample_ID),]
#Data_ted2 <-  Data_ted2[-grep("T4S1",Data_ted2$Sample_ID),]
#Data_ted2 <-  Data_ted2[-grep("T28S2",Data_ted2$Sample_ID),]

vec_comb <- c("all","1","2","3","1_2","1_3","2_3")
output_conf_mat <- list()

for ( xx in c(1:length(vec_comb))) {

slelector <- vec_comb[xx]

Data_ted2 <-  Data_ted[-grep("leaf_wounded",Data_ted$Plant_part),]
Data_ted2 <-  Data_ted2[-grep("whorle",Data_ted2$Plant_part),]

if(slelector == "all") {Data_ted2 <- Data_ted2}
if(slelector == "1") {Data_ted2 <-  Data_ted2[grep(1,Data_ted2$Repetition_measurment),]}
if(slelector == "2") {Data_ted2 <-  Data_ted2[grep(2,Data_ted2$Repetition_measurment),]}
if(slelector == "3") {Data_ted2 <-  Data_ted2[grep(3,Data_ted2$Repetition_measurment),]}
if(slelector == "1_2") {Data_ted2 <-  Data_ted2[-grep(3,Data_ted2$Repetition_measurment),]}
if(slelector == "1_3") {Data_ted2 <-  Data_ted2[-grep(2,Data_ted2$Repetition_measurment),]}
if(slelector == "2_3") {Data_ted2 <-  Data_ted2[-grep(1,Data_ted2$Repetition_measurment),]}
 
 #Data_ted2 <-  Data_ted2[-grep("leaf_systemic",Data_ted2$),]
 
 #Data_ted2 <-  Data_ted2[-grep("leaf_systemic",Data_ted2$Plant_part),]
 #Data_ted2 <-  Data_ted2[grep(1,Data_ted2$Repetition_measurment),]
 #Data_ted2 <-  Data_ted


 vec_date <- as.Date(Data_ted2$Date, format = "%m/%d/%y")
 vec_date_unique <- unique(vec_date)
 

 p <- ggplot(Data_ted2, aes(Treatment,C3H5O.))
 p + geom_boxplot()

 
#### loop sec 
 
 Sample_ID <- unique(Data_ted2$Sample_ID)
 Data_ted2x=NULL
 
 for (i in c(1:length(Sample_ID))) {
   
   
   Data_tedx <- Data_ted2[Data_ted2$Sample_ID == Sample_ID[i],]
   Data_tedx$chrono <- c(1:nrow(Data_tedx))
   Data_tedx <- Data_tedx[,c(1:7,ncol(Data_tedx),(8):(ncol(Data_tedx)-1))]
   
   Data_ted2x <- rbind(Data_ted2x,Data_tedx)
 
 }
 
 
 
 matt_max_bck_full = NULL
 
 for (j in c(1:length(vec_date_unique))) {
 
 
Data_ted1 <- Data_ted2x[vec_date == vec_date_unique[j],]
Data_ted1 <- Data_ted1[Data_ted1$chrono <26,]

matt_max_sp<-   data.frame(aggregate(Data_ted1[,c(6,9:36)],
                             by = list(Data_ted1$Treatment,Data_ted1$Plant), # ,Data_ted1$Sampling_order_within_day
                             FUN = function(x) quantile(x, probs = 0.9))) # function(x) quantile(x, probs = 0.95)

matt_max_bckgr<-   data.frame(aggregate(Data_ted1[,c(6,9:36)],
                                     by = list(Data_ted1$Treatment,Data_ted1$Plant), # ,Data_ted1$Sampling_order_within_day
                                     FUN = function(x) quantile(x, probs = 0.9))) # function(x) quantile(x, probs = 0.95)


matt_max_bck <- matt_max_bckgr[(matt_max_bckgr$Group.1 == "background"),] 
matt_max_inter <- matt_max_sp[!(matt_max_sp$Group.1 == "background"),] 

  
  for ( i in c(4:ncol(matt_max_inter))) { 
    
    matt_max_inter[,i] <- matt_max_inter[,i] - matt_max_bck[1,i]
     #colMeans(apply(background_inter,2,as.numeric))[i] ## apply(apply(background_inter,2,as.numeric),2,max)[i]
  }
  


date <- rep(j,nrow(matt_max_inter))

matt_max2 <- data.frame(date,matt_max_inter)




matt_max_bck_full <- rbind(matt_max_bck_full,matt_max2)

 


}



 
 
 
matt_max_bck_save <- matt_max_bck_full
if(slelector == "all") {matt_max_bck_save2 <- matt_max_bck_save}


p <- ggplot(matt_max_bck_save, aes(Group.1,matt_max_bck_save$C8H8N.))
p + geom_boxplot()

#######################################################################################
  
  
  
  
  matt_max_bck <- matt_max_bck_save 
 #matt_max_bck<-  matt_max_bck[grep(3,matt_max_bck$Repetition_measurment),]
  
  sample <- paste(matt_max_bck$Group.1,matt_max_bck$Group.2)
  
  matt_max_bck<-   aggregate(matt_max_bck[,c(5:32)],
                             by = list(sample), # 
                             FUN = mean) # function(x) quantile(x, probs = 0.95)
  
  
  ########## heat map 
  
  matt_max_bck_stat <- matt_max_bck
  matt_max_bck_stat$treatment <-  rep(1:nrow(matt_max_bck_stat))
  colnames(matt_max_bck_stat)[1] <- "sample"
  matt_max_bck_stat$treatment[grep("healthy",matt_max_bck_stat$sample)] <- "healthy"
  matt_max_bck_stat$treatment[grep("induced",matt_max_bck_stat$sample)] <- "induced"

  #GCMS_delprim <- GCMS_delprim[,c(1,2,16,21)]
  
  
#matt_max_bck[matt_max_bck<0] <- 0
matt_max_class <- matt_max_bck ### 
#???matt_max_class <- matt_max_class[-grep("background",matt_max_class$Group.1),]

train_data_full <- matt_max_class
train_data_full$Group.1 <- matt_max_bck_stat$treatment


train_data <- train_data_full[,2:ncol(train_data_full)]
train_label <- train_data_full$Group.1
train_label[train_label=="healthy"] <- 2
train_label[train_label=="induced"] <- 1
train <- data.frame(train_label,train_data)   
train_label2 <- train_label
train_label2[train_label2=="2"] <- "0"

bstSparse <- xgboost(data = as.matrix(train_data), 
                     label = train_label2, 
                     max.depth = 5, eta = 0.05, 
                     nthread = 10, nrounds = 500, 
                     objective = "binary:logistic",
                     print_every_n = 100)

importance_matrix <- xgb.importance(model = bstSparse)

if(slelector == "all") {plot_imp <- xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")}
if(slelector == "all") {vec_imp_save <- c("date","Group.1","Group.2",importance_matrix$Feature[1:10])}

## auto selection
vec_imp <- c("date","Group.1","Group.2",importance_matrix$Feature[1:10])

matt_max_class_filter <- matt_max_class[,colnames(matt_max_class) %in% vec_imp]
########### gbm  classification 
library(xgboost)



matt_max_class_filter$Group.1 <- matt_max_bck_stat$treatment


valid_tree <- NULL
valid_xgboost <- NULL

for (i in c (1:nrow(matt_max_class_filter))) {
  
  
  
  train_data_full <- matt_max_class_filter[-i,] 
  train_data <- train_data_full[,2:ncol(train_data_full)]
  train_label <- train_data_full$Group.1
  train_label[train_label=="healthy"] <- 2
  train_label[train_label=="induced"] <- 1
  train <- data.frame(train_label,train_data)
  
   tree=rpart(train_label~.,data=train,method="class")
   train_label2 <- train_label
   train_label2[train_label2=="2"] <- "0"
   bstSparse <- xgboost(data = as.matrix(train_data), 
                        label = train_label2, 
                        max.depth = 5, eta = 0.05, 
                        nthread = 10, nrounds = 500, 
                        objective = "binary:logistic",
                        print_every_n = 100)
  
  test_data_full <- matt_max_class_filter[i,]  
  test_data <- test_data_full[,2:ncol(test_data_full)]
  test_label <- test_data_full$Group.1
  test_label[test_label=="healthy"] <- 2
  test_label[test_label=="induced"] <- 1

  xgb.pred = predict(tree,test_data,type="class")
  
  pred <- predict(bstSparse, as.matrix(test_data))
  pred[pred>=0.5] <- 1
  pred[pred<0.5] <- 2
  
  valid_tree <- c(valid_tree,xgb.pred)
  valid_xgboost <- c(valid_xgboost,pred)
  
}


real_data <-  matt_max_class_filter$Group.1
real_data[real_data=="healthy"] <- 2
real_data[real_data=="induced"] <- 1


table(real_data,valid_xgboost)
confi <-caret::confusionMatrix(as.factor(real_data),as.factor(valid_xgboost))


output_conf_mat[[xx]] <- confi

}

data_plot_conf = NULL

for (s in c(1:length(vec_comb))) {

combi <-as.numeric(output_conf_mat[[s]]$overall[c(1,3,4)])
data_plot_conf <- rbind(data_plot_conf,combi)
}
colnames(data_plot_conf) <- c("Accuracy","AccuracyUpper","AccuracyLower")

vec_group <- c("three_run_25s","one_run_25s","one_run_25s","one_run_25s","two_run_25s","two_run_25s","two_run_25s")
data_plot_conf<- data.frame(vec_group,vec_comb,data_plot_conf)

data_plot_conf_int1 <- data_plot_conf[,c(1,2,3)]
colnames(data_plot_conf_int1)[3] <- "Accuracy"
data_plot_conf_int2 <- data_plot_conf[,c(1,2,4)]
colnames(data_plot_conf_int2)[3] <- "Accuracy"
data_plot_conf_int3 <- data_plot_conf[,c(1,2,5)]
colnames(data_plot_conf_int3)[3] <- "Accuracy"
data_plot_conf_stack <-rbind(data_plot_conf_int1,data_plot_conf_int2,data_plot_conf_int3)
##### different expé$
 
data_plot_conf_stack$vec_group <- factor(data_plot_conf_stack$vec_group  , levels=c("one_run_25s", "two_run_25s", "three_run_25s"))

accura_plot<-  ggplot() + 
  geom_boxplot(data=data_plot_conf_stack, mapping=aes(x=vec_group, ymin=Accuracy,group=vec_group), width=0.3, size=1, color="gray50")  + ylim(.25,1) +
  theme_classic() + scale_fill_viridis() +  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) 

plot_imp_tot <- ggplot(plot_imp,aes(x=reorder(Feature,Importance),y=Importance,fill=Importance))+   labs(x = "compounds") + theme_classic() +
  geom_bar(stat="identity") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

conf_table <- output_conf_mat[[1]]$table
colnames(conf_table) <- c("healty","induced")
rownames(conf_table) <- c("healty","induced")



lay <- rbind(c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2),
             c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2))

grid.arrange(grobs = list(accura_plot, plot_imp_tot),ncol=2, layout_matrix = lay)



pdf(file = "G:/My Drive/taf/postdoc neuchatel/ted_turling/ERC/PTR-MS field data/Machine_learning_output1.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 4) # The height of the plot in inches

grid.arrange(grobs = list(accura_plot, plot_imp_tot),ncol=2, layout_matrix = lay)

 
fourfoldplot(conf_table)


#### plot 3d maiz2
cm_xg <- output_conf_mat[[1]]

x = c(0, 200, 200, 0)
y = c(0, 0, 200, 200)
xy = cbind(x,y)
p = Polygon(xy)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
#plot(sps)


#plot(spsample(sps,n=100,type="hexagonal"), pch=3, cex=0.5)
virtual_field <- as.data.frame(spsample(sps,n=75,type="hexagonal"))
height <- sample(runif(nrow(virtual_field), min=5, max=8))
virtual_field$height <- height
virtual_field <- virtual_field[1:sum(cm_xg$table),]


ID  <- c(rep("true_pos",(cm_xg$table[1])),rep("false_neg",(cm_xg$table[2])),rep("false_pos",(cm_xg$table[3])),rep("true_neg",(cm_xg$table[4])))
#undefined <- rep("undefined",(nrow(virtual_field)-length(ID)))
#ID <- c(undefined,ID)
virtual_field$result <- sample(ID)
virtual_field$result <- sample(virtual_field$result )

cols <- virtual_field$result 
cols[cols == "true_pos"] <- "chartreuse1"
cols[cols == "true_neg"] <- "chartreuse4"
cols[cols == "false_pos"] <- "red"
cols[cols == "false_neg"] <- "red"
#cols[cols == "undefined"] <- "orange"




#maize <- read.ply("maize.ply", ShowSpecimen = TRUE,addNormals = TRUE)

maize=read.obj("C:/Users/admin/Downloads/maize1/maize/10439_Corn_Field_v1_max2010_it2.obj",convert.rgl=TRUE)
maize <- maize[1]$`10439_Corn_Field_v1_SG`



sim1 <- grf(40000, grid = "reg", cov.pars = c(1, .25)) 
z1 <- sim1$data*10+20
z <-matrix(z1,200,200)# Exaggerate the relief

x <- 10 * (1:nrow(z))   # 10 meter spacing (S to N)
y <- 10 * (1:ncol(z))   # 10 meter spacing (E to W)
zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1

#colorlut <- brewer.pal(n = 8, name = "Dark2") #terrain.colors(zlen) # height color lookup table ###??? brewer.pal(n = 8, name = "Dark2")

cols2 <- myColorRamp(c("lemonchiffon1","wheat1","burlywood3","darkgoldenrod4"),c(range(z1)[1]:range(z1)[2]))



as.matrix(z) %>%
  #sphere_shade(texture = create_texture("lemonchiffon1","wheat1","khaki3","lightgoldenrod4","peru")) %>%
  height_shade(texture = cols2) %>%
  add_shadow(ray_shade(z, zscale = 3), 0.8) %>%
  add_shadow(ambient_shade(z), 0.7) %>%
  plot_3d(z, zscale = 5, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800))


for (i in 1:nrow(virtual_field)) {
  
  plot3d(translate3d(rotate3d(scale3d(maize,x=0.02*virtual_field$height[i], y=0.02*virtual_field$height[i], z=0.04*virtual_field$height[i])
                              ,90,x=1, y=1, z=1),x= virtual_field$x[i]-100,z=rev(virtual_field$y[i]-110),y=0),col=cols[i],add=T,smooth=0.1,
                                specular = "black")

}




render_snapshot(clear = TRUE)

dev.off()
#render_depth(focus = 0.5, focallength = 30,clear=TRUE)




#######################################################################################
######################################################################################3
######################################################################################
##### multi field 

matt_max_bck <- matt_max_bck_save2 


#matt_max_bckgrp <- matt_max_bck[grep(1,matt_max_bck$date),]


sample <- paste(matt_max_bck$Group.1,matt_max_bck$Group.2)

## data_prep


matt_max_bck_stat <- matt_max_bck[,-3]
matt_max_bck_stat <- matt_max_bck_stat[,-3]
colnames(matt_max_bck_stat)[2] <- "treatment"


vec_imp <- c("date","treatment",importance_matrix$Feature[1:10])
vec_imp_save[2]  <- "treatment"

matt_max_class_filter <- matt_max_bck_stat[,colnames(matt_max_bck_stat) %in% vec_imp_save]

cm_table_full = NULL

vec_loops <- unique(matt_max_class_filter$date)


#matt_max_class_filter[matt_max_class_filter<0] <- 0

for ( i in c(1:length(vec_loops))) {

grp <- vec_loops[i]

train_data_full <- matt_max_class_filter[-grep(grp,matt_max_class_filter$date),] 
#train_data_full <- matt_max_class[tirage,] 


train_data <- train_data_full[,3:ncol(train_data_full)]
train_data <- train_data # /(mean(apply(train_data,2,as.numeric))*mean(apply(train_data,2,as.numeric)))


train_label <- train_data_full$treatment
train_label[train_label=="healthy"] <- 2
train_label[train_label=="induced"] <- 1
train <- data.frame(train_label,train_data)

tree=rpart(train_label~.,data=train,method="class")


train_label[train_label=="2"] <- "0"
bstSparse <- xgboost(data = as.matrix(train_data), label = train_label, max.depth = 5, eta = 0.05, nthread = 5, nrounds = 2000, objective = "binary:logistic")

test_data_full <- matt_max_class_filter[grep(grp,matt_max_class_filter$date),] 
#test_data_full <- matt_max_class[tirage,] 

test_data <- test_data_full[,3:ncol(test_data_full)]
test_label <- test_data_full$treatment
test_label[test_label=="healthy"] <- 2
test_label[test_label=="induced"] <- 1

pred = predict(tree,test_data,type="class")
xgb.pred <- predict(bstSparse, as.matrix(test_data))

xgb.pred[xgb.pred>=0.5] <- 1
xgb.pred[xgb.pred<0.5] <- 2

cm_table  <- data.frame(test_label,xgb.pred)
cm_table$date <- rep(vec_loops[i],nrow(cm_table))
cm_table_full <- data.frame(rbind(cm_table_full,cm_table))

}

cm <- confusionMatrix(as.factor(cm_table_full$test_label),as.factor(cm_table_full$xgb.pred))
tocsv4 <- data.frame(t(c(cm$overall,cm$byClass)))
rownames(tocsv3) <- c("dual_model_FULL")


tot_table2 <- data.frame(test_label,xgb.pred,pred)
tot_table2$sur <- tot_table2$xgb.pred-as.numeric(tot_table2$pred)
tot_table2$sur[tot_table2$sur>0] = NA
tot_table2$sur[tot_table2$sur<0] = NA

cross_result2 <- tot_table2$xgb.pred[!is.na(tot_table2$sur)]
cross_real2 <- tot_table2$test_label[!is.na(tot_table2$sur)]

table(cross_real2,cross_result2)
confusionMatrix(as.factor(cross_real2),as.factor(cross_result2))


cm1 <- cm_table_full[cm_table_full$date == 1,]
cm1x <- caret::confusionMatrix(as.factor(cm1$test_label),as.factor(cm1$xgb.pred))

cm2 <- cm_table_full[cm_table_full$date == 2,]
cm2x <- caret::confusionMatrix(as.factor(cm2$test_label),as.factor(cm2$xgb.pred))


cm3 <- cm_table_full[cm_table_full$date == 3,]
cm3x <- caret::confusionMatrix(as.factor(cm3$test_label),as.factor(cm3$xgb.pred))

cm1x
cm2x
cm3x




combi1 <-as.numeric(cm1x$overall[c(1,3,4)])
combi2 <-as.numeric(cm2x$overall[c(1,3,4)])
combi3 <-as.numeric(cm3x$overall[c(1,3,4)])
data_plot_conf_field <- rbind(combi1,combi2,combi3)


colnames(data_plot_conf_field) <- c("Accuracy","AccuracyUpper","AccuracyLower")

vec_group <- c("field_one","field_two","field_three")
data_plot_conf_field<- data.frame(vec_group,data_plot_conf_field)

data_plot_conf_int1 <- data_plot_conf_field[,c(1,2)]
colnames(data_plot_conf_int1)[2] <- "Accuracy"
data_plot_conf_int2 <- data_plot_conf_field[,c(1,3)]
colnames(data_plot_conf_int2)[2] <- "Accuracy"
data_plot_conf_int3 <- data_plot_conf_field[,c(1,4)]
colnames(data_plot_conf_int3)[2] <- "Accuracy"
data_plot_conf_stack <-rbind(data_plot_conf_int1,data_plot_conf_int2,data_plot_conf_int3)
##### different expé$

data_plot_conf_stack$vec_group <- factor(data_plot_conf_stack$vec_group  , levels=c("field_one", "field_two", "field_three"))



accura_plot_field<-  ggplot() + 
  geom_boxplot(data=data_plot_conf_stack, mapping=aes(x=vec_group, ymin=Accuracy,group=vec_group), width=0.3, size=1, color="gray50")  + ylim(0,1) +
  theme_classic() + scale_fill_viridis() +  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) 

table_compil <- cm1x$table + cm2x$table + cm3x$table 
colnames(table_compil) <- c("healty","induced")
rownames(table_compil) <- c("healty","induced")

pdf(file = "G:/My Drive/taf/postdoc neuchatel/ted_turling/ERC/PTR-MS field data/Machine_learning_output2.pdf",   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

accura_plot_field

fourfoldplot(table_compil)

dev.off()




