# Load libraries
library(vegan)
library(rfPermute)
library(ggplot2)
library(dplyr)
library(ape)
library(ggalt)

# Set working directory
setwd("./data")

# Set seed for reproducibility
set.seed(6)

# Load data

Data_PTRMS_outdoor <- read.csv("Data.csv", sep = ";")
Data_PTRMS_outdoor$Repetition_measurment[Data_PTRMS_outdoor$Treatment=="background"] <- 0

# Convert date to Date format and get unique dates
vec_date <- as.Date(Data_PTRMS_outdoor$Date, format = "%m/%d/%y")
vec_date_unique <- unique(vec_date)

# Get unique sample IDs
Sample_ID <- unique(Data_PTRMS_outdoor$Sample_ID)


# Initialize empty data frame
Data_PTRMS_outdoorx <- NULL

# Build time column for each sample
for (i in c(1:length(Sample_ID))) {
  # Subset data for each sample ID
  Data_PTRMS_outdoorx_inter <- Data_PTRMS_outdoor[Data_PTRMS_outdoor$Sample_ID == Sample_ID[i], ]
  # Create a chronological column
  Data_PTRMS_outdoorx_inter$chrono <- c(1:nrow(Data_PTRMS_outdoorx_inter))
  # Reorder columns
  Data_PTRMS_outdoorx_inter <- Data_PTRMS_outdoorx_inter[, c(1:5, ncol(Data_PTRMS_outdoorx_inter), (6):(ncol(Data_PTRMS_outdoorx_inter) - 1))]
  # Append to the main data frame
  Data_PTRMS_outdoorx <- rbind(Data_PTRMS_outdoorx, Data_PTRMS_outdoorx_inter)
}

#Initialize empty data frame 
matt_quant_bck <- NULL
# Calculate 90th percentile for each sample and background by date and remove background
for (j in c(1:length(vec_date_unique))) {
  # Subset data for each unique date
  Data_PTRMS_outdoorx1 <- Data_PTRMS_outdoorx[vec_date == vec_date_unique[j], ]
  # Ensure measurement time is 25 sec
  Data_PTRMS_outdoorx1 <- Data_PTRMS_outdoorx1[Data_PTRMS_outdoorx1$chrono < 26, ]
  # Calculate 90th percentile for each group
  matt_quant <- data.frame(aggregate(Data_PTRMS_outdoorx1[, 7:34],
    by = list(Data_PTRMS_outdoorx1$Treatment, Data_PTRMS_outdoorx1$Plant,Data_PTRMS_outdoorx1$Repetition_measurment),
    FUN = function(x) quantile(x, probs = 0.9)
  ))



  Data_PTRMS_outdoorx1$Plant
  unique(matt_quant$Group.1)
  # Separate background data
  background <- matt_quant[matt_quant$Group.1 == "background", ]
  background
  # Remove background data from main data
  matt_quant <- matt_quant[!(matt_quant$Group.1 == "background"), ]

  # Copy data for modification
  matt_quant2 <- matt_quant

  matt_quant2 <- subset(matt_quant2, !(Group.2 == 29))
  # Subtract background values
  for (i in c(3:ncol(matt_quant2))) {
    matt_quant2[, i] <- matt_quant2[, i] - background[nrow(background), i]
  }
  # Add date column
  date <- rep(j, nrow(matt_quant2))
  matt_quant2 <- data.frame(date, matt_quant2)

  # Append to the main results data frame
  matt_quant_bck <- rbind(matt_quant_bck, matt_quant2)
}
## loop stops here
colnames(matt_quant_bck)[colnames(matt_quant_bck) == "Group.3"]  <- "Repetition_measurment"


vec_comb <- c("all","1","2","3","1_2","1_3","2_3")
output_conf_mat <- list()

for (xx in c(1:length(vec_comb))) {

slelector <- vec_comb[xx]

if(slelector == "all") {matt_quant_bck_inter <- matt_quant_bck}
if(slelector == "1") {matt_quant_bck_inter <-  matt_quant_bck[grep(1,matt_quant_bck$Repetition_measurment),]}
if(slelector == "2") {matt_quant_bck_inter <-  matt_quant_bck[grep(2,matt_quant_bck$Repetition_measurment),]}
if(slelector == "3") {matt_quant_bck_inter <-  matt_quant_bck[grep(3,matt_quant_bck$Repetition_measurment),]}
if(slelector == "1_2") {matt_quant_bck_inter <-  matt_quant_bck[-grep(3,matt_quant_bck$Repetition_measurment),]}
if(slelector == "1_3") {matt_quant_bck_inter <-  matt_quant_bck[-grep(2,matt_quant_bck$Repetition_measurment),]}
if(slelector == "2_3") {matt_quant_bck_inter <-  matt_quant_bck[-grep(1,matt_quant_bck$Repetition_measurment),]}


######################################################################
######################################################################
######################################################################

  
  matt_max_bck <- matt_quant_bck_inter 
 #matt_max_bck<-  matt_max_bck[grep(3,matt_max_bck$Repetition_measurment),]

  sample_id <- paste(matt_max_bck$Group.1,matt_max_bck$Group.2)
  
  matt_max_bck <- data.frame(sample_id,matt_max_bck)
  
  
  
  ########## heat map 
  
  matt_max_bck_stat <- matt_max_bck

  # Define a mapping of old names to new names
colnames(matt_max_bck_stat)[colnames(matt_max_bck_stat) == "Group.1"]  <- "treatment"


train_data <- matt_max_bck_stat[,6:33]
train_label <- matt_max_bck_stat$treatment
train_label[train_label=="healthy"] <- 0
train_label[train_label=="induced"] <- 1


bstSparse <- xgboost(data = as.matrix(train_data), 
                     label = train_label, 
                     max.depth = 5, eta = 0.05, 
                     nthread = 10, nrounds = 500, 
                     objective = "binary:logistic",
                     print_every_n = 100)

importance_matrix <- xgb.importance(model = bstSparse)

if(slelector == "all") {plot_imp <- xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")}
if(slelector == "all") {vec_imp_save <- c("date","Group.1","Group.2",importance_matrix$Feature[1:10])}

## auto selection
vec_imp <- c("date","treatment",importance_matrix$Feature[1:10])

matt_max_class_filter <- matt_max_bck_stat[,colnames(matt_max_bck_stat) %in% vec_imp]
########### gbm  classification 
library(xgboost)


valid_tree <- NULL
valid_xgboost <- NULL

for (i in c (1:nrow(matt_max_class_filter))) {
  
  
  train_data_full <- matt_max_class_filter[-i,] 
  train_data <- train_data_full[,3:ncol(train_data_full)]
  train_label <- train_data_full$treatment
  train_label[train_label=="healthy"] <- 0
  train_label[train_label=="induced"] <- 1

  bstSparse <- xgboost(data = as.matrix(train_data), 
                        label = train_label, 
                        max.depth = 5, eta = 0.05, 
                        nthread = 10, nrounds = 500, 
                        objective = "binary:logistic",
                        print_every_n = 100)
  
  test_data_full <- matt_max_class_filter[i,]  
  test_data <- test_data_full[,3:ncol(test_data_full)]
  test_label <- test_data_full$treatment
  test_label[test_label=="healthy"] <- 0
  test_label[test_label=="induced"] <- 1


  pred <- predict(bstSparse, as.matrix(test_data))

  pred[pred>=0.5] <- 1
  pred[pred<0.5] <- 0
  
  valid_xgboost <- c(valid_xgboost,pred)
  
}


real_data <-  matt_max_class_filter$treatment
real_data[real_data=="healthy"] <- 0
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
##### different exp�$
 
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
##### different exp�$

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




