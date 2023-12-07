# Load libraries
library(xgboost)
library(caret)
library(ggplot2)
library(dplyr)
library(viridis)

# Set working directory
setwd("./data")

# Set seed for reproducibility
set.seed(6)

# Load data
Data_PTRMS_outdoor <- read.csv("PTR_Outdoor.csv", sep = ",")

# Define combinations of measurement repetitions
vec_comb <- c("all", "1", "2", "3", "1_2", "1_3", "2_3")
output_conf_mat <- list()

for (xx in c(1:length(vec_comb))) {
  selector <- vec_comb[xx]

  # Subset backgrounds and samples based on selector
  Data_PTRMS_outdoor_sp <- Data_PTRMS_outdoor[!(Data_PTRMS_outdoor$Treatment == "background"), ]
  Data_PTRMS_outdoor_back <- Data_PTRMS_outdoor[(Data_PTRMS_outdoor$Treatment == "background"), ]

  if (selector == "all") {
    Data_PTRMS_outdoor_Rep <- Data_PTRMS_outdoor_sp
  }
  if (selector == "all") {
    Data_PTRMS_outdoor_Rep <- rbind(Data_PTRMS_outdoor_Rep, Data_PTRMS_outdoor_back)
  }

  if (selector == "1") {
    Data_PTRMS_outdoor_Rep <- Data_PTRMS_outdoor_sp[grep(1, Data_PTRMS_outdoor_sp$Repetition_measurment), ]
  }
  if (selector == "1") {
    Data_PTRMS_outdoor_rep_bck <- Data_PTRMS_outdoor_back[Data_PTRMS_outdoor_back$back_link %in% Data_PTRMS_outdoor_Rep$back_link, ]
  }
  if (selector == "1") {
    Data_PTRMS_outdoor_Rep <- rbind(Data_PTRMS_outdoor_Rep, Data_PTRMS_outdoor_rep_bck)
  }

  if (selector == "2") {
    Data_PTRMS_outdoor_Rep <- Data_PTRMS_outdoor_sp[grep(2, Data_PTRMS_outdoor_sp$Repetition_measurment), ]
  }
  if (selector == "2") {
    Data_PTRMS_outdoor_rep_bck <- Data_PTRMS_outdoor_back[Data_PTRMS_outdoor_back$back_link %in% Data_PTRMS_outdoor_Rep$back_link, ]
  }
  if (selector == "2") {
    Data_PTRMS_outdoor_Rep <- rbind(Data_PTRMS_outdoor_Rep, Data_PTRMS_outdoor_rep_bck)
  }

  if (selector == "3") {
    Data_PTRMS_outdoor_Rep <- Data_PTRMS_outdoor_sp[grep(3, Data_PTRMS_outdoor_sp$Repetition_measurment), ]
  }
  if (selector == "3") {
    Data_PTRMS_outdoor_rep_bck <- Data_PTRMS_outdoor_back[Data_PTRMS_outdoor_back$back_link %in% Data_PTRMS_outdoor_Rep$back_link, ]
  }
  if (selector == "3") {
    Data_PTRMS_outdoor_Rep <- rbind(Data_PTRMS_outdoor_Rep, Data_PTRMS_outdoor_rep_bck)
  }

  if (selector == "1_2") {
    Data_PTRMS_outdoor_Rep <- Data_PTRMS_outdoor_sp[-grep(3, Data_PTRMS_outdoor_sp$Repetition_measurment), ]
  }
  if (selector == "1_2") {
    Data_PTRMS_outdoor_rep_bck <- Data_PTRMS_outdoor_back[Data_PTRMS_outdoor_back$back_link %in% Data_PTRMS_outdoor_Rep$back_link, ]
  }
  if (selector == "1_2") {
    Data_PTRMS_outdoor_Rep <- rbind(Data_PTRMS_outdoor_Rep, Data_PTRMS_outdoor_rep_bck)
  }

  if (selector == "1_3") {
    Data_PTRMS_outdoor_Rep <- Data_PTRMS_outdoor_sp[-grep(2, Data_PTRMS_outdoor_sp$Repetition_measurment), ]
  }
  if (selector == "1_3") {
    Data_PTRMS_outdoor_rep_bck <- Data_PTRMS_outdoor_back[Data_PTRMS_outdoor_back$back_link %in% Data_PTRMS_outdoor_Rep$back_link, ]
  }
  if (selector == "1_3") {
    Data_PTRMS_outdoor_Rep <- rbind(Data_PTRMS_outdoor_Rep, Data_PTRMS_outdoor_rep_bck)
  }

  if (selector == "2_3") {
    Data_PTRMS_outdoor_Rep <- Data_PTRMS_outdoor_sp[-grep(3, Data_PTRMS_outdoor_sp$Repetition_measurment), ]
  }
  if (selector == "2_3") {
    Data_PTRMS_outdoor_rep_bck <- Data_PTRMS_outdoor_back[Data_PTRMS_outdoor_back$back_link %in% Data_PTRMS_outdoor_Rep$back_link, ]
  }
  if (selector == "2_3") {
    Data_PTRMS_outdoor_Rep <- rbind(Data_PTRMS_outdoor_Rep, Data_PTRMS_outdoor_rep_bck)
  }


  # Convert date to Date format and get unique dates
  vec_date <- as.Date(Data_PTRMS_outdoor_Rep$Date, format = "%m/%d/%y")
  vec_date_unique <- unique(vec_date)

  # Get unique sample IDs
  Sample_ID <- unique(Data_PTRMS_outdoor_Rep$Sample_ID)

  # Initialize empty data frame
  Data_PTRMS_outdoorx <- NULL

  # Build time column for each sample
  for (i in c(1:length(Sample_ID))) {
    # Subset data for each sample ID
    Data_PTRMS_outdoorx_inter <- Data_PTRMS_outdoor_Rep[Data_PTRMS_outdoor_Rep$Sample_ID == Sample_ID[i], ]
    # Create a chronological column
    Data_PTRMS_outdoorx_inter$chrono <- c(1:nrow(Data_PTRMS_outdoorx_inter))
    # Reorder columns
    Data_PTRMS_outdoorx_inter <- Data_PTRMS_outdoorx_inter[, c(1:6, ncol(Data_PTRMS_outdoorx_inter), (7):(ncol(Data_PTRMS_outdoorx_inter) - 1))]
    # Append to the main data frame
    Data_PTRMS_outdoorx <- rbind(Data_PTRMS_outdoorx, Data_PTRMS_outdoorx_inter)
  }


  # Initialize empty data frame
  matt_quant_bck <- NULL
  bck_unique <- unique(Data_PTRMS_outdoorx$back_link)

  # Calculate 90th percentile for each sample and background by date and remove background
  for (j in c(1:length(vec_date_unique))) {
    # Subset data for each unique date
    Data_PTRMS_outdoorx1 <- Data_PTRMS_outdoorx[vec_date == vec_date_unique[j], ]
    # Ensure measurement time is 25 sec
    Data_PTRMS_outdoorx1 <- Data_PTRMS_outdoorx1[Data_PTRMS_outdoorx1$chrono < 26, ]
    # Calculate 90th percentile for each group

    matt_quant <- data.frame(aggregate(Data_PTRMS_outdoorx1[, 8:35],
      by = list(Data_PTRMS_outdoorx1$Treatment, Data_PTRMS_outdoorx1$Plant),
      FUN = function(x) quantile(x, probs = 0.9)
    ))

    # Separate background data
    background <- matt_quant[matt_quant$Group.1 == "background", ]

    # Remove background data from main data
    matt_quant <- matt_quant[!(matt_quant$Group.1 == "background"), ]

    # Copy data for modification
    matt_quant2 <- matt_quant

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




  # First gradient boosting model to select 10th most important variables
  matt_max_bck <- matt_quant_bck

  sample_id <- paste(matt_max_bck$Group.1, matt_max_bck$Group.2)

  matt_max_bck <- data.frame(sample_id, matt_max_bck)

  matt_max_bck_stat <- matt_max_bck

  # Define a mapping of old names to new names
  colnames(matt_max_bck_stat)[colnames(matt_max_bck_stat) == "Group.1"] <- "treatment"

  train_data <- matt_max_bck_stat[, 5:32]
  train_label <- matt_max_bck_stat$treatment
  train_label[train_label == "healthy"] <- 0
  train_label[train_label == "induced"] <- 1


  bstSparse <- xgboost(
    data = as.matrix(train_data),
    label = train_label,
    max.depth = 5, eta = 0.05,
    nthread = 10, nrounds = 1000,
    objective = "binary:logistic",
    print_every_n = 100
  )

  importance_matrix <- xgb.importance(model = bstSparse)

  if (selector == "all") {
    plot_imp <- xgb.plot.importance(importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")
  }
  if (selector == "all") {
    vec_imp_save <- c("date", "Group.1", "Group.2", importance_matrix$Feature[1:10])
  }

  ## auto selection
  vec_imp <- c("date", "treatment", importance_matrix$Feature[1:10])

  matt_max_class_filter <- matt_max_bck_stat[, colnames(matt_max_bck_stat) %in% vec_imp]

  # gbm  classification
  valid_tree <- NULL
  valid_xgboost <- NULL

  for (i in c(1:nrow(matt_max_class_filter))) {
    train_data_full <- matt_max_class_filter[-i, ]
    train_data <- train_data_full[, 3:ncol(train_data_full)]
    train_label <- train_data_full$treatment
    train_label[train_label == "healthy"] <- 0
    train_label[train_label == "induced"] <- 1

    bstSparse <- xgboost(
      data = as.matrix(train_data),
      label = train_label,
      max.depth = 5, eta = 0.05,
      nthread = 10, nrounds = 500,
      objective = "binary:logistic",
      print_every_n = 100
    )

    test_data_full <- matt_max_class_filter[i, ]
    test_data <- test_data_full[, 3:ncol(test_data_full)]
    test_label <- test_data_full$treatment
    test_label[test_label == "healthy"] <- 0
    test_label[test_label == "induced"] <- 1


    pred <- predict(bstSparse, as.matrix(test_data))

    pred[pred >= 0.5] <- 1
    pred[pred < 0.5] <- 0

    valid_xgboost <- c(valid_xgboost, pred)
  }


  real_data <- matt_max_class_filter$treatment
  real_data[real_data == "healthy"] <- 0
  real_data[real_data == "induced"] <- 1


  table(real_data, valid_xgboost)
  confi <- caret::confusionMatrix(as.factor(real_data), as.factor(valid_xgboost))

  output_conf_mat[[xx]] <- confi
}
##### PLOTS

# Accuracy plot
data_plot_conf <- NULL

for (s in c(1:length(vec_comb))) {
  combi <- as.numeric(output_conf_mat[[s]]$overall[c(1, 3, 4)])
  data_plot_conf <- rbind(data_plot_conf, combi)
}
colnames(data_plot_conf) <- c("Accuracy", "AccuracyUpper", "AccuracyLower")

vec_group <- c("three_run_25s", "one_run_25s", "one_run_25s", "one_run_25s", "two_run_25s", "two_run_25s", "two_run_25s")
data_plot_conf <- data.frame(vec_group, vec_comb, data_plot_conf)

data_plot_conf_int1 <- data_plot_conf[, c(1, 2, 3)]
colnames(data_plot_conf_int1)[3] <- "Accuracy"
data_plot_conf_int2 <- data_plot_conf[, c(1, 2, 4)]
colnames(data_plot_conf_int2)[3] <- "Accuracy"
data_plot_conf_int3 <- data_plot_conf[, c(1, 2, 5)]
colnames(data_plot_conf_int3)[3] <- "Accuracy"
data_plot_conf_stack <- rbind(data_plot_conf_int1, data_plot_conf_int2, data_plot_conf_int3)


data_plot_conf_stack$vec_group <- factor(data_plot_conf_stack$vec_group, levels = c("one_run_25s", "two_run_25s", "three_run_25s"))

accura_plot <- ggplot() +
  geom_boxplot(data = data_plot_conf_stack, mapping = aes(x = vec_group, ymin = Accuracy, group = vec_group), width = 0.3, size = 1, color = "gray50") +
  ylim(.25, 1) +
  theme_classic() +
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
ggsave("accuracy.pdf", accura_plot)

# Importance plot
plot_imp_tot <- ggplot(plot_imp, aes(x = reorder(Feature, Importance), y = Importance, fill = Importance)) +
  labs(x = "compounds") +
  theme_classic() +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("Importance.pdf", plot_imp_tot)


# Confusion matrix
conf_table <- output_conf_mat[[1]]$table
colnames(conf_table) <- c("healty", "induced")
rownames(conf_table) <- c("healty", "induced")
fourfoldplot(conf_table)
