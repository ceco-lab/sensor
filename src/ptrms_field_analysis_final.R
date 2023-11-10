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
    by = list(Data_PTRMS_outdoorx1$Treatment, Data_PTRMS_outdoorx1$Plant),
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

# Copy data for statistical analysis
matt_quant_stat <- matt_quant_bck

# Create a sample column
matt_quant_stat$sample <- paste(matt_quant_stat$Group.1, matt_quant_stat$Group.2)
# Rename the second column to "treatment"
colnames(matt_quant_stat)[2] <- "treatment"

# Select columns for NMDS and PERMANOVA analysis (all ion masses)
matt_quant_x <- matt_quant_stat[, 4:(ncol(matt_quant_stat) - 1)]
# Set row names as sample names
row.names(matt_quant_x) <- matt_quant_stat$sample
# Set negative values to zero
matt_quant_x[matt_quant_x < 0] <- 0


# Perform NMDS analysis
nmds <- metaMDS(matt_quant_x, distance = "gower")

# Prepare data for NMDS plot
matt_plot_nmds <- vegan::scores(nmds)[1] %>% cbind(matt_quant_stat)
matt_plot_nmds <- matt_plot_nmds[-40, ]

# Plot NMDS
cols <- c("khaki", "darkred")
nmds_plot <- ggplot(data = matt_plot_nmds, aes(x = sites.NMDS1, y = sites.NMDS2)) +
  xlab(label = "NMDS1") +
  ylab(label = "NMDS2") +
  geom_encircle(data = subset(matt_plot_nmds, treatment == "healthy"), color = cols[2], fill = cols[2], alpha = 0.5, s_shape = 0.8) +
  geom_encircle(data = subset(matt_plot_nmds, treatment == "induced"), color = cols[1], fill = cols[1], alpha = 0.5, s_shape = 0.8) +
  geom_point(data = subset(matt_plot_nmds, treatment == "healthy"), color = cols[2], fill = cols[2]) +
  geom_point(data = subset(matt_plot_nmds, treatment == "induced"), color = cols[1], fill = cols[1]) +
  xlim(-0.6, 0.7) +
  ylim(-0.5, 0.5) +
  theme_bw() +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  annotate("text", x = -0.6, y = 0.5, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0)
nmds_plot

# PERMANOVA analysis
ado <- adonis2(matt_quant_x ~ matt_quant_stat$treatment, permutations = 999, method = "gower")

# Prepare data for Random Forest analysis
# Select input variables
matt_quant_stat_rf <- matt_quant_stat[, 4:(ncol(matt_quant_stat) - 1)]
# Set treatment as output variable
matt_quant_stat_rf$treatment <- as.factor(matt_quant_stat$treatment)

# Perform Random Forest analysis with permutation test
rp <- rfPermute(treatment ~ ., data = matt_quant_stat_rf, na.action = na.omit, ntree = 1000, num.rep = 150)

# Plot variable importance
var_imp <- plotImportance(rp, scale = TRUE, size = 3)

# print output

ado

rp

title <- "PTRMS field result"
grid::grid.text(title, x = (0.5), y = (0.6))
nmds_plot
plotImportance(rp, scale = TRUE, size = 3)
