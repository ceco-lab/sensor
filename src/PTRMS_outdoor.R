# Load libraries
library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)

# Set working directory
setwd("./data")

# Set seed for reproducibility
set.seed(6)

# Load data
Data_PTRMS_outdoor <- read.csv("PTR_Outdoor.csv", sep = ",")

# Convert date to Date format and get unique dates
vec_date <- as.Date(Data_PTRMS_outdoor$Date, format = "%m/%d/%y")
vec_date_unique <- unique(vec_date)

# Get unique sample IDs
Sample_ID <- unique(Data_PTRMS_outdoor$Sample_ID)
colnames(Data_PTRMS_outdoor)
# Initialize empty data frame
Data_PTRMS_outdoorx <- NULL

# Build time column for each sample
for (i in c(1:length(Sample_ID))) {
  # Subset data for each sample ID
  Data_PTRMS_outdoorx_inter <- Data_PTRMS_outdoor[Data_PTRMS_outdoor$Sample_ID == Sample_ID[i], ]
  # Create a chronological column
  Data_PTRMS_outdoorx_inter$chrono <- c(1:nrow(Data_PTRMS_outdoorx_inter))
  # Reorder columns
  Data_PTRMS_outdoorx_inter <- Data_PTRMS_outdoorx_inter[, c(1:6, ncol(Data_PTRMS_outdoorx_inter), (7):(ncol(Data_PTRMS_outdoorx_inter) - 1))]
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
  matt_quant <- data.frame(aggregate(Data_PTRMS_outdoorx1[, 8:ncol(Data_PTRMS_outdoorx1)],
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

# Plot NMDS
cols <-c("darkgreen","#972D15")
nmds_plot <- vegan::scores(nmds)[1] %>% 
cbind(matt_quant_stat)%>%
ggplot(aes(x = sites.NMDS1, y = sites.NMDS2)) +
  xlab(label = "NMDS1") +
  ylab(label = "NMDS2") +
  geom_encircle(data = subset(matt_plot_nmds, treatment == "healthy"), color = cols[1], fill = cols[1], alpha = 0.6, s_shape = 0.8) +
  geom_encircle(data = subset(matt_plot_nmds, treatment == "induced"), color = cols[2], fill = cols[2], alpha = 0.6, s_shape = 0.8) +
  geom_point(data = subset(matt_plot_nmds, treatment == "healthy"), color = cols[1], fill = cols[1]) +
  geom_point(data = subset(matt_plot_nmds, treatment == "induced"), color = cols[2], fill = cols[2]) +
  xlim(-0.5, 0.5) +
  ylim(-0.45, 0.45) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
    theme(legend.position="none",
  aspect.ratio = 1)+
  annotate("text", x = -0.5, y = 0.45, label = paste0("stress: ", format(nmds$stress, digits = 3)), hjust = 0)
ggsave("NMDS_PTR_outdoor.png", dpi = 600,height=12,width=12, units = "cm")

# PERMANOVA analysis
sink("permanova_PTR_outdoor.txt")  # Redirect output to a file
adonis2(matt_quant_x ~ matt_quant_stat$treatment, permutations = 999, method = "gower")
sink()  # Stop redirecting output

#Random Forest analysis
matt_quant_stat_rf <- matt_quant_stat[, 4:(ncol(matt_quant_stat) - 1)]
matt_quant_stat_rf$treatment <- as.factor(matt_quant_stat$treatment)
sink("rf_PTR_outdoor.txt")  
rfPermute(treatment ~ ., data = matt_quant_stat_rf, na.action = na.omit, ntree = 500, num.rep = 150)
sink()  

#Univariate analysis

# Important variables selected
Variables <- c("C8H8N.", "C6H11.", "C7H9.", "C11H19.")

# Create an empty data frame to store the summary data
Summary_sub2 <- data.frame()

# Define a function to plot barplots and calculate Wilcoxon test
Univar <- function(variable) {

  # Extract relevant columns, pivot the data, and calculate mean and standard error
  Summary_sub <- matt_quant_stat %>%
    select(treatment, !!sym(variable)) %>%
    pivot_longer(!c(treatment), names_to = "compounds", values_to = "value") %>%
    group_by(treatment, compounds) %>%
    summarise(mean = mean(value), se = sd(value) / sqrt(n())) %>%
    as.data.frame()

  # Subset data for the specified variable
  X <- Summary_sub %>% filter(compounds == variable)

  # Perform Wilcoxon test
  sink(paste0("Wilcox", variable, ".txt"))  
  cat(paste0(variable, ": ", round((X[X$treatment == "induced", "mean"] - X[X$treatment == "healthy", "mean"]) / abs(X[X$treatment == "healthy", "mean"]) * 100, 0), "% increase\n"))
  print(wilcox.test(unlist(matt_quant_stat[, variable]) ~ treatment, data = matt_quant_stat))
  sink()

# Plot barplots
  gg <- ggplot(X, aes(x = treatment, y = mean, fill = treatment)) +
    geom_bar(position = position_dodge(), color = "black", stat = "identity", linewidth = 0.25, width = 0.85) +
    geom_errorbar(position = position_dodge(0.9), aes(ymin = mean - se, ymax = mean + se),
                  width = 0.3, linewidth = 0.25) +
    ggtitle(variable) +
    scale_fill_manual(values = c("green4", "red3"), labels = c("Healthy", "Induced")) +
    scale_x_discrete(labels = c("Control", "Damaged")) +
    ylab(expression(paste("Ions/s"))) +
    theme_classic() +
    theme(
      axis.line = element_line(size = 0.25, linetype = "solid"),
      axis.ticks = element_line(size = 0.25),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 12, vjust = 3),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.text = element_text(size = 9),
      legend.title = element_blank(),
      legend.key.size = unit(0.3, 'cm'),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )

  ggsave(paste0("Barplot_PTR_outdoor_", variable, ".pdf"), width = 4.25, height = 5, units = "cm", dpi = 600, scale = 2)
  return(Summary_sub)

}
# Apply the function to each variable
Barplots <- lapply(Variables, Univar)
Summary_sub2 <- do.call(rbind, Barplots) 
write.csv(Summary_sub2, "Summary_PTR_outdoor.csv", row.names = FALSE)


# Trace plots

# Subset data for the specified sample IDs and variables
Comparison <- Data_PTRMS_outdoorx|> filter(Sample_ID == "T18S3"|Sample_ID == "C18L3")
Variables <- c("C8H8N.", "C6H11.", "C7H9.", "C11H19.")

# Define a function to plot trace plots
create_trace_plot <- function(variable) {
  plot <- ggplot(Comparison, aes(x = chrono, y = .data[[variable]], group = Sample_ID, colour = Sample_ID)) +
    geom_line(linewidth = 1) +
    scale_x_continuous(expand = c(0, 0)) +
     scale_color_manual(values = c("green4", "red3"), labels = c("Healthy", "Induced")) +
    xlab("Time (s)") +
    ylab("Ions/s") +
    ggtitle(variable) +
    theme(
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 10),
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.position="none"
    )
  ggsave(paste0("Trace_", variable, ".pdf"), width = 5, height = 5, units = "cm", dpi = 600, scale = 2)
  }

# Apply the function to each variable
Traces <- lapply(Variables, create_trace_plot)

