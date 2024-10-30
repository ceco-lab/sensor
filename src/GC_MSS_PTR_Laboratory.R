# Load necessary libraries
library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)
library(car)
library(readxl)
library(writexl)
library(broom)
library(multcomp)
library(sandwich)


# Set working directory
setwd("./data")

# Set seed for reproducibility
set.seed(6)

#Load data----

#Specify the method
Method <- "MSS" #Change to "PTR", "GC" or "MSS"

if (Method == "PTR") {
  Data60 <- read_excel("PTR_Lab.xlsx")
  #Use the average of the 60-second measurements
  Data <- Data60 %>%
    group_by(Sample_ID, Treatment, Variety) %>%
    summarize(across(C3H5O:C15H27O, mean)) %>%
    as.data.frame()
} else if (Method == "GC") {
  Data<- read_excel("GC_Lab.xlsx")
} else if (Method == "MSS") {
  Data <- read_excel("MSS_Lab.xlsx")
}

# Specify the variety
Variety <- "Delprim" # Change to "Delprim" or "Aventicum" 

# Read data based on variety
if (Variety == "Delprim") {
  Data <- Data %>% filter(Variety=="Delprim")%>%droplevels()
} else if (Variety == "Aventicum") {
  Data <- Data %>% filter(Variety=="Aventicum")%>%droplevels()
}

# Prepare data----

# Change characters in factors and drop columns with only zeros
Data <- Data %>% 
  mutate(across(where(is.character), as.factor))%>% 
  select_if(~ !all(. == 0))

# Reordering the levels
order <- c("Control", "C. graminicola", "S. exigua", "S. frugiperda")
Data$Treatment <- factor(Data$Treatment, levels=order)

# Keep only features
X <- Data %>% dplyr::select(-c("Treatment","Variety","Sample_ID"))

# NMDS----

view(X)
str(X)

cols <- c("#0B775E","#E1BD6D", "darkred","darkblue")


if (Method == "GC") {
  nmds <- metaMDS(X, distance = "gower", autotransform = FALSE)
} else {
  nmds <- metaMDS(X, distance = "gower",autotransform=TRUE)
}

nmds_plot <- vegan::scores(nmds,display="site")%>%
  cbind(Data) %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_encircle(aes(group = Treatment, color = Treatment, fill = Treatment), alpha = 0.6, s_shape = 0.8, expand = 0.03) +
  geom_point(aes(color = Treatment)) +
  xlab(label = "NMDS1") +
  ylab(label = "NMDS2") +
  annotate("text", x = -Inf, y = Inf, label = paste0("stress: ", round(nmds$stress, digits = 3)), hjust = -0.1, vjust = 1.1) + 
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2)))

nmds_plot
ggsave(paste0("NMDS_",Method,"_Lab_", Variety, ".pdf"), dpi = 1200)

# PERMANOVA----
Treatment <- Data$Treatment
sink(paste0("Permanova_",Method,"_Lab_", Variety, ".txt"))
adonis2(X ~ Treatment, method = "gower")
sink()

#Random Forest----
X_rf <- X %>% 
  mutate(Treatment = Data$Treatment)

rp <- rfPermute(Treatment ~ ., X_rf,num.rep=500,ntree = 1000)
sink(paste0("RF_",Method,"_Lab_", Variety, ".txt"))
summary(rp)
sink()


### robust ANOVAs and pairwise comparisons - For GC and PTR data only---

if (Method == "GC" || Method == "PTR") {
  if(Method == "GC") {
    Data <- Data %>%
      mutate(Total = rowSums(X))
  } else {
    Data <- Data
  }
# Extracting the order of compounds
X <- Data %>% dplyr::select(-c("Treatment","Variety","Sample_ID"))
compound_order <- colnames(X)

# Get the mean and standard error for each compound by treatment
Summary <- Data %>%
  dplyr::select(-c(Variety,Sample_ID)) %>%
  pivot_longer(!c(Treatment), names_to = "compounds", values_to = "value") %>%
  mutate(compounds = factor(compounds, levels = compound_order)) %>%
  group_by(Treatment, compounds) %>%
  summarise(mean = mean(value), stand.error = sd(value) / sqrt(n())) 

# Cleaner summary
Clean_Summary <- Summary %>%
  filter(Treatment %in% c("Control", "C. graminicola", "S. exigua", "S. frugiperda")) %>%
  mutate(
    mean = ifelse(mean<0.01,format(mean, scientific = TRUE,digits = 3),round(mean,3)), 
    stand.error = ifelse(stand.error<0.01,format(stand.error, scientific = TRUE,digits = 3),round(stand.error,3))) %>% 
  pivot_wider(names_from = "Treatment", values_from = c("mean", "stand.error")) %>%
  arrange(match(compounds, compound_order)) %>%
  unite("Control", starts_with(c("mean_Control", "stand.error_Control")), sep = " ± ") %>%
  unite("C. graminicola", starts_with(c("mean_C. graminicola", "stand.error_C. graminicola")), sep = " ± ") %>%
  unite("S. exigua", starts_with(c("mean_S. exigua", "stand.error_S. exigua")), sep = " ± ") %>%
  unite("S. frugiperda", starts_with(c("mean_S. frugiperda", "stand.error_S. frugiperda")), sep = " ± ") 

# Initialize empty data frames to store results
comp_results_all <- data.frame(compounds = character(), Letters = character(), stringsAsFactors = FALSE)
anova_results <- data.frame(compounds = character(), p.value = numeric(), statistic = numeric(), stringsAsFactors = FALSE)

# Loop over each compound
for (compound in compound_order) {
  
  # Check if any treatment level contains only zeros = no variance
  zero_treatments <- Data %>%
    group_by(Treatment) %>%
    summarize(all_zeros = all(get(compound) == 0)) %>%
    filter(all_zeros == TRUE) %>%
    pull(Treatment)
  
  # Add 1 to half of the zeros in the treatment levels that contain only zeros
  if (length(zero_treatments) > 0) {
    for (treatment in zero_treatments) {
      zero_indices <- which(Data$Treatment == treatment & Data[[compound]] == 0)
      num_to_modify <- ceiling(length(zero_indices) / 2)
      modify_indices <- sample(zero_indices, num_to_modify)
      Data[modify_indices, compound] <- 1
    }
  }
  
  # Create model
  model <- lm(as.formula(paste(compound, " ~ Treatment")), data = Data)
  summary(model)
  
  # Perform ANOVA and extract p-value and statistic with heteroskedasticity-consistent (HC) standard error
  anova_result <- Anova(model, type = "II", vcov = vcovHC) %>%
    tidy() %>%
    filter(term == "Treatment") %>%
    dplyr::select(p.value, statistic) %>%
    mutate(compounds = compound)
  anova_results <- bind_rows(anova_results, anova_result)
  
  # Check if p-value is less than 0.05
  if (anova_result$p.value < 0.05) {
    # Compute emmeans and contrasts only if ANOVA p-value is significant
    emmeans_results <- glht(model, mcp(Treatment = "Tukey"), vcov = vcovHC)
    # Use cld() to obtain letters with fdr correction
    cld_result <- cld(emmeans_results, adjust = "fdr", Letters = letters, sort = FALSE, reversed = FALSE)
    # Extract letters and append to comp_results_all data frame
    Comp_letters <- data.frame(compounds = compound, Letters = cld_result$mcletters$Letters)
    comp_results_all <- bind_rows(comp_results_all, Comp_letters)
  }
}

# Round p-values and statistics
anova_results$p.value <- as.numeric(round((anova_results$p.value), 3))
anova_results$statistic <- as.numeric(round((anova_results$statistic), 2))

# Reshape comp_results_all to have three columns
comp_results_all2 <- comp_results_all %>%
  group_by(compounds) %>%
  mutate(row_id = row_number()) %>%
  pivot_wider(names_from = row_id, values_from = Letters) %>%
  ungroup()

# Merge comp_results_all with Clean_Summary
Clean_Summary_all <- left_join(Clean_Summary, comp_results_all2, by = "compounds") %>%
  mutate(
    Control = paste(Control, `1`),
    `C. graminicola` = paste(`C. graminicola`, `2`),
    `S. exigua` = paste(`S. exigua`, `3`),
    `S. frugiperda` = paste(`S. frugiperda`, `4`)) %>%
  dplyr::select(-`1`, -`2`, -`3`, -`4`) %>%
  mutate(
    Control = str_replace(Control, "NA$", ""),
    `C. graminicola` = str_replace(`C. graminicola`, "NA$", ""),
    `S. exigua`= str_replace(`S. exigua`, "NA$", ""),
    `S. frugiperda` = str_replace(`S. frugiperda`, "NA$", "")
  )

# Merge the anova_results with Clean_Summary_all
Clean_Summary_all2 <- left_join(Clean_Summary_all, anova_results, by = "compounds")

# Format p-values
Clean_Summary_all3 <- Clean_Summary_all2 %>%
  mutate(
    p.values = case_when(
      p.value < 0.001 ~ "< 0.001",
      p.value < 0.01 ~ "< 0.01",
      TRUE ~ as.character(round(p.value, 3)))) %>%
  dplyr::select(-p.value)

# Writing the combined summary to Excel
write_xlsx(Clean_Summary_all3, paste0("Table_",Method,"_Lab_", Variety, ".xlsx"))

} else if (Method == "MSS") {
  print("Method is MSS. Skipping the loop.")
}