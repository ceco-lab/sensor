library(vegan)
library(rfPermute)
library(ggplot2)
library(tidyverse)
library(ggalt)
library(car)
library(openxlsx)
library(broom)
library(emmeans)
library(multcomp)


setwd("./data")

set.seed(6)

#Select plant variety and load corresponding data 

Variety="Delprim"
#Variety="Aventicum"

if (Variety == "Delprim") {
  GCMS <- read.csv("GCMS_Delprim.csv", sep = ",") # Delprim data
} else if (Variety == "Aventicum") {
  GCMS <- read.csv("GCMS_Aventicum.csv", sep = ",") # Aventicum data
}

GCMS$treatment <- factor(GCMS$treatment, levels = c("Control", "CG", "SE", "SF"))

#Keep only features
GCMS_x <- GCMS[, 2:ncol(GCMS)]

#NMDS
cols <-c("#D8B70A","darkgreen","#972D15","royalblue4")

nmds <- metaMDS(GCMS_x, distance = "gower", autotransform = FALSE)
nmds_plot <- vegan::scores(nmds)[1] %>%  
  cbind(GCMS) %>%
  ggplot(aes(x = sites.NMDS1, y = sites.NMDS2)) +
  xlab(label = "NMDS1") +
  ylab(label = "NMDS2") +
  geom_encircle(aes(group = treatment, color = treatment, fill = treatment), alpha = 0.6, s_shape = 0.8,expand=0.03) +
  geom_point(aes(color = treatment)) +
  annotate("text", x = -0.35, y = 0.3, label = paste0("stress: ", round(nmds$stress, digits = 4)), hjust = 0) +
  xlim(-0.4, 0.5) + ylim(-0.3, 0.3) +
    coord_fixed(ratio = 1) +
  theme_bw() + scale_fill_manual(values = cols) + scale_color_manual(values = cols)+
  theme(legend.position="none",
  axis.title=element_blank(),
  aspect.ratio = 1)

nmds_plot
ggsave(paste0("NMDS_GC_", Variety, ".pdf"), width=8, height=8,dpi = 600)

# PERMANOVA 
sink(paste0("permanova_GC_",Variety,".txt"))  
adonis2(GCMS_x ~ GCMS$treatment, permutations = 999, method = "gower")
sink()  

#Random forest 
GCMS_x$treatment <- as.factor(GCMS$treatment)
sink(paste0("rf_GC_",Variety,".txt"))  # Redirect output to a file
rfPermute(treatment ~ ., data = GCMS_x, na.action = na.omit, ntree = 500, num.rep = 150)
sink()  


GCMS$total <- rowSums(GCMS[, !colnames(GCMS) %in% "treatment"])

# Extracting the order of compounds from the original "PTRMS" dataframe
compound_order <- colnames(GCMS)[-c(which(colnames(GCMS) %in% c( "treatment")))]

Summary <- GCMS %>%
    pivot_longer(!c(treatment), names_to = "compounds", values_to = "value") %>%
    group_by(treatment, compounds) %>%
    summarise(mean = mean(value), stand.error = sd(value)/sqrt(n()))  # Using n() for simplicity

# Creating a cleaner summary
Clean_Summary <- Summary %>%
  filter(treatment %in% c("Control", "CG", "SF", "SE")) %>%
  mutate(
    mean = round(mean, 2),  # Round mean to 2 decimal places
    stand.error = round(stand.error, 2)  # Round standard error to 2 decimal places
  ) %>%
  pivot_wider(names_from = "treatment", values_from = c("mean", "stand.error")) %>%
  arrange(match(compounds, compound_order)) %>%
  dplyr::select(
    compounds,
    starts_with("mean_Control"),
    starts_with("stand.error_Control"),
    starts_with("mean_CG"),
    starts_with("stand.error_CG"),
    starts_with("mean_SE"),
    starts_with("stand.error_SE"),
    starts_with("mean_SF"),
    starts_with("stand.error_SF")
  ) %>%
  unite("Control", starts_with(c("mean_Control", "stand.error_Control")), sep = " ± ") %>%
  unite("CG", starts_with(c("mean_CG", "stand.error_CG")), sep = " ± ") %>%
  unite("SE", starts_with(c("mean_SE", "stand.error_SE")), sep = " ± ") %>%
  unite("SF", starts_with(c("mean_SF", "stand.error_SF")), sep = " ± ") %>%
  dplyr::select(compounds, Control, CG, SE, SF)



# Initialize an empty data frame to store comp results
comp_results_all <- data.frame(compounds = character(), Letters = character(), stringsAsFactors = FALSE)
anova_results <- data.frame(compounds = character(), p.value = numeric(), statistic = numeric(), stringsAsFactors = FALSE)


# Loop over each compound
for (compound in compound_order) {
  # Create lm model
  model <- lm(as.formula(paste(compound, " ~ treatment")), data = GCMS)

  # Perform ANOVA and extract p-value
  anova_result <- Anova(model, type = "II") %>%
    tidy() %>%
    filter(term == "treatment") %>%
    dplyr::select(p.value, statistic) %>%
    mutate(compounds = compound)
    anova_results <- bind_rows(anova_results, anova_result)

  # Check if p-value is less than 0.05
  if (anova_result$p.value < 0.05) {
    # Compute emmeans and contrasts only if ANOVA p-value is significant
    emmeans_result <- emmeans(model, ~ treatment)
        # Use cld() to obtain letters with Bonferroni correction
    cld_result <- cld(emmeans_result, adjust="bonferroni",Letters = letters,sort=FALSE,reversed=FALSE)

    # Extract letters and append to comp_results_all data frame
    Comp_letters <- data.frame(
      compounds = compound,
      Letters = as.character(cld_result$.group)
    )
    comp_results_all <- bind_rows(comp_results_all, Comp_letters)
  }
}

anova_results$p.value <- as.numeric(round((anova_results$p.value),3))
anova_results$statistic <- as.numeric(round((anova_results$statistic),2))

# Reshape comp_results_all to have three columns
comp_results_all2 <- comp_results_all %>%
  group_by(compounds) %>%
  mutate(row_id = row_number()) %>%
  pivot_wider(names_from = row_id, values_from = Letters) %>%
  ungroup()

# Merge comp_results_all with Clean_Summary_all
Clean_Summary_all <- left_join(Clean_Summary, comp_results_all2, by = "compounds") %>%
  mutate(
    Control = paste(Control, `1`),
    CG = paste(CG, `2`),
    SE = paste(SE, `3`),
    SF = paste(SF, `4`)
  ) %>%
  dplyr::select(-`1`, -`2`, -`3`, -`4`) %>%
     mutate(
    Control = str_replace(Control, "NA$", ""),
    CG = str_replace(CG, "NA$", ""),
    SE = str_replace(SE, "NA$", ""),
    SF = str_replace(SF, "NA$", "")
  )

# Merge the anovCa_results with Clean_Summary
Clean_Summary_all2 <- left_join(Clean_Summary_all, anova_results, by = "compounds")
Clean_Summary_all3 <- Clean_Summary_all2 %>%
  mutate(
    compounds = paste0(gsub("^X\\.|\\.\\.$", "", compounds)),
    p.values = case_when(
      p.value < 0.001 ~ "< 0.001",
      p.value < 0.01 ~ "< 0.01",
      TRUE ~ as.character(round(p.value, 3))  # Format other p-values to three decimal places
    )
  ) %>%
  dplyr::select(-p.value) 
# Writing the combined summary to Excel
write.xlsx(Clean_Summary_all3, paste0("Table_GC_lm", Variety, ".xlsx"))




