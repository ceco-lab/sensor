#contact: marine.mamin@gmail.com

# Load necessary libraries
library(tidyverse)
library(readxl)
library(writexl)

# Set working directory
setwd("./data")

#Load data
Data1 <- read_excel("PTR_Outdoors.xlsx")

#Prepare data--------

#Ensure 60 first seconds of measurements are kept and take the average
Data <- Data1 %>%
  group_by(Sample_ID) %>%
  slice_head(n = 60) %>%  
  ungroup() %>%
  mutate(across(where(is.character), as.factor),
         Plant = as.factor(Plant)) %>%
  group_by(Sample_ID, Treatment, Distance, Day, Plant) %>%
  summarize(across(C3H5O:C15H27O, mean), .groups = "drop") %>%  # Calculate the mean across specified columns
  as.data.frame()

#Summary and Wilcoxon tests-----

#Create variable combi
X1 <- Data %>%
  dplyr::select(-c("Day","Sample_ID","Plant")) %>%
  mutate(combi = paste0(Treatment,"_",Distance)) %>%
  dplyr::select(!c(Treatment,Distance)) 

#List of compounds
compound_order <- colnames(X1[, -which(colnames(X1) == c("combi"))])

#Get the mean and standard error for each compound by treatment
Summary1 <- X1 %>%
  pivot_longer(!c(combi), names_to = "compounds", values_to = "value") %>%
  mutate(compounds = factor(compounds, levels = compound_order))%>%
  group_by(combi,compounds) %>%
  summarise(mean = mean(value), stand.error = sd(value)/sqrt(n()))

# Creating a cleaner summary
Clean_Summary1 <- Summary1 %>%
  mutate(
    mean_num = mean,
    mean = ifelse(mean<0.01,format(mean, scientific = TRUE,digits = 3),round(mean,3)), 
    stand.error = ifelse(stand.error<0.01,format(stand.error, scientific = TRUE,digits = 3),round(stand.error,3))) %>% 
  pivot_wider(names_from = "combi", values_from = c("mean", "mean_num", "stand.error"))%>%
  mutate(
  Percentage_increase_12cm = paste0(round(((mean_num_Damaged_12cm - mean_num_Control_12cm)/mean_num_Control_12cm) * 100), "%"),
  Percentage_increase_45cm = paste0(round(((mean_num_Damaged_45cm- mean_num_Control_45cm)/mean_num_Control_45cm) * 100), "%"))%>% 
  dplyr::select(-contains("mean_num")) %>% 
  unite("Damaged_12cm_mean", c("mean_Damaged_12cm", "stand.error_Damaged_12cm"), sep = " ± ")%>%
  unite("Damaged_45cm_mean", c("mean_Damaged_45cm", "stand.error_Damaged_45cm"), sep = " ± ") %>%
  unite("Control_12cm_mean", c("mean_Control_12cm", "stand.error_Control_12cm"), sep = " ± ")%>%
  unite("Control_45cm_mean", c("mean_Control_45cm", "stand.error_Control_45cm"), sep = " ± ")

# Prepare the 1-2cm dataframe
D12cm1 <- Data %>%
  dplyr::select(-c("Day", "Sample_ID","Plant")) %>%
  filter(Distance == "12cm") %>%
  dplyr::select(-Distance)%>%as.data.frame()

# Prepare the 4-5cm dataframe
D45cm1 <- Data %>%
  dplyr::select(-c("Day", "Sample_ID","Plant")) %>%
  filter(Distance == "45cm") %>%
  dplyr::select(-Distance)%>%as.data.frame()

# Wilcoxon test 1-2 cm
results12cm1 <- apply(D12cm1[, sapply(D12cm1, is.numeric)], 2, function(column) {
  test_result12cm1 <- wilcox.test(column ~ D12cm1$Treatment)
  return(c(W = test_result12cm1$statistic, p_value = test_result12cm1$p.value))
})
results_12cm1 <- as.data.frame(t(results12cm1))

# Wilcoxon test 4-5 cm
results45cm1 <- apply(D45cm1[, sapply(D45cm1, is.numeric)], 2, function(column) {
  test_result45cm1 <- wilcox.test(column ~ D45cm1$Treatment)
  return(c(W = test_result45cm1$statistic, p_value = test_result45cm1$p.value))
})
results_45cm1 <- as.data.frame(t(results45cm1))

# Format p-values
Results_clean_12cm1 <- results_12cm1%>%
  mutate(
    p_value = case_when(
      p_value < 0.001 ~ "< 0.001",
      p_value < 0.01 ~ "< 0.01",
      TRUE ~ as.character(round(p_value, 3))))  

Results_clean_45cm1 <- results_45cm1%>%
  mutate(
    p_value = case_when(
      p_value < 0.001 ~ "< 0.001",
      p_value < 0.01 ~ "< 0.01",
      TRUE ~ as.character(round(p_value, 3))))  

#Add column compounds from rownames
Results_clean_12cm1$compounds <- rownames(Results_clean_12cm1)
Results_clean_45cm1$compounds <- rownames(Results_clean_45cm1)

#Merge results
final_results1 <- merge(Results_clean_12cm1,Results_clean_45cm1,by = "compounds", suffixes = c("_12cm", "_45cm"))
final_summary1 <- merge(Clean_Summary1, final_results1, by = "compounds")

# Reorder final_summary based on compound_order
final_summary1 <- final_summary1 %>%
  mutate(compound_order = factor(compounds, levels = compound_order)) %>%
  arrange(compound_order) %>%
  dplyr::select(-compound_order)

# Writing the combined summary to Excel
write_xlsx(final_summary1, "Summary_PTR_Outdoors.xlsx")



#Test effect of bagging plants--------

# Set seed for reproducibility
set.seed(10)

#List of compounds
compounds <- Data %>% dplyr::select(where(is.numeric)) %>% colnames()

#Create variable Plantok
Data$Plantok <-paste0(Data$Treatment,"_",Data$Plant,"_",Data$Day)

# Filter out the control plants that were bagged
Bagged<- Data %>%  dplyr::filter(Plantok == "Control_17_1"|
                                  Plantok == "Control_19_1"|
                                  Plantok == "Control_14_1"|
                                  Plantok == "Control_1_1"|
                                  Plantok == "Control_20_1"|
                                  
                                  Plantok == "Control_1_2"|
                                  Plantok == "Control_4_2"|
                                  Plantok == "Control_8_2"|
                                  Plantok == "Control_20_2"|
                                  Plantok == "Control_23_2"|
                                  
                                  Plantok == "Control_23_3"|
                                  Plantok == "Control_29_3"|
                                  Plantok == "Control_13_3"|
                                  Plantok == "Control_19_3"|
                                  Plantok == "Control_11_3" ) %>% 
  mutate(Bag ="present")

filtered_data <- Data %>%
  filter(Treatment == "Control", !(Plantok %in% Bagged$Plantok))

# Filter data for each day
filtered_data_day1 <- filtered_data %>% dplyr::filter(Day == "1")
filtered_data_day2 <- filtered_data %>% filter(Day == "2")
filtered_data_day3 <- filtered_data %>% filter(Day == "3")

# Select 5 random unique Plantok values for each day
random_Plantok_day1 <- sample(unique(filtered_data_day1$Plantok), size = 5)
random_Plantok_day2 <- sample(unique(filtered_data_day2$Plantok), size = 5)
random_Plantok_day3 <- sample(unique(filtered_data_day3$Plantok), size = 5)

# Get all rows corresponding to these Plantok values
random_rows_day1 <- filtered_data_day1 %>% filter(Plantok %in% random_Plantok_day1)
random_rows_day2 <- filtered_data_day2 %>% filter(Plantok %in% random_Plantok_day2)
random_rows_day3 <- filtered_data_day3 %>% filter(Plantok %in% random_Plantok_day3)

# Combine the randomly selected rows for each day
Unbagged <- rbind(random_rows_day1, random_rows_day2, random_rows_day3)

# Add Bag column to Unbagged data
Unbagged <- Unbagged %>% 
  as.data.frame() %>% 
  mutate(Bag = "absent")

#Combine Bagged and Unbagged dataframes
TestBag <- rbind(Bagged,Unbagged)

#Convert Bag to factor
TestBag$Bag <- as.factor(TestBag$Bag)

#Create combination Bag and Distance
X2 <- TestBag %>%
  dplyr::select(all_of(compounds),"Bag","Distance") %>%
  mutate(combi = paste0(Bag,"_",Distance)) %>%
  dplyr::select(!c(Bag,Distance)) 

#List of compounds
compounds_order <- colnames(X2[, -which(colnames(X2) == c("combi"))])

#Get the mean and standard error for each compound by treatment
Summary2 <- X2 %>%
  pivot_longer(!c(combi), names_to = "compounds", values_to = "value") %>%
  mutate(compounds = factor(compounds, levels = compounds_order))%>%
  group_by(combi,compounds) %>%
  summarise(mean = mean(value), stand.error = sd(value, na.rm = TRUE)/sqrt(n()), .groups = 'drop') 

# Creating a cleaner summary
Clean_Summary2 <- Summary2 %>%
  mutate(
    mean_num = mean,
    mean = ifelse(mean<0.01,format(mean, scientific = TRUE,digits = 3),round(mean,3)), 
    stand.error = ifelse(stand.error<0.01,format(stand.error, scientific = TRUE,digits = 3),round(stand.error,3)))%>% 
  pivot_wider(names_from = "combi", values_from = c("mean_num","mean", "stand.error"))%>%
  mutate(
    Percentage_increase_12cm = paste0(round(((mean_num_present_12cm - mean_num_absent_12cm)/mean_num_absent_12cm) * 100), "%"),
    Percentage_increase_45cm = paste0(round(((mean_num_present_45cm- mean_num_absent_45cm)/mean_num_absent_45cm) * 100), "%")) %>% 
  dplyr::select(-contains("mean_num")) %>% 
  unite("Yes_mean_12cm", c("mean_present_12cm", "stand.error_present_12cm"), sep = " ± ")%>%
  unite("Yes_mean_45cm", c("mean_present_45cm", "stand.error_present_45cm"), sep = " ± ")%>%
  unite("No_mean_12cm", c("mean_absent_12cm", "stand.error_absent_12cm"), sep = " ± ")%>%
  unite("No_mean_45cm", c("mean_absent_45cm", "stand.error_absent_45cm"), sep = " ± ")

# Prepare the 1-2cm dataframe
D12cm2<- TestBag%>%
  dplyr::select(c(all_of(compounds),"Bag","Distance")) %>%
  dplyr::filter(Distance =="12cm")%>%
  dplyr::select(!Distance) 

# Prepare the 4-5cm dataframe
D45cm2<- TestBag%>%
  dplyr::select(c(all_of(compounds),"Bag","Distance")) %>%
  dplyr::filter(Distance =="45cm")%>%
  dplyr::select(!Distance) 

# Wilcoxon test 1-2 cm
results12cm2 <- apply(D12cm2[, sapply(D12cm2, is.numeric)], 2, function(column) {
  test_result12cm2 <- wilcox.test(column ~ D12cm2$Bag)
  return(c(W = test_result12cm2$statistic, p_value = test_result12cm2$p.value))
})
results_12cm2 <- as.data.frame(t(results12cm2))

# Wilcoxon test 4-5 cm
results45cm2 <- apply(D45cm2[, sapply(D45cm2, is.numeric)], 2, function(column) {
  test_result45cm2 <- wilcox.test(column ~ D45cm2$Bag)
  return(c(W = test_result45cm2$statistic, p_value = test_result45cm2$p.value))
})
results_45cm2 <- as.data.frame(t(results45cm2))

# Format p-values
Results_clean_12cm2 <- results_12cm2%>%
  mutate(
    p_value = case_when(
      p_value < 0.001 ~ "< 0.001",
      p_value < 0.01 ~ "< 0.01",
      TRUE ~ as.character(round(p_value, 3))))  

Results_clean_45cm2 <- results_45cm2%>%
  mutate(
    p_value = case_when(
      p_value < 0.001 ~ "< 0.001",
      p_value < 0.01 ~ "< 0.01",
      TRUE ~ as.character(round(p_value, 3))))  

#Add column compounds from rownames
Results_clean_12cm2$compounds <- rownames(Results_clean_12cm2)
Results_clean_45cm2$compounds <- rownames(Results_clean_45cm2)

#Merge the results
final_results2 <- merge(Results_clean_12cm2,Results_clean_45cm2, by = "compounds", suffixes = c("_12cm", "_45cm"))

#Merge summary and wilcoxon tests results
final_summary2 <- merge(Clean_Summary2, final_results2, by = "compounds")

# Reorder based on compound_order 
final_summary2 <- final_summary2 %>%
  mutate(compound_order = factor(compounds, levels = compound_order)) %>%
  arrange(compound_order) %>%
  dplyr::select(-compound_order)

# Writing the combined summary to Excel
write_xlsx(final_summary2, "Summary_BagEffect.xlsx")



