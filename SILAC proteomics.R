#set working directory
setwd("C:/Users/mbonf/OneDrive - University of Strathclyde/PhD/Proteomics R/SILAC reanalysis")
getwd()

#Load all necessary libraries
library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(grid)
library(ggthemes)
library(ggrepel)


#Load raw datasets

df_raw1<- read_excel("C:/Users/mbonf/OneDrive - University of Strathclyde/PhD/Proteomics R/Unfiltered datasets/240516_Set1_filtered_0.xlsx", 
                                      sheet = "Set 1 raw")
df_raw2<- read_excel("C:/Users/mbonf/OneDrive - University of Strathclyde/PhD/Proteomics R/Unfiltered datasets/240516_Set2_MB.xlsx", 
                  sheet = "Raw data")
df_raw3<- read_excel("C:/Users/mbonf/OneDrive - University of Strathclyde/PhD/Proteomics R/Unfiltered datasets/240516_Set3_MB.xlsx", 
                  sheet = "Raw data")
df_raw4<- read_excel("C:/Users/mbonf/OneDrive - University of Strathclyde/PhD/Proteomics R/Unfiltered datasets/240516_Set4_MB.xlsx", 
                  sheet = "Raw data")
df_raw5<- read_excel("C:/Users/mbonf/OneDrive - University of Strathclyde/PhD/Proteomics R/Unfiltered datasets/240516_Set5_filtered_0.xlsx", 
                  sheet = "Set 5 raw")

# List of raw data frames
raw_dfs <- list(df_raw1, df_raw2, df_raw3, df_raw4, df_raw5)

# Initialize lists to store processed data frames and thresholds
processed_dfs <- list()
thresholds_MH <- list()
thresholds_ML <- list()
thresholds_HL <- list()
dfs_ML <- list()
dfs_HL <- list()
dfs_MH <- list()

# Process each data frame
for (i in seq_along(raw_dfs)) {
  df <- raw_dfs[[i]] %>%
    drop_na(`Heavy/Light`, `Medium/Heavy`, `Medium/Light`) %>%
    filter(`# Unique Peptides` >= 1) %>%
    mutate(
      log2_Heavy_Light = log2(`Heavy/Light`),
      log2_Medium_Heavy = log2(`Medium/Heavy`),
      log2_Medium_Light = log2(`Medium/Light`)
    )
  
  # Store processed data frames
  processed_dfs[[i]] <- df
  
  # Calculate thresholds
  thr_MH <- (1.96 * sd(df$log2_Medium_Heavy)) + mean(df$log2_Medium_Heavy)
  thr_ML <- (1.96 * sd(df$log2_Medium_Light)) + mean(df$log2_Medium_Light)
  thr_HL <- (1.96 * sd(df$log2_Heavy_Light)) + mean(df$log2_Heavy_Light)
  
  thresholds_MH[[i]] <- thr_MH
  thresholds_ML[[i]] <- thr_ML
  thresholds_HL[[i]] <- thr_HL
  
  # Filter entries below threshold
  dfs_ML[[i]] <- df %>% filter(`Medium/Light` >= thr_ML)
  dfs_HL[[i]] <- df %>% filter(`Heavy/Light` >= thr_HL)
  dfs_MH[[i]] <- df %>% filter(`Medium/Heavy` >= thr_MH)
}

# Merge the datasets
ML <- bind_rows(dfs_ML)
HL <- bind_rows(dfs_HL)
MH <- bind_rows(dfs_MH)

# Define a function to process each dataset
process_dataset <- function(data) {
  # Drop entries not present at least in triplicates and calculate averages
  data %>%
    group_by(Accession) %>%
    filter(n() >= 3) %>%  # Filter groups with at least 3 observations
    distinct(Accession, .keep_all = TRUE) %>%
    group_by(Accession, Description) %>%
    summarize(
      avg_Heavy_Light = mean(`Heavy/Light`, na.rm = TRUE),
      avg_Medium_Light = mean(`Medium/Light`, na.rm = TRUE),
      avg_Medium_Heavy = mean(`Medium/Heavy`, na.rm = TRUE),
      .groups = 'drop'  # Automatically drop grouping structure after summarizing
    )
}

# Apply the function to each dataset
ML_processed <- process_dataset(ML)
HL_processed <- process_dataset(HL)
MH_processed <- process_dataset(MH)

# Write the processed datasets to CSV files
write.csv(ML_processed, "Medium_Light.csv")
write.csv(HL_processed, "Heavy_Light.csv")
write.csv(MH_processed, "Medium_Heavy.csv")

#======= HISTOGRAMS==========
df1HL_his <- ggplot(df1, aes(x = log2_Heavy_Light)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.3, color = "black", fill = "white") +
  stat_function(
    fun = dnorm, # Gaussian density function
    args = list(mean = mean(df1$log2_Heavy_Light), sd = sd(df1$log2_Heavy_Light)), # Mean and standard deviation
    aes(color = "Gaussian Curve")
  ) +
  geom_vline(xintercept = c(-df1_thr_HL, df1_thr_HL), linetype = "dashed", color = "black") +
  labs(
    title = "Histogram the distribution of significant PAR4 interacting proteins", 
    x = "log2 Fold Change",
    y = "Count",
    caption = 'Cut-off lines drawn at 95% confidence intervals') + 
  xlim(c(-10, 10))

df1HL_his + theme_Publication()


#Set themes
theme_Publication <- function(base_size=14, base_family="helvetica") {
  
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.margin = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
