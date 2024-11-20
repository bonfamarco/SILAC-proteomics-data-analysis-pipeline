#set working directory
setwd("C:/Users/mbonf/OneDrive - University of Strathclyde/PhD/Proteomics R/Y157C SILAC Proteomics_2023")
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
df_raw1<- read_excel("Experiment 4 2023_Y157C_RAW_SILAC_RS.xlsx")
df_raw2<- read_excel("Experiment 5 2023_Y157C_SILAC_RAW_RS.xlsx")
df_raw3<- read_excel("Experiment 6 2023_Y157C_SILAC_RAW_RS.xlsx")


#Filter dataset by number of peptides and SILAC ratios and add columns with log2 fold change
df1 <- df_raw1 %>%
  drop_na(`Heavy/Light`, `Medium/Heavy`, `Medium/Light`) %>%
  filter(`# Unique Peptides` >= 1) %>%
  mutate(
    log2_Heavy_Light = log2(`Heavy/Light`),
    log2_Medium_Heavy = log2(`Medium/Heavy`),
    log2_Medium_Light = log2(`Medium/Light`))

df2 <- df_raw2 %>%
  drop_na(`Heavy/Light`, `Medium/Heavy`, `Medium/Light`) %>%
  filter(`# Unique Peptides` >= 1) %>%
  mutate(
    log2_Heavy_Light = log2(`Heavy/Light`),
    log2_Medium_Heavy = log2(`Medium/Heavy`),
    log2_Medium_Light = log2(`Medium/Light`))

df3 <- df_raw3 %>%
  drop_na(`Heavy/Light`, `Medium/Heavy`, `Medium/Light`) %>%
  filter(`# Unique Peptides` >= 1) %>%
  mutate(
    log2_Heavy_Light = log2(`Heavy/Light`),
    log2_Medium_Heavy = log2(`Medium/Heavy`),
    log2_Medium_Light = log2(`Medium/Light`))


#Calculate thresholds
df1_thr_MH = (1.96*sd(df1$log2_Medium_Heavy))+mean(df1$log2_Medium_Heavy)
df1_thr_ML = (1.96*sd(df1$log2_Medium_Light))+mean(df1$log2_Medium_Light)
df1_thr_HL = (1.96*sd(df1$log2_Heavy_Light))+mean(df1$log2_Heavy_Light)

df2_thr_MH = (1.96*sd(df2$log2_Medium_Heavy))+mean(df2$log2_Medium_Heavy)
df2_thr_ML = (1.96*sd(df2$log2_Medium_Light))+mean(df2$log2_Medium_Light)
df2_thr_HL = (1.96*sd(df2$log2_Heavy_Light))+mean(df2$log2_Heavy_Light)

df3_thr_MH = (1.96*sd(df3$log2_Medium_Heavy))+mean(df3$log2_Medium_Heavy)
df3_thr_ML = (1.96*sd(df3$log2_Medium_Light))+mean(df3$log2_Medium_Light)
df3_thr_HL = (1.96*sd(df3$log2_Heavy_Light))+mean(df3$log2_Heavy_Light)

#Filter entries below threshold
df1_ML <- df1 %>%
  filter(`Medium/Light` >= df1_thr_ML)
df1_HL <- df1 %>%
  filter(`Heavy/Light` >= df1_thr_HL)
df1_MH <- df1 %>%
  filter(`Medium/Heavy` >= df1_thr_MH)

df2_ML <- df2 %>%
  filter(`Medium/Light` >= df2_thr_ML)
df2_HL <- df2 %>%
  filter(`Heavy/Light` >= df2_thr_HL)
df2_MH <- df2 %>%
  filter(`Medium/Heavy` >= df2_thr_MH)

df3_ML <- df3 %>%
  filter(`Medium/Light` >= df3_thr_ML)
df3_HL <- df3 %>%
  filter(`Heavy/Light` >= df3_thr_HL)
df3_MH <- df3 %>%
  filter(`Medium/Heavy` >= df3_thr_MH)


#Merge the datasets
ML <- bind_rows(df1_ML, df2_ML, df3_ML)
HL <- bind_rows(df1_HL, df2_HL, df3_HL)
MH <- bind_rows(df1_MH, df2_MH, df3_MH)

# Drop entries that are not present at least in duplicates
# Group the dataframe by the "Accession" column and count the number of observations in each group
# Filter the dataframe to only keep the groups with at least 2 observations (i.e., at least duplicate)
ML_grouped <- ML%>%
  group_by(Accession) %>%
  summarize(obs_count = n())

ML <- ML %>%
  left_join(ML_grouped %>% dplyr::select(Accession, obs_count), by = "Accession")

ML <- ML %>%
  filter(obs_count >= 2)

HL_grouped <- HL%>%
  group_by(Accession) %>%
  summarize(obs_count = n())

HL <- HL %>%
  left_join(HL_grouped %>% dplyr::select(Accession, obs_count), by = "Accession")

HL <- HL %>%
  filter(obs_count >= 2)


MH_grouped <- MH%>%
  group_by(Accession) %>%
  summarize(obs_count = n())

MH <- MH %>%
  left_join(MH_grouped %>% dplyr::select(Accession, obs_count), by = "Accession")

MH <- MH %>%
  filter(obs_count >= 2)

# Calculate the average SILAC ratios and drop any repetition
ML <- ML %>%
  distinct(Accession, .keep_all = TRUE) %>%
  group_by(Accession, Description) %>%
  summarize(
    avg_Heavy_Light = mean(`Heavy/Light`),
    avg_Medium_Light = mean(`Medium/Light`),
    avg_Medium_Heavy = mean(`Medium/Heavy`))

HL <- HL %>%
  distinct(Accession, .keep_all = TRUE) %>%
  group_by(Accession, Description) %>%
  summarize(
    avg_Heavy_Light = mean(`Heavy/Light`),
    avg_Medium_Light = mean(`Medium/Light`),
    avg_Medium_Heavy = mean(`Medium/Heavy`))

MH <- MH %>%
  distinct(Accession, .keep_all = TRUE) %>%
  group_by(Accession, Description) %>%
  summarize(
    avg_Heavy_Light = mean(`Heavy/Light`),
    avg_Medium_Light = mean(`Medium/Light`),
    avg_Medium_Heavy = mean(`Medium/Heavy`))

write.csv(ML, "SILAC_Y157C_Medium_Light.csv")
write.csv(HL, "SILAC_Y157C_Heavy_Light.csv")
write.csv(MH, "SILAC_Y157C_Medium_Heavy.csv")
