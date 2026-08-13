# Student project summary tables v6
# Max Korbmacher, Nov 2025

# Clean up
rm(list = ls(all.names = TRUE))
gc()
path <- "/Users/max/Documents/Local/MS/Student_Projects/"

# 0. Prep ####
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, effectsize, rstatix, openxlsx, janitor, stringr)

# Read and clean names
df <- read_excel(paste0(path,"New_Table.xlsx"))


# --- Summarize data ####
binary_vars <- c(
  "Snus",
  "Kjønn",
  "Røyk",
  "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
  "Alkohol_kode",
  "cutoff_5_depletion_anySession",
  "cutoff_5_depletion_Session2",
  "CRP_below_1"
)

continuous_vars <- c(
  "Alder_(år)",
  "CD19_celler_per_µL",
  "CD19_%",
  "BMI",
  "Total_antall_attakker",
  "Tid_siden_første_kliniske_hendelse_(måneder)",
  "Tid_siden_MS_diagnose_(måneder)",
  "Antall_attakker_siste_året",
  "Minimum_EDSS_score",
  "IgG",
  "IgM",
  "Antall_oligoklonale_bånd"
)

# Make sure continuous vars are numeric
df <- df %>%
  mutate(across(all_of(continuous_vars),
                ~ as.numeric(gsub("[^0-9\\.]", "", gsub(",", ".", .)))),
         across(all_of(binary_vars), as.character))

# Binary summary
binary_summary <- df %>%
  filter(`Måling nr.` == 1) %>%
  select(all_of(binary_vars)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(summary = paste0(value, ": ", n)) %>%
  group_by(variable) %>%
  summarise(summary = paste(summary, collapse = ", "), .groups = "drop")

# Continuous summary
continuous_summary <- df %>%
  filter(`Måling nr.` == 1) %>%
  select(all_of(continuous_vars)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(summary = paste0(round(mean(value, na.rm = TRUE), 2),
                             " ± ", round(sd(value, na.rm = TRUE), 2)), .groups = "drop")

# Combine and save
final_summary <- bind_rows(binary_summary, continuous_summary) %>% arrange(variable)
