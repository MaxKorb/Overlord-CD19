# Student project summary -- v3
# Max Korbmacher, Aug 2025

# Clean up
rm(list = ls(all.names = TRUE))
gc()
path <- "/Users/max/Documents/Local/MS/Student_Projects/"

# 0. Prep ####
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readxl, effectsize, openxlsx, rcompanion, boot, rstatix, effsize)

# Read data
df <- read_excel(paste0(path,"2025-08-22-table.xlsx"))
df$CD19_celler_per_µL = ifelse(df$CD19_celler_per_µL == "",NA,df$CD19_celler_per_µL)
df = df[!df$PasientID %in% c(004, 012, 030),]

# Wrangle data
df$Røyk = as.numeric(df$Røyk)
df$`Alder_(år)` <- as.numeric(df$`Alder_(år)`)

fully_missing_participants <- df %>%
  group_by(PasientID) %>%
  summarise(is_variable_fully_missing = all(!complete.cases(CD19_celler_per_µL))) %>%
  filter(is_variable_fully_missing)

print(fully_missing_participants)
df = df[!df$PasientID %in% fully_missing_participants$PasientID,]

# Categorize into depletion groups ####
df <- df %>%
  mutate(CD19_celler_per_µL = as.numeric(gsub("[^0-9\\.]", "", gsub(",", ".", CD19_celler_per_µL)))) %>%
  group_by(PasientID) %>%
  mutate(
    baseline = first(CD19_celler_per_µL),
    drop5_this  = ifelse(!is.na(CD19_celler_per_µL) & `Måling nr.` > 1 & CD19_celler_per_µL > 5, 0, 1),
    drop5_2  = ifelse(!is.na(CD19_celler_per_µL) & `Måling nr.` == 2 & CD19_celler_per_µL < 6, 1, 0),
    drop5_2  = ifelse(is.na(CD19_celler_per_µL) == T & `Måling nr.` == 2, NA, drop5_2),
    cutoff_5_depletion_anySession   = ifelse(any(drop5_this == 0), 0, 1),
    cutoff_5_depletion_Session2     = ifelse(any(drop5_2 == 1), 1, 0)
  ) %>%
  ungroup()

# --- Diagnostic summary ####
diagnostic_summary <- list(
  patients = df %>%
    group_by(PasientID) %>%
    summarise(
      cutoff_5_anySession   = first(cutoff_5_depletion_anySession),
      cutoff_5_Session2     = first(cutoff_5_depletion_Session2),
      .groups = "drop"
    ) %>%
    summarise(across(starts_with("cutoff"), sum, na.rm = TRUE)),
  
  measurements = df %>%
    summarise(
      n_measures_5_this  = sum(drop5_this, na.rm = TRUE),
      n_measures_5_2     = sum(drop5_2, na.rm = TRUE)
    )
)
diagnostic_summary

# --- Summarize data ####
binary_vars <- c(
  "Snus","Kjønn","Røyk","Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
  "Alkohol_kode","cutoff_5_depletion_anySession","cutoff_5_depletion_Session2","CRP_below_1"
)

continuous_vars <- c(
  "Alder_(år)","CD19_celler_per_µL","CD19_%","BMI","Total_antall_attakker",
  "Tid_siden_første_kliniske_hendelse_(måneder)","Tid_siden_MS_diagnose_(måneder)",
  "Antall_attakker_siste_året","Minimum_EDSS_score","IgG","IgM","Antall_oligoklonale_bånd"
)

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

final_summary <- bind_rows(binary_summary, continuous_summary) %>% arrange(variable)
write.xlsx(final_summary, paste0(path,"Descriptives.xlsx"))
write.csv(df, paste0(path,"New_Table.csv"))
write.xlsx(df, paste0(path,"New_Table.xlsx"))

# --- 1. Statistics ####
df$Alkohol_avhold <- ifelse(df$Alkohol_kode == "Avhold", 0, 1)

depletion_vars <- c("cutoff_5_depletion_anySession","cutoff_5_depletion_Session2")
binary_predictors <- c("Snus","Kjønn","Røyk","Alkohol_avhold","CRP_below_1")

df_patients <- df %>%
  group_by(PasientID) %>%
  summarise(across(all_of(c(depletion_vars,binary_predictors,continuous_vars)), ~ first(.)), .groups = "drop") %>%
  mutate(across(all_of(binary_predictors), as.factor),
         across(all_of(continuous_vars), as.numeric))

# --- Chi-square tests with Cramér's V + CI ---
run_chi_sq <- function(dep_var, pred_var) {
  tbl <- table(df_patients[[dep_var]], df_patients[[pred_var]])
  if (all(dim(tbl) > 1)) {
    test <- chisq.test(tbl)
    V <- cramerV(tbl, bias.correct = TRUE)
    boot_fun <- function(data, i) {
      t <- table(data[i, 1], data[i, 2])
      cramerV(t, bias.correct = TRUE)
    }
    boot_res <- boot(data = as.data.frame(tbl), statistic = boot_fun, R = 1000)
    ci <- boot.ci(boot_res, type = "perc")$percent[4:5]
    data.frame(depletion_var = dep_var, predictor = pred_var,
               statistic = unname(test$statistic), df = unname(test$parameter),
               p_value = test$p.value, cramerV = V,
               ci_lower = ci[1], ci_upper = ci[2])
  } else {
    data.frame(depletion_var = dep_var, predictor = pred_var,
               statistic = NA, df = NA, p_value = NA,
               cramerV = NA, ci_lower = NA, ci_upper = NA)
  }
}

chi_sq_results <- do.call(rbind, lapply(depletion_vars, function(dv) {
  do.call(rbind, lapply(binary_predictors, function(bp) run_chi_sq(dv, bp)))
}))
chi_sq_results$p_FDR <- p.adjust(chi_sq_results$p_value, method = "fdr")

# --- Wilcoxon tests with rank-biserial r + CI + Cohen's d ---
wilcox_row <- function(dep, cont) {
  g <- df_patients[[dep]]
  x <- df_patients[[cont]]
  ok <- !is.na(g) & !is.na(x)
  g <- factor(g[ok]); x <- x[ok]
  
  # Skip invalid cases
  if (length(levels(g)) != 2 || any(table(g) == 0)) return(NULL)
  
  lv <- levels(g)
  wt <- wilcox.test(x ~ g, exact = FALSE)
  
  # Safe effect size
  eff <- tryCatch({
    wilcox_effsize(data.frame(x, g), x ~ g, ci = TRUE)
  }, error = function(e) {
    data.frame(effsize = NA, conf.low = NA, conf.high = NA)
  })
  
  # Safe Cohen's d
  d <- tryCatch({
    cohen.d(x ~ g)$estimate
  }, error = function(e) NA)
  
  data.frame(
    comparison        = paste0(cont, " ~ ", dep),
    depletion_var     = dep,
    continuous_var    = cont,
    group0_level      = lv[1],
    group1_level      = lv[2],
    n_group0          = sum(g == lv[1]),
    n_group1          = sum(g == lv[2]),
    median_IQR_group0 = sprintf("%.2f [%.2f]", median(x[g == lv[1]]), IQR(x[g == lv[1]])),
    median_IQR_group1 = sprintf("%.2f [%.2f]", median(x[g == lv[2]]), IQR(x[g == lv[2]])),
    test_used         = "Wilcoxon rank-sum",
    p_value           = wt$p.value,
    rank_biserial_r   = eff$effsize,
    ci_lower          = eff$conf.low,
    ci_upper          = eff$conf.high,
    cohens_d          = d
  )
}

res_list <- lapply(depletion_vars, function(dep) {
  lapply(continuous_vars, function(cont) wilcox_row(dep, cont))
})
wilcox_results <- bind_rows(unlist(res_list, recursive = FALSE)) %>%
  mutate(p_FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(depletion_var, continuous_var)

# --- Export all results ---
write.xlsx(list(
  Descriptives = final_summary,
  ChiSquare = chi_sq_results,
  Wilcoxon = wilcox_results
), paste0(path,"Stats_results.xlsx"))