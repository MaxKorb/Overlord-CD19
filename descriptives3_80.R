# =====================================================================
# MS B-cell Depletion — 80% Sensitivity Definition
# Max Korbmacher - 2026
# =====================================================================

# --- CLEANUP ---------------------------------------------------------
rm(list = ls(all.names = TRUE))
gc()

path <- "/Users/max/Documents/Local/MS/Student_Projects/"

# --- PACKAGES --------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, readxl, openxlsx,
  effectsize, rstatix, effsize, rcompanion
)

# =====================================================================
# IMPORT + CLEANING
# =====================================================================

df <- read_excel(paste0(path, "2025-08-22-table.xlsx"))

df$CD19_celler_per_µL <- ifelse(df$CD19_celler_per_µL == "", NA, df$CD19_celler_per_µL)

# Remove excluded patients
df = df[!df$PasientID %in% c("004", "012", "030", "048"),]

# Convert variables
df$Røyk <- as.numeric(df$Røyk)
df$`Alder_(år)` <- as.numeric(df$`Alder_(år)`)

# Remove patients with completely missing CD19
fully_missing <- df %>%
  group_by(PasientID) %>%
  summarise(all_missing = all(is.na(CD19_celler_per_µL)))

df <- df[!df$PasientID %in% fully_missing$PasientID[fully_missing$all_missing], ]

# Clean CD19 numeric
df <- df %>%
  mutate(CD19_celler_per_µL =
           as.numeric(gsub("[^0-9\\.]", "", gsub(",", ".", CD19_celler_per_µL))))

# =====================================================================
# 80% DEPLETION DEFINITION (EXCLUDING FIRST VISIT)
# =====================================================================

df <- df %>%
  group_by(PasientID) %>%
  mutate(
    is_first_visit = `Måling nr.` == min(`Måling nr.`, na.rm = TRUE),
    
    drop5_binary = if_else(
      is_first_visit,
      NA_integer_,
      if_else(!is.na(CD19_celler_per_µL) &
                CD19_celler_per_µL <= 5, 1L, 0L)
    ),
    
    n_visits   = sum(!is.na(CD19_celler_per_µL) & !is_first_visit),
    n_depleted = sum(drop5_binary == 1L, na.rm = TRUE),
    
    cutoff_5_depletion_80pct = if_else(
      n_visits > 0 & (n_depleted / n_visits) >= 0.8, 1L, 0L
    )
  ) %>%
  ungroup()

# Export processed dataset
write.xlsx(df, paste0(path, "Table_with_80pct_depletion.xlsx"))

# =====================================================================
# VARIABLE SETUP
# =====================================================================

binary_vars <- c(
  "Snus","Kjønn","Røyk",
  "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
  "Alkohol_kode","CRP_below_1"
)

continuous_vars <- c(
  "Alder_(år)","CD19_celler_per_µL","CD19_%","BMI",
  "Total_antall_attakker",
  "Tid_siden_første_kliniske_hendelse_(måneder)",
  "Tid_siden_MS_diagnose_(måneder)",
  "Antall_attakker_siste_året",
  "Minimum_EDSS_score","IgG","IgM","Antall_oligoklonale_bånd"
)

# Clean continuous
df <- df %>%
  mutate(across(all_of(continuous_vars),
                ~ as.numeric(gsub("[^0-9\\.]", "", gsub(",", ".", .)))))

# Alcohol → binary
df$Alkohol_avhold <- ifelse(df$Alkohol_kode == "Avhold", 0, 1)

binary_predictors <- c("Snus","Kjønn","Røyk","Alkohol_avhold","CRP_below_1")

# =====================================================================
# COLLAPSE TO PATIENT LEVEL
# =====================================================================

df_patients <- df %>%
  group_by(PasientID) %>%
  summarise(across(
    all_of(c("cutoff_5_depletion_80pct", binary_predictors, continuous_vars)),
    ~ first(.)
  ),
  .groups = "drop"
  ) %>%
  mutate(
    across(all_of(binary_predictors), as.factor),
    across(all_of(continuous_vars), as.numeric)
  )

# =====================================================================
# STATISTICS
# =====================================================================

# --- CHI-SQUARE ------------------------------------------------------

run_chi_sq <- function(pred_var) {
  
  g <- df_patients$cutoff_5_depletion_80pct
  x <- df_patients[[pred_var]]
  
  ok <- !is.na(g) & !is.na(x)
  g <- g[ok]
  x <- x[ok]
  
  tbl <- table(g, x)
  
  N_total <- sum(tbl)
  n_group0 <- sum(g == 0)
  n_group1 <- sum(g == 1)
  
  if (all(dim(tbl) > 1)) {
    test <- chisq.test(tbl)
    V <- cramerV(tbl, bias.correct = TRUE)
    
    data.frame(
      predictor = pred_var,
      N_total = N_total,
      n_group0 = n_group0,   # ✅ added
      n_group1 = n_group1,   # ✅ added
      statistic = as.numeric(test$statistic),
      df = as.numeric(test$parameter),
      p_value = test$p.value,
      cramerV = V
    )
  } else {
    data.frame(
      predictor = pred_var,
      N_total = N_total,
      n_group0 = n_group0,
      n_group1 = n_group1,
      statistic = NA,
      df = NA,
      p_value = NA,
      cramerV = NA
    )
  }
}
chi_sq_results <- bind_rows(lapply(binary_predictors, run_chi_sq)) %>%
  mutate(p_FDR = p.adjust(p_value, method = "fdr"))

# ---------------------------------------------------------------------
# WILCOXON
# ---------------------------------------------------------------------

wilcox_row <- function(cont) {
  
  g <- df_patients$cutoff_5_depletion_80pct
  x <- df_patients[[cont]]
  
  ok <- !is.na(g) & !is.na(x)
  g <- factor(g[ok])
  x <- x[ok]
  
  N <- length(x)
  
  if (length(levels(g)) != 2 || any(table(g) == 0)) return(NULL)
  
  lv <- levels(g)
  n0 <- sum(g == lv[1])
  n1 <- sum(g == lv[2])
  
  med0 <- median(x[g == lv[1]])
  iqr0 <- IQR(x[g == lv[1]])
  med1 <- median(x[g == lv[2]])
  iqr1 <- IQR(x[g == lv[2]])
  
  wt <- wilcox.test(x ~ g, exact = FALSE)
  eff <- wilcox_effsize(data.frame(x, g), x ~ g, ci = TRUE)
  d <- tryCatch(cohen.d(x ~ g)$estimate, error = function(e) NA)
  
  data.frame(
    continuous_var = cont,
    N = N,
    n_group0 = n0,
    n_group1 = n1,
    median_IQR_group0 = sprintf("%.2f [%.2f]", med0, iqr0),
    median_IQR_group1 = sprintf("%.2f [%.2f]", med1, iqr1),
    p_value = wt$p.value,
    rank_biserial_r = eff$effsize,
    ci_low = eff$conf.low,
    ci_high = eff$conf.high,
    cohens_d = d
  )
}

wilcox_results <- bind_rows(lapply(continuous_vars, wilcox_row)) %>%
  mutate(p_FDR = p.adjust(p_value, method = "fdr"))

# =====================================================================
# EXPORT
# =====================================================================

write.xlsx(list(
  ChiSquare_80pct = chi_sq_results,
  Wilcoxon_80pct = wilcox_results
), paste0(path, "Stats_results_only_80pct.xlsx"))

# =====================================================================
# END
# =====================================================================