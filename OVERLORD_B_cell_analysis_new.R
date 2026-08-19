# ============================================================
# B-cell depletion analysis
# OVERLORD-MS laboratory sub-study
# ============================================================

# This script contains the R code used for data processing,
# descriptive analyses, statistical analyses, and figures.

# ============================================================
# 1. PACKAGES
# ============================================================

library(dplyr)
library(openxlsx)

# ============================================================
# 2. IMPORT DATA
# ============================================================

# Local path to the analysis dataset
data_path <- "~/Documents/data_R/OVERLORD_B_cell_analysis_data.csv"

ny_tabell_raw <- read.csv2(
  data_path,
  header = TRUE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)


# ============================================================
# 3. DEFINE ANALYSIS COHORT
# ============================================================

# Inclusion criterion:
# At least two available absolute CD19 measurements after baseline

followup_count <- tapply(
  !is.na(ny_tabell_raw$`CD19_celler_per_µL`) &
    ny_tabell_raw$`Måling nr.` > 1,
  ny_tabell_raw$PasientID,
  sum
)

followup_count <- data.frame(
  PasientID = as.numeric(names(followup_count)),
  N_followup_CD19 = as.numeric(followup_count)
)

# Identify participants who meet the inclusion criterion
included_ids <- followup_count$PasientID[
  followup_count$N_followup_CD19 >= 2
]

# Identify excluded participants
excluded_ids <- followup_count$PasientID[
  followup_count$N_followup_CD19 < 2
]

# Create analysis dataset
ny_tabell <- ny_tabell_raw[
  ny_tabell_raw$PasientID %in% included_ids,
]

# Display excluded participant IDs
excluded_ids


# ============================================================
# 4. CHECK ANALYSIS COHORT AND CREATE PATIENT-LEVEL DATASET
# ============================================================

# Check dimensions of raw and analysis datasets
dim(ny_tabell_raw)
dim(ny_tabell)

# Check number of participants before and after exclusion
length(unique(ny_tabell_raw$PasientID))
length(unique(ny_tabell$PasientID))

# Add readable treatment labels to the longitudinal dataset
ny_tabell$Behandling <- ifelse(
  ny_tabell$Treatment == 1,
  "Rituximab",
  ifelse(
    ny_tabell$Treatment == 2,
    "Ocrelizumab",
    NA
  )
)

# Keep one row per participant
pasienter <- ny_tabell[
  !duplicated(ny_tabell$PasientID),
]

# Check number of participants in the analysis cohort
nrow(pasienter)

# Check treatment distribution
table(pasienter$Behandling)

# ============================================================
# 5. DEFINE SUSTAINED DEPLETION
# ============================================================

# Exclude baseline (measurement 1)
oppfolging <- ny_tabell[
  ny_tabell$`Måling nr.` != 1,
]

# Create one sustained depletion status per participant
deplesjon_alle_tidspunkt <- data.frame(
  PasientID = unique(oppfolging$PasientID)
)

deplesjon_alle_tidspunkt$gruppe <- sapply(
  deplesjon_alle_tidspunkt$PasientID,
  function(id) {
    
    x <- oppfolging$`CD19_celler_per_µL`[
      oppfolging$PasientID == id
    ]
    
    if (any(x > 5, na.rm = TRUE)) {
      "Ikke-deplesjon"
    } else {
      "Deplesjon"
    }
  }
)

# Check sustained depletion groups
table(deplesjon_alle_tidspunkt$gruppe)

# ============================================================
# 6. DEFINE EARLY DEPLETION AT MONTH 3
# ============================================================

# Select first follow-up visit
# Month 3 corresponds to measurement 2
forste_oppfolging <- ny_tabell[
  ny_tabell$`Måling nr.` == 2,
]

# Classify early depletion
# Depletion is defined as CD19 <= 5 cells/µL
forste_oppfolging$gruppe <- ifelse(
  is.na(forste_oppfolging$`CD19_celler_per_µL`),
  NA,
  ifelse(
    forste_oppfolging$`CD19_celler_per_µL` <= 5,
    "Deplesjon",
    "Ikke-deplesjon"
  )
)

# Check early depletion groups
table(
  forste_oppfolging$gruppe,
  useNA = "ifany"
)

# ============================================================
# 7. DEFINE 80% SENSITIVITY ANALYSIS
# ============================================================

# Use follow-up measurements only
# Baseline (measurement 1) is excluded
oppfolging_80 <- ny_tabell[
  ny_tabell$`Måling nr.` != 1,
]

# Calculate the proportion of available follow-up measurements
# with CD19 <= 5 cells/µL for each participant
sensitivitet80 <- aggregate(
  `CD19_celler_per_µL` ~ PasientID,
  data = oppfolging_80,
  FUN = function(x) {
    mean(x[!is.na(x)] <= 5)
  }
)

# Rename the calculated proportion
names(sensitivitet80)[2] <- "andel_depletert"

# Classify participants according to the 80% definition
sensitivitet80$Deplesjon80 <- ifelse(
  sensitivitet80$andel_depletert >= 0.80,
  "Depletert",
  "Ikke-depletert"
)

# Check number of participants and depletion groups
nrow(sensitivitet80)

table(sensitivitet80$Deplesjon80)

# ============================================================
# 8. CREATE PATIENT-LEVEL DATASET FOR 80% SENSITIVITY ANALYSIS
# ============================================================

# Merge the 80% depletion classification
# with the patient-level dataset
sens80 <- merge(
  pasienter,
  sensitivitet80,
  by = "PasientID"
)

# Check number of participants
nrow(sens80)

# Check depletion groups
table(sens80$Deplesjon80)

# Check whether any participants are missing
setdiff(
  pasienter$PasientID,
  sens80$PasientID
)

# Check treatment distribution
table(sens80$Behandling)


# ============================================================
# 9. CREATE ANALYSIS VARIABLES
# ============================================================

# Binary alcohol variable:
# 0 = abstinent
# 1 = current alcohol use

pasienter$Alkohol_bruk <- ifelse(
  pasienter$Alkohol == "Avhold",
  0,
  1
)

sens80$Alkohol_bruk <- ifelse(
  sens80$Alkohol == "Avhold",
  0,
  1
)

# Check coding
table(
  pasienter$Alkohol_bruk,
  useNA = "ifany"
)

table(
  sens80$Alkohol_bruk,
  useNA = "ifany"
)

# ============================================================
# 10. SEX AND SMOKING DISTRIBUTION
#     IN THE 80% SENSITIVITY ANALYSIS
# ============================================================

# Sex distribution within depletion groups
prop.table(
  table(
    sens80$Kjønn,
    sens80$Deplesjon80
  ),
  margin = 2
) * 100

# Current smoking distribution within depletion groups
prop.table(
  table(
    sens80$Røyk,
    sens80$Deplesjon80
  ),
  margin = 2
) * 100

# ============================================================
# 11. WILCOXON RANK-SUM TESTS
#     80% SENSITIVITY ANALYSIS
# ============================================================

# Continuous baseline variables included in predictor analyses
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

# Function performing a Wilcoxon rank-sum test
# for one continuous variable
run_wilcoxon <- function(var) {
  
  test_data <- data.frame(
    value = sens80[[var]],
    depletion = sens80$Deplesjon80
  )
  
  # Keep complete observations
  test_data <- test_data[
    complete.cases(test_data),
  ]
  
  test <- wilcox.test(
    value ~ depletion,
    data = test_data,
    exact = FALSE
  )
  
  data.frame(
    Variable = var,
    N = nrow(test_data),
    n_non_depleted = sum(
      test_data$depletion == "Ikke-depletert"
    ),
    n_depleted = sum(
      test_data$depletion == "Depletert"
    ),
    p_value = test$p.value
  )
}

# Run the test for all continuous variables
wilcoxon_results <- do.call(
  rbind,
  lapply(
    continuous_vars,
    run_wilcoxon
  )
)

# Display results
wilcoxon_results

# ============================================================
# 12. CHI-SQUARE TESTS
#     80% SENSITIVITY ANALYSIS
# ============================================================

# Function performing Pearson's chi-square test
# with Yates' continuity correction
run_chisq <- function(var) {
  
  tab <- table(
    sens80[[var]],
    sens80$Deplesjon80
  )
  
  test <- chisq.test(tab)
  
  data.frame(
    Variable = var,
    N = sum(tab),
    n_non_depleted = sum(
      tab[, "Ikke-depletert"]
    ),
    n_depleted = sum(
      tab[, "Depletert"]
    ),
    statistic = unname(
      test$statistic
    ),
    df = unname(
      test$parameter
    ),
    p_value = test$p.value
  )
}

# Categorical baseline variables
categorical_vars <- c(
  "Kjønn",
  "Røyk",
  "Snus",
  "CRP_below_1"
)

# Run tests for all categorical variables
chisq_results <- do.call(
  rbind,
  lapply(
    categorical_vars,
    run_chisq
  )
)

# Display results
chisq_results

# ============================================================
# 13. MULTIPLE-TESTING CORRECTION
#     80% SENSITIVITY ANALYSIS
# ============================================================

# ------------------------------------------------------------
# FDR correction
# ------------------------------------------------------------

# Benjamini-Hochberg FDR correction is applied separately
# to the continuous and categorical predictor analyses

wilcoxon_results$p_FDR <- p.adjust(
  wilcoxon_results$p_value,
  method = "BH"
)

chisq_results$p_FDR <- p.adjust(
  chisq_results$p_value,
  method = "BH"
)


# ------------------------------------------------------------
# Bonferroni correction
# ------------------------------------------------------------

# Bonferroni correction is applied across all
# 16 predictor tests included in the sensitivity analysis

all_p_values <- c(
  setNames(
    wilcoxon_results$p_value,
    wilcoxon_results$Variable
  ),
  setNames(
    chisq_results$p_value,
    chisq_results$Variable
  )
)

p_Bonferroni <- p.adjust(
  all_p_values,
  method = "bonferroni"
)


# ------------------------------------------------------------
# Create combined results table
# ------------------------------------------------------------

multiple_testing_results <- rbind(
  
  data.frame(
    Variable = wilcoxon_results$Variable,
    p_value = wilcoxon_results$p_value,
    p_FDR = wilcoxon_results$p_FDR
  ),
  
  data.frame(
    Variable = chisq_results$Variable,
    p_value = chisq_results$p_value,
    p_FDR = chisq_results$p_FDR
  )
)

# Match Bonferroni-adjusted p-values
# to the corresponding variables
multiple_testing_results$p_Bonferroni <-
  p_Bonferroni[
    match(
      multiple_testing_results$Variable,
      names(p_Bonferroni)
    )
  ]

rownames(
  multiple_testing_results
) <- NULL

# Display results
multiple_testing_results

# ============================================================
# 14. CREATE PATIENT-LEVEL DATASET
#     FOR SUSTAINED DEPLETION
# ============================================================

# Rename sustained depletion variable for clarity
deplesjon_status <- deplesjon_alle_tidspunkt

names(deplesjon_status)[
  names(deplesjon_status) == "gruppe"
] <- "Deplesjon_status"

# Merge sustained depletion status
# with the patient-level dataset
pasienter_sustained <- merge(
  pasienter,
  deplesjon_status,
  by = "PasientID"
)

# Create combined treatment/depletion group
pasienter_sustained$Gruppe <- paste(
  pasienter_sustained$Behandling,
  pasienter_sustained$Deplesjon_status,
  sep = " - "
)

# Binary alcohol variable:
# 0 = abstinent
# 1 = current alcohol use
pasienter_sustained$Alkohol_bruk <- ifelse(
  pasienter_sustained$Alkohol == "Avhold",
  0,
  1
)

# Check group sizes
table(pasienter_sustained$Gruppe)


# ============================================================
# 15. HELPER FUNCTIONS FOR SUPPLEMENTARY TABLE 1
# ============================================================

# Continuous variables: mean (SD)
mean_sd_sustained <- function(var) {
  tapply(
    pasienter_sustained[[var]],
    pasienter_sustained$Gruppe,
    function(x) {
      paste0(
        round(mean(x, na.rm = TRUE), 2),
        " (",
        round(sd(x, na.rm = TRUE), 2),
        ")"
      )
    }
  )
}

# Categorical variables: n (%)
n_pct_sustained <- function(var, value) {
  tapply(
    pasienter_sustained[[var]] == value,
    pasienter_sustained$Gruppe,
    function(x) {
      paste0(
        sum(x, na.rm = TRUE),
        " (",
        round(mean(x, na.rm = TRUE) * 100, 2),
        "%)"
      )
    }
  )
}

# ============================================================
# 16. CREATE SUPPLEMENTARY TABLE 1
# ============================================================

supp_table1 <- rbind(
  "Number of participants" = as.character(
    table(pasienter_sustained$Gruppe)
  ),
  
  "Age, years" =
    mean_sd_sustained("Alder_(år)"),
  
  "Sex, female" =
    n_pct_sustained("Kjønn", "Female"),
  
  "Sex, male" =
    n_pct_sustained("Kjønn", "Male"),
  
  "Body Mass Index, kg/m²" =
    mean_sd_sustained("BMI"),
  
  "Current use of smoked tobacco" =
    n_pct_sustained("Røyk", 1),
  
  "Current use of smokeless tobacco (e.g. Snuff)" =
    n_pct_sustained("Snus", "Ja"),
  
  "Current use of alcohol" =
    n_pct_sustained("Alkohol_bruk", 1),
  
  "Time since MS diagnosis, months" =
    mean_sd_sustained("Tid_siden_MS_diagnose_(måneder)"),
  
  "Time since first neurological symptom, months" =
    mean_sd_sustained(
      "Tid_siden_første_kliniske_hendelse_(måneder)"
    ),
  
  "Minimum EDSS score" =
    mean_sd_sustained("Minimum_EDSS_score"),
  
  "Total relapses since onset" =
    mean_sd_sustained("Total_antall_attakker"),
  
  "Relapses in the past 12 months" =
    mean_sd_sustained("Antall_attakker_siste_året"),
  
  "CSF oligoclonal bands tested" =
    n_pct_sustained(
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja"
    ),
  
  "CSF oligoclonal bands" =
    mean_sd_sustained("Antall_oligoklonale_bånd"),
  
  "IgG, g/L" =
    mean_sd_sustained("IgG"),
  
  "IgM, g/L" =
    mean_sd_sustained("IgM"),
  
  "CD19, % of lymphocytes" =
    mean_sd_sustained("CD19_%"),
  
  "CD19, cells/µL" =
    mean_sd_sustained("CD19_celler_per_µL"),
  
  "CRP < 1 mg/L" =
    n_pct_sustained("CRP_below_1", 1)
)

# Convert to data frame
supp_table1 <- as.data.frame(
  supp_table1
)

# Add variable names as first column
supp_table1 <- cbind(
  Variable = rownames(supp_table1),
  supp_table1
)

rownames(supp_table1) <- NULL

# Display table
supp_table1

# ============================================================
# 17. EXPORT SUPPLEMENTARY TABLE 1
# ============================================================

write.xlsx(
  supp_table1,
  file = "~/Documents/data_R/Supplementary_Table_1.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)

# ============================================================
# 18. CREATE PATIENT-LEVEL DATASET
#     FOR EARLY DEPLETION
# ============================================================

# Create early depletion status from month 3
early_deplesjon_status <- data.frame(
  PasientID = forste_oppfolging$PasientID,
  Deplesjon_status = ifelse(
    is.na(forste_oppfolging$`CD19_celler_per_µL`),
    NA,
    ifelse(
      forste_oppfolging$`CD19_celler_per_µL` <= 5,
      "Depletert",
      "Ikke-depletert"
    )
  )
)

# Exclude participants without an available CD19
# measurement at month 3
early_deplesjon_status <- early_deplesjon_status[
  !is.na(early_deplesjon_status$Deplesjon_status),
]

# Merge early depletion status
# with the patient-level dataset
pasienter_early <- merge(
  pasienter,
  early_deplesjon_status,
  by = "PasientID"
)

# Create combined treatment/depletion group
pasienter_early$Gruppe <- paste(
  pasienter_early$Behandling,
  pasienter_early$Deplesjon_status,
  sep = " - "
)

# Binary alcohol variable:
# 0 = abstinent
# 1 = current alcohol use
pasienter_early$Alkohol_bruk <- ifelse(
  pasienter_early$Alkohol == "Avhold",
  0,
  1
)

# Check number of participants included
nrow(pasienter_early)

# Check early depletion groups
table(pasienter_early$Gruppe)

# ============================================================
# 19. HELPER FUNCTIONS FOR SUPPLEMENTARY TABLE 2
# ============================================================

# Continuous variables: mean (SD)
mean_sd_early <- function(var) {
  tapply(
    pasienter_early[[var]],
    pasienter_early$Gruppe,
    function(x) {
      paste0(
        round(mean(x, na.rm = TRUE), 2),
        " (",
        round(sd(x, na.rm = TRUE), 2),
        ")"
      )
    }
  )
}

# Categorical variables: n (%)
n_pct_early <- function(var, value) {
  tapply(
    pasienter_early[[var]] == value,
    pasienter_early$Gruppe,
    function(x) {
      paste0(
        sum(x, na.rm = TRUE),
        " (",
        round(mean(x, na.rm = TRUE) * 100, 2),
        "%)"
      )
    }
  )
}

# ============================================================
# 20. CREATE SUPPLEMENTARY TABLE 2
# ============================================================

supp_table2 <- rbind(
  "Number of participants" = as.character(
    table(pasienter_early$Gruppe)
  ),
  
  "Age, years" =
    mean_sd_early("Alder_(år)"),
  
  "Sex, female" =
    n_pct_early("Kjønn", "Female"),
  
  "Sex, male" =
    n_pct_early("Kjønn", "Male"),
  
  "Body Mass Index, kg/m²" =
    mean_sd_early("BMI"),
  
  "Current use of smoked tobacco" =
    n_pct_early("Røyk", 1),
  
  "Current use of smokeless tobacco (e.g. Snuff)" =
    n_pct_early("Snus", "Ja"),
  
  "Current use of alcohol" =
    n_pct_early("Alkohol_bruk", 1),
  
  "Time since MS diagnosis, months" =
    mean_sd_early("Tid_siden_MS_diagnose_(måneder)"),
  
  "Time since first neurological symptom, months" =
    mean_sd_early(
      "Tid_siden_første_kliniske_hendelse_(måneder)"
    ),
  
  "Minimum EDSS score" =
    mean_sd_early("Minimum_EDSS_score"),
  
  "Total relapses since onset" =
    mean_sd_early("Total_antall_attakker"),
  
  "Relapses in the past 12 months" =
    mean_sd_early("Antall_attakker_siste_året"),
  
  "CSF oligoclonal bands tested" =
    n_pct_early(
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja"
    ),
  
  "CSF oligoclonal bands" =
    mean_sd_early("Antall_oligoklonale_bånd"),
  
  "IgG, g/L" =
    mean_sd_early("IgG"),
  
  "IgM, g/L" =
    mean_sd_early("IgM"),
  
  "CD19, % of lymphocytes" =
    mean_sd_early("CD19_%"),
  
  "CD19, cells/µL" =
    mean_sd_early("CD19_celler_per_µL"),
  
  "CRP < 1 mg/L" =
    n_pct_early("CRP_below_1", 1)
)

# Convert to data frame
supp_table2 <- as.data.frame(
  supp_table2
)

# Add variable names as first column
supp_table2 <- cbind(
  Variable = rownames(supp_table2),
  supp_table2
)

rownames(supp_table2) <- NULL

# Display table
supp_table2


# ============================================================
# 21. EXPORT SUPPLEMENTARY TABLE 2
# ============================================================

write.xlsx(
  supp_table2,
  file = "~/Documents/data_R/Supplementary_Table_2.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)

# ============================================================
# 22. CD19 MEASUREMENT AVAILABILITY BY VISIT
# ============================================================

# Percentage of participants with an available absolute CD19
# measurement at each scheduled visit

cd19_tilgjengelighet <- tapply(
  !is.na(ny_tabell$`CD19_celler_per_µL`),
  ny_tabell$`Måling nr.`,
  mean
) * 100

# Add visit labels
names(cd19_tilgjengelighet) <- c(
  "Baseline",
  "3 months",
  "6 months",
  "12 months",
  "18 months",
  "24 months",
  "30 months",
  "36 months"
)

# Display availability by visit (%)
round(
  cd19_tilgjengelighet,
  1
)

# ============================================================
# 23. OVERALL CD19 AVAILABILITY AND FOLLOW-UP DURATION
# ============================================================

# ------------------------------------------------------------
# Overall availability of absolute CD19 measurements
# during follow-up (months 3–36)
# ------------------------------------------------------------

overall_cd19_availability <- mean(
  !is.na(oppfolging$`CD19_celler_per_µL`)
) * 100

# Display overall availability (%)
round(
  overall_cd19_availability,
  1
)


# ------------------------------------------------------------
# Follow-up duration based on last available CD19 measurement
# ------------------------------------------------------------

# Keep observations with an available absolute CD19 measurement
cd19_malt <- ny_tabell[
  !is.na(ny_tabell$`CD19_celler_per_µL`),
]

# Identify the last visit with an available CD19 measurement
# for each participant
siste_cd19 <- aggregate(
  `Måling nr.` ~ PasientID,
  data = cd19_malt,
  FUN = max
)

# Convert measurement number to follow-up time in months
followup_months <- c(
  0,
  3,
  6,
  12,
  18,
  24,
  30,
  36
)

siste_cd19$followup_mnd <- followup_months[
  siste_cd19$`Måling nr.`
]

# Mean follow-up duration
mean_followup <- mean(
  siste_cd19$followup_mnd
)

# Standard deviation of follow-up duration
sd_followup <- sd(
  siste_cd19$followup_mnd
)

# Display results
round(
  mean_followup,
  2
)

round(
  sd_followup,
  2
)

# ============================================================
# 24. DEPLETION BY FOLLOW-UP VISIT
# ============================================================

# Calculate number with an available CD19 measurement,
# number depleted (CD19 <= 5 cells/µL),
# and percentage depleted at each follow-up visit

depletion_by_visit <- do.call(
  rbind,
  lapply(
    sort(unique(oppfolging$`Måling nr.`)),
    function(visit) {
      
      x <- oppfolging[
        oppfolging$`Måling nr.` == visit &
          !is.na(oppfolging$`CD19_celler_per_µL`),
      ]
      
      n_measured <- nrow(x)
      
      n_depleted <- sum(
        x$`CD19_celler_per_µL` <= 5
      )
      
      data.frame(
        Measurement = visit,
        N_measured = n_measured,
        N_depleted = n_depleted,
        Percent_depleted = 100 *
          n_depleted /
          n_measured
      )
    }
  )
)

# Add follow-up month labels
depletion_by_visit$Month <- c(
  3,
  6,
  12,
  18,
  24,
  30,
  36
)

# Arrange columns
depletion_by_visit <- depletion_by_visit[
  ,
  c(
    "Month",
    "N_measured",
    "N_depleted",
    "Percent_depleted"
  )
]

# Round percentage to one decimal
depletion_by_visit$Percent_depleted <- round(
  depletion_by_visit$Percent_depleted,
  1
)

# Display results
depletion_by_visit

# ============================================================
# 25. B-CELL REPOPULATION AFTER EARLY DEPLETION
# ============================================================

# Identify participants with early depletion at month 3,
# defined as CD19 <= 5 cells/µL

mnd3 <- ny_tabell[
  ny_tabell$`Måling nr.` == 2,
]

early_depletion_id <- mnd3$PasientID[
  !is.na(mnd3$`CD19_celler_per_µL`) &
    mnd3$`CD19_celler_per_µL` <= 5
]

# Keep measurements after month 3 (months 6–36)
# among participants with early depletion

senere <- ny_tabell[
  ny_tabell$PasientID %in% early_depletion_id &
    ny_tabell$`Måling nr.` > 2,
]

# Reappearance of B cells is defined as at least one
# subsequent CD19 measurement > 5 cells/µL

repopulation_id <- unique(
  senere$PasientID[
    !is.na(senere$`CD19_celler_per_µL`) &
      senere$`CD19_celler_per_µL` > 5
  ]
)

# Number of participants with early depletion
n_early_depleted <- length(
  early_depletion_id
)

# Number with subsequent B-cell reappearance
n_repopulation <- length(
  repopulation_id
)

# Percentage with subsequent B-cell reappearance
pct_repopulation <- 100 *
  n_repopulation /
  n_early_depleted

# Display results
n_early_depleted

n_repopulation

round(
  pct_repopulation,
  1
)

# ============================================================
# 26. NUMBER OF LATER VISITS WITH CD19 > 5
#     AMONG PARTICIPANTS WITH EARLY DEPLETION
# ============================================================

# Count the number of later visits (months 6–36)
# with CD19 > 5 cells/µL for each participant

non_depleted_visits <- aggregate(
  `CD19_celler_per_µL` ~ PasientID,
  data = senere,
  FUN = function(x) {
    sum(
      x > 5,
      na.rm = TRUE
    )
  }
)

names(non_depleted_visits)[2] <-
  "N_visits_above_5"

# Ensure all participants with early depletion are represented
all_early_depleted <- data.frame(
  PasientID = early_depletion_id
)

non_depleted_visits <- merge(
  all_early_depleted,
  non_depleted_visits,
  by = "PasientID",
  all.x = TRUE
)

# Participants without later values > 5 receive a count of 0
non_depleted_visits$N_visits_above_5[
  is.na(non_depleted_visits$N_visits_above_5)
] <- 0

# Distribution of number of later visits > 5
table(
  non_depleted_visits$N_visits_above_5
)

# Number with at least one later value > 5
sum(
  non_depleted_visits$N_visits_above_5 >= 1
)

# Number with at least two later values > 5
sum(
  non_depleted_visits$N_visits_above_5 >= 2
)

# Participants with at least one later CD19 measurement > 5
with_repopulation <- non_depleted_visits[
  non_depleted_visits$N_visits_above_5 >= 1,
]

# Mean, SD and range of number of visits > 5
mean(
  with_repopulation$N_visits_above_5
)

sd(
  with_repopulation$N_visits_above_5
)

range(
  with_repopulation$N_visits_above_5
)

# ============================================================
# 27. LAST AVAILABLE CD19 > 5 AFTER EARLY DEPLETION
# ============================================================

# Keep available measurements after month 3
# among participants with early depletion

later_available <- senere[
  !is.na(senere$`CD19_celler_per_µL`),
]

# Order observations by participant and measurement number
later_available <- later_available[
  order(
    later_available$PasientID,
    later_available$`Måling nr.`
  ),
]

# Extract the last available CD19 measurement
# for each participant
last_later_measurement <- later_available[
  !duplicated(
    later_available$PasientID,
    fromLast = TRUE
  ),
]

# Number of participants represented
nrow(
  last_later_measurement
)

# Number whose last available CD19 measurement
# was > 5 cells/µL
n_last_above_5 <- sum(
  last_later_measurement$`CD19_celler_per_µL` > 5
)

# Percentage of participants with early depletion
pct_last_above_5 <- 100 *
  n_last_above_5 /
  length(early_depletion_id)

# Display results
n_last_above_5

round(
  pct_last_above_5,
  1
)

# ============================================================
# 28. ALL LATER CD19 MEASUREMENTS > 5 AFTER EARLY DEPLETION
# ============================================================

# Determine whether all available later CD19 measurements
# were > 5 cells/µL for each participant

all_later_above_5 <- aggregate(
  `CD19_celler_per_µL` ~ PasientID,
  data = senere,
  FUN = function(x) {
    
    x <- x[
      !is.na(x)
    ]
    
    length(x) > 0 &&
      all(x > 5)
  }
)

names(all_later_above_5)[2] <-
  "All_later_CD19_above_5"

# Number of participants with all later measurements > 5
n_all_later_above_5 <- sum(
  all_later_above_5$All_later_CD19_above_5
)

# Percentage of participants with early depletion
pct_all_later_above_5 <- 100 *
  n_all_later_above_5 /
  length(early_depletion_id)

# Display results
n_all_later_above_5

round(
  pct_all_later_above_5,
  1
)

# ============================================================
# 29. MEDIAN AND IQR OF CD19 BY VISIT
# ============================================================

# Calculate median, Q1, Q3 and IQR for absolute CD19
# at each scheduled measurement

cd19_summary_by_visit <- do.call(
  rbind,
  lapply(
    sort(unique(ny_tabell$`Måling nr.`)),
    function(visit) {
      
      x <- ny_tabell$`CD19_celler_per_µL`[
        ny_tabell$`Måling nr.` == visit
      ]
      
      # Exclude missing CD19 measurements
      x <- x[
        !is.na(x)
      ]
      
      data.frame(
        Measurement = visit,
        N = length(x),
        Median = median(x),
        Q1 = unname(
          quantile(x, 0.25)
        ),
        Q3 = unname(
          quantile(x, 0.75)
        ),
        IQR = IQR(x)
      )
    }
  )
)

# Visit labels
visit_labels <- c(
  "Baseline",
  "3 months",
  "6 months",
  "12 months",
  "18 months",
  "24 months",
  "30 months",
  "36 months"
)

# Add visit labels
cd19_summary_by_visit$Visit <-
  visit_labels[
    cd19_summary_by_visit$Measurement
  ]

# Arrange columns
cd19_summary_by_visit <-
  cd19_summary_by_visit[
    ,
    c(
      "Visit",
      "N",
      "Median",
      "Q1",
      "Q3",
      "IQR"
    )
  ]

# Round numerical summaries to two decimals
cd19_summary_by_visit[
  ,
  c(
    "Median",
    "Q1",
    "Q3",
    "IQR"
  )
] <- round(
  cd19_summary_by_visit[
    ,
    c(
      "Median",
      "Q1",
      "Q3",
      "IQR"
    )
  ],
  2
)

# Remove automatic row names
rownames(
  cd19_summary_by_visit
) <- NULL

# Display results
cd19_summary_by_visit

# ============================================================
# 30. MEDIAN AND IQR OF CD19 BY VISIT AND TREATMENT
# ============================================================

# Calculate CD19 summary for each visit
# and treatment group

cd19_summary_by_treatment <- do.call(
  rbind,
  lapply(
    sort(unique(ny_tabell$`Måling nr.`)),
    function(visit) {
      
      do.call(
        rbind,
        lapply(
          c(
            "Rituximab",
            "Ocrelizumab"
          ),
          function(treatment) {
            
            x <- ny_tabell$`CD19_celler_per_µL`[
              ny_tabell$`Måling nr.` == visit &
                ny_tabell$Behandling == treatment
            ]
            
            # Exclude missing CD19 measurements
            x <- x[
              !is.na(x)
            ]
            
            data.frame(
              Measurement = visit,
              Treatment = treatment,
              N = length(x),
              Median = median(x),
              Q1 = unname(
                quantile(x, 0.25)
              ),
              Q3 = unname(
                quantile(x, 0.75)
              ),
              IQR = IQR(x)
            )
          }
        )
      )
    }
  )
)

# Add visit labels
cd19_summary_by_treatment$Visit <-
  visit_labels[
    cd19_summary_by_treatment$Measurement
  ]

# Arrange columns
cd19_summary_by_treatment <-
  cd19_summary_by_treatment[
    ,
    c(
      "Visit",
      "Treatment",
      "N",
      "Median",
      "Q1",
      "Q3",
      "IQR"
    )
  ]

# Round numerical summaries to two decimals
cd19_summary_by_treatment[
  ,
  c(
    "Median",
    "Q1",
    "Q3",
    "IQR"
  )
] <- round(
  cd19_summary_by_treatment[
    ,
    c(
      "Median",
      "Q1",
      "Q3",
      "IQR"
    )
  ],
  2
)

# Remove automatic row names
rownames(
  cd19_summary_by_treatment
) <- NULL

# Display results
cd19_summary_by_treatment

# ============================================================
# 31. NUMBER OF AVAILABLE CD19 MEASUREMENTS
#     BY VISIT AND TREATMENT
# ============================================================

# Count available absolute CD19 measurements
# for each treatment group at each scheduled visit

cd19_n_by_treatment <- with(
  subset(
    ny_tabell,
    !is.na(`CD19_celler_per_µL`)
  ),
  table(
    Behandling,
    `Måling nr.`
  )
)

# Add visit labels
colnames(
  cd19_n_by_treatment
) <- visit_labels

# Display results
cd19_n_by_treatment

# ============================================================
# 32. TABLE 1. BASELINE CHARACTERISTICS BY TREATMENT
# ============================================================

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

# Mean ± SD for continuous variables
mean_sd_table1 <- function(var, treatment = NULL) {
  
  if (is.null(treatment)) {
    
    x <- pasienter[[var]]
    
  } else {
    
    x <- pasienter[[var]][
      pasienter$Behandling == treatment
    ]
  }
  
  # Exclude missing observations
  x <- x[
    !is.na(x)
  ]
  
  sprintf(
    "%.1f ± %.1f",
    mean(x),
    sd(x)
  )
}


# n (%) for categorical variables
n_pct_table1 <- function(var, value, treatment = NULL) {
  
  if (is.null(treatment)) {
    
    x <- pasienter[[var]]
    
  } else {
    
    x <- pasienter[[var]][
      pasienter$Behandling == treatment
    ]
  }
  
  # Exclude missing observations from denominator
  x <- x[
    !is.na(x)
  ]
  
  n <- sum(
    x == value
  )
  
  pct <- 100 *
    n /
    length(x)
  
  sprintf(
    "%d (%.1f)",
    n,
    pct
  )
}


# ------------------------------------------------------------
# Number of participants
# ------------------------------------------------------------

N_total <- nrow(
  pasienter
)

N_RTX <- sum(
  pasienter$Behandling == "Rituximab"
)

N_OCR <- sum(
  pasienter$Behandling == "Ocrelizumab"
)


# ------------------------------------------------------------
# Create Table 1
# ------------------------------------------------------------

table1 <- data.frame(
  
  Variable = c(
    "Number of participants, n (%)",
    "Age, years, mean (±SD)",
    "Sex, female, n (%)",
    "Sex, male, n (%)",
    "Body Mass Index, kg/m², mean (±SD)",
    "Current use of smoked tobacco, n (%)",
    "Current use of smokeless tobacco (e.g. Snuff), n (%)",
    "Current use of alcohol, n (%)",
    "Time since MS diagnosis, months, mean (±SD)",
    "Time since first neurological symptom, months, mean (±SD)",
    "Minimum EDSS score, mean (±SD)",
    "Total relapses since onset, mean (±SD)",
    "Relapses in the past 12 months, mean (±SD)",
    "CSF oligoclonal bands tested, n (%)",
    "CSF oligoclonal bands, mean (±SD)",
    "IgG, g/L, mean (±SD)",
    "IgM, g/L, mean (±SD)",
    "CD19+ B cells, % of lymphocytes, mean (±SD)",
    "CD19+ B cells, cells/µL, mean (±SD)",
    "CRP < 1 mg/L, n (%)"
  ),
  
  Total_population = c(
    
    sprintf(
      "%d (%.1f)",
      N_total,
      100
    ),
    
    mean_sd_table1(
      "Alder_(år)"
    ),
    
    n_pct_table1(
      "Kjønn",
      "Female"
    ),
    
    n_pct_table1(
      "Kjønn",
      "Male"
    ),
    
    mean_sd_table1(
      "BMI"
    ),
    
    n_pct_table1(
      "Røyk",
      1
    ),
    
    n_pct_table1(
      "Snus",
      "Ja"
    ),
    
    n_pct_table1(
      "Alkohol_bruk",
      1
    ),
    
    mean_sd_table1(
      "Tid_siden_MS_diagnose_(måneder)"
    ),
    
    mean_sd_table1(
      "Tid_siden_første_kliniske_hendelse_(måneder)"
    ),
    
    mean_sd_table1(
      "Minimum_EDSS_score"
    ),
    
    mean_sd_table1(
      "Total_antall_attakker"
    ),
    
    mean_sd_table1(
      "Antall_attakker_siste_året"
    ),
    
    n_pct_table1(
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja"
    ),
    
    mean_sd_table1(
      "Antall_oligoklonale_bånd"
    ),
    
    mean_sd_table1(
      "IgG"
    ),
    
    mean_sd_table1(
      "IgM"
    ),
    
    mean_sd_table1(
      "CD19_%"
    ),
    
    mean_sd_table1(
      "CD19_celler_per_µL"
    ),
    
    n_pct_table1(
      "CRP_below_1",
      1
    )
  ),
  
  Rituximab = c(
    
    sprintf(
      "%d (%.1f)",
      N_RTX,
      100 *
        N_RTX /
        N_total
    ),
    
    mean_sd_table1(
      "Alder_(år)",
      "Rituximab"
    ),
    
    n_pct_table1(
      "Kjønn",
      "Female",
      "Rituximab"
    ),
    
    n_pct_table1(
      "Kjønn",
      "Male",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "BMI",
      "Rituximab"
    ),
    
    n_pct_table1(
      "Røyk",
      1,
      "Rituximab"
    ),
    
    n_pct_table1(
      "Snus",
      "Ja",
      "Rituximab"
    ),
    
    n_pct_table1(
      "Alkohol_bruk",
      1,
      "Rituximab"
    ),
    
    mean_sd_table1(
      "Tid_siden_MS_diagnose_(måneder)",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "Tid_siden_første_kliniske_hendelse_(måneder)",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "Minimum_EDSS_score",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "Total_antall_attakker",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "Antall_attakker_siste_året",
      "Rituximab"
    ),
    
    n_pct_table1(
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "Antall_oligoklonale_bånd",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "IgG",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "IgM",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "CD19_%",
      "Rituximab"
    ),
    
    mean_sd_table1(
      "CD19_celler_per_µL",
      "Rituximab"
    ),
    
    n_pct_table1(
      "CRP_below_1",
      1,
      "Rituximab"
    )
  ),
  
  Ocrelizumab = c(
    
    sprintf(
      "%d (%.1f)",
      N_OCR,
      100 *
        N_OCR /
        N_total
    ),
    
    mean_sd_table1(
      "Alder_(år)",
      "Ocrelizumab"
    ),
    
    n_pct_table1(
      "Kjønn",
      "Female",
      "Ocrelizumab"
    ),
    
    n_pct_table1(
      "Kjønn",
      "Male",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "BMI",
      "Ocrelizumab"
    ),
    
    n_pct_table1(
      "Røyk",
      1,
      "Ocrelizumab"
    ),
    
    n_pct_table1(
      "Snus",
      "Ja",
      "Ocrelizumab"
    ),
    
    n_pct_table1(
      "Alkohol_bruk",
      1,
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "Tid_siden_MS_diagnose_(måneder)",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "Tid_siden_første_kliniske_hendelse_(måneder)",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "Minimum_EDSS_score",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "Total_antall_attakker",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "Antall_attakker_siste_året",
      "Ocrelizumab"
    ),
    
    n_pct_table1(
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "Antall_oligoklonale_bånd",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "IgG",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "IgM",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "CD19_%",
      "Ocrelizumab"
    ),
    
    mean_sd_table1(
      "CD19_celler_per_µL",
      "Ocrelizumab"
    ),
    
    n_pct_table1(
      "CRP_below_1",
      1,
      "Ocrelizumab"
    )
  ),
  
  check.names = FALSE
)

# Display Table 1
table1

# Export Table 1 to Excel
write.xlsx(
  table1,
  file = "~/Documents/data_R/Table_1.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)


