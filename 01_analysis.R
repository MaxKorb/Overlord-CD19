# ============================================================
# OVERLORD-MS B-CELL DEPLETION ANALYSIS
# ============================================================
#
# Laboratory sub-study of the OVERLORD-MS trial
#
# This script reproduces:
#   - analysis cohort definition
#   - descriptive analyses
#   - sustained depletion analysis
#   - early depletion analysis
#   - 80% sensitivity analysis
#   - B-cell reappearance analyses
#   - predictor analyses
#   - multiple-testing corrections
#   - Table 1
#   - Supplementary Tables 1-4
#
# Figures are generated separately in:
# scripts/02_figures.R
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(dplyr)
library(openxlsx)
library(rstatix)
library(rcompanion)


# ============================================================
# 2. IMPORT DATA
# ============================================================

data_path <- file.path(
  "data",
  "OVERLORD_B_cell_analysis_data.csv"
)

ny_tabell_raw <- read.csv2(
  data_path,
  header = TRUE,
  check.names = FALSE,
  fileEncoding = "UTF-8"
)


# Basic import checks
dim(ny_tabell_raw)

length(
  unique(ny_tabell_raw$PasientID)
)

# Create output directory if it does not already exist
if (!dir.exists("output")) {
  dir.create("output")
}


# ============================================================
# 3. DEFINE FINAL ANALYSIS COHORT
# ============================================================

# Inclusion criterion:
# At least two available absolute CD19 measurements after baseline.

followup_count <- ny_tabell_raw %>%
  filter(`Måling nr.` > 1) %>%
  group_by(PasientID) %>%
  summarise(
    N_followup_CD19 = sum(!is.na(`CD19_celler_per_µL`)),
    .groups = "drop"
  )


# Participants meeting the inclusion criterion
included_ids <- followup_count %>%
  filter(N_followup_CD19 >= 2) %>%
  pull(PasientID)


# Participants excluded because of fewer than two
# available post-baseline CD19 measurements
excluded_participants <- followup_count %>%
  filter(N_followup_CD19 < 2)


# Create final longitudinal analysis dataset
ny_tabell <- ny_tabell_raw %>%
  filter(PasientID %in% included_ids)


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

# Number of participants before exclusion
n_distinct(ny_tabell_raw$PasientID)

# Excluded participant(s)
excluded_participants

# Number of participants in final analysis cohort
n_distinct(ny_tabell$PasientID)

# Number of longitudinal rows after exclusion
dim(ny_tabell)

# ============================================================
# 4. CREATE PARTICIPANT-LEVEL DATASET AND TREATMENT LABELS
# ============================================================

# Add readable treatment labels to the longitudinal dataset
ny_tabell <- ny_tabell %>%
  mutate(
    Behandling = case_when(
      Treatment == 1 ~ "Rituximab",
      Treatment == 2 ~ "Ocrelizumab",
      TRUE ~ NA_character_
    )
  )


# Keep one row per participant
pasienter <- ny_tabell %>%
  distinct(
    PasientID,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

# Number of participants
nrow(pasienter)

# Treatment distribution
table(pasienter$Behandling)

# Treatment coding check
table(
  pasienter$Treatment,
  pasienter$Behandling
)

# ============================================================
# ROUNDING FUNCTION FOR REPORTED RESULTS
# ============================================================

# Conventional half-up rounding for reported values.
# Examples:
#   11.35 -> 11.4
#   300.35 -> 300.4
#
# This affects presentation only.
# Unrounded values are retained for all analyses.

round_half_up <- function(
    x,
    digits = 0
) {
  
  factor <- 10^digits
  
  sign(x) *
    floor(
      abs(x) * factor + 0.5
    ) /
    factor
}

# ============================================================
# 5. DEFINE SUSTAINED DEPLETION
# ============================================================

# Sustained depletion is defined as:
# CD19+ B cells <= 5 cells/µL at all measured post-baseline
# time points following the first treatment.

# Keep post-baseline measurements only
followup_data <- ny_tabell %>%
  filter(`Måling nr.` > 1)


# Create one sustained depletion status per participant
sustained_status <- followup_data %>%
  group_by(PasientID) %>%
  summarise(
    Deplesjon_status = if_else(
      any(
        !is.na(`CD19_celler_per_µL`) &
          `CD19_celler_per_µL` > 5
      ),
      "Ikke-deplesjon",
      "Deplesjon"
    ),
    .groups = "drop"
  )


# Merge sustained depletion status with participant-level data
pasienter_sustained <- pasienter %>%
  left_join(
    sustained_status,
    by = "PasientID"
  )


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

# Sustained depletion group sizes
table(
  pasienter_sustained$Deplesjon_status
)

# Number and percentage with sustained depletion
n_sustained <- sum(
  pasienter_sustained$Deplesjon_status == "Deplesjon"
)

pct_sustained <- 100 *
  n_sustained /
  nrow(pasienter_sustained)

n_sustained

round_half_up(
  pct_sustained,
  1
)


# ============================================================
# 6. DEFINE EARLY DEPLETION AT MONTH 3
# ============================================================

# Month 3 corresponds to measurement number 2.
# Early depletion is defined as CD19+ B cells <= 5 cells/µL
# at the first post-baseline measurement.

month3_data <- ny_tabell %>%
  filter(`Måling nr.` == 2)


# Create early depletion status
early_status <- month3_data %>%
  transmute(
    PasientID,
    Deplesjon_status = case_when(
      is.na(`CD19_celler_per_µL`) ~ NA_character_,
      `CD19_celler_per_µL` <= 5 ~ "Depletert",
      `CD19_celler_per_µL` > 5 ~ "Ikke-depletert"
    )
  )


# Keep participants with an available month 3 CD19 measurement
early_status_available <- early_status %>%
  filter(!is.na(Deplesjon_status))


# Merge with participant-level data
pasienter_early <- pasienter %>%
  inner_join(
    early_status_available,
    by = "PasientID"
  )


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

# Number with available CD19 measurement at month 3
n_month3_available <- nrow(pasienter_early)

# Early depletion group sizes
table(
  pasienter_early$Deplesjon_status
)

# Number with early depletion
n_early_depleted <- sum(
  pasienter_early$Deplesjon_status == "Depletert"
)

# Percentage with early depletion
pct_early_depleted <- 100 *
  n_early_depleted /
  n_month3_available


n_month3_available

n_early_depleted

round_half_up(
  pct_early_depleted,
  1
)


# ============================================================
# 7. DEFINE 80% SENSITIVITY DEPLETION
# ============================================================

# Sensitivity definition:
# CD19+ B cells <= 5 cells/µL in at least 80% of
# measured post-baseline time points.

sensitivity80_status <- followup_data %>%
  group_by(PasientID) %>%
  summarise(
    proportion_depleted = mean(
      `CD19_celler_per_µL`[
        !is.na(`CD19_celler_per_µL`)
      ] <= 5
    ),
    Deplesjon80 = if_else(
      proportion_depleted >= 0.80,
      "Depletert",
      "Ikke-depletert"
    ),
    .groups = "drop"
  )


# Merge with participant-level data
sens80 <- pasienter %>%
  left_join(
    sensitivity80_status,
    by = "PasientID"
  )


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

# Number of participants
nrow(sens80)

# Sensitivity-analysis depletion groups
table(
  sens80$Deplesjon80
)

# Percentage depleted according to the 80% definition
n_sens80_depleted <- sum(
  sens80$Deplesjon80 == "Depletert"
)

pct_sens80_depleted <- 100 *
  n_sens80_depleted /
  nrow(sens80)

n_sens80_depleted

round_half_up(
  pct_sens80_depleted,
  1
)

# ============================================================
# 8. CREATE ANALYSIS VARIABLES
# ============================================================

# Binary alcohol variable:
# 0 = abstinent
# 1 = current alcohol use

pasienter <- pasienter %>%
  mutate(
    Alkohol_bruk = if_else(
      Alkohol == "Avhold",
      0,
      1
    )
  )


# Add the same variable to sustained depletion dataset
pasienter_sustained <- pasienter_sustained %>%
  mutate(
    Alkohol_bruk = if_else(
      Alkohol == "Avhold",
      0,
      1
    )
  )


# Add the same variable to early depletion dataset
pasienter_early <- pasienter_early %>%
  mutate(
    Alkohol_bruk = if_else(
      Alkohol == "Avhold",
      0,
      1
    )
  )


# Add the same variable to the 80% sensitivity dataset
sens80 <- sens80 %>%
  mutate(
    Alkohol_bruk = if_else(
      Alkohol == "Avhold",
      0,
      1
    )
  )


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

table(
  pasienter$Alkohol_bruk,
  useNA = "ifany"
)

table(
  pasienter_sustained$Alkohol_bruk,
  useNA = "ifany"
)

table(
  pasienter_early$Alkohol_bruk,
  useNA = "ifany"
)

table(
  sens80$Alkohol_bruk,
  useNA = "ifany"
)

# ============================================================
# 9. FOLLOW-UP DURATION AND DATA COMPLETENESS
# ============================================================

# Scheduled post-baseline visits:
# measurement 2 = month 3
# measurement 3 = month 6
# measurement 4 = month 12
# measurement 5 = month 18
# measurement 6 = month 24
# measurement 7 = month 30
# measurement 8 = month 36

visit_months <- c(
  `2` = 3,
  `3` = 6,
  `4` = 12,
  `5` = 18,
  `6` = 24,
  `7` = 30,
  `8` = 36
)


# ------------------------------------------------------------
# 9.1 Available CD19 measurements at each follow-up visit
# ------------------------------------------------------------

availability_by_visit <- followup_data %>%
  group_by(`Måling nr.`) %>%
  summarise(
    N_available = sum(!is.na(`CD19_celler_per_µL`)),
    .groups = "drop"
  ) %>%
  mutate(
    Month = unname(
      visit_months[
        as.character(`Måling nr.`)
      ]
    ),
    N_scheduled = n_distinct(ny_tabell$PasientID),
    Percent_available = 100 * N_available / N_scheduled
  ) %>%
  select(
    Month,
    N_available,
    N_scheduled,
    Percent_available
  )


availability_by_visit


# ------------------------------------------------------------
# 9.2 Overall availability across all scheduled
#     post-baseline visits
# ------------------------------------------------------------

total_available <- sum(
  !is.na(followup_data$`CD19_celler_per_µL`)
)

total_scheduled <- n_distinct(ny_tabell$PasientID) * 7

overall_availability <- 100 *
  total_available /
  total_scheduled


total_available

total_scheduled

round_half_up(
  overall_availability,
  1
)


# Range of availability across follow-up visits
round_half_up(
  range(
    availability_by_visit$Percent_available
  ),
  1
)


# ------------------------------------------------------------
# 9.3 Availability at month 36
# ------------------------------------------------------------

month36_available <- availability_by_visit %>%
  filter(
    Month == 36
  )

month36_available


# ------------------------------------------------------------
# 9.4 Duration of CD19 follow-up
# ------------------------------------------------------------

last_followup <- followup_data %>%
  filter(
    !is.na(`CD19_celler_per_µL`)
  ) %>%
  mutate(
    Month = unname(
      visit_months[
        as.character(`Måling nr.`)
      ]
    )
  ) %>%
  group_by(PasientID) %>%
  summarise(
    Last_followup_month = max(
      Month
    ),
    .groups = "drop"
  )


mean_followup <- mean(
  last_followup$Last_followup_month
)

sd_followup <- sd(
  last_followup$Last_followup_month
)


round_half_up(
  mean_followup,
  1
)

round_half_up(
  sd_followup,
  2
)

# ============================================================
# 10. BASELINE DEMOGRAPHICS OF THE FINAL ANALYSIS COHORT
# ============================================================

# ------------------------------------------------------------
# 10.1 Sex distribution
# ------------------------------------------------------------

sex_distribution <- table(
  pasienter$Kjønn
)

sex_distribution


# Number and percentage of women
n_female <- sum(
  pasienter$Kjønn == "Female",
  na.rm = TRUE
)

pct_female <- 100 *
  n_female /
  sum(!is.na(pasienter$Kjønn))


n_female

round_half_up(
  pct_female,
  1
)


# ------------------------------------------------------------
# 10.2 Age
# ------------------------------------------------------------

mean_age <- mean(
  pasienter$`Alder_(år)`,
  na.rm = TRUE
)

sd_age <- sd(
  pasienter$`Alder_(år)`,
  na.rm = TRUE
)


round_half_up(
  mean_age,
  1
)

round_half_up(
  sd_age,
  1
)

# ============================================================
# 11. DEPLETION OVER TIME AND B-CELL REAPPEARANCE
# ============================================================

# ------------------------------------------------------------
# 11.1 Depletion at each post-baseline visit
# ------------------------------------------------------------

# Depletion at a given visit is defined as:
# CD19+ B cells <= 5 cells/µL.
#
# Percentages are calculated among participants with an
# available CD19 measurement at that specific visit.

depletion_by_visit <- followup_data %>%
  filter(
    !is.na(`CD19_celler_per_µL`)
  ) %>%
  mutate(
    Month = unname(
      visit_months[
        as.character(`Måling nr.`)
      ]
    )
  ) %>%
  group_by(Month) %>%
  summarise(
    N_measured = n(),
    N_depleted = sum(
      `CD19_celler_per_µL` <= 5
    ),
    Percent_depleted = 100 *
      N_depleted /
      N_measured,
    .groups = "drop"
  )


depletion_by_visit


# ------------------------------------------------------------
# 11.2 Participants with early depletion at month 3
# ------------------------------------------------------------

early_depleted_ids <- pasienter_early %>%
  filter(
    Deplesjon_status == "Depletert"
  ) %>%
  pull(PasientID)


length(
  early_depleted_ids
)


# ------------------------------------------------------------
# 11.3 B-cell reappearance after early depletion
# ------------------------------------------------------------

# Reappearance is defined as at least one subsequent
# CD19+ B cell measurement > 5 cells/µL after month 3
# among participants who were depleted at month 3.

later_measurements_early_depleted <- ny_tabell %>%
  filter(
    PasientID %in% early_depleted_ids,
    `Måling nr.` > 2
  )


reappearance_by_participant <- later_measurements_early_depleted %>%
  group_by(PasientID) %>%
  summarise(
    Reappearance = any(
      !is.na(`CD19_celler_per_µL`) &
        `CD19_celler_per_µL` > 5
    ),
    .groups = "drop"
  )


n_reappearance <- sum(
  reappearance_by_participant$Reappearance
)

pct_reappearance <- 100 *
  n_reappearance /
  length(early_depleted_ids)


n_reappearance

round_half_up(
  pct_reappearance,
  1
)


# ------------------------------------------------------------
# 11.4 Number of non-depleted post-baseline visits
#      in the full analysis cohort
# ------------------------------------------------------------

# Among all participants in the final analysis cohort,
# identify participants with at least one post-baseline
# CD19+ B cell measurement > 5 cells/µL.
#
# For these participants, count the number of post-baseline
# visits with CD19+ B cells > 5 cells/µL.

nondepleted_visit_counts <- followup_data %>%
  filter(
    !is.na(`CD19_celler_per_µL`),
    `CD19_celler_per_µL` > 5
  ) %>%
  group_by(PasientID) %>%
  summarise(
    N_nondepleted_visits = n(),
    .groups = "drop"
  )


# Number of participants with at least one
# non-depleted post-baseline visit
n_any_nondepleted <- nrow(
  nondepleted_visit_counts
)


# Summary of number of non-depleted visits
mean_nondepleted_visits <- mean(
  nondepleted_visit_counts$N_nondepleted_visits
)

sd_nondepleted_visits <- sd(
  nondepleted_visit_counts$N_nondepleted_visits
)

range_nondepleted_visits <- range(
  nondepleted_visit_counts$N_nondepleted_visits
)


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

n_any_nondepleted

round_half_up(
  mean_nondepleted_visits,
  2
)

round_half_up(
  sd_nondepleted_visits,
  2
)

range_nondepleted_visits

# ============================================================
# 12. DEPLETION ACCORDING TO TREATMENT
# ============================================================

# ------------------------------------------------------------
# 12.1 Sustained depletion by treatment
# ------------------------------------------------------------

sustained_by_treatment <- table(
  pasienter_sustained$Behandling,
  pasienter_sustained$Deplesjon_status
)

sustained_by_treatment


# Percent sustained depleted within each treatment group
sustained_by_treatment_percent <- prop.table(
  sustained_by_treatment,
  margin = 1
) * 100

round_half_up(
  sustained_by_treatment_percent,
  1
)


# Pearson's chi-square test with Yates' continuity correction
sustained_treatment_test <- chisq.test(
  sustained_by_treatment,
  correct = TRUE
)

sustained_treatment_test


# Extract test results
sustained_treatment_chisq <- unname(
  sustained_treatment_test$statistic
)

sustained_treatment_df <- unname(
  sustained_treatment_test$parameter
)

sustained_treatment_p <- sustained_treatment_test$p.value


round_half_up(
  sustained_treatment_chisq,
  2
)

sustained_treatment_df

round_half_up(
  sustained_treatment_p,
  3
)


# ------------------------------------------------------------
# 12.2 Early depletion by treatment
# ------------------------------------------------------------

early_by_treatment <- table(
  pasienter_early$Behandling,
  pasienter_early$Deplesjon_status
)

early_by_treatment


# Percent early depleted within each treatment group
early_by_treatment_percent <- prop.table(
  early_by_treatment,
  margin = 1
) * 100

round_half_up(
  early_by_treatment_percent,
  1
)


# Pearson's chi-square test with Yates' continuity correction
early_treatment_test <- chisq.test(
  early_by_treatment,
  correct = TRUE
)

early_treatment_test


# Extract test results
early_treatment_chisq <- unname(
  early_treatment_test$statistic
)

early_treatment_df <- unname(
  early_treatment_test$parameter
)

early_treatment_p <- early_treatment_test$p.value


round_half_up(
  early_treatment_chisq,
  2
)

early_treatment_df

round_half_up(
  early_treatment_p,
  3
)


# ------------------------------------------------------------
# 12.3 Check treatment of the four participants
#      without early depletion
# ------------------------------------------------------------

non_early_depleted <- pasienter_early %>%
  filter(
    Deplesjon_status == "Ikke-depletert"
  )

table(
  non_early_depleted$Behandling
)

# ============================================================
# 13. COURSE OF PARTICIPANTS WITHOUT EARLY DEPLETION
# ============================================================

# Identify the four participants who did not meet the
# depletion threshold at month 3.

non_early_ids <- pasienter_early %>%
  filter(
    Deplesjon_status == "Ikke-depletert"
  ) %>%
  pull(PasientID)


# ------------------------------------------------------------
# 13.1 Show their complete post-baseline CD19 course
# ------------------------------------------------------------

non_early_course <- ny_tabell %>%
  filter(
    PasientID %in% non_early_ids,
    `Måling nr.` >= 2
  ) %>%
  mutate(
    Month = unname(
      visit_months[
        as.character(`Måling nr.`)
      ]
    ),
    Visit_depletion = case_when(
      is.na(`CD19_celler_per_µL`) ~ NA_character_,
      `CD19_celler_per_µL` <= 5 ~ "Depletert",
      `CD19_celler_per_µL` > 5 ~ "Ikke-depletert"
    )
  ) %>%
  select(
    PasientID,
    Month,
    `CD19_celler_per_µL`,
    Visit_depletion
  ) %>%
  arrange(
    PasientID,
    Month
  )


non_early_course


# ------------------------------------------------------------
# 13.2 Check depletion during early follow-up
#     (months 3–18)
# ------------------------------------------------------------

non_early_month3_18 <- non_early_course %>%
  filter(
    Month %in% c(3, 6, 12, 18)
  ) %>%
  group_by(PasientID) %>%
  summarise(
    Any_depletion_month3_18 = any(
      `CD19_celler_per_µL` <= 5,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


non_early_month3_18


# Number who achieved depletion during months 3–18
sum(
  non_early_month3_18$Any_depletion_month3_18
)


# ------------------------------------------------------------
# 13.3 Check depletion at month 24
# ------------------------------------------------------------

non_early_month24 <- non_early_course %>%
  filter(
    Month == 24
  )


non_early_month24


# Number depleted at month 24
sum(
  non_early_month24$`CD19_celler_per_µL` <= 5,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 13.4 Display later course (months 24–36)
# ------------------------------------------------------------

non_early_late_course <- non_early_course %>%
  filter(
    Month %in% c(24, 30, 36)
  )


non_early_late_course

# ============================================================
# 14. MEDIAN CD19+ B-CELL COUNTS OVER TIME BY TREATMENT
# ============================================================

# Calculate median and interquartile range (IQR) of absolute
# CD19+ B-cell counts at each post-baseline visit,
# separately for rituximab and ocrelizumab.

cd19_by_treatment <- followup_data %>%
  filter(
    !is.na(`CD19_celler_per_µL`)
  ) %>%
  mutate(
    Month = unname(
      visit_months[
        as.character(`Måling nr.`)
      ]
    )
  ) %>%
  group_by(
    Behandling,
    Month
  ) %>%
  summarise(
    N = n(),
    Median = median(
      `CD19_celler_per_µL`
    ),
    Q1 = quantile(
      `CD19_celler_per_µL`,
      0.25
    ),
    Q3 = quantile(
      `CD19_celler_per_µL`,
      0.75
    ),
    IQR = IQR(
      `CD19_celler_per_µL`
    ),
    .groups = "drop"
  ) %>%
  arrange(
    Behandling,
    Month
  )


cd19_by_treatment


# ------------------------------------------------------------
# 14.1 Range of treatment-specific median CD19 counts
# ------------------------------------------------------------

median_range <- range(
  cd19_by_treatment$Median,
  na.rm = TRUE
)

median_range


# ------------------------------------------------------------
# 14.2 Clean display table
# ------------------------------------------------------------

cd19_by_treatment_display <- cd19_by_treatment %>%
  mutate(
    Median = round_half_up(
      Median,
      2
    ),
    Q1 = round_half_up(
      Q1,
      2
    ),
    Q3 = round_half_up(
      Q3,
      2
    ),
    IQR = round_half_up(
      IQR,
      2
    )
  )


cd19_by_treatment_display

# ============================================================
# 15. DEFINE VARIABLES FOR BASELINE PREDICTOR ANALYSES
# ============================================================

# The same 17 baseline predictors are tested in:
#   - sustained depletion
#   - 80% sensitivity analysis
#
# Early depletion is described descriptively only,
# because only four participants were non-depleted at month 3.
#
# Multiple-testing correction will later be performed across
# all 17 tests within each depletion definition separately.


# ------------------------------------------------------------
# 15.1 Continuous variables
# ------------------------------------------------------------

continuous_variables <- c(
  "Alder_(år)",
  "BMI",
  "Tid_siden_MS_diagnose_(måneder)",
  "Tid_siden_første_kliniske_hendelse_(måneder)",
  "Minimum_EDSS_score",
  "Total_antall_attakker",
  "Antall_attakker_siste_året",
  "Antall_oligoklonale_bånd",
  "IgG",
  "IgM",
  "CD19_%",
  "CD19_celler_per_µL"
)


# ------------------------------------------------------------
# 15.2 Categorical variables
# ------------------------------------------------------------

categorical_variables <- c(
  "Kjønn",
  "Røyk",
  "Snus",
  "Alkohol_bruk",
  "CRP_below_1"
)


# ------------------------------------------------------------
# 15.3 Checks
# ------------------------------------------------------------

length(
  continuous_variables
)

length(
  categorical_variables
)

length(
  c(
    continuous_variables,
    categorical_variables
  )
)


# Check that all variables exist in the participant dataset
missing_analysis_variables <- setdiff(
  c(
    continuous_variables,
    categorical_variables
  ),
  names(pasienter)
)

missing_analysis_variables

# ============================================================
# 16. FUNCTIONS FOR BASELINE PREDICTOR ANALYSES
# ============================================================

# The same 17 baseline predictors are tested in:
#   - sustained depletion
#   - 80% sensitivity analysis
#
# Early depletion is described descriptively only,
# because only four participants were non-depleted at month 3.
#
# Continuous variables:
#   Wilcoxon rank-sum test
#   median [IQR]
#   Wilcoxon effect size r (Z / sqrt(N))
#   with 95% bootstrap confidence interval
#
# Categorical variables:
#   Pearson's chi-square test with Yates' continuity correction
#   bias-corrected Cramer's V


# Set seed to make bootstrap confidence intervals reproducible.
# This does not affect the Wilcoxon tests or other analyses.
set.seed(2026)


# ------------------------------------------------------------
# 16.1 Continuous predictors
# ------------------------------------------------------------

run_continuous_predictor <- function(
    data,
    variable,
    group_variable,
    nondepleted_label,
    depleted_label
) {
  
  # Keep complete observations
  analysis_data <- data.frame(
    value = data[[variable]],
    group = data[[group_variable]]
  ) %>%
    filter(
      !is.na(value),
      !is.na(group)
    )
  
  
  # Explicit group order:
  # group 0 = non-depletion
  # group 1 = depletion
  analysis_data$group <- factor(
    analysis_data$group,
    levels = c(
      nondepleted_label,
      depleted_label
    )
  )
  
  
  # Group-specific values
  x_non_depleted <- analysis_data$value[
    analysis_data$group == nondepleted_label
  ]
  
  x_depleted <- analysis_data$value[
    analysis_data$group == depleted_label
  ]
  
  
  # Wilcoxon rank-sum test
  wilcox_result <- wilcox.test(
    value ~ group,
    data = analysis_data,
    exact = FALSE
  )
  
  
  # Wilcoxon effect size:
  # r = Z / sqrt(N)
  # 95% bootstrap confidence interval
  effect_result <- rstatix::wilcox_effsize(
    analysis_data,
    value ~ group,
    ci = TRUE
  )
  
  
  # Create result row
  data.frame(
    Variable = variable,
    
    N = nrow(
      analysis_data
    ),
    
    n_non_depleted = length(
      x_non_depleted
    ),
    
    n_depleted = length(
      x_depleted
    ),
    
    median_non_depleted = median(
      x_non_depleted
    ),
    
    Q1_non_depleted = quantile(
      x_non_depleted,
      0.25,
      names = FALSE
    ),
    
    Q3_non_depleted = quantile(
      x_non_depleted,
      0.75,
      names = FALSE
    ),
    
    median_depleted = median(
      x_depleted
    ),
    
    Q1_depleted = quantile(
      x_depleted,
      0.25,
      names = FALSE
    ),
    
    Q3_depleted = quantile(
      x_depleted,
      0.75,
      names = FALSE
    ),
    
    p_value = wilcox_result$p.value,
    
    effect_size_r = unname(
      effect_result$effsize
    ),
    
    effect_size_ci_low = unname(
      effect_result$conf.low
    ),
    
    effect_size_ci_high = unname(
      effect_result$conf.high
    )
  )
}


# ------------------------------------------------------------
# 16.2 Categorical predictors
# ------------------------------------------------------------

run_categorical_predictor <- function(
    data,
    variable,
    group_variable,
    nondepleted_label,
    depleted_label
) {
  
  # Keep complete observations
  analysis_data <- data.frame(
    value = data[[variable]],
    group = data[[group_variable]]
  ) %>%
    filter(
      !is.na(value),
      !is.na(group)
    )
  
  
  # Explicit group order
  analysis_data$group <- factor(
    analysis_data$group,
    levels = c(
      nondepleted_label,
      depleted_label
    )
  )
  
  
  # Contingency table
  tab <- table(
    analysis_data$value,
    analysis_data$group
  )
  
  
  # Pearson's chi-square test
  # Yates' continuity correction is applied for 2 x 2 tables
  chisq_result <- chisq.test(
    tab,
    correct = TRUE
  )
  
  
  # Bias-corrected Cramer's V
  cramer_v_result <- rcompanion::cramerV(
    tab,
    bias.correct = TRUE
  )
  
  
  # Create result row
  data.frame(
    Variable = variable,
    
    N = sum(
      tab
    ),
    
    n_non_depleted = sum(
      analysis_data$group == nondepleted_label
    ),
    
    n_depleted = sum(
      analysis_data$group == depleted_label
    ),
    
    statistic = unname(
      chisq_result$statistic
    ),
    
    df = unname(
      chisq_result$parameter
    ),
    
    p_value = chisq_result$p.value,
    
    cramer_v = unname(
      cramer_v_result
    )
  )
}


# ------------------------------------------------------------
# 16.3 Sustained depletion predictor analyses
# ------------------------------------------------------------

sustained_continuous <- do.call(
  rbind,
  lapply(
    continuous_variables,
    function(variable) {
      
      run_continuous_predictor(
        data = pasienter_sustained,
        variable = variable,
        group_variable = "Deplesjon_status",
        nondepleted_label = "Ikke-deplesjon",
        depleted_label = "Deplesjon"
      )
      
    }
  )
)


sustained_categorical <- do.call(
  rbind,
  lapply(
    categorical_variables,
    function(variable) {
      
      run_categorical_predictor(
        data = pasienter_sustained,
        variable = variable,
        group_variable = "Deplesjon_status",
        nondepleted_label = "Ikke-deplesjon",
        depleted_label = "Deplesjon"
      )
      
    }
  )
)


# ------------------------------------------------------------
# 16.4 Checks
# ------------------------------------------------------------

sustained_continuous

sustained_categorical


# Check number of predictor tests
nrow(
  sustained_continuous
)

nrow(
  sustained_categorical
)

nrow(
  sustained_continuous
) +
  nrow(
    sustained_categorical
  )

# ============================================================
# 17. 80% SENSITIVITY PREDICTOR ANALYSES
# ============================================================

# Predictor analyses for early depletion are not performed
# because only four participants were classified as
# non-depleted at month 3. Early depletion is therefore
# described descriptively only.


# ------------------------------------------------------------
# 17.1 Continuous predictors
# ------------------------------------------------------------

sensitivity80_continuous <- do.call(
  rbind,
  lapply(
    continuous_variables,
    function(variable) {
      
      run_continuous_predictor(
        data = sens80,
        variable = variable,
        group_variable = "Deplesjon80",
        nondepleted_label = "Ikke-depletert",
        depleted_label = "Depletert"
      )
      
    }
  )
)


# ------------------------------------------------------------
# 17.2 Categorical predictors
# ------------------------------------------------------------

sensitivity80_categorical <- do.call(
  rbind,
  lapply(
    categorical_variables,
    function(variable) {
      
      run_categorical_predictor(
        data = sens80,
        variable = variable,
        group_variable = "Deplesjon80",
        nondepleted_label = "Ikke-depletert",
        depleted_label = "Depletert"
      )
      
    }
  )
)


# ------------------------------------------------------------
# 17.3 Checks
# ------------------------------------------------------------

sensitivity80_continuous

sensitivity80_categorical


# Confirm 17 predictor tests
nrow(
  sensitivity80_continuous
) +
  nrow(
    sensitivity80_categorical
  )

# ============================================================
# 18. MULTIPLE-TESTING CORRECTION
# ============================================================

# Multiple-testing correction is performed separately for:
#   1. Sustained depletion
#   2. 80% sensitivity analysis
#
# Within each analysis, all 17 baseline predictor tests
# (12 continuous + 5 categorical) are treated as one family.
#
# Two correction methods are applied:
#   - Benjamini-Hochberg false discovery rate (FDR)
#   - Bonferroni correction


# ------------------------------------------------------------
# 18.1 Sustained depletion
# ------------------------------------------------------------

# Combine p-values from continuous and categorical predictors
sustained_all_pvalues <- bind_rows(
  sustained_continuous %>%
    transmute(
      Variable = Variable,
      Variable_type = "Continuous",
      p_value = p_value
    ),
  
  sustained_categorical %>%
    transmute(
      Variable = Variable,
      Variable_type = "Categorical",
      p_value = p_value
    )
)


# Apply corrections across all 17 tests
sustained_all_pvalues <- sustained_all_pvalues %>%
  mutate(
    p_FDR = p.adjust(
      p_value,
      method = "BH"
    ),
    
    p_Bonferroni = p.adjust(
      p_value,
      method = "bonferroni"
    )
  )


# ------------------------------------------------------------
# 18.2 80% sensitivity analysis
# ------------------------------------------------------------

# Combine p-values from continuous and categorical predictors
sensitivity80_all_pvalues <- bind_rows(
  sensitivity80_continuous %>%
    transmute(
      Variable = Variable,
      Variable_type = "Continuous",
      p_value = p_value
    ),
  
  sensitivity80_categorical %>%
    transmute(
      Variable = Variable,
      Variable_type = "Categorical",
      p_value = p_value
    )
)


# Apply corrections across all 17 tests
sensitivity80_all_pvalues <- sensitivity80_all_pvalues %>%
  mutate(
    p_FDR = p.adjust(
      p_value,
      method = "BH"
    ),
    
    p_Bonferroni = p.adjust(
      p_value,
      method = "bonferroni"
    )
  )


# ------------------------------------------------------------
# 18.3 Checks
# ------------------------------------------------------------

# Confirm that each correction is based on 17 tests
nrow(
  sustained_all_pvalues
)

nrow(
  sensitivity80_all_pvalues
)


# Display sustained depletion results
sustained_all_pvalues %>%
  arrange(
    p_value
  )


# Display 80% sensitivity results
sensitivity80_all_pvalues %>%
  arrange(
    p_value
  )


# ------------------------------------------------------------
# 18.4 Significant results after correction
# ------------------------------------------------------------

# Sustained depletion:
# significant after FDR correction
sustained_all_pvalues %>%
  filter(
    p_FDR < 0.05
  )


# Sustained depletion:
# significant after Bonferroni correction
sustained_all_pvalues %>%
  filter(
    p_Bonferroni < 0.05
  )


# 80% sensitivity:
# significant after FDR correction
sensitivity80_all_pvalues %>%
  filter(
    p_FDR < 0.05
  )


# 80% sensitivity:
# significant after Bonferroni correction
sensitivity80_all_pvalues %>%
  filter(
    p_Bonferroni < 0.05
  )

# ============================================================
# 19. TREATMENT COMPARISON IN THE 80% SENSITIVITY ANALYSIS
# ============================================================

# Compare depletion status according to the 80% sensitivity
# definition between rituximab and ocrelizumab.


# ------------------------------------------------------------
# 19.1 Contingency table
# ------------------------------------------------------------

sensitivity80_by_treatment <- table(
  sens80$Behandling,
  sens80$Deplesjon80
)

sensitivity80_by_treatment


# ------------------------------------------------------------
# 19.2 Percentages within treatment groups
# ------------------------------------------------------------

sensitivity80_by_treatment_percent <- prop.table(
  sensitivity80_by_treatment,
  margin = 1
) * 100

round_half_up(
  sensitivity80_by_treatment_percent,
  1
)


# ------------------------------------------------------------
# 19.3 Pearson's chi-square test
# ------------------------------------------------------------

sensitivity80_treatment_test <- chisq.test(
  sensitivity80_by_treatment,
  correct = TRUE
)

sensitivity80_treatment_test


# Extract test results
sensitivity80_treatment_chisq <- unname(
  sensitivity80_treatment_test$statistic
)

sensitivity80_treatment_df <- unname(
  sensitivity80_treatment_test$parameter
)

sensitivity80_treatment_p <- sensitivity80_treatment_test$p.value


# ------------------------------------------------------------
# 19.4 Cramer's V
# ------------------------------------------------------------

sensitivity80_treatment_cramer_v <- unname(
  rcompanion::cramerV(
    sensitivity80_by_treatment,
    bias.correct = TRUE
  )
)


# ------------------------------------------------------------
# 19.5 Checks
# ------------------------------------------------------------

sensitivity80_by_treatment

round_half_up(
  sensitivity80_by_treatment_percent,
  1
)

round_half_up(
  sensitivity80_treatment_chisq,
  2
)

sensitivity80_treatment_df

round_half_up(
  sensitivity80_treatment_p,
  3
)

round_half_up(
  sensitivity80_treatment_cramer_v,
  2
)


# ============================================================
# 20. TABLE 1: PARTICIPANT CHARACTERISTICS BY TREATMENT
# ============================================================

# Table 1 presents baseline characteristics for:
#   - Total population
#   - Rituximab
#   - Ocrelizumab
#
# Continuous variables are summarised as mean ± SD.
# Categorical variables are summarised as n (%).


# ------------------------------------------------------------
# 20.1 Helper functions
# ------------------------------------------------------------

mean_sd_table1 <- function(
    data,
    variable
) {
  
  x <- data[[variable]]
  x <- x[!is.na(x)]
  
  mean_value <- round_half_up(
    mean(x),
    digits = 1
  )
  
  sd_value <- round_half_up(
    sd(x),
    digits = 1
  )
  
  sprintf(
    "%.1f ± %.1f",
    mean_value,
    sd_value
  )
}


n_pct_table1 <- function(
    data,
    variable,
    value
) {
  
  x <- data[[variable]]
  x <- x[!is.na(x)]
  
  n <- sum(
    x == value
  )
  
  pct <- round_half_up(
    100 * n / length(x),
    digits = 1
  )
  
  sprintf(
    "%d (%.1f)",
    n,
    pct
  )
}


# ------------------------------------------------------------
# 20.2 Treatment-specific datasets
# ------------------------------------------------------------

pasienter_rituximab <- pasienter %>%
  filter(
    Behandling == "Rituximab"
  )


pasienter_ocrelizumab <- pasienter %>%
  filter(
    Behandling == "Ocrelizumab"
  )


# ------------------------------------------------------------
# 20.3 Checks
# ------------------------------------------------------------

nrow(
  pasienter
)

nrow(
  pasienter_rituximab
)

nrow(
  pasienter_ocrelizumab
)


# ------------------------------------------------------------
# 20.4 Create Table 1
# ------------------------------------------------------------

table1 <- data.frame(
  
  Characteristic = c(
    
    "Number of participants, n (%)",
    
    "Demographics and Lifestyle Characteristics",
    
    "Age, years, mean (±SD)",
    
    "Sex, female, n (%)",
    
    "Sex, male, n (%)",
    
    "Body Mass Index, kg/m², mean (±SD)",
    
    "Current use of smoked tobacco, n (%)",
    
    "Current use of smokeless tobacco (e.g. snuff), n (%)",
    
    "Current use of alcohol, n (%)",
    
    "Baseline Disease Characteristics",
    
    "Time since MS diagnosis, months, mean (±SD)",
    
    "Time since first neurological symptom, months, mean (±SD)",
    
    "Minimum EDSS score, mean (±SD)",
    
    "Total relapses since onset, mean (±SD)",
    
    "Relapses in the past 12 months, mean (±SD)",
    
    "CSF oligoclonal bands tested, n (%)",
    
    "CSF oligoclonal bands, mean (±SD)",
    
    "Baseline Laboratory Values",
    
    "IgG, g/L, mean (±SD)",
    
    "IgM, g/L, mean (±SD)",
    
    "CD19+ B cells, % of lymphocytes, mean (±SD)",
    
    "CD19+ B cells, cells/µL, mean (±SD)",
    
    "CRP < 1 mg/L, n (%)"
  ),
  
  
  Total_population = c(
    
    sprintf(
      "%d (100.0)",
      nrow(pasienter)
    ),
    
    "",
    
    mean_sd_table1(
      pasienter,
      "Alder_(år)"
    ),
    
    n_pct_table1(
      pasienter,
      "Kjønn",
      "Female"
    ),
    
    n_pct_table1(
      pasienter,
      "Kjønn",
      "Male"
    ),
    
    mean_sd_table1(
      pasienter,
      "BMI"
    ),
    
    n_pct_table1(
      pasienter,
      "Røyk",
      1
    ),
    
    n_pct_table1(
      pasienter,
      "Snus",
      "Ja"
    ),
    
    n_pct_table1(
      pasienter,
      "Alkohol_bruk",
      1
    ),
    
    "",
    
    mean_sd_table1(
      pasienter,
      "Tid_siden_MS_diagnose_(måneder)"
    ),
    
    mean_sd_table1(
      pasienter,
      "Tid_siden_første_kliniske_hendelse_(måneder)"
    ),
    
    mean_sd_table1(
      pasienter,
      "Minimum_EDSS_score"
    ),
    
    mean_sd_table1(
      pasienter,
      "Total_antall_attakker"
    ),
    
    mean_sd_table1(
      pasienter,
      "Antall_attakker_siste_året"
    ),
    
    n_pct_table1(
      pasienter,
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja"
    ),
    
    mean_sd_table1(
      pasienter,
      "Antall_oligoklonale_bånd"
    ),
    
    "",
    
    mean_sd_table1(
      pasienter,
      "IgG"
    ),
    
    mean_sd_table1(
      pasienter,
      "IgM"
    ),
    
    mean_sd_table1(
      pasienter,
      "CD19_%"
    ),
    
    mean_sd_table1(
      pasienter,
      "CD19_celler_per_µL"
    ),
    
    n_pct_table1(
      pasienter,
      "CRP_below_1",
      1
    )
  ),
  
  
  Rituximab = c(
    
    sprintf(
      "%d (%.1f)",
      nrow(pasienter_rituximab),
      round_half_up(
        100 *
          nrow(pasienter_rituximab) /
          nrow(pasienter),
        digits = 1
      )
    ),
    
    "",
    
    mean_sd_table1(
      pasienter_rituximab,
      "Alder_(år)"
    ),
    
    n_pct_table1(
      pasienter_rituximab,
      "Kjønn",
      "Female"
    ),
    
    n_pct_table1(
      pasienter_rituximab,
      "Kjønn",
      "Male"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "BMI"
    ),
    
    n_pct_table1(
      pasienter_rituximab,
      "Røyk",
      1
    ),
    
    n_pct_table1(
      pasienter_rituximab,
      "Snus",
      "Ja"
    ),
    
    n_pct_table1(
      pasienter_rituximab,
      "Alkohol_bruk",
      1
    ),
    
    "",
    
    mean_sd_table1(
      pasienter_rituximab,
      "Tid_siden_MS_diagnose_(måneder)"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "Tid_siden_første_kliniske_hendelse_(måneder)"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "Minimum_EDSS_score"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "Total_antall_attakker"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "Antall_attakker_siste_året"
    ),
    
    n_pct_table1(
      pasienter_rituximab,
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "Antall_oligoklonale_bånd"
    ),
    
    "",
    
    mean_sd_table1(
      pasienter_rituximab,
      "IgG"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "IgM"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "CD19_%"
    ),
    
    mean_sd_table1(
      pasienter_rituximab,
      "CD19_celler_per_µL"
    ),
    
    n_pct_table1(
      pasienter_rituximab,
      "CRP_below_1",
      1
    )
  ),
  
  
  Ocrelizumab = c(
    
    sprintf(
      "%d (%.1f)",
      nrow(pasienter_ocrelizumab),
      round_half_up(
        100 *
          nrow(pasienter_ocrelizumab) /
          nrow(pasienter),
        digits = 1
      )
    ),
    
    "",
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "Alder_(år)"
    ),
    
    n_pct_table1(
      pasienter_ocrelizumab,
      "Kjønn",
      "Female"
    ),
    
    n_pct_table1(
      pasienter_ocrelizumab,
      "Kjønn",
      "Male"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "BMI"
    ),
    
    n_pct_table1(
      pasienter_ocrelizumab,
      "Røyk",
      1
    ),
    
    n_pct_table1(
      pasienter_ocrelizumab,
      "Snus",
      "Ja"
    ),
    
    n_pct_table1(
      pasienter_ocrelizumab,
      "Alkohol_bruk",
      1
    ),
    
    "",
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "Tid_siden_MS_diagnose_(måneder)"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "Tid_siden_første_kliniske_hendelse_(måneder)"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "Minimum_EDSS_score"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "Total_antall_attakker"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "Antall_attakker_siste_året"
    ),
    
    n_pct_table1(
      pasienter_ocrelizumab,
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "Antall_oligoklonale_bånd"
    ),
    
    "",
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "IgG"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "IgM"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "CD19_%"
    ),
    
    mean_sd_table1(
      pasienter_ocrelizumab,
      "CD19_celler_per_µL"
    ),
    
    n_pct_table1(
      pasienter_ocrelizumab,
      "CRP_below_1",
      1
    )
  ),
  
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 20.5 Display Table 1
# ------------------------------------------------------------

table1


# ------------------------------------------------------------
# 20.6 Export Table 1
# ------------------------------------------------------------

# Export as CSV
write.csv(
  table1,
  file = "output/Table_1_participant_characteristics.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# Export as Excel
openxlsx::write.xlsx(
  table1,
  file = "output/Table_1_participant_characteristics.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)


# Check that both files were created
file.exists(
  "output/Table_1_participant_characteristics.csv"
)

file.exists(
  "output/Table_1_participant_characteristics.xlsx"
)


# ============================================================
# 21. SUPPLEMENTARY TABLE 1:
# BASELINE CHARACTERISTICS BY TREATMENT AND
# SUSTAINED DEPLETION STATUS
# ============================================================

# Supplementary Table 1 presents baseline characteristics
# stratified by treatment and sustained depletion status.
#
# Continuous variables are summarised as mean ± SD.
# Categorical variables are summarised as n (%).


# ------------------------------------------------------------
# 21.1 Create the four analysis groups
# ------------------------------------------------------------

rtx_nondepleted <- pasienter_sustained %>%
  filter(
    Behandling == "Rituximab",
    Deplesjon_status == "Ikke-deplesjon"
  )


rtx_depleted <- pasienter_sustained %>%
  filter(
    Behandling == "Rituximab",
    Deplesjon_status == "Deplesjon"
  )


ocr_nondepleted <- pasienter_sustained %>%
  filter(
    Behandling == "Ocrelizumab",
    Deplesjon_status == "Ikke-deplesjon"
  )


ocr_depleted <- pasienter_sustained %>%
  filter(
    Behandling == "Ocrelizumab",
    Deplesjon_status == "Deplesjon"
  )


# ------------------------------------------------------------
# 21.2 Check group sizes
# ------------------------------------------------------------

nrow(rtx_nondepleted)
nrow(rtx_depleted)
nrow(ocr_nondepleted)
nrow(ocr_depleted)


# Expected:
# Rituximab non-depletion = 36
# Rituximab depletion = 34
# Ocrelizumab non-depletion = 15
# Ocrelizumab depletion = 15


# ------------------------------------------------------------
# 21.3 Helper function for one Supplementary Table 1 column
# ------------------------------------------------------------

create_supp1_column <- function(data) {
  
  c(
    
    # Number of participants
    as.character(
      nrow(data)
    ),
    
    # Section heading
    "",
    
    # Demographics and lifestyle
    mean_sd_table1(
      data,
      "Alder_(år)"
    ),
    
    n_pct_table1(
      data,
      "Kjønn",
      "Female"
    ),
    
    n_pct_table1(
      data,
      "Kjønn",
      "Male"
    ),
    
    mean_sd_table1(
      data,
      "BMI"
    ),
    
    n_pct_table1(
      data,
      "Røyk",
      1
    ),
    
    n_pct_table1(
      data,
      "Snus",
      "Ja"
    ),
    
    n_pct_table1(
      data,
      "Alkohol_bruk",
      1
    ),
    
    # Section heading
    "",
    
    # Baseline disease characteristics
    mean_sd_table1(
      data,
      "Tid_siden_MS_diagnose_(måneder)"
    ),
    
    mean_sd_table1(
      data,
      "Tid_siden_første_kliniske_hendelse_(måneder)"
    ),
    
    mean_sd_table1(
      data,
      "Minimum_EDSS_score"
    ),
    
    mean_sd_table1(
      data,
      "Total_antall_attakker"
    ),
    
    mean_sd_table1(
      data,
      "Antall_attakker_siste_året"
    ),
    
    n_pct_table1(
      data,
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja"
    ),
    
    mean_sd_table1(
      data,
      "Antall_oligoklonale_bånd"
    ),
    
    # Section heading
    "",
    
    # Baseline laboratory values
    mean_sd_table1(
      data,
      "IgG"
    ),
    
    mean_sd_table1(
      data,
      "IgM"
    ),
    
    mean_sd_table1(
      data,
      "CD19_%"
    ),
    
    mean_sd_table1(
      data,
      "CD19_celler_per_µL"
    ),
    
    n_pct_table1(
      data,
      "CRP_below_1",
      1
    )
  )
}


# ------------------------------------------------------------
# 21.4 Create Supplementary Table 1
# ------------------------------------------------------------

supplementary_table1 <- data.frame(
  
  Characteristic = c(
    
    "Number of participants, n",
    
    "Demographics and Lifestyle Characteristics",
    
    "Age, years, mean (±SD)",
    
    "Sex, female, n (%)",
    
    "Sex, male, n (%)",
    
    "Body Mass Index, kg/m², mean (±SD)",
    
    "Current use of smoked tobacco, n (%)",
    
    "Current use of smokeless tobacco (e.g. snuff), n (%)",
    
    "Current use of alcohol, n (%)",
    
    "Baseline Disease Characteristics",
    
    "Time since MS diagnosis, months, mean (±SD)",
    
    "Time since first neurological symptom, months, mean (±SD)",
    
    "Minimum EDSS score, mean (±SD)",
    
    "Total relapses since onset, mean (±SD)",
    
    "Relapses in the past 12 months, mean (±SD)",
    
    "CSF oligoclonal bands tested, n (%)",
    
    "CSF oligoclonal bands, mean (±SD)",
    
    "Baseline Laboratory Values",
    
    "IgG, g/L, mean (±SD)",
    
    "IgM, g/L, mean (±SD)",
    
    "CD19, % of lymphocytes, mean (±SD)",
    
    "CD19, cells/µL, mean (±SD)",
    
    "CRP < 1 mg/L, n (%)"
  ),
  
  
  `Rituximab - Non-depletion` =
    create_supp1_column(
      rtx_nondepleted
    ),
  
  `Rituximab - Depletion` =
    create_supp1_column(
      rtx_depleted
    ),
  
  `Ocrelizumab - Non-depletion` =
    create_supp1_column(
      ocr_nondepleted
    ),
  
  `Ocrelizumab - Depletion` =
    create_supp1_column(
      ocr_depleted
    ),
  
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 21.5 Display Supplementary Table 1
# ------------------------------------------------------------

supplementary_table1


# ------------------------------------------------------------
# 21.6 Export Supplementary Table 1
# ------------------------------------------------------------

# CSV
write.csv(
  supplementary_table1,
  file = "output/Supplementary_Table_1.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# Excel
openxlsx::write.xlsx(
  supplementary_table1,
  file = "output/Supplementary_Table_1.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)


# ------------------------------------------------------------
# 21.7 Check exports
# ------------------------------------------------------------

file.exists(
  "output/Supplementary_Table_1.csv"
)

file.exists(
  "output/Supplementary_Table_1.xlsx"
)

# ============================================================
# 22. SUPPLEMENTARY TABLE 2:
# BASELINE CHARACTERISTICS BY TREATMENT AND
# EARLY DEPLETION STATUS
# ============================================================

# Supplementary Table 2 presents baseline characteristics
# stratified by treatment and early depletion status.
#
# Continuous variables are summarised as mean ± SD.
# Categorical variables are summarised as n (%).
#
# No inferential predictor analyses are performed for early
# depletion because only four participants were non-depleted
# at month 3.


# ------------------------------------------------------------
# 22.1 Create the four analysis groups
# ------------------------------------------------------------

rtx_early_nondepleted <- pasienter_early %>%
  filter(
    Behandling == "Rituximab",
    Deplesjon_status == "Ikke-depletert"
  )


rtx_early_depleted <- pasienter_early %>%
  filter(
    Behandling == "Rituximab",
    Deplesjon_status == "Depletert"
  )


ocr_early_nondepleted <- pasienter_early %>%
  filter(
    Behandling == "Ocrelizumab",
    Deplesjon_status == "Ikke-depletert"
  )


ocr_early_depleted <- pasienter_early %>%
  filter(
    Behandling == "Ocrelizumab",
    Deplesjon_status == "Depletert"
  )


# ------------------------------------------------------------
# 22.2 Check group sizes
# ------------------------------------------------------------

nrow(rtx_early_nondepleted)
nrow(rtx_early_depleted)
nrow(ocr_early_nondepleted)
nrow(ocr_early_depleted)

# Expected:
# Rituximab non-depletion = 4
# Rituximab depletion = 63
# Ocrelizumab non-depletion = 0
# Ocrelizumab depletion = 27


# ------------------------------------------------------------
# 22.3 Helper functions that also handle empty groups
# ------------------------------------------------------------

mean_sd_supp2 <- function(
    data,
    variable
) {
  
  if (nrow(data) == 0) {
    return("—")
  }
  
  x <- data[[variable]]
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return("—")
  }
  
  mean_value <- round_half_up(
    mean(x),
    digits = 1
  )
  
  sd_value <- round_half_up(
    sd(x),
    digits = 1
  )
  
  sprintf(
    "%.1f ± %.1f",
    mean_value,
    sd_value
  )
}


n_pct_supp2 <- function(
    data,
    variable,
    value
) {
  
  if (nrow(data) == 0) {
    return("—")
  }
  
  x <- data[[variable]]
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return("—")
  }
  
  n <- sum(
    x == value
  )
  
  pct <- round_half_up(
    100 * n / length(x),
    digits = 1
  )
  
  sprintf(
    "%d (%.1f)",
    n,
    pct
  )
}


# ------------------------------------------------------------
# 22.4 Helper function for one Supplementary Table 2 column
# ------------------------------------------------------------

create_supp2_column <- function(data) {
  
  c(
    
    # Number of participants
    ifelse(
      nrow(data) == 0,
      "0",
      as.character(
        nrow(data)
      )
    ),
    
    # Section heading
    "",
    
    # Demographics and lifestyle
    mean_sd_supp2(
      data,
      "Alder_(år)"
    ),
    
    n_pct_supp2(
      data,
      "Kjønn",
      "Female"
    ),
    
    n_pct_supp2(
      data,
      "Kjønn",
      "Male"
    ),
    
    mean_sd_supp2(
      data,
      "BMI"
    ),
    
    n_pct_supp2(
      data,
      "Røyk",
      1
    ),
    
    n_pct_supp2(
      data,
      "Snus",
      "Ja"
    ),
    
    n_pct_supp2(
      data,
      "Alkohol_bruk",
      1
    ),
    
    # Section heading
    "",
    
    # Baseline disease characteristics
    mean_sd_supp2(
      data,
      "Tid_siden_MS_diagnose_(måneder)"
    ),
    
    mean_sd_supp2(
      data,
      "Tid_siden_første_kliniske_hendelse_(måneder)"
    ),
    
    mean_sd_supp2(
      data,
      "Minimum_EDSS_score"
    ),
    
    mean_sd_supp2(
      data,
      "Total_antall_attakker"
    ),
    
    mean_sd_supp2(
      data,
      "Antall_attakker_siste_året"
    ),
    
    n_pct_supp2(
      data,
      "Cerebrospinalvæske_analysert_for_oligoklonale_bånd",
      "Ja"
    ),
    
    mean_sd_supp2(
      data,
      "Antall_oligoklonale_bånd"
    ),
    
    # Section heading
    "",
    
    # Baseline laboratory values
    mean_sd_supp2(
      data,
      "IgG"
    ),
    
    mean_sd_supp2(
      data,
      "IgM"
    ),
    
    mean_sd_supp2(
      data,
      "CD19_%"
    ),
    
    mean_sd_supp2(
      data,
      "CD19_celler_per_µL"
    ),
    
    n_pct_supp2(
      data,
      "CRP_below_1",
      1
    )
  )
}


# ------------------------------------------------------------
# 22.5 Create Supplementary Table 2
# ------------------------------------------------------------

supplementary_table2 <- data.frame(
  
  Characteristic = c(
    
    "Number of participants, n",
    
    "Demographics and Lifestyle Characteristics",
    
    "Age, years, mean (±SD)",
    
    "Sex, female, n (%)",
    
    "Sex, male, n (%)",
    
    "Body Mass Index, kg/m², mean (±SD)",
    
    "Current use of smoked tobacco, n (%)",
    
    "Current use of smokeless tobacco (e.g. snuff), n (%)",
    
    "Current use of alcohol, n (%)",
    
    "Baseline Disease Characteristics",
    
    "Time since MS diagnosis, months, mean (±SD)",
    
    "Time since first neurological symptom, months, mean (±SD)",
    
    "Minimum EDSS score, mean (±SD)",
    
    "Total relapses since onset, mean (±SD)",
    
    "Relapses in the past 12 months, mean (±SD)",
    
    "CSF oligoclonal bands tested, n (%)",
    
    "CSF oligoclonal bands, mean (±SD)",
    
    "Baseline Laboratory Values",
    
    "IgG, g/L, mean (±SD)",
    
    "IgM, g/L, mean (±SD)",
    
    "CD19+ B cells, % of lymphocytes, mean (±SD)",
    
    "CD19+ B cells, cells/µL, mean (±SD)",
    
    "CRP < 1 mg/L, n (%)"
  ),
  
  
  `Rituximab - Non-depletion` =
    create_supp2_column(
      rtx_early_nondepleted
    ),
  
  `Rituximab - Depletion` =
    create_supp2_column(
      rtx_early_depleted
    ),
  
  `Ocrelizumab - Non-depletion` =
    create_supp2_column(
      ocr_early_nondepleted
    ),
  
  `Ocrelizumab - Depletion` =
    create_supp2_column(
      ocr_early_depleted
    ),
  
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 22.6 Display Supplementary Table 2
# ------------------------------------------------------------

supplementary_table2


# ------------------------------------------------------------
# 22.7 Export Supplementary Table 2
# ------------------------------------------------------------

write.csv(
  supplementary_table2,
  file = "output/Supplementary_Table_2.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


openxlsx::write.xlsx(
  supplementary_table2,
  file = "output/Supplementary_Table_2.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)


# ------------------------------------------------------------
# 22.8 Check exports
# ------------------------------------------------------------

file.exists(
  "output/Supplementary_Table_2.csv"
)

file.exists(
  "output/Supplementary_Table_2.xlsx"
)


# ============================================================
# 23. SUPPLEMENTARY TABLE 3:
# PREDICTORS OF SUSTAINED DEPLETION
# ============================================================

# This table combines:
#   - descriptive statistics
#   - Wilcoxon rank-sum tests for continuous predictors
#   - Pearson's chi-square tests for categorical predictors
#   - Wilcoxon effect size r (Z / sqrt(N)) with 95% bootstrap CI
#   - bias-corrected Cramer's V
#   - raw p-values
#   - FDR-adjusted p-values
#   - Bonferroni-adjusted p-values
#
# Multiple-testing correction is based on the 17 predictor
# tests within the sustained depletion analysis.


# ------------------------------------------------------------
# 23.1 Formatting functions
# ------------------------------------------------------------

format_p <- function(x) {
  
  if (is.na(x)) {
    return("—")
  }
  
  if (x < 0.001) {
    return("<0.001")
  }
  
  x <- round_half_up(
    x,
    digits = 3
  )
  
  sprintf(
    "%.3f",
    x
  )
}


format_number <- function(
    x,
    digits = 2
) {
  
  if (is.na(x)) {
    return("—")
  }
  
  x <- round_half_up(
    x,
    digits = digits
  )
  
  sprintf(
    paste0(
      "%.",
      digits,
      "f"
    ),
    x
  )
}


format_median_iqr <- function(
    median,
    q1,
    q3
) {
  
  median <- round_half_up(
    median,
    digits = 1
  )
  
  q1 <- round_half_up(
    q1,
    digits = 1
  )
  
  q3 <- round_half_up(
    q3,
    digits = 1
  )
  
  sprintf(
    "%.1f [%.1f-%.1f]",
    median,
    q1,
    q3
  )
}


format_effect_r <- function(
    r,
    ci_low,
    ci_high
) {
  
  r <- round_half_up(
    r,
    digits = 2
  )
  
  ci_low <- round_half_up(
    ci_low,
    digits = 2
  )
  
  ci_high <- round_half_up(
    ci_high,
    digits = 2
  )
  
  sprintf(
    "%.2f (%.2f-%.2f)",
    r,
    ci_low,
    ci_high
  )
}


# ------------------------------------------------------------
# 23.2 Add corrected p-values to analysis results
# ------------------------------------------------------------

sustained_continuous_final <- sustained_continuous %>%
  left_join(
    sustained_all_pvalues %>%
      select(
        Variable,
        p_FDR,
        p_Bonferroni
      ),
    by = "Variable"
  )


sustained_categorical_final <- sustained_categorical %>%
  left_join(
    sustained_all_pvalues %>%
      select(
        Variable,
        p_FDR,
        p_Bonferroni
      ),
    by = "Variable"
  )


# ------------------------------------------------------------
# 23.3 Helper functions for table rows
# ------------------------------------------------------------

blank_supp3_row <- function(label) {
  
  data.frame(
    Variable = label,
    Statistical_test = "",
    `Non-depletion, n (%)` = "",
    `Depletion, n (%)` = "",
    `Non-depletion, median [IQR]` = "",
    `Depletion, median [IQR]` = "",
    Test_statistic = "",
    Effect_size = "",
    df = "",
    p_value = "",
    FDR_adjusted_p_value = "",
    Bonferroni_adjusted_p_value = "",
    check.names = FALSE
  )
}


continuous_supp3_row <- function(
    variable,
    label
) {
  
  x <- sustained_continuous_final %>%
    filter(
      Variable == variable
    )
  
  nondep_pct <- round_half_up(
    100 * x$n_non_depleted / 51,
    digits = 1
  )
  
  dep_pct <- round_half_up(
    100 * x$n_depleted / 49,
    digits = 1
  )
  
  
  data.frame(
    Variable = label,
    
    Statistical_test =
      "Wilcoxon rank-sum test",
    
    `Non-depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        x$n_non_depleted,
        nondep_pct
      ),
    
    `Depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        x$n_depleted,
        dep_pct
      ),
    
    `Non-depletion, median [IQR]` =
      format_median_iqr(
        x$median_non_depleted,
        x$Q1_non_depleted,
        x$Q3_non_depleted
      ),
    
    `Depletion, median [IQR]` =
      format_median_iqr(
        x$median_depleted,
        x$Q1_depleted,
        x$Q3_depleted
      ),
    
    Test_statistic = "—",
    
    Effect_size =
      format_effect_r(
        x$effect_size_r,
        x$effect_size_ci_low,
        x$effect_size_ci_high
      ),
    
    df = "—",
    
    p_value =
      format_p(
        x$p_value
      ),
    
    FDR_adjusted_p_value =
      format_p(
        x$p_FDR
      ),
    
    Bonferroni_adjusted_p_value =
      format_p(
        x$p_Bonferroni
      ),
    
    check.names = FALSE
  )
}


categorical_supp3_row <- function(
    variable,
    label
) {
  
  x <- sustained_categorical_final %>%
    filter(
      Variable == variable
    )
  
  nondep_pct <- round_half_up(
    100 * x$n_non_depleted / 51,
    digits = 1
  )
  
  dep_pct <- round_half_up(
    100 * x$n_depleted / 49,
    digits = 1
  )
  
  
  data.frame(
    Variable = label,
    
    Statistical_test =
      "Pearson's χ² test",
    
    `Non-depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        x$n_non_depleted,
        nondep_pct
      ),
    
    `Depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        x$n_depleted,
        dep_pct
      ),
    
    `Non-depletion, median [IQR]` = "—",
    
    `Depletion, median [IQR]` = "—",
    
    Test_statistic =
      format_number(
        x$statistic,
        2
      ),
    
    Effect_size =
      format_number(
        x$cramer_v,
        2
      ),
    
    df =
      as.character(
        x$df
      ),
    
    p_value =
      format_p(
        x$p_value
      ),
    
    FDR_adjusted_p_value =
      format_p(
        x$p_FDR
      ),
    
    Bonferroni_adjusted_p_value =
      format_p(
        x$p_Bonferroni
      ),
    
    check.names = FALSE
  )
}


categorical_level_supp3_row <- function(
    label,
    variable,
    value
) {
  
  nondep <- pasienter_sustained %>%
    filter(
      Deplesjon_status == "Ikke-deplesjon"
    )
  
  dep <- pasienter_sustained %>%
    filter(
      Deplesjon_status == "Deplesjon"
    )
  
  
  x0 <- nondep[[variable]]
  x1 <- dep[[variable]]
  
  denominator0 <- sum(
    !is.na(x0)
  )
  
  denominator1 <- sum(
    !is.na(x1)
  )
  
  n0 <- sum(
    x0 == value,
    na.rm = TRUE
  )
  
  n1 <- sum(
    x1 == value,
    na.rm = TRUE
  )
  
  
  pct0 <- round_half_up(
    100 * n0 / denominator0,
    digits = 1
  )
  
  pct1 <- round_half_up(
    100 * n1 / denominator1,
    digits = 1
  )
  
  
  data.frame(
    Variable = label,
    Statistical_test = "—",
    
    `Non-depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        n0,
        pct0
      ),
    
    `Depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        n1,
        pct1
      ),
    
    `Non-depletion, median [IQR]` = "—",
    `Depletion, median [IQR]` = "—",
    Test_statistic = "—",
    Effect_size = "—",
    df = "—",
    p_value = "—",
    FDR_adjusted_p_value = "—",
    Bonferroni_adjusted_p_value = "—",
    
    check.names = FALSE
  )
}


# ------------------------------------------------------------
# 23.4 Treatment rows
# ------------------------------------------------------------

sustained_treatment_cramer_v <- unname(
  rcompanion::cramerV(
    sustained_by_treatment,
    bias.correct = TRUE
  )
)


treatment_row_supp3 <- data.frame(
  
  Variable = "Treatment",
  
  Statistical_test =
    "Pearson's χ² test",
  
  `Non-depletion, n (%)` =
    "51 (100.0)",
  
  `Depletion, n (%)` =
    "49 (100.0)",
  
  `Non-depletion, median [IQR]` = "—",
  
  `Depletion, median [IQR]` = "—",
  
  Test_statistic =
    format_number(
      sustained_treatment_chisq,
      2
    ),
  
  Effect_size =
    format_number(
      sustained_treatment_cramer_v,
      2
    ),
  
  df =
    as.character(
      sustained_treatment_df
    ),
  
  p_value =
    format_p(
      sustained_treatment_p
    ),
  
  FDR_adjusted_p_value = "—",
  
  Bonferroni_adjusted_p_value = "—",
  
  check.names = FALSE
)


treatment_level_supp3 <- function(
    label,
    value
) {
  
  nondep <- pasienter_sustained %>%
    filter(
      Deplesjon_status == "Ikke-deplesjon"
    )
  
  dep <- pasienter_sustained %>%
    filter(
      Deplesjon_status == "Deplesjon"
    )
  
  
  n0 <- sum(
    nondep$Behandling == value
  )
  
  n1 <- sum(
    dep$Behandling == value
  )
  
  
  pct0 <- round_half_up(
    100 * n0 / nrow(nondep),
    digits = 1
  )
  
  pct1 <- round_half_up(
    100 * n1 / nrow(dep),
    digits = 1
  )
  
  
  data.frame(
    Variable = label,
    Statistical_test = "—",
    
    `Non-depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        n0,
        pct0
      ),
    
    `Depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        n1,
        pct1
      ),
    
    `Non-depletion, median [IQR]` = "—",
    `Depletion, median [IQR]` = "—",
    Test_statistic = "—",
    Effect_size = "—",
    df = "—",
    p_value = "—",
    FDR_adjusted_p_value = "—",
    Bonferroni_adjusted_p_value = "—",
    
    check.names = FALSE
  )
}


# ------------------------------------------------------------
# 23.5 Create Supplementary Table 3
# ------------------------------------------------------------

supplementary_table3 <- bind_rows(
  
  # Number of participants
  data.frame(
    Variable = "Number of participants",
    Statistical_test = "—",
    `Non-depletion, n (%)` = "51",
    `Depletion, n (%)` = "49",
    `Non-depletion, median [IQR]` = "—",
    `Depletion, median [IQR]` = "—",
    Test_statistic = "—",
    Effect_size = "—",
    df = "—",
    p_value = "—",
    FDR_adjusted_p_value = "—",
    Bonferroni_adjusted_p_value = "—",
    check.names = FALSE
  ),
  
  
  # Treatment
  treatment_row_supp3,
  
  treatment_level_supp3(
    "Rituximab",
    "Rituximab"
  ),
  
  treatment_level_supp3(
    "Ocrelizumab",
    "Ocrelizumab"
  ),
  
  
  # Demographics and lifestyle
  blank_supp3_row(
    "Demographics and Lifestyle Characteristics"
  ),
  
  continuous_supp3_row(
    "Alder_(år)",
    "Age, years"
  ),
  
  categorical_supp3_row(
    "Kjønn",
    "Sex"
  ),
  
  categorical_level_supp3_row(
    "Female",
    "Kjønn",
    "Female"
  ),
  
  categorical_level_supp3_row(
    "Male",
    "Kjønn",
    "Male"
  ),
  
  continuous_supp3_row(
    "BMI",
    "Body Mass Index, kg/m²"
  ),
  
  categorical_supp3_row(
    "Røyk",
    "Smoking"
  ),
  
  categorical_level_supp3_row(
    "Yes",
    "Røyk",
    1
  ),
  
  categorical_level_supp3_row(
    "No",
    "Røyk",
    0
  ),
  
  categorical_supp3_row(
    "Snus",
    "Smokeless tobacco use (e.g. snuff)"
  ),
  
  categorical_level_supp3_row(
    "Yes",
    "Snus",
    "Ja"
  ),
  
  categorical_level_supp3_row(
    "No",
    "Snus",
    "Nei"
  ),
  
  categorical_supp3_row(
    "Alkohol_bruk",
    "Alcohol use"
  ),
  
  categorical_level_supp3_row(
    "Yes",
    "Alkohol_bruk",
    1
  ),
  
  categorical_level_supp3_row(
    "No",
    "Alkohol_bruk",
    0
  ),
  
  
  # Baseline disease characteristics
  blank_supp3_row(
    "Baseline Disease Characteristics"
  ),
  
  continuous_supp3_row(
    "Tid_siden_MS_diagnose_(måneder)",
    "Time since MS diagnosis, months"
  ),
  
  continuous_supp3_row(
    "Tid_siden_første_kliniske_hendelse_(måneder)",
    "Time since first neurological symptom, months"
  ),
  
  continuous_supp3_row(
    "Minimum_EDSS_score",
    "Minimum EDSS score"
  ),
  
  continuous_supp3_row(
    "Total_antall_attakker",
    "Total relapses since onset"
  ),
  
  continuous_supp3_row(
    "Antall_attakker_siste_året",
    "Relapses in the past 12 months"
  ),
  
  continuous_supp3_row(
    "Antall_oligoklonale_bånd",
    "CSF oligoclonal bands"
  ),
  
  
  # Baseline laboratory values
  blank_supp3_row(
    "Baseline Laboratory Values"
  ),
  
  continuous_supp3_row(
    "IgG",
    "IgG, g/L"
  ),
  
  continuous_supp3_row(
    "IgM",
    "IgM, g/L"
  ),
  
  continuous_supp3_row(
    "CD19_%",
    "CD19+ B cells, % of lymphocytes"
  ),
  
  continuous_supp3_row(
    "CD19_celler_per_µL",
    "CD19+ B cells, cells/µL"
  ),
  
  categorical_supp3_row(
    "CRP_below_1",
    "CRP < 1 mg/L"
  ),
  
  categorical_level_supp3_row(
    "Yes",
    "CRP_below_1",
    1
  ),
  
  categorical_level_supp3_row(
    "No",
    "CRP_below_1",
    0
  )
)


# ------------------------------------------------------------
# 23.6 Display Supplementary Table 3
# ------------------------------------------------------------

supplementary_table3


# ------------------------------------------------------------
# 23.7 Export Supplementary Table 3
# ------------------------------------------------------------

write.csv(
  supplementary_table3,
  file = "output/Supplementary_Table_3.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


openxlsx::write.xlsx(
  supplementary_table3,
  file = "output/Supplementary_Table_3.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)


# ------------------------------------------------------------
# 23.8 Check exports
# ------------------------------------------------------------

file.exists(
  "output/Supplementary_Table_3.csv"
)

file.exists(
  "output/Supplementary_Table_3.xlsx"
)

# ============================================================
# 24. SUPPLEMENTARY TABLE 4:
# PREDICTORS OF DEPLETION IN THE 80% SENSITIVITY ANALYSIS
# ============================================================

# This table combines:
#   - descriptive statistics
#   - Wilcoxon rank-sum tests for continuous predictors
#   - Pearson's chi-square tests for categorical predictors
#   - Wilcoxon effect size r (Z / sqrt(N)) with 95% bootstrap CI
#   - bias-corrected Cramer's V
#   - raw p-values
#   - FDR-adjusted p-values
#   - Bonferroni-adjusted p-values
#
# Multiple-testing correction is based on the 17 predictor
# tests within the 80% sensitivity analysis.


# ------------------------------------------------------------
# 24.1 Add corrected p-values to analysis results
# ------------------------------------------------------------

sensitivity80_continuous_final <- sensitivity80_continuous %>%
  left_join(
    sensitivity80_all_pvalues %>%
      select(
        Variable,
        p_FDR,
        p_Bonferroni
      ),
    by = "Variable"
  )


sensitivity80_categorical_final <- sensitivity80_categorical %>%
  left_join(
    sensitivity80_all_pvalues %>%
      select(
        Variable,
        p_FDR,
        p_Bonferroni
      ),
    by = "Variable"
  )


# ------------------------------------------------------------
# 24.2 Helper functions for table rows
# ------------------------------------------------------------

blank_supp4_row <- function(label) {
  
  data.frame(
    Variable = label,
    Statistical_test = "",
    `Non-depletion, n (%)` = "",
    `Depletion, n (%)` = "",
    `Non-depletion, median [IQR]` = "",
    `Depletion, median [IQR]` = "",
    Test_statistic = "",
    Effect_size = "",
    df = "",
    p_value = "",
    FDR_adjusted_p_value = "",
    Bonferroni_adjusted_p_value = "",
    check.names = FALSE
  )
}


continuous_supp4_row <- function(
    variable,
    label
) {
  
  x <- sensitivity80_continuous_final %>%
    filter(
      Variable == variable
    )
  
  nondep_pct <- round_half_up(
    100 * x$n_non_depleted / 33,
    digits = 1
  )
  
  dep_pct <- round_half_up(
    100 * x$n_depleted / 67,
    digits = 1
  )
  
  
  data.frame(
    Variable = label,
    
    Statistical_test =
      "Wilcoxon rank-sum test",
    
    `Non-depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        x$n_non_depleted,
        nondep_pct
      ),
    
    `Depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        x$n_depleted,
        dep_pct
      ),
    
    `Non-depletion, median [IQR]` =
      format_median_iqr(
        x$median_non_depleted,
        x$Q1_non_depleted,
        x$Q3_non_depleted
      ),
    
    `Depletion, median [IQR]` =
      format_median_iqr(
        x$median_depleted,
        x$Q1_depleted,
        x$Q3_depleted
      ),
    
    Test_statistic = "—",
    
    Effect_size =
      format_effect_r(
        x$effect_size_r,
        x$effect_size_ci_low,
        x$effect_size_ci_high
      ),
    
    df = "—",
    
    p_value =
      format_p(
        x$p_value
      ),
    
    FDR_adjusted_p_value =
      format_p(
        x$p_FDR
      ),
    
    Bonferroni_adjusted_p_value =
      format_p(
        x$p_Bonferroni
      ),
    
    check.names = FALSE
  )
}


categorical_supp4_row <- function(
    variable,
    label
) {
  
  x <- sensitivity80_categorical_final %>%
    filter(
      Variable == variable
    )
  
  nondep_pct <- round_half_up(
    100 * x$n_non_depleted / 33,
    digits = 1
  )
  
  dep_pct <- round_half_up(
    100 * x$n_depleted / 67,
    digits = 1
  )
  
  
  data.frame(
    Variable = label,
    
    Statistical_test =
      "Pearson's χ² test",
    
    `Non-depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        x$n_non_depleted,
        nondep_pct
      ),
    
    `Depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        x$n_depleted,
        dep_pct
      ),
    
    `Non-depletion, median [IQR]` = "—",
    
    `Depletion, median [IQR]` = "—",
    
    Test_statistic =
      format_number(
        x$statistic,
        2
      ),
    
    Effect_size =
      format_number(
        x$cramer_v,
        2
      ),
    
    df =
      as.character(
        x$df
      ),
    
    p_value =
      format_p(
        x$p_value
      ),
    
    FDR_adjusted_p_value =
      format_p(
        x$p_FDR
      ),
    
    Bonferroni_adjusted_p_value =
      format_p(
        x$p_Bonferroni
      ),
    
    check.names = FALSE
  )
}


categorical_level_supp4_row <- function(
    label,
    variable,
    value
) {
  
  nondep <- sens80 %>%
    filter(
      Deplesjon80 == "Ikke-depletert"
    )
  
  dep <- sens80 %>%
    filter(
      Deplesjon80 == "Depletert"
    )
  
  
  x0 <- nondep[[variable]]
  x1 <- dep[[variable]]
  
  denominator0 <- sum(
    !is.na(x0)
  )
  
  denominator1 <- sum(
    !is.na(x1)
  )
  
  n0 <- sum(
    x0 == value,
    na.rm = TRUE
  )
  
  n1 <- sum(
    x1 == value,
    na.rm = TRUE
  )
  
  
  pct0 <- round_half_up(
    100 * n0 / denominator0,
    digits = 1
  )
  
  pct1 <- round_half_up(
    100 * n1 / denominator1,
    digits = 1
  )
  
  
  data.frame(
    Variable = label,
    Statistical_test = "—",
    
    `Non-depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        n0,
        pct0
      ),
    
    `Depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        n1,
        pct1
      ),
    
    `Non-depletion, median [IQR]` = "—",
    `Depletion, median [IQR]` = "—",
    Test_statistic = "—",
    Effect_size = "—",
    df = "—",
    p_value = "—",
    FDR_adjusted_p_value = "—",
    Bonferroni_adjusted_p_value = "—",
    
    check.names = FALSE
  )
}


# ------------------------------------------------------------
# 24.3 Treatment rows
# ------------------------------------------------------------

treatment_row_supp4 <- data.frame(
  
  Variable = "Treatment",
  
  Statistical_test =
    "Pearson's χ² test",
  
  `Non-depletion, n (%)` =
    "33 (100.0)",
  
  `Depletion, n (%)` =
    "67 (100.0)",
  
  `Non-depletion, median [IQR]` = "—",
  
  `Depletion, median [IQR]` = "—",
  
  Test_statistic =
    format_number(
      sensitivity80_treatment_chisq,
      2
    ),
  
  Effect_size =
    format_number(
      sensitivity80_treatment_cramer_v,
      2
    ),
  
  df =
    as.character(
      sensitivity80_treatment_df
    ),
  
  p_value =
    format_p(
      sensitivity80_treatment_p
    ),
  
  FDR_adjusted_p_value = "—",
  
  Bonferroni_adjusted_p_value = "—",
  
  check.names = FALSE
)


treatment_level_supp4 <- function(
    label,
    value
) {
  
  nondep <- sens80 %>%
    filter(
      Deplesjon80 == "Ikke-depletert"
    )
  
  dep <- sens80 %>%
    filter(
      Deplesjon80 == "Depletert"
    )
  
  
  n0 <- sum(
    nondep$Behandling == value
  )
  
  n1 <- sum(
    dep$Behandling == value
  )
  
  
  pct0 <- round_half_up(
    100 * n0 / nrow(nondep),
    digits = 1
  )
  
  pct1 <- round_half_up(
    100 * n1 / nrow(dep),
    digits = 1
  )
  
  
  data.frame(
    Variable = label,
    Statistical_test = "—",
    
    `Non-depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        n0,
        pct0
      ),
    
    `Depletion, n (%)` =
      sprintf(
        "%d (%.1f)",
        n1,
        pct1
      ),
    
    `Non-depletion, median [IQR]` = "—",
    `Depletion, median [IQR]` = "—",
    Test_statistic = "—",
    Effect_size = "—",
    df = "—",
    p_value = "—",
    FDR_adjusted_p_value = "—",
    Bonferroni_adjusted_p_value = "—",
    
    check.names = FALSE
  )
}


# ------------------------------------------------------------
# 24.4 Create Supplementary Table 4
# ------------------------------------------------------------

supplementary_table4 <- bind_rows(
  
  # Number of participants
  data.frame(
    Variable = "Number of participants",
    Statistical_test = "—",
    `Non-depletion, n (%)` = "33",
    `Depletion, n (%)` = "67",
    `Non-depletion, median [IQR]` = "—",
    `Depletion, median [IQR]` = "—",
    Test_statistic = "—",
    Effect_size = "—",
    df = "—",
    p_value = "—",
    FDR_adjusted_p_value = "—",
    Bonferroni_adjusted_p_value = "—",
    check.names = FALSE
  ),
  
  
  # Treatment
  treatment_row_supp4,
  
  treatment_level_supp4(
    "Rituximab",
    "Rituximab"
  ),
  
  treatment_level_supp4(
    "Ocrelizumab",
    "Ocrelizumab"
  ),
  
  
  # Demographics and lifestyle
  blank_supp4_row(
    "Demographics and Lifestyle Characteristics"
  ),
  
  continuous_supp4_row(
    "Alder_(år)",
    "Age, years"
  ),
  
  categorical_supp4_row(
    "Kjønn",
    "Sex"
  ),
  
  categorical_level_supp4_row(
    "Female",
    "Kjønn",
    "Female"
  ),
  
  categorical_level_supp4_row(
    "Male",
    "Kjønn",
    "Male"
  ),
  
  continuous_supp4_row(
    "BMI",
    "Body Mass Index, kg/m²"
  ),
  
  categorical_supp4_row(
    "Røyk",
    "Smoking"
  ),
  
  categorical_level_supp4_row(
    "Yes",
    "Røyk",
    1
  ),
  
  categorical_level_supp4_row(
    "No",
    "Røyk",
    0
  ),
  
  categorical_supp4_row(
    "Snus",
    "Smokeless tobacco use (e.g. snuff)"
  ),
  
  categorical_level_supp4_row(
    "Yes",
    "Snus",
    "Ja"
  ),
  
  categorical_level_supp4_row(
    "No",
    "Snus",
    "Nei"
  ),
  
  categorical_supp4_row(
    "Alkohol_bruk",
    "Alcohol use"
  ),
  
  categorical_level_supp4_row(
    "Yes",
    "Alkohol_bruk",
    1
  ),
  
  categorical_level_supp4_row(
    "No",
    "Alkohol_bruk",
    0
  ),
  
  
  # Baseline disease characteristics
  blank_supp4_row(
    "Baseline Disease Characteristics"
  ),
  
  continuous_supp4_row(
    "Tid_siden_MS_diagnose_(måneder)",
    "Time since MS diagnosis, months"
  ),
  
  continuous_supp4_row(
    "Tid_siden_første_kliniske_hendelse_(måneder)",
    "Time since first neurological symptom, months"
  ),
  
  continuous_supp4_row(
    "Minimum_EDSS_score",
    "Minimum EDSS score"
  ),
  
  continuous_supp4_row(
    "Total_antall_attakker",
    "Total relapses since onset"
  ),
  
  continuous_supp4_row(
    "Antall_attakker_siste_året",
    "Relapses in the past 12 months"
  ),
  
  continuous_supp4_row(
    "Antall_oligoklonale_bånd",
    "CSF oligoclonal bands"
  ),
  
  
  # Baseline laboratory values
  blank_supp4_row(
    "Baseline Laboratory Values"
  ),
  
  continuous_supp4_row(
    "IgG",
    "IgG, g/L"
  ),
  
  continuous_supp4_row(
    "IgM",
    "IgM, g/L"
  ),
  
  continuous_supp4_row(
    "CD19_%",
    "CD19+ B cells, % of lymphocytes"
  ),
  
  continuous_supp4_row(
    "CD19_celler_per_µL",
    "CD19+ B cells, cells/µL"
  ),
  
  categorical_supp4_row(
    "CRP_below_1",
    "CRP < 1 mg/L"
  ),
  
  categorical_level_supp4_row(
    "Yes",
    "CRP_below_1",
    1
  ),
  
  categorical_level_supp4_row(
    "No",
    "CRP_below_1",
    0
  )
)


# ------------------------------------------------------------
# 24.5 Display Supplementary Table 4
# ------------------------------------------------------------

supplementary_table4


# ------------------------------------------------------------
# 24.6 Export Supplementary Table 4
# ------------------------------------------------------------

write.csv(
  supplementary_table4,
  file = "output/Supplementary_Table_4.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


openxlsx::write.xlsx(
  supplementary_table4,
  file = "output/Supplementary_Table_4.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)


# ------------------------------------------------------------
# 24.7 Check exports
# ------------------------------------------------------------

file.exists(
  "output/Supplementary_Table_4.csv"
)

file.exists(
  "output/Supplementary_Table_4.xlsx"
)

# ============================================================
# 25. REPRODUCIBILITY INFORMATION
# ============================================================

sessionInfo()