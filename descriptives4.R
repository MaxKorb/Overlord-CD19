# ============================================================
# OVERLORD-MS B-cell depletion analysis
# Updated after 36-month data update
# Max Korbmacher
# ============================================================

rm(list = ls(all.names = TRUE))
gc()

# ============================================================
# PACKAGES
# ============================================================

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse,
  readxl,
  openxlsx,
  rstatix,
  effectsize,
  rcompanion,
  boot,
  effsize
)

# ============================================================
# IMPORT DATA
# ============================================================

path <- "/Users/max/Documents/Local/MS/Student_Projects/"

df_raw <- read_excel(
  paste0(
    path,
    "OVERLORD_B_cell_analysis_data.xlsx"
  )
)

# ============================================================
# DATA CLEANING
# ============================================================

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

for(v in continuous_vars){
  
  df_raw[[v]] <- as.numeric(
    gsub(
      "[^0-9\\.-]",
      "",
      gsub(
        ",",
        ".",
        as.character(df_raw[[v]])
      )
    )
  )
  
}

df_raw$Røyk <- as.numeric(
  df_raw$Røyk
)

# ============================================================
# DEFINE ANALYSIS COHORT
# At least two follow-up CD19 measurements
# ============================================================

followup_count <- df_raw %>%
  group_by(
    PasientID
  ) %>%
  summarise(
    n_followup_cd19 =
      sum(
        !is.na(
          CD19_celler_per_µL
        ) &
          `Måling nr.` > 1
      ),
    .groups = "drop"
  )

included_ids <- followup_count %>%
  filter(
    n_followup_cd19 >= 2
  ) %>%
  pull(
    PasientID
  )

excluded_ids <- followup_count %>%
  filter(
    n_followup_cd19 < 2
  ) %>%
  pull(
    PasientID
  )

cat(
  "\nExcluded participants:\n"
)
print(
  excluded_ids
)

df <- df_raw %>%
  filter(
    PasientID %in% included_ids
  )

# ============================================================
# TREATMENT LABELS
# ============================================================

df <- df %>%
  mutate(
    Behandling =
      case_when(
        Treatment == 1 ~ "Rituximab",
        Treatment == 2 ~ "Ocrelizumab",
        TRUE ~ NA_character_
      )
  )

# ============================================================
# PATIENT-LEVEL DATASET
# ============================================================

patients <- df %>%
  arrange(
    PasientID,
    `Måling nr.`
  ) %>%
  group_by(
    PasientID
  ) %>%
  slice(1) %>%
  ungroup()

patients$Alkohol_bruk <- ifelse(
  patients$Alkohol == "Avhold",
  0,
  1
)

# ============================================================
# SUSTAINED DEPLETION
# ============================================================

sustained <- df %>%
  filter(
    `Måling nr.` > 1
  ) %>%
  group_by(
    PasientID
  ) %>%
  summarise(
    sustained_depletion =
      ifelse(
        any(
          CD19_celler_per_µL > 5,
          na.rm = TRUE
        ),
        0,
        1
      ),
    .groups = "drop"
  )

# ============================================================
# EARLY DEPLETION (MONTH 3)
# ============================================================

early <- df %>%
  filter(
    `Måling nr.` == 2
  ) %>%
  transmute(
    PasientID,
    early_depletion =
      case_when(
        is.na(
          CD19_celler_per_µL
        ) ~ NA_real_,
        CD19_celler_per_µL <= 5 ~ 1,
        TRUE ~ 0
      )
  )

# ============================================================
# 80% DEPLETION ANALYSIS
# ============================================================

sens80 <- df %>%
  filter(
    `Måling nr.` > 1
  ) %>%
  group_by(
    PasientID
  ) %>%
  summarise(
    prop_depleted =
      mean(
        CD19_celler_per_µL[
          !is.na(
            CD19_celler_per_µL
          )
        ] <= 5
      ),
    depletion80 =
      ifelse(
        prop_depleted >= 0.80,
        1,
        0
      ),
    .groups = "drop"
  )

patients <- patients %>%
  left_join(
    sustained,
    by = "PasientID"
  ) %>%
  left_join(
    early,
    by = "PasientID"
  ) %>%
  left_join(
    sens80,
    by = "PasientID"
  )

# ============================================================
# CHECKS
# ============================================================

cat(
  "\nAnalysis cohort size:\n"
)
print(
  nrow(patients)
)

cat(
  "\nSustained depletion:\n"
)
print(
  table(
    patients$sustained_depletion
  )
)

cat(
  "\n80% depletion:\n"
)
print(
  table(
    patients$depletion80
  )
)

# ============================================================
# ANALYSIS VARIABLES
# ============================================================

categorical_vars <- c(
  "Kjønn",
  "Røyk",
  "Snus",
  "CRP_below_1",
  "Alkohol_bruk"
)

outcomes <- c(
  "sustained_depletion",
  "early_depletion",
  "depletion80"
)

median_iqr <- function(x){
  
  sprintf(
    "%.2f [%.2f-%.2f]",
    median(x, na.rm = TRUE),
    quantile(x, .25, na.rm = TRUE),
    quantile(x, .75, na.rm = TRUE)
  )
  
}

n_pct <- function(x, value){
  
  n <- sum(
    x == value,
    na.rm = TRUE
  )
  
  pct <- 100 * n / sum(!is.na(x))
  
  sprintf(
    "%d (%.1f%%)",
    n,
    pct
  )
  
}
# ============================================================
# DESCRIPTIVES
# ============================================================

make_descriptives <- function(
    outcome_var,
    group_labels = c("0","1")
){
  
  dat <- patients %>%
    filter(
      !is.na(.data[[outcome_var]])
    )
  
  continuous_table <- bind_rows(
    
    lapply(
      continuous_vars,
      
      function(v){
        
        g0 <- dat[[v]][dat[[outcome_var]] == 0]
        g1 <- dat[[v]][dat[[outcome_var]] == 1]
        
        data.frame(
          Variable = v,
          Group0 = median_iqr(g0),
          Group1 = median_iqr(g1)
        )
        
      }
      
    )
    
  )
  
  categorical_table <- bind_rows(
    
    lapply(
      
      categorical_vars,
      
      function(v){
        
        vals <- sort(unique(na.omit(dat[[v]])))
        
        bind_rows(
          
          lapply(
            
            vals,
            
            function(val){
              
              data.frame(
                Variable =
                  paste0(
                    v,
                    ": ",
                    val
                  ),
                
                Group0 =
                  n_pct(
                    dat[[v]][
                      dat[[outcome_var]] == 0
                    ],
                    val
                  ),
                
                Group1 =
                  n_pct(
                    dat[[v]][
                      dat[[outcome_var]] == 1
                    ],
                    val
                  )
              )
              
            }
            
          )
          
        )
        
      }
      
    )
    
  )
  
  bind_rows(
    continuous_table,
    categorical_table
  )
  
}

desc_sustained <- make_descriptives(
  "sustained_depletion"
)

n0 <- sum(
  patients$sustained_depletion == 0,
  na.rm = TRUE
)

n1 <- sum(
  patients$sustained_depletion == 1,
  na.rm = TRUE
)

desc_sustained <- rbind(
  data.frame(
    Variable = sprintf(
      "=== SUSTAINED DEPLETION (0=%d, 1=%d) ===",
      n0,
      n1
    ),
    Group0 = "",
    Group1 = ""
  ),
  desc_sustained
)

desc_early <- make_descriptives(
  "early_depletion"
)

n0 <- sum(
  patients$early_depletion == 0,
  na.rm = TRUE
)

n1 <- sum(
  patients$early_depletion == 1,
  na.rm = TRUE
)

desc_early <- rbind(
  data.frame(
    Variable = sprintf(
      "=== EARLY DEPLETION (0=%d, 1=%d) ===",
      n0,
      n1
    ),
    Group0 = "",
    Group1 = ""
  ),
  desc_early
)

desc_80 <- make_descriptives(
  "depletion80"
)

n0 <- sum(
  patients$depletion80 == 0,
  na.rm = TRUE
)

n1 <- sum(
  patients$depletion80 == 1,
  na.rm = TRUE
)

desc_80 <- rbind(
  data.frame(
    Variable = sprintf(
      "=== 80%% DEPLETION (0=%d, 1=%d) ===",
      n0,
      n1
    ),
    Group0 = "",
    Group1 = ""
  ),
  desc_80
)
descriptives <- bind_rows(
  desc_sustained,
  desc_early,
  desc_80
)


# ============================================================
# CHI-SQUARE FUNCTION
# ============================================================

run_chi_sq <- function(
    dep_var,
    pred_var
){
  
  dat <- patients %>%
    select(
      all_of(
        c(
          dep_var,
          pred_var
        )
      )
    ) %>%
    drop_na()
  
  tab <- table(
    dat[[dep_var]],
    dat[[pred_var]]
  )
  
  if(any(dim(tab) < 2)){
    return(NULL)
  }
  
  test <- chisq.test(
    tab
  )
  
  V <- cramerV(
    tab,
    bias.correct = TRUE
  )
  
  data.frame(
    outcome = dep_var,
    predictor = pred_var,
    statistic =
      unname(
        test$statistic
      ),
    df =
      unname(
        test$parameter
      ),
    p_value =
      test$p.value,
    cramer_v =
      V
  )
  
}

# ============================================================
# CHI-SQUARE ANALYSES
# ============================================================

chi_results <- bind_rows(
  
  lapply(
    outcomes,
    
    function(dep){
      
      bind_rows(
        
        lapply(
          
          categorical_vars,
          
          function(pred){
            
            run_chi_sq(
              dep,
              pred
            )
            
          }
          
        )
        
      )
      
    }
    
  )
  
)

chi_results$p_FDR <- p.adjust(
  chi_results$p_value,
  method = "BH"
)

chi_results <- chi_results %>%
  arrange(
    outcome,
    predictor
  )

# ============================================================
# WILCOXON FUNCTION
# ============================================================

run_wilcox <- function(
    dep_var,
    pred_var
){
  
  dat <- patients %>%
    select(
      all_of(
        c(
          dep_var,
          pred_var
        )
      )
    ) %>%
    drop_na()
  
  names(dat) <- c(
    "group",
    "value"
  )
  
  dat$group <- factor(
    dat$group
  )
  
  if(
    length(
      unique(dat$group)
    ) != 2
  ){
    return(NULL)
  }
  
  x0 <- dat$value[
    dat$group ==
      levels(dat$group)[1]
  ]
  
  x1 <- dat$value[
    dat$group ==
      levels(dat$group)[2]
  ]
  
  wt <- wilcox.test(
    value ~ group,
    data = dat,
    exact = FALSE
  )
  
  rb <- tryCatch(
    
    wilcox_effsize(
      dat,
      value ~ group,
      ci = TRUE
    ),
    
    error = function(e){
      NULL
    }
    
  )
  
  d_res <- tryCatch(
    
    effectsize::cohens_d(
      value ~ group,
      data = dat,
      ci = 0.95
    ),
    
    error = function(e){
      NULL
    }
    
  )
  
  data.frame(
    
    outcome = dep_var,
    predictor = pred_var,
    
    n_group0 = length(x0),
    n_group1 = length(x1),
    
    median_IQR_group0 =
      sprintf(
        "%.2f [%.2f-%.2f]",
        median(x0),
        quantile(x0, 0.25),
        quantile(x0, 0.75)
      ),
    
    median_IQR_group1 =
      sprintf(
        "%.2f [%.2f-%.2f]",
        median(x1),
        quantile(x1, 0.25),
        quantile(x1, 0.75)
      ),
    
    p_value =
      wt$p.value,
    
    rank_biserial_r =
      ifelse(
        is.null(rb),
        NA,
        rb$effsize
      ),
    
    rank_biserial_ci_low =
      ifelse(
        is.null(rb),
        NA,
        rb$conf.low
      ),
    
    rank_biserial_ci_high =
      ifelse(
        is.null(rb),
        NA,
        rb$conf.high
      ),
    
    cohens_d =
      ifelse(
        is.null(d_res),
        NA,
        d_res$Cohens_d
      ),
    
    cohens_d_ci_low =
      ifelse(
        is.null(d_res),
        NA,
        d_res$CI_low
      ),
    
    cohens_d_ci_high =
      ifelse(
        is.null(d_res),
        NA,
        d_res$CI_high
      )
    
  )
  
}

# ============================================================
# WILCOXON ANALYSES
# ============================================================

wilcox_results <- bind_rows(
  
  lapply(
    outcomes,
    
    function(dep){
      
      bind_rows(
        
        lapply(
          
          continuous_vars,
          
          function(pred){
            
            run_wilcox(
              dep,
              pred
            )
            
          }
          
        )
        
      )
      
    }
    
  )
  
)

wilcox_results$p_FDR <- p.adjust(
  wilcox_results$p_value,
  method = "BH"
)

wilcox_results <- wilcox_results %>%
  arrange(
    outcome,
    predictor
  )

# ============================================================
# EXPORT RESULTS
# ============================================================

write.xlsx(
  
  list(
    
    Descriptives =
      descriptives,
    
    ChiSquare =
      chi_results,
    
    Wilcoxon =
      wilcox_results
    
  ),
  
  file = paste0(
    path,
    "OVERLORD_updated_results.xlsx"
  ),
  
  overwrite = TRUE
  
)

# ============================================================
# FINAL CHECKS
# ============================================================

cat(
  "\nWorkbook written:\n"
)

cat(
  paste0(
    path,
    "OVERLORD_updated_results.xlsx"
  )
)

cat(
  "\n\nExpected cohort:\n"
)

cat(
  "n = 100\n"
)

cat(
  "Sustained depletion = 49\n"
)

cat(
  "No sustained depletion = 51\n"
)

# ============================================================
# END
# ============================================================