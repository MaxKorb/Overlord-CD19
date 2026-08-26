# ============================================================
# OVERLORD-MS B-CELL DEPLETION FIGURES
# ============================================================
#
# Laboratory sub-study of the OVERLORD-MS trial
#
# This script reproduces:
#   - Figure 1A: Study design, treatment regimen, and sampling schedule
#   - Figure 2: CD19+ B-cell counts over time
#   - Figure 3: Proportion of participants with B-cell depletion over time
#   - Figure 4: Median CD19+ B-cell counts over time by treatment
#   - Supplementary Figure 1: Longitudinal trajectories among participants
#     not achieving sustained depletion
#   - Supplementary Figure 2: Longitudinal trajectories among participants
#     not achieving early depletion
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(dplyr)
library(ggplot2)


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


# ============================================================
# 3. CREATE FIGURE OUTPUT DIRECTORY
# ============================================================

figure_output_dir <- file.path(
  "output",
  "figures"
)

if (!dir.exists(figure_output_dir)) {
  dir.create(
    figure_output_dir,
    recursive = TRUE
  )
}


# ============================================================
# 4. DEFINE FINAL ANALYSIS COHORT
# ============================================================

# Inclusion criterion:
# At least two available absolute CD19 measurements after baseline.

followup_count <- ny_tabell_raw %>%
  filter(`Måling nr.` > 1) %>%
  group_by(PasientID) %>%
  summarise(
    N_followup_CD19 = sum(
      !is.na(`CD19_celler_per_µL`)
    ),
    .groups = "drop"
  )


included_ids <- followup_count %>%
  filter(
    N_followup_CD19 >= 2
  ) %>%
  pull(PasientID)


excluded_participants <- followup_count %>%
  filter(
    N_followup_CD19 < 2
  )


ny_tabell <- ny_tabell_raw %>%
  filter(
    PasientID %in% included_ids
  )


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

n_distinct(
  ny_tabell_raw$PasientID
)

excluded_participants

n_distinct(
  ny_tabell$PasientID
)

dim(
  ny_tabell
)


# ============================================================
# 5. CREATE PARTICIPANT-LEVEL DATASET AND TREATMENT LABELS
# ============================================================

ny_tabell <- ny_tabell %>%
  mutate(
    Behandling = case_when(
      Treatment == 1 ~ "Rituximab",
      Treatment == 2 ~ "Ocrelizumab",
      TRUE ~ NA_character_
    )
  )


pasienter <- ny_tabell %>%
  distinct(
    PasientID,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

nrow(
  pasienter
)

table(
  pasienter$Behandling
)


# ============================================================
# 6. DEFINE SUSTAINED DEPLETION
# ============================================================

followup_data <- ny_tabell %>%
  filter(
    `Måling nr.` > 1
  )


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


pasienter_sustained <- pasienter %>%
  left_join(
    sustained_status,
    by = "PasientID"
  )


table(
  pasienter_sustained$Deplesjon_status
)


# ============================================================
# 7. DEFINE EARLY DEPLETION AT MONTH 3
# ============================================================

month3_data <- ny_tabell %>%
  filter(
    `Måling nr.` == 2
  )


early_status <- month3_data %>%
  transmute(
    PasientID,
    Deplesjon_status = case_when(
      is.na(`CD19_celler_per_µL`) ~ NA_character_,
      `CD19_celler_per_µL` <= 5 ~ "Depletert",
      `CD19_celler_per_µL` > 5 ~ "Ikke-depletert"
    )
  )


early_status_available <- early_status %>%
  filter(
    !is.na(Deplesjon_status)
  )


pasienter_early <- pasienter %>%
  inner_join(
    early_status_available,
    by = "PasientID"
  )


table(
  pasienter_early$Deplesjon_status
)

# ============================================================
# 8. FIGURE 1A. STUDY DESIGN, TREATMENT REGIMEN,
#    AND SAMPLING SCHEDULE
# ============================================================


# ------------------------------------------------------------
# 8.1 Time points and treatment schedules
# ------------------------------------------------------------

visit_times <- c(
  0, 3, 6, 12, 18, 24, 30, 36
)

ocrelizumab_times <- c(
  0, 6, 12, 18, 24, 30, 36
)

rituximab_times <- c(
  0, 6, 12, 18, 24, 30, 36
)


# ------------------------------------------------------------
# 8.2 Colors
# ------------------------------------------------------------

ocrelizumab_col <- "#2166AC"
rituximab_col <- "darkorange"
text_col <- "black"


# ------------------------------------------------------------
# 8.3 Create Figure 1A
# ------------------------------------------------------------

p1a <- ggplot() +
  
  # Vertical guide lines
  geom_segment(
    aes(
      x = visit_times,
      xend = visit_times,
      y = 0.65,
      yend = 3.25
    ),
    color = "grey90",
    linewidth = 0.6
  ) +
  
  # Ocrelizumab treatment line
  geom_segment(
    aes(
      x = 0,
      xend = 36,
      y = 3,
      yend = 3
    ),
    color = ocrelizumab_col,
    linewidth = 1.2
  ) +
  
  geom_point(
    aes(
      x = ocrelizumab_times,
      y = 3
    ),
    shape = 21,
    size = 3.0,
    stroke = 0.7,
    fill = ocrelizumab_col,
    color = ocrelizumab_col
  ) +
  
  # Rituximab treatment line
  geom_segment(
    aes(
      x = 0,
      xend = 36,
      y = 2,
      yend = 2
    ),
    color = rituximab_col,
    linewidth = 1.2
  ) +
  
  geom_point(
    aes(
      x = rituximab_times,
      y = 2
    ),
    shape = 21,
    size = 3.0,
    stroke = 0.7,
    fill = rituximab_col,
    color = rituximab_col
  ) +
  
  # Sampling schedule
  geom_segment(
    aes(
      x = 0,
      xend = 36,
      y = 1,
      yend = 1
    ),
    color = "black",
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(
      x = visit_times,
      y = 1
    ),
    shape = 21,
    size = 2.7,
    stroke = 0.7,
    fill = "white",
    color = "grey30"
  ) +
  
  # Left-side labels
  annotate(
    "text",
    x = -1.0,
    y = 3,
    label = "Ocrelizumab",
    hjust = 1,
    family = "serif",
    size = 4.5,
    color = text_col
  ) +
  
  annotate(
    "text",
    x = -1.0,
    y = 2,
    label = "Rituximab",
    hjust = 1,
    family = "serif",
    size = 4.5,
    color = text_col
  ) +
  
  annotate(
    "text",
    x = -1.0,
    y = 1.17,
    label = "Visits:",
    hjust = 1,
    family = "serif",
    size = 4.0,
    color = text_col
  ) +
  
  annotate(
    "text",
    x = -1.0,
    y = 0.98,
    label = "CD19+ B cells",
    hjust = 1,
    family = "serif",
    size = 4.0,
    color = text_col
  ) +
  
  annotate(
    "text",
    x = -1.0,
    y = 0.78,
    label = "and routine blood tests",
    hjust = 1,
    family = "serif",
    size = 4.0,
    color = text_col
  ) +
  
  # Dose annotations
  annotate(
    "text",
    x = 36.5,
    y = 3,
    label = "600 mg every 6 months",
    hjust = 0,
    family = "serif",
    size = 3.6,
    color = text_col
  ) +
  
  annotate(
    "text",
    x = 36.5,
    y = 2.12,
    label = "1000 mg at baseline",
    hjust = 0,
    family = "serif",
    size = 3.6,
    color = text_col
  ) +
  
  annotate(
    "text",
    x = 36.5,
    y = 1.88,
    label = "500 mg every 6 months",
    hjust = 0,
    family = "serif",
    size = 3.6,
    color = text_col
  ) +
  
  # Legend
  annotate(
    "text",
    x = 14.3,
    y = 3.90,
    label = "Treatment:",
    family = "serif",
    size = 4.0,
    color = text_col
  ) +
  
  geom_segment(
    aes(
      x = 19.3,
      xend = 20.7,
      y = 3.90,
      yend = 3.90
    ),
    color = ocrelizumab_col,
    linewidth = 1.2
  ) +
  
  geom_point(
    aes(
      x = 20.0,
      y = 3.90
    ),
    shape = 21,
    size = 2.7,
    stroke = 0.7,
    fill = ocrelizumab_col,
    color = ocrelizumab_col
  ) +
  
  annotate(
    "text",
    x = 21.2,
    y = 3.90,
    label = "Ocrelizumab",
    hjust = 0,
    family = "serif",
    size = 3.8,
    color = text_col
  ) +
  
  geom_segment(
    aes(
      x = 27.3,
      xend = 28.7,
      y = 3.90,
      yend = 3.90
    ),
    color = rituximab_col,
    linewidth = 1.2
  ) +
  
  geom_point(
    aes(
      x = 28.0,
      y = 3.90
    ),
    shape = 21,
    size = 2.7,
    stroke = 0.7,
    fill = rituximab_col,
    color = rituximab_col
  ) +
  
  annotate(
    "text",
    x = 29.2,
    y = 3.90,
    label = "Rituximab",
    hjust = 0,
    family = "serif",
    size = 3.8,
    color = text_col
  ) +
  
  scale_x_continuous(
    breaks = visit_times,
    labels = visit_times
  ) +
  
  coord_cartesian(
    xlim = c(-7, 43),
    ylim = c(0.45, 4.10),
    clip = "off"
  ) +
  
  labs(
    x = "Time (months)",
    y = NULL
  ) +
  
  theme_minimal(
    base_family = "serif",
    base_size = 14
  ) +
  
  theme(
    text = element_text(
      family = "serif",
      color = text_col
    ),
    axis.title.x = element_text(
      family = "serif",
      size = 14,
      color = text_col,
      margin = margin(t = 10)
    ),
    axis.text.x = element_text(
      family = "serif",
      size = 11,
      color = text_col
    ),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 50,
      b = 20,
      l = 60
    )
  )


# ------------------------------------------------------------
# 8.4 Display Figure 1A
# ------------------------------------------------------------

print(p1a)


# ------------------------------------------------------------
# 8.5 Export Figure 1A
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_output_dir,
    "Figure_1A_Study_Design.pdf"
  ),
  plot = p1a,
  width = 10,
  height = 4.5,
  units = "in"
)


ggsave(
  filename = file.path(
    figure_output_dir,
    "Figure_1A_Study_Design.tiff"
  ),
  plot = p1a,
  width = 10,
  height = 4.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 8.6 Confirm exports
# ------------------------------------------------------------

file.exists(
  file.path(
    figure_output_dir,
    "Figure_1A_Study_Design.pdf"
  )
)

file.exists(
  file.path(
    figure_output_dir,
    "Figure_1A_Study_Design.tiff"
  )
)

# ============================================================
# 9. FIGURE 2. CD19+ B-CELL COUNTS OVER TIME
# ============================================================


# ------------------------------------------------------------
# 9.1 Prepare data
# ------------------------------------------------------------

figure2_data <- ny_tabell %>%
  mutate(
    Time = case_when(
      `Måling nr.` == 1 ~ 0,
      `Måling nr.` == 2 ~ 3,
      `Måling nr.` == 3 ~ 6,
      `Måling nr.` == 4 ~ 12,
      `Måling nr.` == 5 ~ 18,
      `Måling nr.` == 6 ~ 24,
      `Måling nr.` == 7 ~ 30,
      `Måling nr.` == 8 ~ 36,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(Time),
    !is.na(`CD19_celler_per_µL`)
  )


# ------------------------------------------------------------
# 9.2 Check maximum CD19 value
# ------------------------------------------------------------

max(
  figure2_data$`CD19_celler_per_µL`,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 9.3 Create Figure 2
# ------------------------------------------------------------

p2 <- ggplot(
  figure2_data,
  aes(
    x = Time,
    y = `CD19_celler_per_µL`
  )
) +
  
  # Boxplots
  geom_boxplot(
    aes(
      group = factor(Time)
    ),
    width = 2.2,
    fill = "#FFC266",
    color = "darkorange",
    linewidth = 0.8,
    outlier.shape = NA
  ) +
  
  # Individual participant measurements
  geom_jitter(
    width = 0.6,
    height = 0,
    shape = 21,
    size = 2.0,
    alpha = 0.70,
    fill = "#7FAFD4",
    color = "#2166AC",
    stroke = 0.9
  ) +
  
  # Depletion threshold
  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey30"
  ) +
  
  scale_x_continuous(
    breaks = c(
      0, 3, 6, 12,
      18, 24, 30, 36
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(
      0,
      800,
      by = 100
    )
  ) +
  
  coord_cartesian(
    xlim = c(-1.5, 37.5),
    ylim = c(0, 800),
    expand = FALSE
  ) +
  
  labs(
    x = "Time (months)",
    y = "CD19+ B cells, cells/µL"
  ) +
  
  theme_minimal(
    base_family = "serif",
    base_size = 14
  ) +
  
  theme(
    text = element_text(
      family = "serif",
      color = "black"
    ),
    axis.title = element_text(
      family = "serif",
      color = "black",
      size = 14
    ),
    axis.text = element_text(
      family = "serif",
      color = "black",
      size = 11
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    axis.ticks = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 30,
      b = 25,
      l = 30
    )
  )


# ------------------------------------------------------------
# 9.4 Display Figure 2
# ------------------------------------------------------------

print(p2)


# ------------------------------------------------------------
# 9.5 Export Figure 2
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_output_dir,
    "Figure_2_CD19_B_cells.pdf"
  ),
  plot = p2,
  width = 8.5,
  height = 8,
  units = "in"
)


ggsave(
  filename = file.path(
    figure_output_dir,
    "Figure_2_CD19_B_cells.tiff"
  ),
  plot = p2,
  width = 8.5,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 9.6 Confirm exports
# ------------------------------------------------------------

file.exists(
  file.path(
    figure_output_dir,
    "Figure_2_CD19_B_cells.pdf"
  )
)

file.exists(
  file.path(
    figure_output_dir,
    "Figure_2_CD19_B_cells.tiff"
  )
)

# ============================================================
# 10. FIGURE 3. PROPORTION OF PARTICIPANTS WITH B-CELL
#     DEPLETION OVER TIME ACCORDING TO TREATMENT
# ============================================================


# ------------------------------------------------------------
# 10.1 Prepare data
# Baseline is excluded.
# Depletion = CD19+ B cells <= 5 cells/µL.
# ------------------------------------------------------------

figure3_data <- ny_tabell %>%
  filter(
    `Måling nr.` > 1,
    !is.na(`CD19_celler_per_µL`)
  ) %>%
  mutate(
    Treatment_group = case_when(
      Treatment == 1 ~ "Rituximab",
      Treatment == 2 ~ "Ocrelizumab",
      TRUE ~ NA_character_
    ),
    Time = case_when(
      `Måling nr.` == 2 ~ "3",
      `Måling nr.` == 3 ~ "6",
      `Måling nr.` == 4 ~ "12",
      `Måling nr.` == 5 ~ "18",
      `Måling nr.` == 6 ~ "24",
      `Måling nr.` == 7 ~ "30",
      `Måling nr.` == 8 ~ "36",
      TRUE ~ NA_character_
    ),
    Time = factor(
      Time,
      levels = c(
        "3", "6", "12",
        "18", "24", "30", "36"
      )
    ),
    Depleted = `CD19_celler_per_µL` <= 5
  ) %>%
  filter(
    !is.na(Time),
    !is.na(Treatment_group)
  )


# ------------------------------------------------------------
# 10.2 Exact binomial 95% confidence interval
# ------------------------------------------------------------

get_binom_ci <- function(
    x,
    n
) {
  
  ci <- binom.test(
    x = x,
    n = n,
    conf.level = 0.95
  )$conf.int
  
  data.frame(
    lower = ci[1] * 100,
    upper = ci[2] * 100
  )
}


# ------------------------------------------------------------
# 10.3 Treatment-specific proportions
# ------------------------------------------------------------

figure3_treatment <- figure3_data %>%
  group_by(
    Treatment_group,
    Time
  ) %>%
  summarise(
    N_measured = n(),
    N_depleted = sum(Depleted),
    Percent_depleted =
      100 * N_depleted / N_measured,
    .groups = "drop"
  )


# Add exact 95% confidence intervals
treatment_ci <- do.call(
  rbind,
  lapply(
    seq_len(
      nrow(figure3_treatment)
    ),
    function(i) {
      
      get_binom_ci(
        figure3_treatment$N_depleted[i],
        figure3_treatment$N_measured[i]
      )
    }
  )
)


figure3_treatment <- bind_cols(
  figure3_treatment,
  treatment_ci
)


# ------------------------------------------------------------
# 10.4 Overall proportions
# ------------------------------------------------------------

figure3_overall <- figure3_data %>%
  group_by(Time) %>%
  summarise(
    N_measured = n(),
    N_depleted = sum(Depleted),
    Percent_depleted =
      100 * N_depleted / N_measured,
    .groups = "drop"
  ) %>%
  mutate(
    Treatment_group = "Overall"
  )


# ------------------------------------------------------------
# 10.5 Checks
# ------------------------------------------------------------

figure3_treatment

figure3_overall


# ------------------------------------------------------------
# 10.6 X-axis labels
# Month on first line and number measured below
# ------------------------------------------------------------

x_labels <- setNames(
  paste0(
    figure3_overall$Time,
    "\n",
    "n = ",
    figure3_overall$N_measured
  ),
  figure3_overall$Time
)


# ------------------------------------------------------------
# 10.7 Colors
# ------------------------------------------------------------

group_colors <- c(
  "Ocrelizumab" = "#2166AC",
  "Rituximab" = "darkorange",
  "Overall" = "black"
)


# ------------------------------------------------------------
# 10.8 Create Figure 3
# ------------------------------------------------------------

p3 <- ggplot() +
  
  # Exact 95% confidence intervals
  geom_errorbar(
    data = figure3_treatment,
    aes(
      x = Time,
      ymin = lower,
      ymax = upper,
      color = Treatment_group
    ),
    width = 0.12,
    linewidth = 0.45,
    alpha = 0.80
  ) +
  
  # Treatment-specific lines
  geom_line(
    data = figure3_treatment,
    aes(
      x = Time,
      y = Percent_depleted,
      color = Treatment_group,
      group = Treatment_group
    ),
    linewidth = 1.0
  ) +
  
  # Treatment-specific points
  geom_point(
    data = figure3_treatment,
    aes(
      x = Time,
      y = Percent_depleted,
      color = Treatment_group
    ),
    size = 2.6
  ) +
  
  # Overall line
  geom_line(
    data = figure3_overall,
    aes(
      x = Time,
      y = Percent_depleted,
      color = Treatment_group,
      group = 1
    ),
    linewidth = 0.9
  ) +
  
  # Overall points
  geom_point(
    data = figure3_overall,
    aes(
      x = Time,
      y = Percent_depleted,
      color = Treatment_group
    ),
    shape = 21,
    size = 2.8,
    stroke = 0.8,
    fill = "white"
  ) +
  
  scale_color_manual(
    values = group_colors,
    breaks = c(
      "Ocrelizumab",
      "Rituximab",
      "Overall"
    )
  ) +
  
  scale_x_discrete(
    labels = x_labels
  ) +
  
  scale_y_continuous(
    breaks = seq(
      0,
      100,
      by = 20
    )
  ) +
  
  coord_cartesian(
    ylim = c(0, 105),
    expand = FALSE
  ) +
  
  labs(
    x = "Time (months)",
    y = "Participants with depletion (%)",
    color = NULL
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = c(1.0, 1.0, 0.9),
        shape = c(16, 16, 21),
        size = c(2.6, 2.6, 2.8),
        fill = c(NA, NA, "white")
      )
    )
  ) +
  
  theme_minimal(
    base_family = "serif",
    base_size = 14
  ) +
  
  theme(
    text = element_text(
      family = "serif",
      color = "black"
    ),
    axis.title = element_text(
      family = "serif",
      color = "black",
      size = 14
    ),
    axis.text = element_text(
      family = "serif",
      color = "black",
      size = 11
    ),
    axis.text.x = element_text(
      lineheight = 1.3,
      margin = margin(t = 6)
    ),
    axis.title.x = element_text(
      margin = margin(t = 10)
    ),
    legend.position = "top",
    legend.text = element_text(
      family = "serif",
      color = "black",
      size = 12
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    axis.ticks = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 25,
      b = 25,
      l = 25
    )
  )


# ------------------------------------------------------------
# 10.9 Display Figure 3
# ------------------------------------------------------------

print(p3)


# ------------------------------------------------------------
# 10.10 Export Figure 3
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_output_dir,
    "Figure_3_Depletion_Over_Time.pdf"
  ),
  plot = p3,
  width = 8.5,
  height = 6.5,
  units = "in"
)


ggsave(
  filename = file.path(
    figure_output_dir,
    "Figure_3_Depletion_Over_Time.tiff"
  ),
  plot = p3,
  width = 8.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 10.11 Confirm exports
# ------------------------------------------------------------

file.exists(
  file.path(
    figure_output_dir,
    "Figure_3_Depletion_Over_Time.pdf"
  )
)

file.exists(
  file.path(
    figure_output_dir,
    "Figure_3_Depletion_Over_Time.tiff"
  )
)

# ============================================================
# 11. FIGURE 4. MEDIAN CD19+ B-CELL COUNTS OVER TIME
#     BY TREATMENT GROUP
# ============================================================


# ------------------------------------------------------------
# 11.1 Prepare data
# Baseline is excluded.
# ------------------------------------------------------------

figure4_data <- ny_tabell %>%
  filter(
    `Måling nr.` > 1,
    !is.na(`CD19_celler_per_µL`)
  ) %>%
  mutate(
    Treatment_group = case_when(
      Treatment == 1 ~ "Rituximab",
      Treatment == 2 ~ "Ocrelizumab",
      TRUE ~ NA_character_
    ),
    Time = case_when(
      `Måling nr.` == 2 ~ 3,
      `Måling nr.` == 3 ~ 6,
      `Måling nr.` == 4 ~ 12,
      `Måling nr.` == 5 ~ 18,
      `Måling nr.` == 6 ~ 24,
      `Måling nr.` == 7 ~ 30,
      `Måling nr.` == 8 ~ 36,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(Time),
    !is.na(Treatment_group)
  ) %>%
  group_by(
    Treatment_group,
    Time
  ) %>%
  summarise(
    N = n(),
    Median = median(
      `CD19_celler_per_µL`,
      na.rm = TRUE
    ),
    Q1 = quantile(
      `CD19_celler_per_µL`,
      0.25,
      na.rm = TRUE
    ),
    Q3 = quantile(
      `CD19_celler_per_µL`,
      0.75,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 11.2 Check calculated values
# ------------------------------------------------------------

figure4_data


# ------------------------------------------------------------
# 11.3 Create equally spaced visit positions
# ------------------------------------------------------------

figure4_data <- figure4_data %>%
  mutate(
    Visit_position = case_when(
      Time == 3 ~ 1,
      Time == 6 ~ 2,
      Time == 12 ~ 3,
      Time == 18 ~ 4,
      Time == 24 ~ 5,
      Time == 30 ~ 6,
      Time == 36 ~ 7
    ),
    x_position = case_when(
      Treatment_group == "Ocrelizumab" ~ Visit_position - 0.08,
      Treatment_group == "Rituximab" ~ Visit_position + 0.08
    )
  )


# ------------------------------------------------------------
# 11.4 Check plotting positions
# ------------------------------------------------------------

figure4_data %>%
  select(
    Treatment_group,
    Time,
    Visit_position,
    x_position,
    Median,
    Q1,
    Q3
  )


# ------------------------------------------------------------
# 11.5 Colors
# ------------------------------------------------------------

ocrelizumab_col <- "#2166AC"
rituximab_col <- "darkorange"


# ------------------------------------------------------------
# 11.6 Create Figure 4
# ------------------------------------------------------------

p4 <- ggplot(
  figure4_data,
  aes(
    x = x_position,
    y = Median,
    color = Treatment_group,
    group = Treatment_group
  )
) +
  
  # IQR bars
  geom_errorbar(
    aes(
      ymin = Q1,
      ymax = Q3
    ),
    width = 0.10,
    linewidth = 0.5,
    alpha = 0.90
  ) +
  
  # Depletion threshold
  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey30"
  ) +
  
  # Median lines
  geom_line(
    linewidth = 1.0
  ) +
  
  # Median points
  geom_point(
    size = 2.6
  ) +
  
  # Colors
  scale_color_manual(
    values = c(
      "Ocrelizumab" = ocrelizumab_col,
      "Rituximab" = rituximab_col
    ),
    breaks = c(
      "Ocrelizumab",
      "Rituximab"
    )
  ) +
  
  # X-axis
  scale_x_continuous(
    breaks = 1:7,
    labels = c(
      "3",
      "6",
      "12",
      "18",
      "24",
      "30",
      "36"
    )
  ) +
  
  # Y-axis
  scale_y_continuous(
    breaks = seq(
      0,
      25,
      by = 5
    )
  ) +
  
  coord_cartesian(
    xlim = c(0.75, 7.25),
    ylim = c(-0.6, 25),
    expand = FALSE
  ) +
  
  labs(
    x = "Time (months)",
    y = "Median CD19+ B cells, cells/µL",
    color = NULL
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 1.0,
        size = 2.6
      )
    )
  ) +
  
  theme_minimal(
    base_family = "serif",
    base_size = 14
  ) +
  
  theme(
    text = element_text(
      family = "serif",
      color = "black"
    ),
    axis.title = element_text(
      family = "serif",
      color = "black",
      size = 14
    ),
    axis.text = element_text(
      family = "serif",
      color = "black",
      size = 11
    ),
    legend.position = "top",
    legend.text = element_text(
      family = "serif",
      color = "black",
      size = 12
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    axis.ticks = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 25,
      b = 20,
      l = 25
    )
  )


# ------------------------------------------------------------
# 11.7 Display Figure 4
# ------------------------------------------------------------

print(p4)


# ------------------------------------------------------------
# 11.8 Export Figure 4
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_output_dir,
    "Figure_4_CD19_Median_IQR_by_Treatment.pdf"
  ),
  plot = p4,
  width = 8.5,
  height = 6.5,
  units = "in"
)


ggsave(
  filename = file.path(
    figure_output_dir,
    "Figure_4_CD19_Median_IQR_by_Treatment.tiff"
  ),
  plot = p4,
  width = 8.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 11.9 Confirm exports
# ------------------------------------------------------------

file.exists(
  file.path(
    figure_output_dir,
    "Figure_4_CD19_Median_IQR_by_Treatment.pdf"
  )
)

file.exists(
  file.path(
    figure_output_dir,
    "Figure_4_CD19_Median_IQR_by_Treatment.tiff"
  )
)

# ============================================================
# 12. SUPPLEMENTARY FIGURE 1.
#     LONGITUDINAL TRAJECTORIES OF CD19+ B-CELL COUNTS
#     AMONG PARTICIPANTS NOT ACHIEVING SUSTAINED DEPLETION
# ============================================================


# ------------------------------------------------------------
# 12.1 Identify participants not achieving sustained depletion
# ------------------------------------------------------------

non_sustained_ids <- pasienter_sustained %>%
  filter(
    Deplesjon_status == "Ikke-deplesjon"
  ) %>%
  pull(PasientID)


# ------------------------------------------------------------
# 12.2 Prepare longitudinal data
# ------------------------------------------------------------

supp_fig1_data <- ny_tabell %>%
  filter(
    PasientID %in% non_sustained_ids,
    !is.na(`CD19_celler_per_µL`)
  ) %>%
  mutate(
    Time = case_when(
      `Måling nr.` == 1 ~ 0,
      `Måling nr.` == 2 ~ 3,
      `Måling nr.` == 3 ~ 6,
      `Måling nr.` == 4 ~ 12,
      `Måling nr.` == 5 ~ 18,
      `Måling nr.` == 6 ~ 24,
      `Måling nr.` == 7 ~ 30,
      `Måling nr.` == 8 ~ 36,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(Time)
  )


# ------------------------------------------------------------
# 12.3 Calculate mean and SD at each time point
# ------------------------------------------------------------

supp_fig1_summary <- supp_fig1_data %>%
  group_by(Time) %>%
  summarise(
    N = n(),
    
    Mean = mean(
      `CD19_celler_per_µL`,
      na.rm = TRUE
    ),
    
    SD = sd(
      `CD19_celler_per_µL`,
      na.rm = TRUE
    ),
    
    Lower = pmax(
      Mean - SD,
      0
    ),
    
    Upper = Mean + SD,
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 12.4 Checks
# ------------------------------------------------------------

supp_fig1_summary

length(
  unique(supp_fig1_data$PasientID)
)

max(
  supp_fig1_data$`CD19_celler_per_µL`,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 12.5 Colors
# ------------------------------------------------------------

individual_col <- "#2166AC"
mean_col <- "darkorange"
mean_fill <- "#FFC266"


# ------------------------------------------------------------
# 12.6 Create Supplementary Figure 1
# ------------------------------------------------------------

supp_fig1 <- ggplot() +
  
  # Mean ± SD shaded area
  geom_ribbon(
    data = supp_fig1_summary,
    aes(
      x = Time,
      ymin = Lower,
      ymax = Upper
    ),
    fill = mean_fill,
    alpha = 0.35
  ) +
  
  # Individual participant trajectories
  geom_line(
    data = supp_fig1_data,
    aes(
      x = Time,
      y = `CD19_celler_per_µL`,
      group = PasientID
    ),
    color = individual_col,
    linewidth = 0.55,
    alpha = 0.55
  ) +
  
  # Individual participant measurements
  geom_point(
    data = supp_fig1_data,
    aes(
      x = Time,
      y = `CD19_celler_per_µL`
    ),
    color = individual_col,
    size = 1.5,
    alpha = 0.65
  ) +
  
  # Cohort mean
  geom_line(
    data = supp_fig1_summary,
    aes(
      x = Time,
      y = Mean
    ),
    color = mean_col,
    linewidth = 1.2
  ) +
  
  # Depletion threshold
  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey30"
  ) +
  
  scale_x_continuous(
    breaks = c(
      0, 3, 6, 12,
      18, 24, 30, 36
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(
      0,
      800,
      by = 100
    )
  ) +
  
  coord_cartesian(
    xlim = c(-1, 37),
    ylim = c(0, 800),
    expand = FALSE
  ) +
  
  labs(
    x = "Time (months)",
    y = "CD19+ B cells, cells/µL"
  ) +
  
  theme_minimal(
    base_family = "serif",
    base_size = 14
  ) +
  
  theme(
    text = element_text(
      family = "serif",
      color = "black"
    ),
    axis.title = element_text(
      family = "serif",
      color = "black",
      size = 14
    ),
    axis.text = element_text(
      family = "serif",
      color = "black",
      size = 11
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    axis.ticks = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 25,
      b = 20,
      l = 25
    )
  )


# ------------------------------------------------------------
# 12.7 Display Supplementary Figure 1
# ------------------------------------------------------------

print(supp_fig1)


# ------------------------------------------------------------
# 12.8 Export Supplementary Figure 1
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_output_dir,
    "Supplementary_Figure_1_Non_Sustained_Depletion.pdf"
  ),
  plot = supp_fig1,
  width = 8.5,
  height = 8,
  units = "in"
)


ggsave(
  filename = file.path(
    figure_output_dir,
    "Supplementary_Figure_1_Non_Sustained_Depletion.tiff"
  ),
  plot = supp_fig1,
  width = 8.5,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 12.9 Confirm exports
# ------------------------------------------------------------

file.exists(
  file.path(
    figure_output_dir,
    "Supplementary_Figure_1_Non_Sustained_Depletion.pdf"
  )
)

file.exists(
  file.path(
    figure_output_dir,
    "Supplementary_Figure_1_Non_Sustained_Depletion.tiff"
  )
)

# ============================================================
# 13. SUPPLEMENTARY FIGURE 2.
#     LONGITUDINAL TRAJECTORIES OF CD19+ B-CELL COUNTS
#     AMONG PARTICIPANTS NOT ACHIEVING EARLY DEPLETION
# ============================================================


# ------------------------------------------------------------
# 13.1 Identify participants not achieving early depletion
# ------------------------------------------------------------

non_early_ids <- pasienter_early %>%
  filter(
    Deplesjon_status == "Ikke-depletert"
  ) %>%
  pull(PasientID)


# ------------------------------------------------------------
# 13.2 Prepare longitudinal data
# ------------------------------------------------------------

supp_fig2_data <- ny_tabell %>%
  filter(
    PasientID %in% non_early_ids,
    !is.na(`CD19_celler_per_µL`)
  ) %>%
  mutate(
    Time = case_when(
      `Måling nr.` == 1 ~ 0,
      `Måling nr.` == 2 ~ 3,
      `Måling nr.` == 3 ~ 6,
      `Måling nr.` == 4 ~ 12,
      `Måling nr.` == 5 ~ 18,
      `Måling nr.` == 6 ~ 24,
      `Måling nr.` == 7 ~ 30,
      `Måling nr.` == 8 ~ 36,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(Time)
  )


# ------------------------------------------------------------
# 13.3 Calculate mean and SD at each time point
# ------------------------------------------------------------

supp_fig2_summary <- supp_fig2_data %>%
  group_by(Time) %>%
  summarise(
    N = n(),
    
    Mean = mean(
      `CD19_celler_per_µL`,
      na.rm = TRUE
    ),
    
    SD = sd(
      `CD19_celler_per_µL`,
      na.rm = TRUE
    ),
    
    Lower = pmax(
      Mean - SD,
      0
    ),
    
    Upper = Mean + SD,
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 13.4 Checks
# ------------------------------------------------------------

supp_fig2_summary

length(
  unique(supp_fig2_data$PasientID)
)

max(
  supp_fig2_data$`CD19_celler_per_µL`,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 13.5 Colors
# ------------------------------------------------------------

individual_col <- "#2166AC"
mean_col <- "darkorange"
mean_fill <- "#FFC266"


# ------------------------------------------------------------
# 13.6 Create Supplementary Figure 2
# ------------------------------------------------------------

supp_fig2 <- ggplot() +
  
  # Mean ± SD shaded area
  geom_ribbon(
    data = supp_fig2_summary,
    aes(
      x = Time,
      ymin = Lower,
      ymax = Upper
    ),
    fill = mean_fill,
    alpha = 0.35
  ) +
  
  # Individual participant trajectories
  geom_line(
    data = supp_fig2_data,
    aes(
      x = Time,
      y = `CD19_celler_per_µL`,
      group = PasientID
    ),
    color = individual_col,
    linewidth = 0.55,
    alpha = 0.55
  ) +
  
  # Individual participant measurements
  geom_point(
    data = supp_fig2_data,
    aes(
      x = Time,
      y = `CD19_celler_per_µL`
    ),
    color = individual_col,
    size = 1.5,
    alpha = 0.65
  ) +
  
  # Cohort mean
  geom_line(
    data = supp_fig2_summary,
    aes(
      x = Time,
      y = Mean
    ),
    color = mean_col,
    linewidth = 1.2
  ) +
  
  # Depletion threshold
  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey30"
  ) +
  
  scale_x_continuous(
    breaks = c(
      0, 3, 6, 12,
      18, 24, 30, 36
    )
  ) +
  
  scale_y_continuous(
    breaks = seq(
      0,
      500,
      by = 100
    )
  ) +
  
  coord_cartesian(
    xlim = c(-1, 37),
    ylim = c(0, 500),
    expand = FALSE
  ) +
  
  labs(
    x = "Time (months)",
    y = "CD19+ B cells, cells/µL"
  ) +
  
  theme_minimal(
    base_family = "serif",
    base_size = 14
  ) +
  
  theme(
    text = element_text(
      family = "serif",
      color = "black"
    ),
    axis.title = element_text(
      family = "serif",
      color = "black",
      size = 14
    ),
    axis.text = element_text(
      family = "serif",
      color = "black",
      size = 11
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(
      color = "grey85",
      linewidth = 0.4
    ),
    axis.ticks = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 25,
      b = 20,
      l = 25
    )
  )


# ------------------------------------------------------------
# 13.7 Display Supplementary Figure 2
# ------------------------------------------------------------

print(supp_fig2)


# ------------------------------------------------------------
# 13.8 Export Supplementary Figure 2
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figure_output_dir,
    "Supplementary_Figure_2_Early_Non_Depletion.pdf"
  ),
  plot = supp_fig2,
  width = 8.5,
  height = 8,
  units = "in"
)


ggsave(
  filename = file.path(
    figure_output_dir,
    "Supplementary_Figure_2_Early_Non_Depletion.tiff"
  ),
  plot = supp_fig2,
  width = 8.5,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


# ------------------------------------------------------------
# 13.9 Confirm exports
# ------------------------------------------------------------

file.exists(
  file.path(
    figure_output_dir,
    "Supplementary_Figure_2_Early_Non_Depletion.pdf"
  )
)

file.exists(
  file.path(
    figure_output_dir,
    "Supplementary_Figure_2_Early_Non_Depletion.tiff"
  )
)
