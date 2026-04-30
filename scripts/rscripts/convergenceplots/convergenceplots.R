#install.packages(c("tidyverse", "shiny", "leaflet", "DT", "hash"))
library(shiny)
library(DT)
library(tidyverse)
library(dplyr)

#---------------------------------------------------------------#

IBMblue = "#648FFF"
IBMpurple = "#785EDC"
IBMpink = "#DC267F"
IBMorange = "#FE6100"
IBMyellow = "#FFB000"

#---------------------------------------------------------------#

#Read in data
magdata <- read_csv('convergence_expex4_exp5_mag_rundate2026-03-23.csv')
sagdata <- read_csv('convergence_expex4_exp29_sag_rundate2026-03-23.csv')
pbcgdata <- read_csv('convergence_expex4_exp21_cg_rundate2026-04-29.csv')

#---------------------------------------------------------------#

df_all <- bind_rows(magdata, sagdata, pbcgdata)
df_all <- df_all %>%
  mutate(method = recode(method, mag = "MAG"))
df_all <- df_all %>%
  mutate(method = recode(method, sag = "SAG"))
df_all <- df_all %>%
  mutate(method = recode(method, cg = "PBCG"))

#---------------------------------------------------------------#

df_long <- df_all %>%
  pivot_longer(
    cols = c(lowerbound, upperbound),
    names_to = "bound",
    values_to = "value"
  ) %>%
  mutate(
    bound = recode(bound,
                        lowerbound = "lower",
                        upperbound = "upper")
  )

df <- df_long %>%
  select(method, iteration, bound, value)

#---------------------------------------------------------------#

ggplot(df, aes(x = iteration, y = value,
               color = method,
               linetype = bound)) +
  geom_line(size = 1) +
  
  scale_x_log10() +
  
  scale_color_manual(values = c(
    "MAG" = IBMorange,   # orange
    "PBCG" = IBMblue,  # dark blue
    "SAG" = IBMpink    # pink
  )) +
  
  scale_linetype_manual(values = c(
    "upper" = "solid",
    "lower" = "dashed"
  )) +
  
  labs(
    x = "Iteration (log scale)",
    y = "Objective value",
    color = NULL,
    linetype = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

#---------------------------------------------------------------#

ggplot(df_all, aes(x = iteration)) +
  
  # --- PATHS (shaded areas) ---
  geom_area(aes(y = path_count, fill = method),
            alpha = 0.35, position = "identity") +
  
  # --- ARCS (lines) ---
  geom_line(aes(y = arc_count, color = method),
            size = 1.2) +
  
  # log scale on x-axis
  scale_x_log10() +
  
  # colors (match your figure)
  scale_color_manual(values = c(
    "MAG"  = IBMorange,
    "PBCG" = IBMblue,
    "SAG"  = IBMpink
  )) +
  
  scale_fill_manual(values = c(
    "MAG"  = IBMorange,
    "PBCG" = IBMblue,
    "SAG"  = IBMpink
  )) +
  
  labs(
    x = "Iteration (log scale)",
    y = "Count",
    color = NULL,
    fill = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

