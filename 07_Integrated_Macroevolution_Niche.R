# ==============================================================================
# Script Name: 07_Integrated_Macroevolution_Niche.R
# Purpose: Final visualization pipeline for Manuscript.
#          Focused analysis for Brachiopoda (Taxonomic, Niche, Paleolatitude).
# ==============================================================================

rm(list = ls())
library(tidyverse)
library(ggridges)
library(here)
library(stringr)
library(patchwork)

setwd(here::here())
theme_set(theme_classic(base_size = 14))

# ---- 1. Data Setup & Filter for Brachiopoda ----
# Load Brachiopoda analysis data directly
f <- "Brachiopoda_analysis_data.rds"

if(file.exists(f)) {
  occurrence_data <- readRDS(f) %>% 
    mutate(Clade = "Brachiopoda") %>%
    mutate(mid_age = (as.numeric(max_ma) + as.numeric(min_ma)) / 2)
} else {
  stop("File 'Brachiopoda_analysis_data.rds' not found. Please ensure the file exists in the working directory.")
}

# Combine and restrict to Brachiopoda only
occurrence_data <- bind_rows(occ_list) %>%
  filter(Clade == "Brachiopoda") %>%
  mutate(mid_age = (as.numeric(max_ma) + as.numeric(min_ma)) / 2)

# Define geological periods for niche & latitude figures
period_boundaries <- tribble(
  ~Period,         ~min_ma,  ~max_ma,
  "Quaternary",    0,        2.58,
  "Neogene",       2.58,     23.03,
  "Paleogene",     23.03,    66.0,
  "Cretaceous",    66.0,     145.0,
  "Jurassic",      145.0,    201.4,
  "Triassic",      201.4,    251.9
) %>%
  mutate(Period = factor(Period, levels = rev(.$Period)))

# Assign Period to data
assign_period <- function(mid_age) {
  case_when(
    mid_age >= 0   & mid_age < 2.58  ~ "Quaternary",
    mid_age >= 2.58  & mid_age < 23.03 ~ "Neogene",
    mid_age >= 23.03 & mid_age < 66.0  ~ "Paleogene",
    mid_age >= 66.0  & mid_age < 145.0 ~ "Cretaceous",
    mid_age >= 145.0 & mid_age < 201.4 ~ "Jurassic",
    mid_age >= 201.4 & mid_age <= 251.9 ~ "Triassic",
    TRUE ~ NA_character_
  )
}
occurrence_data <- occurrence_data %>%
  mutate(Period = assign_period(mid_age)) %>%
  filter(!is.na(Period)) %>%
  mutate(Period = factor(Period, levels = c("Triassic","Jurassic","Cretaceous","Paleogene","Neogene","Quaternary")))

# ---- 2. Fig 1: Taxonomic Scale Dependence (Species, Genus, Family Richness) ----

bin_breaks <- seq(0, 260, by = 10)
bin_labels <- paste0(bin_breaks[-length(bin_breaks)], "-", bin_breaks[-1])

tax_data <- occurrence_data %>%
  filter(!is.na(genus)) %>%
  mutate(
    bin = cut(mid_age, breaks = bin_breaks, labels = bin_labels, right = FALSE),
    species = if_else(str_detect(accepted_name, " "), accepted_name, NA_character_),
    family = if_else(!is.na(family), family, NA_character_)
  ) %>%
  filter(!is.na(bin))

tax_richness <- tax_data %>%
  group_by(bin) %>%
  summarise(
    n_species = n_distinct(species, na.rm = TRUE),
    n_genus   = n_distinct(genus, na.rm = TRUE),
    n_family  = n_distinct(family, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("n_"), names_to = "Level", values_to = "Richness") %>%
  mutate(
    Level = recode(Level, n_species = "Species", n_genus = "Genus", n_family = "Family"),
    Level = factor(Level, levels = c("Species", "Genus", "Family")),
    bin_mid = as.numeric(str_extract(bin, "^[0-9]+")) + 5
  )

p_richness <- ggplot(tax_richness, aes(x = bin_mid, y = Richness, linetype = Level)) +
  geom_line(linewidth = 1.0, color = "#E41A1C") +
  scale_x_reverse(limits = c(252, 0), breaks = seq(250, 0, by = -50)) +
  labs(x = "Time (Ma)", y = "Raw richness", linetype = "Level") +
  theme(legend.position = "top", legend.box = "vertical")

ratio_data <- tax_data %>%
  group_by(bin) %>%
  summarise(
    n_species = n_distinct(species, na.rm = TRUE),
    n_genus   = n_distinct(genus, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ratio = n_species / n_genus,
    bin_mid = as.numeric(str_extract(bin, "^[0-9]+")) + 5
  )

p_ratio <- ggplot(ratio_data, aes(x = bin_mid, y = ratio)) +
  geom_line(linewidth = 1.0, color = "#E41A1C") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_x_reverse(limits = c(252, 0), breaks = seq(250, 0, by = -50)) +
  labs(x = "Time (Ma)", y = "Species/Genus ratio")

p_taxonomic <- p_richness / p_ratio + plot_layout(heights = c(2, 1))

# ---- 3. Fig 2: Ecological Niche Shift (REVISED) ----
classify_env <- function(env) {
  env <- tolower(env)
  case_when(
    # 1. Deep Marine / Basinal Niche (The modern refuge)
    str_detect(env, "deep|slope|basin|fan|abyss|bathyal") ~ "Deep Marine",
    
    # 2. Carbonate / Reef Niche (The Paleozoic stronghold)
    str_detect(env, "carbonate|reef|buildup|bioherm|shoal|lagoonal|perireef|backreef|subreef") ~ "Carbonate/Reef",
    
    # 3. Shallow / Open Shelf Niche (Broadest, normal marine category)
    # Includes generic 'marine', 'subtidal', 'offshore', 'shoreface', 'coastal'
    str_detect(env, "marine|subtidal|offshore|shelf|ramp|shore|coastal|peritidal") ~ "Shallow/Open Shelf",
    
    # Exclude any unmapped artifacts (though prior filtering should have caught non-marine)
    TRUE ~ "Other"
  )
}

occ_eco <- occurrence_data %>%
  mutate(eco_niche = classify_env(environment)) %>%
  # Exclude 'Other' just in case, ensuring we only plot the big three
  filter(eco_niche != "Other") %>% 
  # Set factor levels to ensure consistent color ordering in the plot
  mutate(eco_niche = factor(eco_niche, levels = c("Shallow/Open Shelf", "Carbonate/Reef", "Deep Marine")))

niche_summary <- occ_eco %>%
  count(Period, eco_niche) %>%
  group_by(Period) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_niche <- ggplot(niche_summary, aes(x = Period, y = prop, fill = eco_niche)) +
  geom_col(position = "fill", color = "white", linewidth = 0.2) +
  # Using a distinct 3-color palette suited for these environments
  scale_fill_manual(values = c("Shallow/Open Shelf" = "#66c2a5", 
                               "Carbonate/Reef"     = "#fc8d62", 
                               "Deep Marine"        = "#8da0cb"), 
                    name = "Ecological Niche") +
  labs(
    title = "Brachiopoda Ecological Niche Shift (Triassic to Quaternary)",
    y = "Proportion",
    x = ""
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 4. Fig 3: Paleolatitudinal Distribution ----
lat_data <- occurrence_data %>%
  filter(!is.na(paleolat), !is.na(Period))

p_lat <- ggplot(lat_data, aes(x = paleolat, y = Period)) +
  geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01,
                      quantile_lines = TRUE, quantiles = 2, fill = "#E41A1C") +
  labs(
    title = "Brachiopoda Paleolatitudinal Distribution by Period",
    x = "Paleolatitude (°)",
    y = ""
  )

# ---- 5. Save Outputs ----
if (!dir.exists("Figures")) dir.create("Figures")

ggsave("Figures/Fig3_Taxonomic_Scale.jpg", p_taxonomic, width = 8, height = 7, dpi = 600, device = "jpeg")
ggsave("Figures/Fig4_Niche.jpg", p_niche, width = 8, height = 5, dpi = 600, device = "jpeg")
ggsave("Figures/Fig5_Latitude.jpg", p_lat, width = 8, height = 6, dpi = 600, device = "jpeg")

message("Figures saved successfully.")
