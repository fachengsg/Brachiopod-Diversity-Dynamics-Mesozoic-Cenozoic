# ==============================================================================
# Script Name: 07_Integrated_Macroevolution_Niche.R
# Purpose: Reproducible visualization pipeline for brachiopod post-Paleozoic
#          evolution. Generates Figures 3-6 (taxonomic scale, ecological niches,
#          paleolatitudinal distribution, and paleogeographic maps) plus
#          supplementary family-level summaries.
# ==============================================================================

# Clear the workspace. Comment out this line if you want to keep existing objects.
rm(list = ls())

# ---- 1. Environment Setup & Packages ----
required_packages <- c("dplyr", "readr", "stringr", "tidyr", "ggplot2",
                       "ggridges", "patchwork", "sf", "cowplot")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(sf)
library(cowplot)

# --- Working directory: choose one of the two options below ---
# Option A: Set a fixed directory (e.g., "D:/PBDB_Project" on Windows).
# Option B: Leave custom_dir as NULL or "" to automatically use the
#           folder that contains this script (works in RStudio).
custom_dir <- "D:/PBDB_Project"   # <-- Change this to the local project folder, or set to NULL

if (!is.null(custom_dir) && nchar(custom_dir) > 0) {
  if (!dir.exists(custom_dir)) {
    dir.create(custom_dir, recursive = TRUE)
    cat("Folder created:", custom_dir, "\n")
  }
  setwd(custom_dir)
} else {
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
    if (nchar(script_path) > 0) {
      setwd(dirname(script_path))
    }
  }
}
cat("Working directory:", getwd(), "\n")

# Global ggplot2 theme
theme_set(theme_classic(base_size = 14))

# ---- 2. Load and prepare data ----
if (!file.exists("Brachiopoda_analysis_data.rds")) {
  stop("Input file 'Brachiopoda_analysis_data.rds' not found. Please run Scripts 01-02 first.")
}

occurrence_data <- readRDS("Brachiopoda_analysis_data.rds") %>%
  mutate(
    Clade = "Brachiopoda",
    mid_age = (as.numeric(max_ma) + as.numeric(min_ma)) / 2
  )

# Check that paleocoordinates exist (required for Figures 5 and 6)
if (!all(c("paleolat", "paleolng") %in% names(occurrence_data))) {
  stop("Columns 'paleolat' and 'paleolng' are required. Please reconstruct paleocoordinates before running this script.")
}

# ---- 3. Assign geological periods ----
assign_period <- function(mid_age) {
  case_when(
    mid_age >= 0      & mid_age < 2.58  ~ "Quaternary",
    mid_age >= 2.58   & mid_age < 23.03 ~ "Neogene",
    mid_age >= 23.03  & mid_age < 66.0  ~ "Paleogene",
    mid_age >= 66.0   & mid_age < 145.0 ~ "Cretaceous",
    mid_age >= 145.0  & mid_age < 201.4 ~ "Jurassic",
    mid_age >= 201.4  & mid_age <= 251.9 ~ "Triassic",
    TRUE ~ NA_character_
  )
}

occurrence_data <- occurrence_data %>%
  mutate(Period = assign_period(mid_age)) %>%
  filter(!is.na(Period)) %>%
  mutate(Period = factor(Period, levels = c("Triassic", "Jurassic", "Cretaceous",
                                            "Paleogene", "Neogene", "Quaternary")))

# ---- 4. Fig 3: Taxonomic scale dependence ----
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

# Richness across hierarchical levels
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

# Species-to-genus ratio
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

# Combine into Fig 3
p_taxonomic <- p_richness / p_ratio + plot_layout(heights = c(2, 1))

# ---- 5. Fig 4: Ecological niche shift ----
classify_env <- function(env) {
  env <- tolower(env)
  case_when(
    str_detect(env, "deep|slope|basin|fan|abyss|bathyal") ~ "Deep Marine",
    str_detect(env, "carbonate|reef|buildup|bioherm|shoal|lagoonal|perireef|backreef|subreef") ~ "Carbonate/Reef",
    str_detect(env, "marine|subtidal|offshore|shelf|ramp|shore|coastal|peritidal") ~ "Shallow/Open Shelf",
    TRUE ~ "Other"
  )
}

occ_eco <- occurrence_data %>%
  mutate(eco_niche = classify_env(environment)) %>%
  filter(eco_niche != "Other") %>%
  mutate(eco_niche = factor(eco_niche, levels = c("Shallow/Open Shelf", "Carbonate/Reef", "Deep Marine")))

niche_summary <- occ_eco %>%
  count(Period, eco_niche) %>%
  group_by(Period) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_niche <- ggplot(niche_summary, aes(x = Period, y = prop, fill = eco_niche)) +
  geom_col(position = "fill", color = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("Shallow/Open Shelf" = "#66c2a5",
               "Carbonate/Reef"     = "#fc8d62",
               "Deep Marine"        = "#8da0cb"),
    name = "Ecological Niche"
  ) +
  labs(
    title = "Brachiopoda Ecological Niche Shift (Triassic to Quaternary)",
    y = "Proportion", x = ""
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 6. Fig 5: Paleolatitudinal distribution ----
lat_data <- occurrence_data %>%
  filter(!is.na(paleolat), !is.na(Period))

p_lat <- ggplot(lat_data, aes(x = paleolat, y = Period)) +
  geom_density_ridges(
    alpha = 0.7, scale = 1.2, rel_min_height = 0.01,
    quantile_lines = TRUE, quantiles = 2, fill = "#E41A1C"
  ) +
  scale_x_continuous(
    labels = function(x) ifelse(x < 0, paste0("\u2212", abs(x)), as.character(x))
  ) +
  labs(
    title = "Brachiopoda Paleolatitudinal Distribution by Period",
    x = "Paleolatitude (°)", y = ""
  )

# ---- 7. Fig 6: Paleogeographic maps ----
# Requires PALEOMAP paleocoastlines (v7, Scotese & Wright).
# Download from https://zenodo.org/records/5469129 and extract the 'CS' folder
# into './paleocoastlines/CS/' within the working directory.
coast_dir <- file.path(getwd(), "paleocoastlines", "CS")

if (!dir.exists(coast_dir)) {
  stop("Paleocoastlines not found at '", coast_dir, 
       "'. Please download the PALEOMAP paleocoastlines (v7) and place the 'CS' folder there.")
}

extract_age <- function(filename) {
  as.numeric(str_extract(filename, "\\d+(\\.\\d+)?"))
}

all_coast_files <- list.files(coast_dir, pattern = "\\.shp$", full.names = TRUE)
coast_ages <- sapply(basename(all_coast_files), extract_age)

get_coastline <- function(target_age) {
  idx <- which.min(abs(coast_ages - target_age))
  cat(sprintf("Age %d Ma: using %s (age %.1f Ma)\n",
              target_age, basename(all_coast_files[idx]), coast_ages[idx]))
  st_read(all_coast_files[idx], quiet = TRUE)
}

interval_ages <- c(Triassic   = 230,
                   Jurassic   = 170,
                   Cretaceous = 100,
                   Paleogene  = 50,
                   Neogene    = 15)

plot_period_map <- function(data, period_name, age_ma) {
  pts <- data %>%
    filter(Period == period_name, !is.na(paleolat), !is.na(paleolng)) %>%
    st_as_sf(coords = c("paleolng", "paleolat"), crs = 4326)
  
  coast_sf <- get_coastline(age_ma)
  
  ggplot() +
    geom_sf(data = coast_sf, fill = "grey95", color = "grey40", linewidth = 0.2) +
    geom_sf(data = pts, color = "#E41A1C", alpha = 0.5, size = 0.8, shape = 16) +
    coord_sf(xlim = c(-180, 180), ylim = c(-80, 80), expand = FALSE) +
    labs(title = paste0(period_name, " (", age_ma, " Ma)")) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_line(color = "grey80", linewidth = 0.2),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.margin  = margin(0, 0, 0, 0, "mm"),
      axis.title   = element_blank(),
      axis.text    = element_text(size = 8)
    )
}

plot_list <- lapply(names(interval_ages), function(p) {
  plot_period_map(occurrence_data, p, interval_ages[p])
})
names(plot_list) <- names(interval_ages)

p_paleomaps <- wrap_plots(plot_list, ncol = 3, nrow = 2) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag   = element_text(size = 14, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, "mm"))

# ---- 8. Save all figures ----
if (!dir.exists("Figures")) dir.create("Figures")

ggsave("Figures/Fig3_Taxonomic_Scale.jpg", p_taxonomic,
       width = 8, height = 7, dpi = 600)
ggsave("Figures/Fig4_Niche.jpg", p_niche,
       width = 8, height = 5, dpi = 600)
ggsave("Figures/Fig5_Latitude.jpg", p_lat,
       width = 8, height = 6, dpi = 600)
ggsave("Figures/Fig6_Paleogeographic_Maps.jpg", p_paleomaps,
       width = 18, height = 7, dpi = 600)

cat("Processing complete. All figures successfully saved in the 'Figures' directory.\n")

# ---- 9. Supplementary family-level summaries ----
# Frequency of families by period
family_period_counts <- occurrence_data %>%
  filter(!is.na(family)) %>%
  count(Period, family) %>%
  arrange(Period, desc(n))

top5_per_period <- family_period_counts %>%
  group_by(Period) %>%
  slice_max(order_by = n, n = 5)

cat("\n--- Top 5 families by period ---\n")
print(top5_per_period, n = 50)

# Overall most abundant families
family_total <- occurrence_data %>%
  filter(!is.na(family)) %>%
  count(family, sort = TRUE)

cat("\n--- Top 20 most abundant families overall ---\n")
print(head(family_total, 20))

# Families surviving into the Paleogene and/or Neogene
ceno_families <- occurrence_data %>%
  filter(Period %in% c("Paleogene", "Neogene"), !is.na(family)) %>%
  distinct(family) %>%
  pull(family)

cat("\n--- Families present in Paleogene and/or Neogene ---\n")
cat(paste(ceno_families, collapse = ", "), "\n")

ceno_genera <- occurrence_data %>%
  filter(Period == "Neogene", !is.na(genus)) %>%
  distinct(genus) %>%
  pull(genus)

cat("\n--- Genera present in the Neogene ---\n")
cat(paste(ceno_genera, collapse = ", "), "\n")

# Mesozoic-only vs. Cenozoic-survivor families
meso_families <- occurrence_data %>%
  filter(Period %in% c("Triassic", "Jurassic", "Cretaceous"), !is.na(family)) %>%
  distinct(family) %>%
  pull(family)

extinct_families <- setdiff(meso_families, ceno_families)
cat("\n--- Families present in Mesozoic but ABSENT in Paleogene/Neogene (extinct) ---\n")
cat(paste(extinct_families, collapse = ", "), "\n")

# Environment breakdown for selected top families
top5_families <- c("Rhynchonellidae", "Tetrarhynchiidae", "Zeilleriidae", 
                   "Lobothyrididae", "Cyclothyrididae")

family_env <- occ_eco %>%
  filter(family %in% top5_families) %>%
  count(family, eco_niche) %>%
  group_by(family) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

print(family_env)
