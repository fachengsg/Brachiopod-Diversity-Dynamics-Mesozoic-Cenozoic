# ==============================================================================
# Script Name: 08_Combined_Main_Plots.R
# Purpose: Combine Brachiopoda and Bivalvia diversity curves into side-by-side
#          panels using identical plot logic and aesthetics as Script 04/05.
#          Echinodermata excluded due to sparse data.
#          System labels placed above plot area, below facet strips.
#          Paleogene–Neogene merged to avoid overlap.
# ==============================================================================

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(divDyn)

# ---- 1. Load plot data (Brachiopoda & Bivalvia only) ----
brach_bin   <- readRDS("Plot_Data/brach_bin_095.rds")
biv_bin     <- readRDS("Plot_Data/bivalvia_bin_095.rds")
brach_stage <- readRDS("Plot_Data/brach_stage_095.rds")
biv_stage   <- readRDS("Plot_Data/bivalvia_stage_095.rds")

bin_data <- bind_rows(brach_bin, biv_bin) %>%
  mutate(Clade = factor(Clade, levels = c("Brachiopoda", "Bivalvia")))
stage_data <- bind_rows(brach_stage, biv_stage) %>%
  mutate(Clade = factor(Clade, levels = c("Brachiopoda", "Bivalvia")))

# ---- 2. Aesthetics (aligned with Script 04 & 05) ----
data(stages)
stages_ph <- stages %>%
  filter(stg >= 250 | system %in% c("Triassic", "Jurassic", "Cretaceous",
                                    "Paleogene", "Neogene", "Quaternary"))

sys_colours <- c("Triassic"   = "#812B92", "Jurassic"   = "#34B2C9",
                 "Cretaceous" = "#6E9E44", "Paleogene"  = "#FD9A52",
                 "Neogene"    = "#F3E13C", "Quaternary" = "#F9F9D9")

region_colours <- c("Global (excl. China)" = "#2c7bb6", "China" = "#d7191c")
linetype_all   <- c("Standardised" = "solid", "Raw" = "dashed", "Missing" = "dotted")

theme_common <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = -10, b = 0), 
    legend.spacing.x = unit(0.2, "cm"),     
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    axis.title = element_text(face = "plain", size = 11),
    axis.text  = element_text(size = 10),
    plot.margin = margin(t = 16, r = 8, b = 8, l = 8),
    strip.text = element_text(margin = margin(t = 10, b = 5), face = "bold", size = 12)
  )

sys_rect <- stages_ph %>%
  group_by(system) %>%
  summarise(xmin = min(top), xmax = max(bottom), .groups = "drop") %>%
  filter(system %in% names(sys_colours))

system_labels <- sys_rect %>%
  mutate(mid = (xmin + xmax) / 2, label = system)

# ---- 3. Missing segment helpers (Customized thresholds for clades) ----
make_segments_bin <- function(df) {
  df %>%
    filter(!is.na(value)) %>%
    group_by(Clade, region, type) %>%
    arrange(age) %>%
    mutate(next_age = lead(age), next_val = lead(value)) %>%
    filter(!is.na(next_age), abs(next_age - age) > 10) %>%
    filter(
      (Clade == "Brachiopoda" & abs(next_age - age) <= 35) |
        (Clade == "Bivalvia" & abs(next_age - age) <= 85)
    ) %>%
    ungroup() %>%
    select(Clade, region, type, x = age, xend = next_age, y = value, yend = next_val)
}

stage_order <- stages_ph %>%
  filter(stage %in% unique(stage_data$stage)) %>%
  pull(stage)

make_segments_stage <- function(df) {
  df %>%
    filter(!is.na(value)) %>%
    left_join(data.frame(stage = stage_order, idx = seq_along(stage_order)),
              by = "stage") %>%
    group_by(Clade, region, type) %>%
    arrange(idx) %>%
    mutate(next_idx  = lead(idx),
           next_mid  = lead(mid_age),
           next_val  = lead(value)) %>%
    filter(!is.na(next_mid), abs(next_mid - mid_age) > 0) %>%
    filter(
      (Clade == "Brachiopoda" & abs(next_mid - mid_age) <= 35) |
        (Clade == "Bivalvia" & abs(next_mid - mid_age) <= 85)
    ) %>%
    ungroup() %>%
    select(Clade, region, type, x = mid_age, xend = next_mid, y = value, yend = next_val)
}

seg_bin   <- make_segments_bin(bin_data)
seg_stage <- make_segments_stage(stage_data)

# ---- 4. Combined 10-Myr bin plot ----
p_bin <- ggplot() +
  geom_line(data = bin_data,
            aes(x = age, y = value, color = region, linetype = type),
            linewidth = 1.0, na.rm = TRUE) +
  geom_point(data = bin_data,
             aes(x = age, y = value, color = region, shape = type),
             size = 2.5, na.rm = TRUE) +
  geom_segment(data = seg_bin,
               aes(x = x, xend = xend, y = y, yend = yend,
                   color = region, linetype = "Missing"),
               linewidth = 0.7) +
  facet_wrap(~ Clade, ncol = 2) +
  scale_x_reverse(breaks = seq(0, 250, by = 50), name = "Age (Ma)") +
  scale_y_continuous(name = "Genus Richness",
                     expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_manual(values = region_colours, guide = guide_legend(title = NULL, nrow = 1)) +
  scale_linetype_manual(values = linetype_all, guide = guide_legend(title = NULL, nrow = 1)) +
  scale_shape_manual(values = c("Standardised" = 16, "Raw" = 17, "Missing" = 1), guide = guide_legend(title = NULL, nrow = 1)) +
  labs(title = NULL) + 
  theme_common

# ---- 5. Combined stage-level plot (system labels above plot, below strips) ----
p_stage <- ggplot() +
  geom_rect(data = sys_rect,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = system),
            alpha = 0.15, inherit.aes = FALSE) +
  scale_fill_manual(values = sys_colours, guide = "none") +
  # System labels placed above plot area but below facet strips
  geom_text(data = system_labels,
            aes(x = mid, y = Inf, label = label),
            inherit.aes = FALSE,
            colour = "grey20", fontface = "bold", size = 2.4, vjust = 2.2) +
  coord_cartesian(clip = "off") +
  geom_line(data = stage_data,
            aes(x = mid_age, y = value, color = region, linetype = type),
            linewidth = 1.0, na.rm = TRUE) +
  geom_point(data = stage_data,
             aes(x = mid_age, y = value, color = region, shape = type),
             size = 2.5, na.rm = TRUE) +
  geom_segment(data = seg_stage,
               aes(x = x, xend = xend, y = y, yend = yend,
                   color = region, linetype = "Missing"),
               linewidth = 0.7) +
  facet_wrap(~ Clade, ncol = 2) +
  scale_x_reverse(breaks = seq(0, 250, by = 50), name = "Age (Ma)") +
  scale_y_continuous(name = "Genus Richness",
                     expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_manual(values = region_colours, guide = guide_legend(title = NULL, nrow = 1)) +
  scale_linetype_manual(values = linetype_all, guide = guide_legend(title = NULL, nrow = 1)) +
  scale_shape_manual(values = c("Standardised" = 16, "Raw" = 17, "Missing" = 1), guide = guide_legend(title = NULL, nrow = 1)) +
  labs(title = NULL) + 
  theme_common +
  theme(
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10),  # enough space for labels
    axis.text.x.top = element_blank()
  )

# ---- 6. Save outputs ----
if (!dir.exists("Main_Figures")) dir.create("Main_Figures")

ggsave("Main_Figures/Fig_Main_Bin_Combined.jpg", p_bin, width = 10, height = 5, dpi = 600)
ggsave("Main_Figures/Fig_Main_Stage_Combined.jpg", p_stage, width = 10, height = 5, dpi = 600)

message("Combined main figures saved (JPG, 600 dpi).")