# ==============================================================================
# Script Name: 04_Brachiopoda_Diversity_Analysis.R
# Purpose: Calculate and plot coverage-based standardized diversity (q=0) for 
#          Brachiopoda. Includes filtered data (ecological guild) and a total 
#          unfiltered global baseline.
# Features: Unified logic mirroring Script 05. Exports statistical results as RDS 
#           for final Excel consolidation.
# ==============================================================================
rm(list = ls())

# ---- 1. Environment Setup & Packages ----
library(tidyverse)
library(iNEXT)
library(divDyn)
library(ggplot2)

# Note: Set the working directory to the local project folder before running.
# setwd("C:/path/to/your/project")
# setwd("D:/PBDB_Project")   # original line kept as reference

if (!file.exists("Brachiopoda_analysis_data_Final.rds")) {
  stop("Cleaned dataset not found. Please run Script 03 first.")
}
analysis_data <- readRDS("Brachiopoda_analysis_data_Final.rds")

# Subset occurrences for China and for global data with/without China
china_data <- analysis_data %>% filter((!is.na(cc) & cc == "CN") | (!is.na(country) & country == "China"))
global_data <- analysis_data
global_ex_china_data <- analysis_data %>% filter(!((!is.na(cc) & cc == "CN") | (!is.na(country) & country == "China")))

# ---- 2. Helper Functions & Global Aesthetics ----
#' Helper to reshape iNEXT diversity estimates (wide format, one column per diversity order)
pivot_est <- function(est) {
  est %>% as.data.frame() %>% select(Assemblage, Order.q, qD) %>%
    pivot_wider(names_from = Order.q, values_from = qD, names_prefix = "q")
}

#' Create an abundance matrix (genera × time bin) from occurrence data
make_abundance_matrix <- function(data, bin_col) {
  data %>%
    filter(!is.na(genus) & !is.na(.data[[bin_col]])) %>%
    group_by(.data[[bin_col]], genus) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = all_of(bin_col), values_from = count, values_fill = 0) %>%
    column_to_rownames("genus") %>%
    as.matrix()
}

#' List of three datasets with pre-built abundance matrices: Global (incl. China),
#' Global (excl. China), and China only.
dataset_list <- list(
  "Global_incl" = list(matrix = make_abundance_matrix(global_data, "early_interval"), data = global_data),
  "Global_excl" = list(matrix = make_abundance_matrix(global_ex_china_data, "early_interval"), data = global_ex_china_data),
  "China"       = list(matrix = make_abundance_matrix(china_data, "early_interval"), data = china_data)
)

data(stages)
stages_ph <- stages %>%
  filter(stg >= 250 | system %in% c("Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary")) %>%
  # Remove informal stages like "Early Jurassic", "Middle Triassic", etc.
  filter(!grepl("Early|Middle|Late", stage, ignore.case = TRUE))

stage_order_all <- stages_ph$stage
sys_colours <- c("Triassic" = "#812B92", "Jurassic" = "#34B2C9", "Cretaceous" = "#6E9E44", "Paleogene" = "#FD9A52", "Neogene" = "#F3E13C", "Quaternary" = "#F9F9D9")

# Colour palette for regions (includes the unfiltered baseline colour)
region_colours <- c("Global (excl. China)" = "#2c7bb6", "Global (incl. China)" = "#2c7bb6", "China" = "#d7191c", "Global Total (Unfiltered)" = "#4daf4a")
linetype_all <- c("Standardised" = "solid", "Raw" = "dashed", "Missing" = "dotted")

# Base ggplot theme used for all diversity plots
# Longer legend key width makes linetype differences clearly visible
theme_common <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(1.2, "cm"),   # Longer lines in legend to distinguish solid/dashed
    plot.title = element_text(face = "bold", size = 13, hjust = 0),
    axis.title = element_text(face = "plain", size = 11),
    axis.text = element_text(size = 10),
    plot.margin = margin(t = 16, r = 8, b = 8, l = 8)
  )

sys_rect <- stages_ph %>% group_by(system) %>% summarise(xmin = min(top), xmax = max(bottom), .groups = "drop") %>% filter(system %in% names(sys_colours))
system_labels <- sys_rect %>% mutate(mid = (xmin + xmax) / 2, label = system)

#' Generic ggplot builder for stage- and bin-level richness curves.
#' All legend titles are suppressed (name = NULL). Linetype and shape are merged
#' because they share the same set of levels (Standardized, Raw). Missing segments
#' use a fixed dotted line, not shown in the legend.
build_plot <- function(df, segments, region_cols, title_text, x_var, is_stage = FALSE) {
  active_regions <- unique(df$region)
  active_lines   <- unique(df$type)
  if (nrow(segments) > 0) active_lines <- c(active_lines, "Missing")
  plot_colors <- region_cols[names(region_cols) %in% active_regions]
  plot_lines  <- linetype_all[names(linetype_all) %in% active_lines]
  if (length(plot_colors) == 0) plot_colors <- region_cols
  if (length(plot_lines) == 0) plot_lines <- linetype_all
  p <- ggplot()
  if (is_stage) {
    p <- p + geom_rect(data = sys_rect, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = system), alpha = 0.15, inherit.aes = FALSE) + 
      scale_fill_manual(values = sys_colours, guide = "none") +
      geom_text(data = system_labels, aes(x = mid, y = Inf, label = label), inherit.aes = FALSE,
                colour = "grey20", fontface = "bold", size = 3.8, vjust = -1.0) +
      coord_cartesian(clip = "off")
  }
  p <- p + 
    geom_line(data = df, aes(x = !!sym(x_var), y = value, color = region, linetype = type), linewidth = 1.0, na.rm = TRUE) +
    geom_point(data = df, aes(x = !!sym(x_var), y = value, color = region, shape = type), size = 2.5, na.rm = TRUE) +
    scale_x_reverse(breaks = seq(0, 250, by = 50), name = "Age (Ma)") +
    scale_y_continuous(name = "Genus Richness") +
    scale_color_manual(values = plot_colors, name = NULL) +   # Remove legend title "region"
    scale_linetype_manual(values = plot_lines, name = NULL) +
    scale_shape_manual(values = c("Standardised" = 16, "Raw" = 17, "Missing" = 1), name = NULL) +
    guides(
      color = guide_legend(nrow = 1)
    ) +
    theme_common +
    theme(
      plot.margin = margin(t = 60, r = 10, b = 10, l = 10),
      plot.title = element_text(face = "bold", size = 13, hjust = 0, vjust = 8),
      axis.text.x.top = element_text(angle = 0, hjust = 0.5, size = 8, face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = -10, b = 0), 
      legend.spacing.x = unit(0.2, "cm")      
    )
  if (nrow(segments) > 0) {
    p <- p + geom_segment(data = segments,
                          aes(x = x, xend = xend, y = y, yend = yend,
                              color = region, linetype = "Missing", shape = "Missing"),
                          linewidth = 0.7)
  }
  return(p)
}

# ==============================================================================
# PART A: STAGE-LEVEL ANALYSIS (ECOLOGICALLY FILTERED)
# ==============================================================================
stage_results <- list()
for (set_name in names(dataset_list)) {
  mat <- dataset_list[[set_name]]$matrix
  iNEXT_out <- iNEXT(mat, q = c(0, 1, 2), datatype = "abundance", endpoint = 1.5, knots = 20, nboot = 30, conf = 0.95)
  div_095 <- estimateD(mat, q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = 0.95)
  
  raw <- data.frame(stage = colnames(mat), raw_richness = colSums(mat > 0))
  datainfo <- DataInfo(mat, datatype = "abundance")
  raw <- raw %>%
    left_join(datainfo %>% select(Assemblage, SC), by = c("stage" = "Assemblage")) %>%
    mutate(raw_richness = ifelse(SC < 0.6, NA, raw_richness)) %>%
    select(-SC)
  
  stage_results[[set_name]] <- list(iNEXT_out = iNEXT_out, div_095 = div_095, raw = raw, mat = mat)
}

# Use the union of Global_incl and Global_excl stages, restricted to valid stages (non-informal)
common_stages_all <- union(colnames(stage_results$Global_incl$mat),
                           colnames(stage_results$Global_excl$mat))
common_stages_all <- common_stages_all[common_stages_all %in% stage_order_all]

# Coverage data (some stages may be NA for China)
sc_list <- lapply(names(dataset_list), function(set_name) {
  avail <- common_stages_all
  sc_vec <- DataInfo(stage_results[[set_name]]$mat, datatype = "abundance")$SC
  names(sc_vec) <- DataInfo(stage_results[[set_name]]$mat, datatype = "abundance")$Assemblage
  sc_vec[avail]
})
names(sc_list) <- names(dataset_list)

stage_sc_df <- data.frame(Assemblage = common_stages_all) %>%
  left_join(DataInfo(stage_results$Global_incl$mat, datatype = "abundance") %>% select(Assemblage, Global_incl_SC = SC), by = "Assemblage") %>%
  left_join(DataInfo(stage_results$Global_excl$mat, datatype = "abundance") %>% select(Assemblage, Global_excl_SC = SC), by = "Assemblage") %>%
  left_join(DataInfo(stage_results$China$mat, datatype = "abundance") %>% select(Assemblage, China_SC = SC), by = "Assemblage")

# Keep stages where any region has SC >= 0.6 (NA allowed)
keep_stage <- (sc_list$China >= 0.6) | (sc_list$Global_excl >= 0.6) | (sc_list$Global_incl >= 0.6)
keep_stage[is.na(keep_stage)] <- FALSE
valid_stages_all <- common_stages_all[keep_stage]

# Helper to safely get SC values for valid stages that exist in each matrix
get_sc_for_stages <- function(mat, stages) {
  exist <- intersect(stages, colnames(mat))
  if (length(exist) == 0) return(NULL)
  DataInfo(mat[, exist, drop = FALSE], datatype = "abundance")$SC
}

sc_values_list <- lapply(stage_results, function(x) get_sc_for_stages(x$mat, valid_stages_all))
sc_values <- unlist(sc_values_list)
min_cov_stage <- min(sc_values[sc_values >= 0.6], na.rm = TRUE)

for (set_name in names(dataset_list)) {
  stages <- intersect(valid_stages_all, colnames(stage_results[[set_name]]$mat))
  if (length(stages) > 0) {
    stage_results[[set_name]]$div_min <- estimateD(
      stage_results[[set_name]]$mat[, stages, drop = FALSE],
      q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = min_cov_stage
    )
  } else {
    stage_results[[set_name]]$div_min <- NULL
  }
}

stage_order <- stage_order_all[stage_order_all %in% common_stages_all]
stage_info <- stages_ph %>% select(stage, top, bottom, system) %>% mutate(mid_age = (top + bottom) / 2) %>% filter(stage %in% stage_order)

build_stage_long <- function(global_set, label_global) {
  build_ts <- function(div_metric) {
    data.frame(stage = stage_order) %>%
      left_join(pivot_est(stage_results[[global_set]][[div_metric]]) %>% rename(global_q0 = q0), by = c("stage" = "Assemblage")) %>%
      left_join(pivot_est(stage_results$China[[div_metric]]) %>% rename(china_q0 = q0), by = c("stage" = "Assemblage")) %>%
      left_join(stage_info, by = "stage") %>%
      left_join(stage_results[[global_set]]$raw %>% rename(global_raw = raw_richness), by = "stage") %>%
      left_join(stage_results$China$raw %>% rename(china_raw = raw_richness), by = "stage") %>%
      pivot_longer(c(global_q0, china_q0, global_raw, china_raw), names_to = c("region", "type"), names_pattern = "(global|china)_(.*)") %>%
      mutate(region = ifelse(region == "global", label_global, "China"), type = ifelse(type == "q0", "Standardised", "Raw"))
  }
  list(plot_data_095 = build_ts("div_095"), plot_data_min = build_ts("div_min"))
}

stage_plot_excl <- build_stage_long("Global_excl", "Global (excl. China)")
stage_plot_incl <- build_stage_long("Global_incl", "Global (incl. China)")

make_stage_segments <- function(plot_data, stage_order_vec = stage_order) {
  plot_data %>% 
    filter(!is.na(value)) %>% 
    left_join(data.frame(stage = stage_order_vec, idx = seq_along(stage_order_vec)), by = "stage") %>%
    arrange(region, type, idx) %>% 
    group_by(region, type) %>% 
    mutate(next_idx = lead(idx), next_mid = lead(mid_age), next_val = lead(value)) %>%
    filter(!is.na(next_idx), next_idx - idx > 1) %>% 
    filter(abs(next_mid - mid_age) <= 35) %>% 
    ungroup() %>% 
    select(region, type, x = mid_age, xend = next_mid, y = value, yend = next_val)
}

seg_stage_095_excl <- make_stage_segments(stage_plot_excl$plot_data_095)
seg_stage_min_excl <- make_stage_segments(stage_plot_excl$plot_data_min)
seg_stage_095_incl <- make_stage_segments(stage_plot_incl$plot_data_095)
seg_stage_min_incl <- make_stage_segments(stage_plot_incl$plot_data_min)

# ==============================================================================
# PART B: 10-MYR BIN-LEVEL ANALYSIS (ECOLOGICALLY FILTERED)
# ==============================================================================
bin_breaks <- c(seq(0, 250, by = 10), 260)
bin_labels <- paste0(seq(0, 250, by = 10), "-", seq(10, 260, by = 10))

all_bins_sorted <- bin_labels[order(bin_breaks[-length(bin_breaks)], decreasing = TRUE)]
bin_start_age_sorted <- as.numeric(sub("-.*", "", all_bins_sorted))

for (set_name in names(dataset_list)) {
  dataset_list[[set_name]]$data <- dataset_list[[set_name]]$data %>% mutate(mid_age = (as.numeric(max_ma) + as.numeric(min_ma)) / 2, bin = cut(mid_age, breaks = bin_breaks, labels = bin_labels, right = FALSE)) %>% filter(!is.na(bin))
}

# Build abundance matrices for each region using bins (original full matrices)
bin_matrices_full <- lapply(dataset_list, function(x) make_abundance_matrix(x$data, "bin"))

# Use the union of Global_incl and Global_excl bins
common_bins <- union(colnames(bin_matrices_full$Global_incl),
                     colnames(bin_matrices_full$Global_excl))
common_bins <- common_bins[common_bins %in% all_bins_sorted]
common_bins <- common_bins[order(as.numeric(sub("-.*", "", common_bins)))]

# For each dataset, keep only the common bins that actually exist in its matrix
bin_matrices <- lapply(bin_matrices_full, function(m) {
  keep <- intersect(common_bins, colnames(m))
  m[, keep, drop = FALSE]
})

# Now run iNEXT and estimateD on these safe matrices
bin_results <- list()
for (set_name in names(dataset_list)) {
  mat <- bin_matrices[[set_name]]
  bin_results[[set_name]] <- list(
    iNEXT_out = iNEXT(mat, q = c(0, 1, 2), datatype = "abundance", endpoint = 1.5, knots = 20, nboot = 30, conf = 0.95),
    div_095 = estimateD(mat, q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = 0.95),
    mat = mat
  )
}

# Coverage data (some bins may be NA for China)
sc_bin_list <- lapply(bin_matrices, function(m) {
  sc <- DataInfo(m, datatype = "abundance")$SC
  names(sc) <- colnames(m)
  sc[common_bins]  # returns NA for missing bins
})
names(sc_bin_list) <- names(dataset_list)

bin_sc_df <- data.frame(Assemblage = common_bins) %>%
  left_join(DataInfo(bin_matrices$Global_incl, datatype = "abundance") %>% select(Assemblage, Global_incl_SC = SC), by = "Assemblage") %>%
  left_join(DataInfo(bin_matrices$Global_excl, datatype = "abundance") %>% select(Assemblage, Global_excl_SC = SC), by = "Assemblage") %>%
  left_join(DataInfo(bin_matrices$China, datatype = "abundance") %>% select(Assemblage, China_SC = SC), by = "Assemblage")

keep_bin <- (sc_bin_list$China >= 0.6) | (sc_bin_list$Global_excl >= 0.6) | (sc_bin_list$Global_incl >= 0.6)
keep_bin[is.na(keep_bin)] <- FALSE
valid_bins <- common_bins[keep_bin]

# Helper to safely get SC values for valid bins that exist in each matrix
get_sc_for_bins <- function(mat, bins) {
  exist <- intersect(bins, colnames(mat))
  if (length(exist) == 0) return(NULL)
  DataInfo(mat[, exist, drop = FALSE], datatype = "abundance")$SC
}

sc_bin_values_list <- lapply(bin_results, function(x) get_sc_for_bins(x$mat, valid_bins))
sc_bin_values <- unlist(sc_bin_values_list)
min_cov_bin <- min(sc_bin_values[sc_bin_values >= 0.6], na.rm = TRUE)

for (set_name in names(dataset_list)) {
  bins <- intersect(valid_bins, colnames(bin_results[[set_name]]$mat))
  if (length(bins) > 0) {
    bin_results[[set_name]]$div_min <- estimateD(
      bin_results[[set_name]]$mat[, bins, drop = FALSE],
      q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = min_cov_bin
    )
  } else {
    bin_results[[set_name]]$div_min <- NULL
  }
}

# Raw bin data with low coverage trimming
raw_bin_list <- lapply(names(dataset_list), function(set_name) {
  x <- dataset_list[[set_name]]
  raw_df <- x$data %>% filter(!is.na(genus)) %>% group_by(bin) %>% summarise(raw = n_distinct(genus), .groups = "drop")
  mat <- bin_matrices[[set_name]]
  datainfo <- DataInfo(mat, datatype = "abundance")
  sc_map <- data.frame(bin = colnames(mat), SC = datainfo$SC)
  raw_df <- raw_df %>% left_join(sc_map, by = "bin") %>% mutate(raw = ifelse(SC < 0.6, NA, raw)) %>% select(-SC)
  return(raw_df)
})
names(raw_bin_list) <- names(dataset_list)

build_bin_long <- function(global_set, label_global) {
  raw_df <- data.frame(bin = all_bins_sorted, age = bin_start_age_sorted) %>% left_join(raw_bin_list[[global_set]], by = "bin") %>% rename(global_raw = raw) %>% left_join(raw_bin_list$China, by = "bin") %>% rename(china_raw = raw)
  build_std <- function(div_metric) {
    data.frame(bin = all_bins_sorted, age = bin_start_age_sorted,
               global_q0 = bin_results[[global_set]][[div_metric]]$qD[match(all_bins_sorted, bin_results[[global_set]][[div_metric]]$Assemblage)],
               china_q0  = bin_results$China[[div_metric]]$qD[match(all_bins_sorted, bin_results$China[[div_metric]]$Assemblage)]) %>%
      pivot_longer(c(global_q0, china_q0), names_to = "region", values_to = "value") %>% mutate(region = ifelse(region == "global_q0", label_global, "China"), type = "Standardised")
  }
  raw_long <- raw_df %>% pivot_longer(c(global_raw, china_raw), names_to = "region", values_to = "value") %>% mutate(region = ifelse(region == "global_raw", label_global, "China"), type = "Raw")
  bind_rows(list(C095 = bind_rows(build_std("div_095"), raw_long), mincov = bind_rows(build_std("div_min"), raw_long)), .id = "coverage_type")
}

bin_plot_excl <- build_bin_long("Global_excl", "Global (excl. China)")
bin_plot_incl <- build_bin_long("Global_incl", "Global (incl. China)")

make_dashed_segments_bin <- function(df) {
  df %>% filter(!is.na(value)) %>% group_by(region, type) %>% arrange(age) %>% 
    mutate(next_age = lead(age), next_val = lead(value)) %>%
    filter(!is.na(next_age), abs(next_age - age) > 10) %>% 
    filter(abs(next_age - age) <= 35) %>% ungroup() %>% 
    select(region, type, x = age, xend = next_age, y = value, yend = next_val)
}

seg_bin_095_excl <- make_dashed_segments_bin(filter(bin_plot_excl, coverage_type == "C095"))
seg_bin_min_excl <- make_dashed_segments_bin(filter(bin_plot_excl, coverage_type == "mincov"))
seg_bin_095_incl <- make_dashed_segments_bin(filter(bin_plot_incl, coverage_type == "C095"))
seg_bin_min_incl <- make_dashed_segments_bin(filter(bin_plot_incl, coverage_type == "mincov"))

# ==============================================================================
# PART C: DATA VISUALIZATION & DATA EXPORT
# ==============================================================================
main_plots <- list(
  stage_095 = build_plot(stage_plot_excl$plot_data_095, seg_stage_095_excl, region_colours, "Stage-level genus richness (q = 0, C = 0.95)", "mid_age", TRUE),
  stage_min = build_plot(stage_plot_excl$plot_data_min, seg_stage_min_excl, region_colours, paste0("Stage-level genus richness (q = 0, C = ", round(min_cov_stage, 3), ")"), "mid_age", TRUE),
  bin_095   = build_plot(filter(bin_plot_excl, coverage_type == "C095"), seg_bin_095_excl, region_colours, "10-Myr bin genus richness (q = 0, C = 0.95)", "age"),
  bin_min   = build_plot(filter(bin_plot_excl, coverage_type == "mincov"), seg_bin_min_excl, region_colours, paste0("10-Myr bin genus richness (q = 0, C = ", round(min_cov_bin, 3), ")"), "age")
)
supp_plots <- list(
  stage_095 = build_plot(stage_plot_incl$plot_data_095, seg_stage_095_incl, region_colours, "Stage-level genus richness (q = 0, C = 0.95)", "mid_age", TRUE),
  stage_min = build_plot(stage_plot_incl$plot_data_min, seg_stage_min_incl, region_colours, paste0("Stage-level genus richness (q = 0, C = ", round(min_cov_stage, 3), ")"), "mid_age", TRUE),
  bin_095   = build_plot(filter(bin_plot_incl, coverage_type == "C095"), seg_bin_095_incl, region_colours, "10-Myr bin genus richness (q = 0, C = 0.95)", "age"),
  bin_min   = build_plot(filter(bin_plot_incl, coverage_type == "mincov"), seg_bin_min_incl, region_colours, paste0("10-Myr bin genus richness (q = 0, C = ", round(min_cov_bin, 3), ")"), "age")
)

if (!dir.exists("Main_Figures")) dir.create("Main_Figures")
if (!dir.exists("Supplementary_Figures")) dir.create("Supplementary_Figures")
for (fig_name in names(main_plots)) ggsave(paste0("Main_Figures/Brachiopoda_Main_", fig_name, ".jpg"), plot = main_plots[[fig_name]], width = 8, height = 5, dpi = 600)
for (fig_name in names(supp_plots)) ggsave(paste0("Supplementary_Figures/Brachiopoda_Supp_", fig_name, ".jpg"), plot = supp_plots[[fig_name]], width = 8, height = 5, dpi = 600)

dir.create("Plot_Data", showWarnings = FALSE, recursive = TRUE)
bin_095_brach <- bin_plot_excl %>% filter(coverage_type == "C095") %>% mutate(Clade = "Brachiopoda")
stage_095_brach <- stage_plot_excl$plot_data_095 %>% mutate(Clade = "Brachiopoda")
saveRDS(bin_095_brach, "Plot_Data/brach_bin_095.rds")
saveRDS(stage_095_brach, "Plot_Data/brach_stage_095.rds")

# ==============================================================================
# PART D: TOTAL GLOBAL BRACHIOPODA (UNFILTERED BASELINE)
# ==============================================================================
cat("\n--- Running Unfiltered Global Brachiopoda Analysis ---\n")

if (file.exists("Brachiopoda_analysis_data.rds")) {
  unfiltered_data <- readRDS("Brachiopoda_analysis_data.rds") %>%
    mutate(
      mid_age = (as.numeric(max_ma) + as.numeric(min_ma)) / 2,
      bin = cut(mid_age, breaks = bin_breaks, labels = bin_labels, right = FALSE)
    ) %>%
    filter(!is.na(bin))
  
  mat_stage_unf <- make_abundance_matrix(unfiltered_data, "early_interval")
  mat_bin_unf   <- make_abundance_matrix(unfiltered_data, "bin")
  
  stage_order_unf <- stages_ph %>%
    filter(stage %in% colnames(mat_stage_unf)) %>%
    arrange(desc((top + bottom) / 2)) %>%
    pull(stage)
  
  stage_info_unf <- stages_ph %>%
    select(stage, top, bottom, system) %>%
    mutate(mid_age = (top + bottom) / 2) %>%
    filter(stage %in% stage_order_unf)
  
  div_stage_unf <- estimateD(mat_stage_unf, q = 0, datatype = "abundance", base = "coverage", level = 0.95)
  raw_stage_unf <- data.frame(stage = colnames(mat_stage_unf), raw = colSums(mat_stage_unf > 0))
  datainfo_stage_unf <- DataInfo(mat_stage_unf, datatype = "abundance")
  raw_stage_unf <- raw_stage_unf %>%
    left_join(datainfo_stage_unf %>% select(Assemblage, SC), by = c("stage" = "Assemblage")) %>%
    mutate(raw = ifelse(SC < 0.6, NA, raw))
  
  df_stage_unf <- data.frame(stage = stage_order_unf) %>%
    left_join(div_stage_unf %>% select(Assemblage, qD), by = c("stage" = "Assemblage")) %>%
    left_join(raw_stage_unf %>% select(stage, raw), by = "stage") %>%
    left_join(stage_info_unf, by = "stage") %>%
    pivot_longer(c(qD, raw), names_to = "type", values_to = "value") %>%
    mutate(
      type = ifelse(type == "qD", "Standardised", "Raw"),
      region = "Global Total (Unfiltered)"
    )
  
  seg_stage_unf <- make_stage_segments(df_stage_unf, stage_order_vec = stage_order_unf)
  
  p_stage_unf <- build_plot(
    df_stage_unf, seg_stage_unf, region_colours, 
    "Global Brachiopoda (Unfiltered) Stage-level Genus Richness", 
    "mid_age", TRUE
  ) +
    guides(colour = "none",
           linetype = guide_legend(override.aes = list(colour = "#4daf4a")))
  
  # Bin-level unfiltered
  bin_order_unf <- bin_start_age_sorted[all_bins_sorted %in% colnames(mat_bin_unf)]
  bin_labels_unf <- all_bins_sorted[all_bins_sorted %in% colnames(mat_bin_unf)]
  
  div_bin_unf <- estimateD(mat_bin_unf, q = 0, datatype = "abundance", base = "coverage", level = 0.95)
  raw_bin_unf <- data.frame(bin = colnames(mat_bin_unf), raw = colSums(mat_bin_unf > 0))
  datainfo_bin_unf <- DataInfo(mat_bin_unf, datatype = "abundance")
  raw_bin_unf <- raw_bin_unf %>%
    left_join(datainfo_bin_unf %>% select(Assemblage, SC), by = c("bin" = "Assemblage")) %>%
    mutate(raw = ifelse(SC < 0.6, NA, raw))
  
  df_bin_unf <- data.frame(bin = bin_labels_unf, age = bin_order_unf) %>%
    left_join(div_bin_unf %>% select(Assemblage, qD), by = c("bin" = "Assemblage")) %>%
    left_join(raw_bin_unf %>% select(bin, raw), by = "bin") %>%
    pivot_longer(c(qD, raw), names_to = "type", values_to = "value") %>%
    mutate(
      type = ifelse(type == "qD", "Standardised", "Raw"),
      region = "Global Total (Unfiltered)"
    )
  
  seg_bin_unf <- make_dashed_segments_bin(df_bin_unf)
  
  p_bin_unf <- build_plot(
    df_bin_unf, seg_bin_unf, region_colours, 
    "Global Brachiopoda (Unfiltered) 10-Myr Bin Genus Richness", 
    "age", FALSE
  ) +
    guides(colour = "none",
           linetype = guide_legend(override.aes = list(colour = "#4daf4a")))
  
  ggsave("Main_Figures/Brachiopoda_Unfiltered_Stage.jpg", plot = p_stage_unf, width = 8, height = 5, dpi = 600)
  ggsave("Main_Figures/Brachiopoda_Unfiltered_Bin.jpg", plot = p_bin_unf, width = 8, height = 5, dpi = 600)
  cat("Unfiltered analysis complete. Baseline plots exported.\n")
} else {
  warning("Raw 'Brachiopoda_analysis_data.rds' not found. Skipping unfiltered baseline analysis.")
}

export_results <- list(
  stage_results = stage_results,
  bin_results = bin_results,
  stage_sc_df = stage_sc_df,
  bin_sc_df = bin_sc_df
)
saveRDS(export_results, "Results_Brachiopoda.rds")
cat("\nBrachiopoda pipeline complete. Plots and statistical RDS files exported.\n")

# ==============================================================================
# PART E: TAXONOMIC HIERARCHY DATA EXPORT
# ==============================================================================
cat("\n--- Running Taxonomic Hierarchy Data Extraction ---\n")
if (exists("unfiltered_data")) {
  make_taxa_matrix <- function(data, time_col, tax_col) {
    data %>%
      filter(!is.na(.data[[tax_col]]) & !is.na(.data[[time_col]])) %>%
      group_by(.data[[time_col]], .data[[tax_col]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      pivot_wider(names_from = all_of(time_col), values_from = count, values_fill = 0) %>%
      column_to_rownames(tax_col) %>%
      as.matrix()
  }
  
  calc_taxa_div <- function(data, time_col, tax_col, tax_name) {
    raw_counts <- data %>%
      filter(!is.na(.data[[tax_col]]) & !is.na(.data[[time_col]])) %>%
      group_by(.data[[time_col]]) %>%
      summarise(raw = n_distinct(.data[[tax_col]]), .groups = "drop") %>%
      rename(Time = .data[[time_col]])
    mat <- make_taxa_matrix(data, time_col, tax_col)
    datainfo <- DataInfo(mat, datatype = "abundance")
    suppressMessages(est <- estimateD(mat, q = 0, datatype = "abundance", base = "coverage", level = 0.95))
    std_counts <- est %>% select(Assemblage, std = qD) %>% rename(Time = Assemblage)
    res <- raw_counts %>% left_join(std_counts, by = "Time") %>%
      left_join(datainfo %>% select(Time = Assemblage, SC), by = "Time") %>%
      mutate(raw = ifelse(SC < 0.6, NA, raw), std = ifelse(SC < 0.6, NA, std)) %>%
      select(Time, raw, std)
    colnames(res)[2:3] <- paste0(tax_name, "_", c("Raw", "Std"))
    res
  }
  
  stage_species <- calc_taxa_div(unfiltered_data, "early_interval", "accepted_name", "Species")
  stage_genus   <- calc_taxa_div(unfiltered_data, "early_interval", "genus", "Genus")
  stage_family  <- calc_taxa_div(unfiltered_data, "early_interval", "family", "Family")
  
  stage_summary <- stage_species %>% full_join(stage_genus, by = "Time") %>% full_join(stage_family, by = "Time")
  
  bin_species <- calc_taxa_div(unfiltered_data, "bin", "accepted_name", "Species")
  bin_genus   <- calc_taxa_div(unfiltered_data, "bin", "genus", "Genus")
  bin_family  <- calc_taxa_div(unfiltered_data, "bin", "family", "Family")
  bin_summary <- bin_species %>% full_join(bin_genus, by = "Time") %>% full_join(bin_family, by = "Time")
  
  stage_summary_sorted <- stage_summary %>%
    filter(Time %in% stage_info$stage) %>%
    left_join(stage_info %>% select(stage, mid_age), by = c("Time" = "stage")) %>%
    arrange(desc(mid_age)) %>% relocate(mid_age, .after = Time)
  
  bin_summary_sorted <- bin_summary %>%
    mutate(start_age = as.numeric(sub("-.*", "", Time))) %>%
    arrange(desc(start_age)) %>% select(-start_age)
  
  write_csv(stage_summary_sorted, "Unfiltered_Brachiopoda_Stage_Summary.csv")
  write_csv(bin_summary_sorted, "Unfiltered_Brachiopoda_Bin_Summary.csv")
  
  cat("\n--- Taxonomic hierarchy data exported ---\n")
} else {
  warning("Object 'unfiltered_data' not found. Ensure PART D was executed successfully.")
}
# ==============================================================================
# END OF SCRIPT 04
# ==============================================================================
