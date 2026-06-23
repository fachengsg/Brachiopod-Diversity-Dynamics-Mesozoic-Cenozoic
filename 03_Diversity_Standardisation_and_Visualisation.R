# ==============================================================================
# Script Name: 03_Diversity_Standardisation_and_Visualisation.R
# Purpose: Calculate and plot coverage-based standardised diversity (q=0) for 
#          Brachiopoda across three datasets (Global incl. China, Global excl. 
#          China, China) using both stage-level and 10-Myr binning.
# Outputs: 8 PDF figures comparing standardized and raw richness.
# ==============================================================================

library(tidyverse)
library(iNEXT)
library(divDyn)
library(ggplot2)

# ==============================================================================
# 0. Load unified data and partition datasets
# ==============================================================================
# Build subsets directly from the cleaned dataset (Script 02) to ensure reproducibility.

if (!file.exists("Brachiopoda_analysis_data.rds")) {
  stop("Cleaned dataset not found. Please run Script 02 first.")
}

analysis_data <- readRDS("Brachiopoda_analysis_data.rds")

# Partition datasets (PBDB uses 'cc'="CN", GBDB uses 'country'="China")
china_data <- analysis_data %>% 
  filter((!is.na(cc) & cc == "CN") | (!is.na(country) & country == "China"))

global_data <- analysis_data

global_ex_china_data <- analysis_data %>% 
  filter(!((!is.na(cc) & cc == "CN") | (!is.na(country) & country == "China")))

# Function: Create Stage-level abundance matrices (Genus x Stage)
make_stage_matrix <- function(data) {
  data %>%
    filter(!is.na(accepted_genus) & !is.na(early_interval)) %>%
    group_by(early_interval, accepted_genus) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = early_interval, values_from = count, values_fill = 0) %>%
    column_to_rownames("accepted_genus") %>%
    as.matrix()
}

# Generate stage matrices
global_incl_genus <- make_stage_matrix(global_data)
global_excl_genus <- make_stage_matrix(global_ex_china_data)
china_genus       <- make_stage_matrix(china_data)

# Define datasets for loop processing
dataset_list <- list(
  "Global_incl" = list(matrix = global_incl_genus, data = global_data),
  "Global_excl" = list(matrix = global_excl_genus, data = global_ex_china_data),
  "China"       = list(matrix = china_genus,       data = china_data)
)

cat("Data partitioned successfully.\n",
    "Global:", nrow(global_data), "records\n",
    "Global (excl. China):", nrow(global_ex_china_data), "records\n",
    "China:", nrow(china_data), "records\n\n")

# ==============================================================================
# 1. Helper Functions & Plot Aesthetics
# ==============================================================================

pivot_est <- function(est) {
  est %>%
    as.data.frame() %>%
    select(Assemblage, Order.q, qD) %>%
    pivot_wider(names_from = Order.q, values_from = qD, names_prefix = "q")
}

make_abundance_matrix <- function(data, taxon_col, bin_col = "bin") {
  data %>%
    filter(!is.na(.data[[taxon_col]])) %>%
    group_by(.data[[bin_col]], .data[[taxon_col]]) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = all_of(bin_col), values_from = count, values_fill = 0) %>%
    column_to_rownames(taxon_col) %>%
    as.matrix()
}

# Define chronological stage order and international system colours
data(stages)
stages_ph <- stages %>%
  filter(stg >= 250 | system %in% c("Triassic", "Jurassic", "Cretaceous",
                                    "Paleogene", "Neogene", "Quaternary"))
stage_order_all <- stages_ph$stage

sys_colours <- c(
  "Triassic"   = "#812B92", "Jurassic"  = "#34B2C9", "Cretaceous" = "#6E9E44",
  "Paleogene"  = "#FD9A52", "Neogene"   = "#F3E13C", "Quaternary" = "#F9F9D9"
)

# ==============================================================================
# PART A: Stage-Level Analysis
# ==============================================================================
stage_results <- list()

# Calculate iNEXT and base diversity (C=0.95)
for (set_name in names(dataset_list)) {
  mat <- dataset_list[[set_name]]$matrix
  
  iNEXT_out <- iNEXT(mat, q = c(0, 1, 2), datatype = "abundance",
                     endpoint = 1.5, knots = 20, nboot = 30, conf = 0.95)
  div_095 <- estimateD(mat, q = c(0, 1, 2), datatype = "abundance",
                       base = "coverage", level = 0.95)
  
  raw <- data.frame(stage = colnames(mat), raw_richness = colSums(mat > 0))
  stage_results[[set_name]] <- list(iNEXT_out = iNEXT_out, div_095 = div_095, raw = raw, mat = mat)
}

# Isolate common stages across datasets
common_stages_all <- Reduce(intersect, lapply(stage_results, function(x) colnames(x$mat)))

sc_list <- lapply(names(dataset_list), function(set_name) {
  info <- DataInfo(stage_results[[set_name]]$mat, datatype = "abundance")
  info$SC[match(common_stages_all, info$Assemblage)]
})
names(sc_list) <- names(dataset_list)

# ---- DIAGNOSTIC: Inspect Sample Coverage (SC) across all Stages ----
cat("\n--- STAGE-LEVEL SAMPLE COVERAGE DIAGNOSTIC ---\n")

# Combine stage SCs into a single dataframe
stage_sc_df <- data.frame(Assemblage = common_stages_all) %>%
  left_join(DataInfo(stage_results$Global_incl$mat, datatype = "abundance") %>% select(Assemblage, Global_incl_SC = SC), by = "Assemblage") %>%
  left_join(DataInfo(stage_results$Global_excl$mat, datatype = "abundance") %>% select(Assemblage, Global_excl_SC = SC), by = "Assemblage") %>%
  left_join(DataInfo(stage_results$China$mat, datatype = "abundance") %>% select(Assemblage, China_SC = SC), by = "Assemblage")

# Reorder chronologically based on international stages
stage_sc_df <- stage_sc_df %>%
  filter(Assemblage %in% stage_order_all) %>%
  mutate(Assemblage = factor(Assemblage, levels = stage_order_all)) %>%
  arrange(Assemblage)

# Safely print top 30 stages
print(head(stage_sc_df, 30))

# Isolate stages falling below the 0.6 threshold
cat("\n--- STAGES FAILING THE 0.6 THRESHOLD ---\n")
print(stage_sc_df %>% filter(Global_excl_SC < 0.6 | China_SC < 0.6))

# ---- Filter & Interpolate ----
# Filter out stages where Sample Coverage (SC) < 0.6 for either Global_excl or China
keep_stage <- sc_list$China >= 0.6 & sc_list$Global_excl >= 0.6
valid_stages_all <- common_stages_all[keep_stage]
cat("\nStages removed (SC < 0.6):", paste(common_stages_all[!keep_stage], collapse = ", "), "\n")

# Compute diversity at minimal shared coverage across valid stages
sc_valid <- sapply(names(dataset_list), function(set_name) {
  mat_valid <- stage_results[[set_name]]$mat[, valid_stages_all, drop = FALSE]
  DataInfo(mat_valid, datatype = "abundance")$SC
})
min_cov_stage <- min(unlist(sc_valid), na.rm = TRUE)
cat("Minimal shared coverage (stage-level):", round(min_cov_stage, 4), "\n\n")

for (set_name in names(dataset_list)) {
  mat_valid <- stage_results[[set_name]]$mat[, valid_stages_all, drop = FALSE]
  stage_results[[set_name]]$div_min <- estimateD(mat_valid, q = c(0, 1, 2), datatype = "abundance",
                                                 base = "coverage", level = min_cov_stage)
}

# ---- Plotting Preparation (Stage-Level) ----
stage_order <- stage_order_all[stage_order_all %in% common_stages_all]
stage_info <- stages_ph %>%
  select(stage, top, bottom, system) %>%
  mutate(mid_age = (top + bottom) / 2) %>%
  filter(stage %in% stage_order)

build_stage_long <- function(global_set, label_global) {
  build_ts <- function(div_metric) {
    data.frame(stage = stage_order) %>%
      left_join(pivot_est(stage_results[[global_set]][[div_metric]]) %>% rename(global_q0 = q0), by = c("stage" = "Assemblage")) %>%
      left_join(pivot_est(stage_results$China[[div_metric]]) %>% rename(china_q0 = q0), by = c("stage" = "Assemblage")) %>%
      left_join(stage_info, by = "stage") %>%
      left_join(stage_results[[global_set]]$raw %>% rename(global_raw = raw_richness), by = "stage") %>%
      left_join(stage_results$China$raw %>% rename(china_raw = raw_richness), by = "stage") %>%
      pivot_longer(c(global_q0, china_q0, global_raw, china_raw), names_to = c("region", "type"), names_pattern = "(global|china)_(.*)") %>%
      mutate(region = ifelse(region == "global", label_global, "China"),
             type   = ifelse(type == "q0", "Standardised", "Raw"))
  }
  list(plot_data_095 = build_ts("div_095"), plot_data_min = build_ts("div_min"))
}

stage_plot_excl <- build_stage_long("Global_excl", "Global (excl. China)")
stage_plot_incl <- build_stage_long("Global_incl", "Global (incl. China)")

# Generate segments for missing intervals
make_stage_segments <- function(plot_data) {
  plot_data %>%
    filter(type == "Standardised", !is.na(value)) %>%
    left_join(data.frame(stage = stage_order, idx = seq_along(stage_order)), by = "stage") %>%
    arrange(region, idx) %>%
    group_by(region) %>%
    mutate(next_idx = lead(idx), next_mid = lead(mid_age), next_val = lead(value)) %>%
    filter(!is.na(next_idx), next_idx - idx > 1) %>%
    ungroup() %>%
    select(region, x = mid_age, xend = next_mid, y = value, yend = next_val)
}

seg_stage_095_excl <- make_stage_segments(stage_plot_excl$plot_data_095)
seg_stage_min_excl <- make_stage_segments(stage_plot_excl$plot_data_min)
seg_stage_095_incl <- make_stage_segments(stage_plot_incl$plot_data_095)
seg_stage_min_incl <- make_stage_segments(stage_plot_incl$plot_data_min)

# ==============================================================================
# PART B: 10-Myr Bin-Level Analysis
# ==============================================================================

# Ensure valid bin column exists
bin_breaks <- c(seq(0, 250, by = 10), 260)
bin_labels <- paste0(seq(0, 250, by = 10), "-", seq(10, 260, by = 10))

for (set_name in names(dataset_list)) {
  if (!"bin" %in% names(dataset_list[[set_name]]$data)) {
    dataset_list[[set_name]]$data <- dataset_list[[set_name]]$data %>%
      mutate(
        mid_age = (as.numeric(max_ma) + as.numeric(min_ma)) / 2,
        bin = cut(mid_age, breaks = bin_breaks, labels = bin_labels, right = FALSE)
      ) %>%
      filter(!is.na(bin))
  }
}

bin_matrices <- lapply(dataset_list, function(x) make_abundance_matrix(x$data, "accepted_genus"))

common_bins <- Reduce(intersect, lapply(bin_matrices, colnames))
common_bins <- common_bins[order(as.numeric(sub("-.*", "", common_bins)))] # youngest first
bin_matrices <- lapply(bin_matrices, function(m) m[, common_bins, drop = FALSE])

# iNEXT standardisation for bins
bin_results <- list()
for (set_name in names(dataset_list)) {
  mat <- bin_matrices[[set_name]]
  bin_results[[set_name]] <- list(
    iNEXT_out = iNEXT(mat, q = c(0, 1, 2), datatype = "abundance", endpoint = 1.5, knots = 20, nboot = 30, conf = 0.95),
    div_095   = estimateD(mat, q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = 0.95),
    mat       = mat
  )
}

# Determine valid bins and calculate minimal shared coverage
sc_bin_list <- lapply(bin_matrices, function(m) DataInfo(m, datatype = "abundance")$SC)
names(sc_bin_list) <- names(dataset_list)

# ---- DIAGNOSTIC: Inspect Sample Coverage (SC) across all 10-Myr Bins ----
cat("\n--- BIN-LEVEL SAMPLE COVERAGE DIAGNOSTIC ---\n")

bin_sc_df <- data.frame(Assemblage = common_bins) %>%
  left_join(DataInfo(bin_matrices$Global_incl, datatype = "abundance") %>% select(Assemblage, Global_incl_SC = SC), by = "Assemblage") %>%
  left_join(DataInfo(bin_matrices$Global_excl, datatype = "abundance") %>% select(Assemblage, Global_excl_SC = SC), by = "Assemblage") %>%
  left_join(DataInfo(bin_matrices$China, datatype = "abundance") %>% select(Assemblage, China_SC = SC), by = "Assemblage")

# Reorder chronologically
bin_sc_df <- bin_sc_df %>%
  mutate(start_age = as.numeric(sub("-.*", "", Assemblage))) %>%
  arrange(desc(start_age)) %>%
  select(-start_age)

print(head(bin_sc_df, 30))

cat("\n--- BINS FAILING THE 0.6 THRESHOLD ---\n")
print(bin_sc_df %>% filter(Global_excl_SC < 0.6 | China_SC < 0.6))

# ---- Filter & Interpolate ----
keep_bin <- sc_bin_list$China >= 0.6 & sc_bin_list$Global_excl >= 0.6
valid_bins <- common_bins[keep_bin]
cat("\nBins removed (SC < 0.6):", paste(common_bins[!keep_bin], collapse = ", "), "\n")

sc_bin_valid <- sapply(names(dataset_list), function(set_name) {
  DataInfo(bin_matrices[[set_name]][, valid_bins, drop = FALSE], datatype = "abundance")$SC
})
min_cov_bin <- min(unlist(sc_bin_valid), na.rm = TRUE)
cat("Minimal shared coverage (bin-level):", round(min_cov_bin, 4), "\n\n")

for (set_name in names(dataset_list)) {
  bin_results[[set_name]]$div_min <- estimateD(bin_matrices[[set_name]][, valid_bins, drop = FALSE], 
                                               q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = min_cov_bin)
}

# ---- Plotting Preparation (Bin-Level) ----
all_bins_sorted <- bin_labels[order(bin_breaks[-length(bin_breaks)], decreasing = TRUE)]
bin_start_age_sorted <- as.numeric(sub("-.*", "", all_bins_sorted))

raw_bin_list <- lapply(dataset_list, function(x) {
  x$data %>% 
    filter(!is.na(accepted_genus)) %>% 
    group_by(bin) %>% 
    summarise(raw = n_distinct(accepted_genus), .groups = "drop")
})

build_bin_long <- function(global_set, label_global) {
  raw_df <- data.frame(bin = all_bins_sorted, age = bin_start_age_sorted) %>%
    left_join(raw_bin_list[[global_set]], by = "bin") %>% rename(global_raw = raw) %>%
    left_join(raw_bin_list$China, by = "bin") %>% rename(china_raw = raw)
  
  build_std <- function(div_metric) {
    data.frame(bin = all_bins_sorted, age = bin_start_age_sorted,
               global_q0 = bin_results[[global_set]][[div_metric]]$qD[match(all_bins_sorted, bin_results[[global_set]][[div_metric]]$Assemblage)],
               china_q0  = bin_results$China[[div_metric]]$qD[match(all_bins_sorted, bin_results$China[[div_metric]]$Assemblage)]) %>%
      pivot_longer(c(global_q0, china_q0), names_to = "region", values_to = "value") %>%
      mutate(region = ifelse(region == "global_q0", label_global, "China"), type = "Standardised")
  }
  
  raw_long <- raw_df %>%
    pivot_longer(c(global_raw, china_raw), names_to = "region", values_to = "value") %>%
    mutate(region = ifelse(region == "global_raw", label_global, "China"), type = "Raw")
  
  bind_rows(list(C095 = bind_rows(build_std("div_095"), raw_long), 
                 mincov = bind_rows(build_std("div_min"), raw_long)), .id = "coverage_type")
}

bin_plot_excl <- build_bin_long("Global_excl", "Global (excl. China)")
bin_plot_incl <- build_bin_long("Global_incl", "Global (incl. China)")

make_dashed_segments_bin <- function(df) {
  df %>%
    filter(type == "Standardised", !is.na(value)) %>%
    group_by(region) %>% arrange(age) %>%
    mutate(next_age = lead(age), next_val = lead(value)) %>%
    filter(!is.na(next_age), abs(next_age - age) > 10) %>%
    ungroup() %>%
    select(region, x = age, xend = next_age, y = value, yend = next_val)
}

seg_bin_095_excl <- make_dashed_segments_bin(filter(bin_plot_excl, coverage_type == "C095"))
seg_bin_min_excl <- make_dashed_segments_bin(filter(bin_plot_excl, coverage_type == "mincov"))
seg_bin_095_incl <- make_dashed_segments_bin(filter(bin_plot_incl, coverage_type == "C095"))
seg_bin_min_incl <- make_dashed_segments_bin(filter(bin_plot_incl, coverage_type == "mincov"))

# ==============================================================================
# PART C: Data Visualisation Export
# ==============================================================================

theme_common <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(face = "bold"))

region_colours <- c("Global (excl. China)" = "#2c7bb6", "Global (incl. China)" = "#2c7bb6", "China" = "#d7191c")
linetype_all <- c("Standardised" = "solid", "Raw" = "dashed", "Missing" = "dotted")

sys_rect <- stages_ph %>%
  group_by(system) %>% summarise(xmin = min(top), xmax = max(bottom), .groups = "drop") %>%
  filter(system %in% names(sys_colours))

# ---- Base plot helper ----
build_plot <- function(df, segments, region_cols, title_text, x_var, is_stage = FALSE) {
  
  # Match active attributes to prevent empty scale warnings in ggplot
  active_regions <- unique(df$region)
  active_lines   <- unique(df$type)
  if (nrow(segments) > 0) active_lines <- c(active_lines, "Missing")
  
  plot_colors <- region_cols[names(region_cols) %in% active_regions]
  plot_lines  <- linetype_all[names(linetype_all) %in% active_lines]
  
  p <- ggplot()
  
  # 1. Draw background rectangles at the bottom layer
  if (is_stage) {
    p <- p + geom_rect(data = sys_rect, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = system), alpha = 0.15, inherit.aes = FALSE) +
      scale_fill_manual(values = sys_colours, guide = "none")
  }
  
  # 2. Add main diversity lines
  p <- p + geom_line(data = df, aes(x = !!sym(x_var), y = value, color = region, linetype = type), linewidth = 1.0, na.rm = TRUE) +
    scale_x_reverse(breaks = seq(0, 250, by = 50), name = "Age (Ma)") +
    scale_y_continuous(name = "Number of genera") +
    scale_color_manual(values = plot_colors) +
    scale_linetype_manual(values = plot_lines) +
    labs(title = title_text, linetype = "Type") +
    theme_common
  
  # 3. Add dashed segments for missing data intervals
  if (nrow(segments) > 0) {
    p <- p + geom_segment(data = segments, aes(x = x, xend = xend, y = y, yend = yend, color = region, linetype = "Missing"), linewidth = 0.7)
  }
  return(p)
}

# ---- 1. Generate Main Text Figures (Global excl. China) ----
main_plots <- list(
  stage_095 = build_plot(stage_plot_excl$plot_data_095, seg_stage_095_excl, region_colours, 
                         "Stage-level genus richness (q = 0, C = 0.95)", "mid_age", is_stage = TRUE),
  stage_min = build_plot(stage_plot_excl$plot_data_min, seg_stage_min_excl, region_colours, 
                         paste0("Stage-level genus richness (q = 0, C = ", round(min_cov_stage, 3), ")"), "mid_age", is_stage = TRUE),
  bin_095   = build_plot(filter(bin_plot_excl, coverage_type == "C095"), seg_bin_095_excl, region_colours, 
                         "10-Myr bin genus richness (q = 0, C = 0.95)", "age"),
  bin_min   = build_plot(filter(bin_plot_excl, coverage_type == "mincov"), seg_bin_min_excl, region_colours, 
                         paste0("10-Myr bin genus richness (q = 0, C = ", round(min_cov_bin, 3), ")"), "age")
)

# ---- 2. Generate Supplementary Figures (Global incl. China) ----
supp_plots <- list(
  stage_095 = build_plot(stage_plot_incl$plot_data_095, seg_stage_095_incl, region_colours, 
                         "Stage-level genus richness (q = 0, C = 0.95)", "mid_age", is_stage = TRUE),
  stage_min = build_plot(stage_plot_incl$plot_data_min, seg_stage_min_incl, region_colours, 
                         paste0("Stage-level genus richness (q = 0, C = ", round(min_cov_stage, 3), ")"), "mid_age", is_stage = TRUE),
  bin_095   = build_plot(filter(bin_plot_incl, coverage_type == "C095"), seg_bin_095_incl, region_colours, 
                         "10-Myr bin genus richness (q = 0, C = 0.95)", "age"),
  bin_min   = build_plot(filter(bin_plot_incl, coverage_type == "mincov"), seg_bin_min_incl, region_colours, 
                         paste0("10-Myr bin genus richness (q = 0, C = ", round(min_cov_bin, 3), ")"), "age")
)

# ==============================================================================
# Export & Preview Figures
# ==============================================================================

# ---- Export to PDF (Uncomment to save files) ----
# walk2(names(main_plots), main_plots, ~ ggsave(paste0("Fig_Main_", .x, ".pdf"), plot = .y, width = 8, height = 5))
# walk2(names(supp_plots), supp_plots, ~ ggsave(paste0("Fig_Supp_", .x, ".pdf"), plot = .y, width = 8, height = 5))

# ---- Specific Individual Previews ----
# Main Figures:
print(main_plots$stage_095)
print(main_plots$stage_min)
print(main_plots$bin_095)
print(main_plots$bin_min)

# Supplementary Figures:
print(supp_plots$stage_095)
print(supp_plots$stage_min)
print(supp_plots$bin_095)
print(supp_plots$bin_min)
