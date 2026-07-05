# ==============================================================================
# Script Name: 05_Clade_Specific_Evolution.R
# Purpose: Separately calculate and plot coverage-based standardised diversity 
#          (q=0) for Bivalvia and Echinodermata.
# Features: Function-based exact duplicate of Script 04 logic. Exports results
#           as RDS per clade for final Excel consolidation.
# ==============================================================================

rm(list = ls())

# ---- 1. Environment Setup & Packages ----
library(tidyverse)
library(iNEXT)
library(divDyn)
library(ggplot2)

setwd("D:/PBDB_Project")

if (!file.exists("Comparative_Bivalvia_Echinodermata.rds")) {
  stop("Dataset not found. Please run Script 03 first.")
}
comp_data <- readRDS("Comparative_Bivalvia_Echinodermata.rds")

bivalvia_df <- comp_data %>% filter(class == "Bivalvia")
echino_df   <- comp_data %>% filter(phylum == "Echinodermata")

# ---- 2. Helper Functions & Global Aesthetics ----
theme_common <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#f0f0f0", linewidth = 0.5),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    
    axis.title = element_text(
      face = "plain",
      size = 11
    ),
    
    axis.text = element_text(size = 10),
    
    # Added for geological system labels alignment
    plot.margin = margin(
      t = 16,
      r = 8,
      b = 8,
      l = 8
    )
  )

data(stages)
stages_ph <- stages %>% filter(stg >= 250 | system %in% c("Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary"))
stage_order_all <- stages_ph$stage

sys_colours <- c("Triassic" = "#812B92", "Jurassic" = "#34B2C9", "Cretaceous" = "#6E9E44", "Paleogene" = "#FD9A52", "Neogene" = "#F3E13C", "Quaternary" = "#F9F9D9")
region_colours <- c("Global (excl. China)" = "#2c7bb6", "Global (incl. China)" = "#2c7bb6", "China" = "#d7191c", "Global Total (Unfiltered)" = "#4daf4a")
linetype_all <- c("Standardised" = "solid", "Raw" = "dashed", "Missing" = "dotted")

sys_rect <- stages_ph %>% group_by(system) %>% summarise(xmin = min(top), xmax = max(bottom), .groups = "drop") %>% filter(system %in% names(sys_colours))

# Geological system labels for secondary x-axis
system_labels <- sys_rect %>%
  mutate(
    mid = (xmin + xmax) / 2,
    label = system
  )

bin_breaks <- c(seq(0, 250, by = 10), 260)
bin_labels <- paste0(seq(0, 250, by = 10), "-", seq(10, 260, by = 10))
all_bins_sorted <- bin_labels[order(bin_breaks[-length(bin_breaks)], decreasing = TRUE)]
bin_start_age_sorted <- as.numeric(sub("-.*", "", all_bins_sorted))

pivot_est <- function(est) {
  est %>% as.data.frame() %>% select(Assemblage, Order.q, qD) %>% pivot_wider(names_from = Order.q, values_from = qD, names_prefix = "q")
}

make_abundance_matrix <- function(data, bin_col) {
  data %>% filter(!is.na(genus) & !is.na(.data[[bin_col]])) %>% group_by(.data[[bin_col]], genus) %>%
    summarise(count = n(), .groups = "drop") %>% pivot_wider(names_from = all_of(bin_col), values_from = count, values_fill = 0) %>%
    column_to_rownames("genus") %>% as.matrix()
}

# Added geom_point to show data points and merged shape legend with linetype legend
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
    scale_color_manual(values = plot_colors) + 
    scale_linetype_manual(values = plot_lines, name = "Data Profile") +
    scale_shape_manual(values = c("Standardised" = 16, "Raw" = 17, "Missing" = 1), name = "Data Profile") +
    guides(
      color = guide_legend(title = NULL, nrow = 1),
      linetype = guide_legend(title = NULL, nrow = 1),
      shape = guide_legend(title = NULL, nrow = 1)
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
    p <- p + geom_segment(
      data = segments, 
      aes(x = x, xend = xend, y = y, yend = yend, color = region, linetype = "Missing"), 
      linewidth = 0.7
    )
  }
  return(p)
}
# ==============================================================================
# PART A: CORE ANALYTICAL PIPELINE (THE MASTER FUNCTION)
# ==============================================================================

analyze_clade_evolution <- function(clade_name, clade_data) {
  cat(paste0("\n========================================\n", "STARTING ANALYSIS FOR: ", clade_name, "\n========================================\n"))
  
  china_data <- clade_data %>% filter((!is.na(cc) & cc == "CN") | (!is.na(country) & country == "China"))
  global_ex_china_data <- clade_data %>% filter(!((!is.na(cc) & cc == "CN") | (!is.na(country) & country == "China")))
  global_incl_data <- clade_data 
  
  if (nrow(china_data) == 0) {
    warning(paste("No regional China records found for clade:", clade_name, "Bypassing analysis."))
    return(NULL)
  }
  
  dataset_list <- list(
    "Global_incl" = list(matrix = make_abundance_matrix(global_incl_data, "early_interval"), data = global_incl_data),
    "Global_excl" = list(matrix = make_abundance_matrix(global_ex_china_data, "early_interval"), data = global_ex_china_data),
    "China"       = list(matrix = make_abundance_matrix(china_data, "early_interval"), data = china_data)
  )
  
  # ---------------- STAGE-LEVEL ----------------
  stage_results <- list()
  for (set_name in names(dataset_list)) {
    mat <- dataset_list[[set_name]]$matrix
    if (ncol(mat) == 0 || nrow(mat) == 0) stop(paste("Empty abundance matrix generated for", set_name, "in", clade_name))
    
    iNEXT_out <- iNEXT(mat, q = c(0, 1, 2), datatype = "abundance", endpoint = 1.5, knots = 20, nboot = 30, conf = 0.95)
    div_095 <- estimateD(mat, q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = 0.95)
    
    # Directly compute SC inside loop to trim raw tails for low sampling intervals
    raw <- data.frame(stage = colnames(mat), raw_richness = colSums(mat > 0))
    datainfo <- DataInfo(mat, datatype = "abundance")
    raw <- raw %>%
      left_join(datainfo %>% select(Assemblage, SC), by = c("stage" = "Assemblage")) %>%
      mutate(raw_richness = ifelse(SC < 0.6, NA, raw_richness)) %>%
      select(-SC)
    
    stage_results[[set_name]] <- list(iNEXT_out = iNEXT_out, div_095 = div_095, raw = raw, mat = mat)
  }
  
  common_stages_all <- Reduce(intersect, lapply(stage_results, function(x) colnames(x$mat)))
  sc_list <- lapply(names(dataset_list), function(set_name) DataInfo(stage_results[[set_name]]$mat, datatype = "abundance")$SC[match(common_stages_all, DataInfo(stage_results[[set_name]]$mat, datatype = "abundance")$Assemblage)])
  names(sc_list) <- names(dataset_list)
  
  stage_sc_df <- data.frame(Assemblage = common_stages_all) %>%
    left_join(DataInfo(stage_results$Global_incl$mat, datatype = "abundance") %>% select(Assemblage, Global_incl_SC = SC), by = "Assemblage") %>%
    left_join(DataInfo(stage_results$Global_excl$mat, datatype = "abundance") %>% select(Assemblage, Global_excl_SC = SC), by = "Assemblage") %>%
    left_join(DataInfo(stage_results$China$mat, datatype = "abundance") %>% select(Assemblage, China_SC = SC), by = "Assemblage")
  
  keep_stage <- sc_list$China >= 0.6 | sc_list$Global_excl >= 0.6
  valid_stages_all <- common_stages_all[keep_stage]
  
  sc_valid_vec <- unlist(sapply(names(dataset_list), function(set_name) DataInfo(stage_results[[set_name]]$mat[, valid_stages_all, drop = FALSE], datatype = "abundance")$SC))
  min_cov_stage <- min(sc_valid_vec[sc_valid_vec >= 0.6], na.rm = TRUE)
  
  for (set_name in names(dataset_list)) {
    stage_results[[set_name]]$div_min <- estimateD(stage_results[[set_name]]$mat[, valid_stages_all, drop = FALSE], q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = min_cov_stage)
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
  
  # Grouped by type (Raw & Standardised) and added max gap limitation to prevent 50Ma jump
  make_stage_segments <- function(plot_data) {
    plot_data %>% 
      filter(!is.na(value)) %>% 
      left_join(data.frame(stage = stage_order, idx = seq_along(stage_order)), by = "stage") %>%
      arrange(region, type, idx) %>% 
      group_by(region, type) %>% 
      mutate(next_idx = lead(idx), next_mid = lead(mid_age), next_val = lead(value)) %>%
      filter(!is.na(next_idx), next_idx - idx > 1) %>% 
      filter(abs(next_mid - mid_age) <= 85) %>% 
      ungroup() %>% 
      select(region, type, x = mid_age, xend = next_mid, y = value, yend = next_val)
  }
  
  seg_stage_095_excl <- make_stage_segments(stage_plot_excl$plot_data_095)
  seg_stage_min_excl <- make_stage_segments(stage_plot_excl$plot_data_min)
  seg_stage_095_incl <- make_stage_segments(stage_plot_incl$plot_data_095)
  seg_stage_min_incl <- make_stage_segments(stage_plot_incl$plot_data_min)
  
  # ---------------- BIN-LEVEL ----------------
  for (set_name in names(dataset_list)) {
    dataset_list[[set_name]]$data <- dataset_list[[set_name]]$data %>% mutate(mid_age = (as.numeric(max_ma) + as.numeric(min_ma)) / 2, bin = cut(mid_age, breaks = bin_breaks, labels = bin_labels, right = FALSE)) %>% filter(!is.na(bin))
  }
  
  bin_matrices <- lapply(dataset_list, function(x) make_abundance_matrix(x$data, "bin"))
  common_bins <- Reduce(intersect, lapply(bin_matrices, colnames))
  common_bins <- common_bins[order(as.numeric(sub("-.*", "", common_bins)))]
  bin_matrices <- lapply(bin_matrices, function(m) m[, common_bins, drop = FALSE])
  
  bin_results <- list()
  for (set_name in names(dataset_list)) {
    mat <- bin_matrices[[set_name]]
    bin_results[[set_name]] <- list(iNEXT_out = iNEXT(mat, q = c(0, 1, 2), datatype = "abundance", endpoint = 1.5, knots = 20, nboot = 30, conf = 0.95),
                                    div_095 = estimateD(mat, q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = 0.95), mat = mat)
  }
  
  sc_bin_list <- lapply(bin_matrices, function(m) DataInfo(m, datatype = "abundance")$SC)
  names(sc_bin_list) <- names(dataset_list)
  
  bin_sc_df <- data.frame(Assemblage = common_bins) %>%
    left_join(DataInfo(bin_matrices$Global_incl, datatype = "abundance") %>% select(Assemblage, Global_incl_SC = SC), by = "Assemblage") %>%
    left_join(DataInfo(bin_matrices$Global_excl, datatype = "abundance") %>% select(Assemblage, Global_excl_SC = SC), by = "Assemblage") %>%
    left_join(DataInfo(bin_matrices$China, datatype = "abundance") %>% select(Assemblage, China_SC = SC), by = "Assemblage")
  
  keep_bin <- sc_bin_list$China >= 0.6 | sc_bin_list$Global_excl >= 0.6
  valid_bins <- common_bins[keep_bin]
  
  sc_bin_valid_vec <- unlist(sapply(names(dataset_list), function(set_name) DataInfo(bin_matrices[[set_name]][, valid_bins, drop = FALSE], datatype = "abundance")$SC))
  min_cov_bin <- min(sc_bin_valid_vec[sc_bin_valid_vec >= 0.6], na.rm = TRUE)
  
  for (set_name in names(dataset_list)) {
    bin_results[[set_name]]$div_min <- estimateD(bin_matrices[[set_name]][, valid_bins, drop = FALSE], q = c(0, 1, 2), datatype = "abundance", base = "coverage", level = min_cov_bin)
  }
  
  # Robust raw data processing for bins (Directly use DataInfo to trim low coverage tails)
  raw_bin_list <- lapply(names(dataset_list), function(set_name) {
    x <- dataset_list[[set_name]]
    raw_df <- x$data %>% filter(!is.na(genus)) %>% group_by(bin) %>% summarise(raw = n_distinct(genus), .groups = "drop")
    mat <- bin_matrices[[set_name]]
    datainfo <- DataInfo(mat, datatype = "abundance")
    sc_map <- data.frame(bin = colnames(mat), SC = datainfo$SC)
    
    raw_df <- raw_df %>% 
      left_join(sc_map, by = "bin") %>%
      mutate(raw = ifelse(SC < 0.6, NA, raw)) %>%
      select(-SC)
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
  
  # Grouped by type (Raw & Standardised) and added max gap limitation to prevent 50Ma jump
  make_dashed_segments_bin <- function(df) {
    df %>% 
      filter(!is.na(value)) %>% 
      group_by(region, type) %>% 
      arrange(age) %>% 
      mutate(next_age = lead(age), next_val = lead(value)) %>%
      filter(!is.na(next_age), abs(next_age - age) > 10) %>% 
      filter(abs(next_age - age) <= 85) %>% 
      ungroup() %>% 
      select(region, type, x = age, xend = next_age, y = value, yend = next_val)
  }
  
  seg_bin_095_excl <- make_dashed_segments_bin(filter(bin_plot_excl, coverage_type == "C095"))
  seg_bin_min_excl <- make_dashed_segments_bin(filter(bin_plot_excl, coverage_type == "mincov"))
  seg_bin_095_incl <- make_dashed_segments_bin(filter(bin_plot_incl, coverage_type == "C095"))
  seg_bin_min_incl <- make_dashed_segments_bin(filter(bin_plot_incl, coverage_type == "mincov"))
  
  # ---------------- EXPORT PLOT DATA FOR COMBINED FIGURES ----------------
  dir.create("Plot_Data", showWarnings = FALSE, recursive = TRUE)
  bin_095_data <- bin_plot_excl %>% filter(coverage_type == "C095") %>% mutate(Clade = clade_name)
  stage_095_data <- stage_plot_excl$plot_data_095 %>% mutate(Clade = clade_name)
  saveRDS(bin_095_data, paste0("Plot_Data/", tolower(clade_name), "_bin_095.rds"))
  saveRDS(stage_095_data, paste0("Plot_Data/", tolower(clade_name), "_stage_095.rds"))
  
  # ---------------- RENDER PLOTS & EXPORT ----------------
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
  
  for (name in names(main_plots)) ggsave(paste0("Main_Figures/", clade_name, "_Main_", name, ".jpg"), plot = main_plots[[name]], width = 8, height = 5, dpi = 600)
  for (name in names(supp_plots)) ggsave(paste0("Supplementary_Figures/", clade_name, "_Supp_", name, ".jpg"), plot = supp_plots[[name]], width = 8, height = 5, dpi = 600)
  
  # Export raw statistical blocks as RDS
  export_results <- list(
    stage_results = stage_results,
    bin_results = bin_results,
    stage_sc_df = stage_sc_df,
    bin_sc_df = bin_sc_df
  )
  saveRDS(export_results, paste0("Results_", clade_name, ".rds"))
  
  cat(paste("Successfully saved plots and statistical RDS for:", clade_name, "\n"))
}

# ---- Run analyses ----
if (nrow(bivalvia_df) > 0) analyze_clade_evolution(clade_name = "Bivalvia", clade_data = bivalvia_df)
if (nrow(echino_df) > 0) analyze_clade_evolution(clade_name = "Echinodermata", clade_data = echino_df)

cat("\nAll clade-specific evolution analyses completed.\n")

