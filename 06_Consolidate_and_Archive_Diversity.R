# ==============================================================================
# Script Name: 06_Consolidate_and_Archive_Diversity.R
# Purpose: Consolidate standardized diversity results, apply quality filtering
#          (SC >= 0.6, consistent with Scripts 04 and 05), and archive the
#          workspace for Script 07 (Ecological Analysis).
# ==============================================================================

# Clear the workspace. Comment out this line if you want to keep existing objects.
rm(list = ls())

# ---- 1. Environment Setup & Packages ----
required_packages <- c("dplyr", "purrr", "tidyr", "openxlsx")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(purrr)
library(tidyr)
library(openxlsx)

# --- Working directory: choose one of the two options below ---
# Option A: Set a fixed directory (e.g., "D:/PBDB_Project" on Windows).
# Option B: Leave custom_dir as NULL or "" to automatically use the
#           folder that contains this script (works in RStudio).
custom_dir <- "D:/PBDB_Project"   # <-- Change this to the preferred directory, or set to NULL

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

# ---- 2. Load Results from Scripts 04 and 05 ----
clades <- c("Brachiopoda", "Bivalvia", "Echinodermata")
all_results <- list()

for (clade in clades) {
  file_path <- paste0("Results_", clade, ".rds")
  if (file.exists(file_path)) {
    all_results[[clade]] <- readRDS(file_path)
    cat("Loaded:", file_path, "\n")
  } else {
    warning("File not found: ", file_path)
  }
}

if (length(all_results) == 0) {
  stop("No result files found. Please run Scripts 04 and 05 first.")
}

# ---- 3. Quality Filtering: Sample Coverage (SC) >= 0.6 ----
SC_THRESHOLD <- 0.6

# Helper function to extract and validate diversity estimates for one time level.
# Applies the rule established in Scripts 04 and 05: 
# Retains the time bin/stage if ANY region (Global_incl, Global_excl, China) 
# achieves a Sample Coverage (SC) >= 0.6.
process_results <- function(level_name, sc_df_name, results_list_name) {
  map_dfr(names(all_results), function(clade) {
    diag_df <- all_results[[clade]][[sc_df_name]]       # SC data frame
    res_list <- all_results[[clade]][[results_list_name]] # Diversity lists per region
    
    if (is.null(diag_df) || is.null(res_list)) return(NULL)
    
    map_dfr(names(res_list), function(region) {
      # Combine C = 0.95 and C = min estimates for this region
      div_df <- bind_rows(
        res_list[[region]]$div_095 %>% mutate(SC_Level = "C=0.95"),
        res_list[[region]]$div_min %>% mutate(SC_Level = "C=min")
      ) %>% mutate(Region = region)
      
      # Merge with sample coverage and flag records below threshold
      div_df %>%
        left_join(diag_df, by = "Assemblage") %>%
        mutate(
          Clade = clade,
          Time_Format = level_name,
          # A stage/bin is marked as valid if ANY of the three regions has SC >= 0.6.
          # NA values (missing data for a region) are treated as not meeting the threshold.
          Status = ifelse(
            (Global_incl_SC >= SC_THRESHOLD & !is.na(Global_incl_SC)) |
              (Global_excl_SC >= SC_THRESHOLD & !is.na(Global_excl_SC)) |
              (China_SC       >= SC_THRESHOLD & !is.na(China_SC)),
            "Valid (Kept)",
            "Excluded (< 0.6 SC)"
          )
        )
    })
  })
}

# Consolidate stage-level and 10-Myr bin-level results
final_raw_data <- bind_rows(
  process_results("Stage", "stage_sc_df", "stage_results"),
  process_results("10-Myr Bin", "bin_sc_df", "bin_results")
)

# ---- 4. Create Publication-Ready Workbook ----
# Add a formatted value +/- margin column
publication_table <- final_raw_data %>%
  mutate(
    Margin = round((qD.UCL - qD.LCL) / 2, 1),
    Value_String = paste0(round(qD, 1), " ± ", Margin)
  ) %>%
  pivot_wider(names_from = Order.q, values_from = Value_String, names_prefix = "q=") %>%
  select(Clade, Time_Format, Region, SC_Level, Assemblage, `q=0`, `q=1`, `q=2`, Status)

wb <- createWorkbook()
addWorksheet(wb, "Table1_Publication_Diversity")
addWorksheet(wb, "Table2_Raw_Diversity_Diagnostic")
writeData(wb, "Table1_Publication_Diversity", publication_table)
writeData(wb, "Table2_Raw_Diversity_Diagnostic", final_raw_data)
saveWorkbook(wb, "Publication_Diversity_Data.xlsx", overwrite = TRUE)

# ---- 5. Archive for Script 07 (Ecological Analysis) ----
saveRDS(final_raw_data, "Archive_Consolidated_Diversity_Results.rds")

cat("\nProcessing complete. Workbook exported and results archived for Script 07.\n")
