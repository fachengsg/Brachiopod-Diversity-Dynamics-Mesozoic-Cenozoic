# ==============================================================================
# Script Name: 06_Consolidate_and_Archive_Diversity.R
# Purpose: Consolidate standardized results, perform quality filtering (SC >= 0.6),
#          and archive workspace for Script 07 (Ecological Analysis).
# Authorship: [Your Name/Project Team]
# ==============================================================================

rm(list = ls())
library(tidyverse)
library(openxlsx)

# Set working directory to project root
setwd("D:/PBDB_Project")

# 1. Load Data
clades <- c("Brachiopoda", "Bivalvia", "Echinodermata")
all_results <- list()

for (clade in clades) {
  file_path <- paste0("Results_", clade, ".rds")
  if (file.exists(file_path)) {
    all_results[[clade]] <- readRDS(file_path)
  }
}

# 2. Coverage-based Filtering Logic (Strict Quality Control)
# Define SC threshold for valid data points
SC_THRESHOLD <- 0.6 

# Function to extract and validate results
process_results <- function(level_name, sc_df_name, results_list_name) {
  map_dfr(names(all_results), function(clade) {
    # Extract diagnostic data
    diag_df <- all_results[[clade]][[sc_df_name]]
    # Extract diversity calculations
    res_list <- all_results[[clade]][[results_list_name]]
    
    if (is.null(diag_df) || is.null(res_list)) return(NULL)
    
    # Merge diagnostic info into diversity results
    map_dfr(names(res_list), function(region) {
      div_df <- bind_rows(
        res_list[[region]]$div_095 %>% mutate(SC_Level = "C=0.95"),
        res_list[[region]]$div_min %>% mutate(SC_Level = "C=min")
      ) %>% mutate(Region = region)
      
      # Join with diagnostics to flag unreliable data
      div_df %>%
        left_join(diag_df, by = c("Assemblage")) %>%
        mutate(
          Clade = clade,
          Time_Format = level_name,
          Status = ifelse(Global_incl_SC >= SC_THRESHOLD | China_SC >= SC_THRESHOLD, 
                          "Valid (Kept)", "Excluded (< 0.6 SC)")
        )
    })
  })
}

# Consolidate all levels
final_raw_data <- bind_rows(
  process_results("Stage", "stage_sc_df", "stage_results"),
  process_results("10-Myr Bin", "bin_sc_df", "bin_results")
)

# 3. Export to Publication-Ready Workbook
# Format: Value ± Margin (Symmetric Error)
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

# 4. Archiving for Script 07 (Ecology)
# This saves the filtered data object so Script 07 doesn't need to re-run iNEXT
saveRDS(final_raw_data, "Archive_Consolidated_Diversity_Results.rds")

cat("Script 06 complete: Workbook exported and results archived for Script 07.\n")

