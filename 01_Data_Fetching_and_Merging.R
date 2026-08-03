# ==============================================================================
# Script Name: 01_Data_Fetching_and_Merging.R
# Purpose: Fetch Mesozoic–Cenozoic brachiopod, bivalve, and echinoderm
#          occurrences from the Paleobiology Database (PBDB) and integrate
#          local data from the Geobiodiversity Database (GBDB).
# Data availability:
#   - PBDB data are fetched directly from the PBDB API.
#   - GBDB data are archived on Zenodo (https://doi.org/10.5281/zenodo.21765787)
#     and will be automatically downloaded if not present in the working directory.
# ==============================================================================

# Clear the workspace. Comment out this line if you want to keep existing objects.
rm(list = ls())

# ---- 1. Environment Setup ----
# Install missing packages
required_packages <- c("httr", "readr", "dplyr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(httr)
library(readr)
library(dplyr)
library(tidyr)

# --- Working directory: choose one of the two options below ---
# Option A: Set a fixed directory (e.g., "D:/PBDB_Project" on Windows).
# Option B: Leave custom_dir as NULL or "" to automatically use the
#           folder that contains this script (works in RStudio).
custom_dir <- "D:/PBDB_Project"   # <-- Change this to your preferred directory, or set to NULL

if (!is.null(custom_dir) && nchar(custom_dir) > 0) {
  # Use the user-specified directory
  if (!dir.exists(custom_dir)) {
    dir.create(custom_dir, recursive = TRUE)
    cat("Folder created:", custom_dir, "\n")
  }
  setwd(custom_dir)
} else {
  # Automatically use the script's location (requires RStudio)
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
    if (nchar(script_path) > 0) {
      setwd(dirname(script_path))
    }
  }
}
cat("Working directory:", getwd(), "\n")

# ---- 2. PBDB Download Function ----
# Downloads occurrences in pages, with a pause between requests to be polite.
fetch_pbdb_simple <- function(taxon_group, interval, outfile, limit = 5000, sleep = 0.5) {
  base_url <- "https://paleobiodb.org/data1.2/occs/list.csv"
  all_pages <- list()
  offset <- 0
  page <- 1
  total_so_far <- 0  # running count of records downloaded
  
  cat("Starting download:", taxon_group, "-", interval, "\n")
  
  repeat {
    cat("  Page", page, "(offset", offset, ") ... ")
    
    resp <- tryCatch(
      GET(base_url,
          query = list(base_name = taxon_group,
                       interval  = interval,
                       show      = "full",
                       limit     = limit,
                       offset    = offset)),
      error = function(e) {
        cat("\nNetwork error:", e$message, "\n")
        return(NULL)
      }
    )
    
    if (is.null(resp) || status_code(resp) != 200) {
      warning("Failed to retrieve page at offset ", offset)
      break
    }
    
    csv_text <- content(resp, as = "text", encoding = "UTF-8")
    if (nchar(csv_text) == 0) break
    
    # Read all columns as character to avoid type conflicts
    df <- read_csv(csv_text, col_types = cols(.default = col_character()),
                   show_col_types = FALSE, progress = FALSE)
    
    n_rows <- nrow(df)
    if (n_rows == 0) break
    
    all_pages[[length(all_pages) + 1]] <- df
    total_so_far <- total_so_far + n_rows
    cat(n_rows, "records (cumulative:", total_so_far, ")\n")
    
    if (n_rows < limit) break
    offset <- offset + limit
    page <- page + 1
    Sys.sleep(sleep)
  }
  
  if (length(all_pages) == 0) {
    cat("No records found.\n")
    return(invisible(NULL))
  }
  
  final_data <- bind_rows(all_pages)
  write_csv(final_data, outfile)
  cat("Download completed. Total records saved:", nrow(final_data), "\n\n")
  return(final_data)
}

# ---- 3. Download PBDB Data ----
# Define target groups and geological periods
taxa <- c("Brachiopoda", "Bivalvia", "Echinodermata")
periods <- c("Triassic", "Jurassic", "Cretaceous", "Paleogene", "Neogene", "Quaternary")

# The download loop is commented out to prevent accidental large re-downloads.
# Uncomment the lines below when you need to fetch fresh data.
# for (t in taxa) {
#   for (p in periods) {
#     outfile_name <- paste0(p, "_", t, ".csv")
#     if (!file.exists(outfile_name)) {
#       fetch_pbdb_simple(taxon_group = t, interval = p, outfile = outfile_name)
#     } else {
#       cat("File already exists, skipping:", outfile_name, "\n")
#     }
#   }
# }

# ---- 4. Merge PBDB Data ----
cat("\n--- Merging PBDB Data ---\n")
tasks <- expand.grid(phylum = taxa, period = periods, stringsAsFactors = FALSE)
tasks$file <- paste0(tasks$period, "_", tasks$phylum, ".csv")

# Check for missing files
missing <- tasks$file[!file.exists(tasks$file)]
if (length(missing) > 0) {
  warning("Missing files: ", paste(missing, collapse = ", "))
  tasks <- tasks[file.exists(tasks$file), ]
}

# Read each file, add the period column, then bind all rows together
pbdb_combined <- lapply(seq_len(nrow(tasks)), function(i) {
  df <- read_csv(tasks$file[i], col_types = cols(.default = col_character()), progress = FALSE)
  df$period <- tasks$period[i]
  return(df)
}) %>% bind_rows()

# Make sure the expected column 'phylum' exists
stopifnot("phylum" %in% colnames(pbdb_combined))

cat("Merged PBDB records:", nrow(pbdb_combined), "\n")

# Basic summary
cat("\n--- PBDB Records by Phylum ---\n")
pbdb_combined %>% count(phylum, sort = TRUE) %>% print()

cat("\n--- Cross-table: PBDB Phylum × Period ---\n")
pbdb_combined %>%
  count(phylum, period) %>%
  pivot_wider(names_from = period, values_from = n, values_fill = list(n = 0)) %>%
  print()

# Save merged PBDB dataset
write_csv(pbdb_combined, "PBDB_All_Mesozoic_Cenozoic.csv")
saveRDS(pbdb_combined, "PBDB_All_Mesozoic_Cenozoic.rds")
cat("Merged PBDB data saved.\n")

# ---- 5. Integrate GBDB Data ----
cat("\n--- Processing GBDB Data ---\n")
# GBDB data are archived on Zenodo (DOI: 10.5281/zenodo.21765787).
# The script will attempt to download them if they are not already present locally.
gbdb_files <- c("Bivalvia 0-252 GBDB.csv",
                "Echinodermata 0-252 GBDB.csv",
                "Brachiopoda 0-252 GBDB.csv")

# Download any missing GBDB files from Zenodo
for (f in gbdb_files) {
  if (!file.exists(f)) {
    cat("Downloading", f, "from Zenodo...\n")
    url <- paste0("https://zenodo.org/record/21765787/files/", URLencode(f))
    tryCatch({
      download.file(url, destfile = f, mode = "wb")
    }, error = function(e) {
      warning("Failed to download ", f, ": ", e$message)
    })
  }
}

if (all(file.exists(gbdb_files))) {
  gbdb_raw <- lapply(gbdb_files, function(f) {
    read_csv(f, col_types = cols(.default = col_character()), show_col_types = FALSE)
  }) %>% bind_rows()
  
  target_phyla <- c("Brachiopoda", "Echinodermata", "Mollusca")
  
  gbdb_clean <- gbdb_raw %>%
    filter(
      !is.na(phylum),
      phylum %in% target_phyla,
      period %in% periods
    ) %>%
    # Keep only Bivalvia within Mollusca
    filter(
      (phylum == "Mollusca" & class == "Bivalvia") | phylum != "Mollusca"
    )
  
  cat("GBDB records after filtering:", nrow(gbdb_clean), "\n")
  
  cat("\n--- Cross-table: GBDB Phylum × Period ---\n")
  gbdb_clean %>%
    count(phylum, period) %>%
    pivot_wider(names_from = period, values_from = n, values_fill = list(n = 0)) %>%
    print()
  
  write_csv(gbdb_clean, "GBDB_Mesozoic_Cenozoic_clean.csv")
  saveRDS(gbdb_clean, "GBDB_Mesozoic_Cenozoic_clean.rds")
  cat("Filtered GBDB data saved.\n")
} else {
  warning("GBDB files not found and could not be downloaded. Skipping GBDB integration.")
}
