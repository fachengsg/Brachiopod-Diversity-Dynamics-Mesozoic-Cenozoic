# ==============================================================================
# Script Name: 03_Ecological_and_Comparative_Processing.R
# Purpose: 
#   1. Harmonize ecological vocabulary between PBDB and GBDB.
#   2. Implement Genus-level Trait Imputation (Taxonomic Inheritance) to rescue 
#      occurrences with NA ecological traits.
#   3. Process and export clean Brachiopoda data.
#   4. Extract, harmonize, and filter comparable Bivalvia & Echinodermata data.
# ==============================================================================

rm(list = ls())

# ---- 1. Environment Setup & Packages ----
library(tidyverse)

setwd("D:/PBDB_Project") 

cat("--- Starting Ecological and Comparative Processing with Trait Imputation ---\n")

# ---- 2. Universal Harmonization Function ----
harmonize_ecology <- function(df) {
  df %>%
    mutate(
      diet = case_when(
        diet == "Filter/suspension" ~ "suspension feeder",
        TRUE ~ diet
      ),
      motility = case_when(
        motility %in% c("Non-motile/attached", "Non-motile/unattached", 
                        "stationary", "stationary, attached", 
                        "stationary, attached, epibiont") ~ "stationary",
        TRUE ~ motility
      ),
      life_habit = case_when(
        life_habit %in% c("Surficial", "Erect", "epifaunal", 
                          "low-level epifaunal", "epifaunal, gregarious") ~ "epifaunal",
        life_habit %in% c("Shallow_infaunal", "Semi-infaunal") ~ "infaunal",
        TRUE ~ life_habit
      )
    )
}

# ==============================================================================
# PART A: BRACHIOPODA PROCESSING & SAFE ECO-PURIFICATION
# ==============================================================================
cat("\n--- Processing Brachiopoda Data ---\n")

if (!file.exists("Brachiopoda_analysis_data.rds")) stop("Raw dataset not found.")
brach_raw <- readRDS("Brachiopoda_analysis_data.rds")
cat("Brachiopoda raw records loaded:", nrow(brach_raw), "\n")

# 1. Harmonize vocabulary
brach_harmonized <- harmonize_ecology(brach_raw)

# 2. SAFE PURIFICATION: Instead of a generic whitelist, we explicitly
#    REMOVE rows that are confirmed to be non-comparable niches (like Infaunal Lingula)
brach_data_final <- brach_harmonized %>%
  filter(
    # Keep if it is our target diet, OR if the diet info is missing (NA)
    (diet == "suspension feeder" | is.na(diet)),
    
    # Keep if it is stationary, OR if motility info is missing (NA)
    (motility == "stationary" | is.na(motility)),
    
    # CRITICAL: Explicitly drop entries that are confirmed Infaunal (like Lingula)
    # If life_habit is NA, we safely keep it because >95% of brachiopods are epifaunal
    (life_habit == "epifaunal" | is.na(life_habit))
  )

cat("Brachiopoda records retained after safe eco-purification:", nrow(brach_data_final), "\n")

saveRDS(brach_data_final, "Brachiopoda_analysis_data_Final.rds")
write_csv(brach_data_final, "Brachiopoda_analysis_data_Final.csv")

# ==============================================================================
# PART B: COMPARATIVE TAXA (BIVALVIA & ECHINODERMATA) PROCESSING
# ==============================================================================
cat("\n--- Processing Comparative Taxa (Bivalvia & Echinodermata) ---\n")

pbdb_raw <- readRDS("PBDB_All_Mesozoic_Cenozoic.rds")
gbdb_raw <- readRDS("GBDB_Mesozoic_Cenozoic_clean.rds")

# Extract PBDB
pbdb_comp <- pbdb_raw %>%
  filter(class == "Bivalvia" | phylum == "Echinodermata") %>%
  mutate(source = "PBDB", country = NA_character_) %>%
  select(occurrence_no, collection_no, accepted_name, genus, family, order, class, phylum,
         early_interval, max_ma, min_ma, period, lng, lat, paleolng, paleolat,
         environment, motility, life_habit, diet, source, 
         cc = `cc...35`, country) 

# Extract GBDB
gbdb_comp <- gbdb_raw %>%
  filter(class == "Bivalvia" | phylum == "Echinodermata") %>%
  rename(
    occurrence_no = occurrence_id,
    collection_no = collection_id,
    accepted_name = identified_name,
    genus = pbdb_genus,             
    environment = paleoenvironment,
    motility = mobility,            
    diet = feeding_mode,            
    life_habit = tiering            
  ) %>%
  mutate(source = "GBDB", cc = NA_character_) %>%
  select(occurrence_no, collection_no, accepted_name, genus, family, order, class, phylum,
         early_interval, max_ma, min_ma, period, lng, lat, paleolng, paleolat,
         environment, motility, life_habit, diet, source, 
         cc, country)

# Combine and Harmonize
comp_harmonized <- bind_rows(pbdb_comp, gbdb_comp) %>%
  mutate(across(c(max_ma, min_ma, lng, lat, paleolng, paleolat), as.numeric)) %>%
  harmonize_ecology()

# ------------------------------------------------------------------------------
# Genus-Level Taxonomic Imputation for Comparative Data
# ------------------------------------------------------------------------------

# 1. Identify Target Genera: Which genera EVER exhibit the exact target ecology?
target_comp_genera <- comp_harmonized %>%
  filter(
    diet == "suspension feeder",
    motility == "stationary",
    life_habit == "epifaunal"
  ) %>%
  drop_na(genus) %>%
  pull(genus) %>%
  unique()

# 2. Rescue Occurrences: Extract ALL records of these target genera
comparable_data_rescued <- comp_harmonized %>%
  filter(genus %in% target_comp_genera)

cat("\n--- Final Comparable Dataset Summary (Post-Imputation) ---\n")
cat("Total comparative Genera captured:", length(target_comp_genera), "\n")
print(comparable_data_rescued %>% count(phylum, class, source))

saveRDS(comparable_data_rescued, "Comparative_Bivalvia_Echinodermata.rds")
write_csv(comparable_data_rescued, "Comparative_Bivalvia_Echinodermata.csv")
cat("\nPipeline finished! Trait imputation successful. Database saved.\n")
