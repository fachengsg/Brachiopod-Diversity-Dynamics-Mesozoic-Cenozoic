# ==============================================================================
# Script Name: 02_Data_Cleaning_and_Unification.R
# Purpose: Clean taxonomy, resolve synonyms, unify PBDB and GBDB,
#          assign stratigraphic stages, deduplicate, and filter environments.
# ==============================================================================

# Clear the workspace. Comment out this line if you want to keep existing objects.
rm(list = ls())

# ---- 0. Setup & Load Data ----
# Install missing packages
required_packages <- c("dplyr", "stringr", "stringdist", "readr", "divDyn")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(dplyr)
library(stringr)
library(stringdist)
library(readr)
library(divDyn)

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

# Load raw merged data (produced by Script 01)
pbdb <- readRDS("PBDB_All_Mesozoic_Cenozoic.rds")
gbdb <- readRDS("GBDB_Mesozoic_Cenozoic_clean.rds")

# ---- 1. Filter Brachiopoda and Select Essential Columns ----
pbdb_brach <- pbdb %>%
  filter(phylum == "Brachiopoda") %>%
  select(
    occurrence_no, collection_no,
    accepted_name, accepted_attr, 
    genus, family, order, class, phylum,
    early_interval, max_ma, min_ma, period,
    lng, lat, paleolng, paleolat,
    cc = `cc...35`,                     # PBDB column name may contain dots; adjust if needed
    environment, motility, life_habit, diet
  )

gbdb_brach <- gbdb %>%
  filter(phylum == "Brachiopoda") %>%
  select(
    occurrence_id, collection_id,
    identified_name,
    genus = pbdb_genus, family, order, class, phylum,
    early_interval, max_ma, min_ma, period,
    lng, lat, paleolng, paleolat,
    country, province,
    formation, paleoenvironment,
    coll_lower_depth, coll_upper_depth,
    feeding_mode, mobility, tiering
  )

# ---- 2. Remove Open Nomenclature ----
# Remove uncertain taxonomic identifications (e.g., cf., aff., sp., ?)
patterns <- c("cf\\.", "aff\\.", "\\?\\ ", " sp\\.", " indet\\.", " ex gr\\.",
              " sensu lato", " spp\\.", " informal", "\\?")

clean_names <- function(df, name_col) {
  df %>%
    filter(
      !grepl(paste(patterns, collapse = "|"), !!sym(name_col), ignore.case = TRUE, useBytes = TRUE),
      !grepl("^\\?+", !!sym(name_col), useBytes = TRUE)
    )
}

pbdb_brach <- clean_names(pbdb_brach, "accepted_name")
gbdb_brach <- clean_names(gbdb_brach, "identified_name")

# ---- 3. Standardize Taxonomy and Capitalization ----
pbdb_brach <- pbdb_brach %>% filter(!is.na(genus) | !is.na(family))
gbdb_brach <- gbdb_brach %>% filter(!is.na(genus) | !is.na(family))

std_taxon <- function(x) str_to_title(str_trim(x))

pbdb_brach <- pbdb_brach %>%
  mutate(genus = std_taxon(genus), family = std_taxon(family),
         accepted_genus = str_extract(accepted_name, "^[A-Za-z]+"))

gbdb_brach <- gbdb_brach %>%
  mutate(genus = std_taxon(genus), family = std_taxon(family))

# ---- 3.5 Compute raw (pre-standardization) taxon counts ----
raw_pbdb_counts <- pbdb_brach %>%
  summarise(
    source = "PBDB",
    occurrences = n(),
    genera = n_distinct(accepted_genus, na.rm = TRUE),
    species = n_distinct(accepted_name[str_detect(accepted_name, " ")], na.rm = TRUE),
    families = n_distinct(family, na.rm = TRUE)
  )

raw_gbdb_counts <- gbdb_brach %>%
  summarise(
    source = "GBDB",
    occurrences = n(),
    genera = n_distinct(genus, na.rm = TRUE),
    species = n_distinct(identified_name[str_detect(identified_name, " ")], na.rm = TRUE),
    families = n_distinct(family, na.rm = TRUE)
  )

raw_combined_summary <- bind_rows(raw_pbdb_counts, raw_gbdb_counts) %>%
  bind_rows(
    tibble(
      source = "Total",
      occurrences = sum(.$occurrences),
      genera = n_distinct(c(pbdb_brach$accepted_genus, gbdb_brach$genus), na.rm = TRUE),
      species = n_distinct(
        c(pbdb_brach$accepted_name[str_detect(pbdb_brach$accepted_name, " ")],
          gbdb_brach$identified_name[str_detect(gbdb_brach$identified_name, " ")]),
        na.rm = TRUE),
      families = n_distinct(c(pbdb_brach$family, gbdb_brach$family), na.rm = TRUE)
    )
  )

cat("\n--- Raw brachiopod data (before synonym merging and deduplication) ---\n")
print(raw_combined_summary)

# ---- 4. Fuzzy Matching on Accepted Genera (PBDB only) ----
detect_similar <- function(vec, threshold = 0.1) {
  unique_v <- unique(na.omit(vec))
  if (length(unique_v) < 2) return(data.frame())
  
  dist_mat <- stringdistmatrix(unique_v, unique_v, method = "jw", p = 0.1, nthread = 1)
  pairs <- which(dist_mat <= threshold & dist_mat > 0, arr.ind = TRUE)
  pairs <- pairs[pairs[,1] < pairs[,2], , drop = FALSE] 
  
  if (nrow(pairs) == 0) return(data.frame())
  
  data.frame(
    genus1   = unique_v[pairs[,1]],
    genus2   = unique_v[pairs[,2]],
    distance = round(dist_mat[pairs], 4)
  ) %>% arrange(distance)
}

similar_accepted <- detect_similar(pbdb_brach$accepted_genus, threshold = 0.1)
cat("\nSimilar accepted genera pairs found:", nrow(similar_accepted), "\n")

# ---- 5. Classify Synonym Pairs and Align GBDB Data ----
synonym_map <- data.frame(target = character(), synonym = character())
pairs_labeled <- data.frame()

if (nrow(similar_accepted) > 0) {
  genus_status <- pbdb_brach %>%
    filter(!is.na(accepted_genus), !is.na(accepted_attr)) %>%
    count(accepted_genus, accepted_attr) %>%
    group_by(accepted_genus) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(accepted_genus, status = accepted_attr)
  
  pairs_labeled <- similar_accepted %>%
    left_join(genus_status, by = c("genus1" = "accepted_genus")) %>% rename(status1 = status) %>%
    left_join(genus_status, by = c("genus2" = "accepted_genus")) %>% rename(status2 = status) %>%
    mutate(
      category = case_when(
        status1 == "synonym" & status2 == "synonym" ~ "both_synonym",
        (status1 == "synonym" & status2 == "valid") | 
          (status1 == "valid"   & status2 == "synonym") ~ "needs_merge",
        status1 == "valid" & status2 == "valid" ~ "both_valid",
        TRUE ~ "other"
      )
    )
  
  to_merge <- pairs_labeled %>% filter(category == "needs_merge")
  if (nrow(to_merge) > 0) {
    synonym_map <- to_merge %>%
      mutate(
        target  = ifelse(status1 == "valid", genus1, genus2),
        synonym = ifelse(status1 == "valid", genus2, genus1)
      ) %>%
      select(target, synonym)
    
    pbdb_brach <- pbdb_brach %>%
      left_join(synonym_map, by = c("accepted_genus" = "synonym")) %>%
      mutate(accepted_genus = ifelse(!is.na(target), target, accepted_genus)) %>%
      select(-target)
    
    gbdb_brach <- gbdb_brach %>%
      left_join(synonym_map, by = c("genus" = "synonym")) %>%
      mutate(genus = ifelse(!is.na(target), target, genus)) %>%
      select(-target)
    
    cat("Merged", nrow(to_merge), "synonym pairs into valid genera across databases.\n")
  }
}

# ---- 5.5 Detect Similar Families ----
unique_families <- unique(na.omit(pbdb_brach$family))
fam_similar <- if(length(unique_families) >= 2) detect_similar(unique_families, 0.1) else data.frame()

# ---- 6. Unify PBDB and GBDB into Analysis Table ----
pbdb_unified <- pbdb_brach %>%
  mutate(source = "PBDB") %>%
  select(occurrence_no, collection_no, accepted_name, accepted_genus, genus, family, order, class, phylum,
         early_interval, max_ma, min_ma, period, lng, lat, paleolng, paleolat,
         environment, motility, life_habit, diet, cc, source)

gbdb_unified <- gbdb_brach %>%
  rename(
    occurrence_no = occurrence_id,
    collection_no = collection_id,
    accepted_name = identified_name,
    environment = paleoenvironment,         
    motility = mobility,            
    diet = feeding_mode,            
    life_habit = tiering            
  ) %>%
  mutate(
    accepted_genus = genus, 
    source = "GBDB",
    cc = NA_character_
  ) %>%
  select(occurrence_no, collection_no, accepted_name, accepted_genus, genus, family, order, class, phylum,
         early_interval, max_ma, min_ma, period, lng, lat, paleolng, paleolat,
         environment, motility, life_habit, diet, cc, source, 
         country, province, formation, coll_lower_depth, coll_upper_depth)

analysis_data <- bind_rows(pbdb_unified, gbdb_unified) %>%
  mutate(across(c(max_ma, min_ma, lng, lat, paleolng, paleolat, coll_lower_depth, coll_upper_depth), as.numeric))

# ---- 6.5 Count taxa after merging but before deduplication ----
cat("\n--- After standardization and merging (before deduplication) ---\n")
merged_counts <- analysis_data %>%
  summarise(
    occurrences = n(),
    genera = n_distinct(accepted_genus, na.rm = TRUE),
    species = n_distinct(accepted_name[str_detect(accepted_name, " ")], na.rm = TRUE),
    families = n_distinct(family, na.rm = TRUE)
  )
print(merged_counts)

# ---- 7. Strict Stage Assignment (Using divDyn) ----
analysis_data$mid_age <- (analysis_data$max_ma + analysis_data$min_ma) / 2

data(stages)
stages_ph <- stages %>%
  filter(stg >= 250 | system %in% c("Triassic", "Jurassic", "Cretaceous",
                                    "Paleogene", "Neogene", "Quaternary"))

assign_stage_vec <- function(age, stages_df) {
  sapply(age, function(x) {
    if (is.na(x)) return(NA_character_)
    idx <- which(x >= stages_df$top & x <= stages_df$bottom)[1]
    if (length(idx) == 0) return(NA_character_)
    stages_df$stage[idx]
  })
}

analysis_data$stage <- assign_stage_vec(analysis_data$mid_age, stages_ph)
analysis_data <- filter(analysis_data, !is.na(stage))

# ---- 7.1 Count after stage assignment but before deduplication ----
cat("\n--- After stage assignment (before deduplication) ---\n")
pre_dedup_counts <- analysis_data %>%
  summarise(
    occurrences = n(),
    genera = n_distinct(accepted_genus, na.rm = TRUE),
    species = n_distinct(accepted_name[str_detect(accepted_name, " ")], na.rm = TRUE),
    families = n_distinct(family, na.rm = TRUE)
  )
print(pre_dedup_counts)

# ---- 8. Identify and Remove Cross-Database Duplicates ----
cat("\n--- Cross-database duplicate diagnostic ---\n")

diag_data <- analysis_data %>%
  mutate(
    coord_lat = if_else(!is.na(lat) & !is.na(lng), round(lat, 1), NA_real_),
    coord_lng = if_else(!is.na(lat) & !is.na(lng), round(lng, 1), NA_real_)
  )

match_summary <- diag_data %>%
  group_by(accepted_genus, stage, coord_lat, coord_lng) %>%
  summarise(
    n_records = n(),
    sources   = paste(unique(source), collapse = "+"),
    .groups   = "drop"
  )

mixed_groups <- match_summary %>% filter(n_records > 1, grepl("PBDB", sources) & grepl("GBDB", sources))
cat("Match groups (genus + stage + 0.1° coord.) containing both PBDB and GBDB records:", nrow(mixed_groups), "\n")

if (nrow(mixed_groups) > 0) {
  gbdb_to_remove <- diag_data %>%
    inner_join(mixed_groups %>% select(accepted_genus, stage, coord_lat, coord_lng),
               by = c("accepted_genus", "stage", "coord_lat", "coord_lng")) %>%
    filter(source == "GBDB")
  
  analysis_data <- analysis_data %>% anti_join(gbdb_to_remove, by = c("occurrence_no", "source"))
  cat("Removed", nrow(gbdb_to_remove), "GBDB duplicate records (PBDB given priority).\n")
}

analysis_data <- distinct(analysis_data)

# ---- 8.5 Strict Environment Filtering ----
# Remove records from non-marine or highly suspicious transitional settings
cat("\n--- Applying strict marine environment filter ---\n")
non_marine_terms <- "terrestrial|fluvial|lacustrine|floodplain|channel|transitional|paralic|volcanic|delta|estuary"

analysis_data <- analysis_data %>%
  filter(!str_detect(str_to_lower(environment), non_marine_terms))

cat("Non-marine and transitional records removed. Proceeding with analysis.\n")

# ---- 8.6 Environment Data Audit ----
cat("\n--- Running Environment Data Audit (Post-filter) ---\n")

env_audit <- analysis_data %>%
  select(source, environment) %>%
  filter(!is.na(environment)) %>%
  mutate(environment = str_to_lower(environment))

env_summary <- env_audit %>%
  separate_rows(environment, sep = ",\\s*|;\\s*") %>%
  group_by(source, environment) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(source, desc(count))

write_csv(env_summary, "Audit_Environment_Terminology.csv")
cat("Audit complete. Please check 'Audit_Environment_Terminology.csv' for data quality.\n")

# ---- 9. Export Outputs ----
saveRDS(analysis_data, "Brachiopoda_analysis_data.rds")
write_csv(analysis_data, "Brachiopoda_analysis_data.csv")

if (nrow(pairs_labeled) > 0) write_csv(pairs_labeled, "similar_genus_pairs_all.csv")
if (nrow(fam_similar) > 0) write_csv(fam_similar, "similar_family_pairs.csv")

# ---- 9.1 Final Dataset Summary (after deduplication and environment filtering) ----
cat("\n--- Final Brachiopod Dataset (after deduplication and marine filter) ---\n")
final_counts <- analysis_data %>%
  summarise(
    total_occurrences = n(),
    unique_genera     = n_distinct(accepted_genus, na.rm = TRUE),
    unique_species    = n_distinct(
      accepted_name[str_detect(accepted_name, " ")], na.rm = TRUE),
    unique_families   = n_distinct(family, na.rm = TRUE)
  )
print(final_counts)

cat("\nTotal records:", nrow(analysis_data),
    "(PBDB:", sum(analysis_data$source == "PBDB"),
    "| GBDB:", sum(analysis_data$source == "GBDB"), ")\n")
