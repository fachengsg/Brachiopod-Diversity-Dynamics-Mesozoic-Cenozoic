# ==============================================================================
# Script Name: 02_Data_Cleaning_and_Unification.R
# Purpose: Clean taxonomy (remove open nomenclature), resolve synonyms via 
#          fuzzy matching, unify PBDB & GBDB datasets, assign strict 
#          stratigraphic stages, and perform cross-database deduplication.
# ==============================================================================

rm(list = ls())

# ---- 0. Setup & Load Data ----
library(dplyr)
library(stringr)
library(stringdist)
library(readr)
library(divDyn) # Added for strict stage assignment

# Set working directory
setwd("D:/PBDB_Project")

# Load raw merged data (imported as character from previous step)
pbdb <- readRDS("PBDB_All_Mesozoic_Cenozoic.rds")
gbdb <- readRDS("GBDB_Mesozoic_Cenozoic_clean.rds")

# ---- 1. Filter Brachiopoda & Select Essential Columns ----
pbdb_brach <- pbdb %>%
  filter(phylum == "Brachiopoda") %>%
  select(
    occurrence_no, collection_no,
    accepted_name, accepted_attr, 
    genus, family, order, class, phylum,
    early_interval, max_ma, min_ma, period,
    lng, lat, paleolng, paleolat,
    cc = `cc...35`,
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

# ---- 3. Standardise Taxonomy & Capitalisation ----
pbdb_brach <- pbdb_brach %>% filter(!is.na(genus) | !is.na(family))
gbdb_brach <- gbdb_brach %>% filter(!is.na(genus) | !is.na(family))

std_taxon <- function(x) str_to_title(str_trim(x))

pbdb_brach <- pbdb_brach %>%
  mutate(genus = std_taxon(genus), family = std_taxon(family),
         accepted_genus = str_extract(accepted_name, "^[A-Za-z]+"))

gbdb_brach <- gbdb_brach %>%
  mutate(genus = std_taxon(genus), family = std_taxon(family))

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

# ---- 5. Classify Synonym Pairs & Align GBDB ----
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
    diet = feeding_mode       
  ) %>%
  mutate(
    accepted_genus = genus, 
    source = "GBDB",
    life_habit = NA_character_,
    cc = NA_character_
  ) %>%
  select(occurrence_no, collection_no, accepted_name, accepted_genus, genus, family, order, class, phylum,
         early_interval, max_ma, min_ma, period, lng, lat, paleolng, paleolat,
         environment, motility, life_habit, diet, cc, source, 
         country, province, formation, coll_lower_depth, coll_upper_depth, tiering)

analysis_data <- bind_rows(pbdb_unified, gbdb_unified) %>%
  mutate(across(c(max_ma, min_ma, lng, lat, paleolng, paleolat, coll_lower_depth, coll_upper_depth), as.numeric))

# ---- 7. Strict Stage Assignment (Using divDyn) ----
# Calculate midpoint age to assign robust international stages
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

# ---- 8. Cross-Database Duplicate Diagnostic & Removal ----
# Detects overlapping records between PBDB and GBDB based on genus, stage, and spatial coordinates.
cat("\n--- Cross-database duplicate diagnostic ---\n")

diag_data <- analysis_data %>%
  mutate(
    coord_lat = if_else(!is.na(paleolat) & !is.na(paleolng), round(paleolat, 1),
                        if_else(!is.na(lat) & !is.na(lng), round(lat, 1), NA_real_)),
    coord_lng = if_else(!is.na(paleolat) & !is.na(paleolng), round(paleolng, 1),
                        if_else(!is.na(lat) & !is.na(lng), round(lng, 1), NA_real_))
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
  # Identify GBDB duplicates that overlap with PBDB
  gbdb_to_remove <- diag_data %>%
    inner_join(mixed_groups %>% select(accepted_genus, stage, coord_lat, coord_lng),
               by = c("accepted_genus", "stage", "coord_lat", "coord_lng")) %>%
    filter(source == "GBDB")
  
  # Remove the overlapping GBDB records (giving priority to PBDB)
  analysis_data <- analysis_data %>% anti_join(gbdb_to_remove, by = c("occurrence_no", "source"))
  cat("Removed", nrow(gbdb_to_remove), "GBDB duplicate records (PBDB given priority).\n")
} else {
  cat("No cross-database duplicates detected. All records retained.\n")
}

# Final safety check: remove exact duplicate rows
analysis_data <- distinct(analysis_data)

# ---- 9. Export Outputs ----
saveRDS(analysis_data, "Brachiopoda_analysis_data.rds")
write_csv(analysis_data, "Brachiopoda_analysis_data.csv")

if (nrow(pairs_labeled) > 0) write_csv(pairs_labeled, "similar_genus_pairs_all.csv")
if (nrow(fam_similar) > 0) write_csv(fam_similar, "similar_family_pairs.csv")

cat("\n--- Final Dataset Summary ---\n")
cat("Total records:", nrow(analysis_data), "(PBDB:", sum(analysis_data$source == "PBDB"), "| GBDB:", sum(analysis_data$source == "GBDB"), ")\n")

analysis_data %>%
  summarise(
    unique_genera = n_distinct(accepted_genus, na.rm = TRUE),
    unique_families = n_distinct(family, na.rm = TRUE),
    unique_species = n_distinct(accepted_name, na.rm = TRUE)
  ) %>%
  print()

