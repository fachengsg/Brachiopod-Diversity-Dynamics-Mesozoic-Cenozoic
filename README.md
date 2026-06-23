# Brachiopod Diversity Dynamics

This repository contains the R code pipeline used to analyze the taxonomic diversity of Brachiopoda across the Mesozoic and Cenozoic. It integrates fossil occurrence data from the Paleobiology Database (PBDB) and the Geobiodiversity Database (GBDB).

## Pipeline Overview

The analysis is divided into three sequential scripts:

1. **`01_Data_Fetching_and_Merging.R`**: Downloads raw occurrence data for Brachiopoda from the PBDB and formats locally stored GBDB data.
2. **`02_Data_Cleaning_and_Unification.R`**: Removes open nomenclature, resolves synonyms using fuzzy matching, aligns geographic coordinates, assigns strict international stratigraphic stages, and removes cross-database duplicates.
3. **`03_Diversity_Standardisation_and_Visualisation.R`**: Uses coverage-based rarefaction and extrapolation (`iNEXT`) to calculate standardized diversity ($q=0$) at both stage-level and 10-Myr binning. Generates visual comparisons of global vs. regional (China) trends.

## Prerequisites
To run this code, you will need R installed along with the following packages:
* `tidyverse`
* `iNEXT`
* `divDyn`
* `stringdist`

## Usage
Run the scripts in numerical order. Script 02 generates a cleaned dataset (`Brachiopoda_analysis_data.rds`) which is directly imported by Script 03.
