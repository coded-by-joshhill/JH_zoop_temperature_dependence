# JH_zoop_temperature_dependence

A repository for 'A global synthesis of the temperature dependence of feeding, growth, and metabolism across marine zooplankton' - data cleaning and analysis estimating Q10.

The rate process data for this project consists of respiration, ingestion and clearance, and growth. Data was compiled by previous scientific articles.

# Folder structure

## Data

Here we save data files such as RDS and CSV.

- CSV files are the original, raw data files of rate observations

- RDS files are the cleaned data, model data frames and model estimates, which are used throughout the scripts

## R

This folder contains all .R scripts that are used to for data cleaning and data analysis. Explanations of each script are outlined below:

- **0_Helpers:** a source file that stores custom functions and is used throughout data cleaning and data analyses.

- **1.1_cleaning_feeding_data:** This script harmonises and cleans the raw rate process (ingestion and clearance) data so that all rates are in a common unit, and all carbon weights are in the same unit (mg). This script also converts feeding rates to mass-specific rates in preparation for estimating temperature dependence (using Q10).

  - First, we assign unique identifier to each row specific to the dataset.

  - Subset the data by taxa (genus or species depending on the record), extract AphiaIDs using the World Register of Marine Species *(worrms)* package. Add taxon classifications using AphiaIDs.

  - Create customised groupings for all zooplankton groups following similar naming conventions to Ikeda (2014) and Kiørboe and Hirst (2014).

  - Harmonise and clean the data using conversion functions in the 0_Helpers.R source file.

  - Convert harmonised rates to carbon mass-specific rate using functions in the 0_Helper.R source file.

  - Finally, save the data as an RDS file for later use.

- **1.2_cleaning_growth_data:** This script follows script 1.1 format and workflow to harmonise and clean growth rate data.

- **1.3_cleaning_respiration_data:** This script follows script 1.1 format and workflow to harmonise and clean respiration rate data.

- **1.4_cleaning_excretion_data:** This script follows script 1.1 format and workflow to harmonise and clean excretion rate data.

- **2.0_taxaCoverage_suppInfo:** This script generates supplementary information and figures. Such as taxonomic coverage across all rate datasets, primary author table etc.

- **3.0_overallZ_model:** data analysis for zooplankton overall. Estimates of Q10 are generated here.

- **3.1_overallZ_model_dryMass:** a version of data analysis for zooplankton overall, but as dryMass-specific rates.

- **4_sizeGrps_model:** data analysis for size groups. Estimates of Q10 are generated here.

- **5_funcGrps_model:** data analysis for functional groups. Estimates of Q10 are generated here.

- **6_zoopGrps_model:** data analysis for various zooplankton groups. Estimates of Q10 are generated here.

- **7_ratio_estimates:** a script that combines all data produced in previous scripts and then estimates the ratios of the rates (i.e., ecological efficiencies of zooplankton), and associated temperature dependence of the ratios.

## Output

A folder to save any RAW outputs such as preliminary figures and tables etc.

## Final figures

A folder containing our publication-ready figures.
