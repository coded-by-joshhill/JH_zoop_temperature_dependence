# JH_zoop_temperature_dependence

A repository for munging zooplankton rate process data and estimating Q10.

The rate process data for this project consists of respiration, ingestion and clearance, and growth. Data was compiled by previous meta-analyses and extracted from more recent scientific articles.

# Folder structure

## Data

Here we save data files such as RDS and/or small csv files. Other data used in the analysis are provided via web links.

## R

This folder holds all scripts that are used to complete data cleaning and analysis. Explanations of each script are as follows:

-   0_Helpers: This script is a home for custom functions. Some of these functions help clean the data and estimate Q10s.

-   1_Scripts: These scripts harmoise and clean the original data so that all rate processes are in a common unit and all carbon weights are in the same unit. These scripts also convert rate processes to mass-specific rates for estimating Q10s.

-   2_Script: Here we use the lat, lon of each recorded rate to create a map to illustrate the global distribution of rates which are used to estimate Q10 values and historic Q10 values.

-   3_Script: This script summarises the compiled historic Q10 values for zooplankton and reports them in a table.

-   4_Script: This script is a work in progress... here we read in the feeding data and prepare it for analysis. We then generate a model to estimate Q10s based on the population-level slopes of each zooplankton group. We then plot those Q10s and some supplementary material Arrhenius plots.

## Output

Here we save any core outputs such as preliminary figures etc.

## Quarto

A folder for any Quarto documents generated from the scripts.

## 
