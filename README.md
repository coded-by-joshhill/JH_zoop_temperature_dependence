# JH_zoop_temperature_dependence

A repository for 'A global synthesis of the temperature dependence of feeding, growth, and metabolism across marine zooplankton' - data cleaning and analysis estimating Q10.

The rate process data for this project consists of respiration, ingestion and clearance, and growth. Data was compiled by previous scientific articles.

# Folder structure

## Data

Here we save data files such as RDS. Other data used in the analysis are provided via web links.

+-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| Files                 | Description                                                                                                                               |
+=======================+===========================================================================================================================================+
| x_y_z_data.rds &\     | These data files are cleaned rate process and Q10 data.                                                                                   |
| historicQ10_dat.rds   |                                                                                                                                           |
+-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| xyz_mdat.rds          | These files are modelling data frames that are used for plotting.                                                                         |
+-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| Q10_estimates_xyz.rds | These files are Q10 estimate output generated using bootstrap models.                                                                     |
+-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| Q10_summary_xyz.rds   | These files are Q10 summaries which include Q10 estimates and confidence intervals which are estimated using the bootstrap models output. |
+-----------------------+-------------------------------------------------------------------------------------------------------------------------------------------+

## R

This folder contains all .R scripts that are used to perform data cleaning and analysis. Explanations of each script are outlined below:

-   **0_Helpers:** This script is a source file which stores custom functions used throughout data cleaning and analyses phases.

-   **1.1_cleaning_feeding_data:** This script harmonises and cleans the raw rate process (ingestion and clearance) data so that all rates are in a common unit, and all carbon weights are in the same unit (mg). This script also converts feeding rates to mass-specific rates in preparation for estimating temperature dependence (using Q~10~).

    -   First, we assign unique identifier to each row specific to the dataset.

    -   Subset the data by taxa (genus or species depending on the record), extract AphiaIDs using the World Register of Marine Species *(worrms)* package. Add taxon classifications using AphiaIDs.

    -   Create customised groupings for all zooplankton groups following similar naming conventions to Ikeda (2014) and Kiørboe and Hirst (2014).

    -   Harmonise and clean the data using conversion functions in the 0_Helpers.R source file.

    -   Convert harmonised rates to carbon mass-specific rate using functions in the 0_Helper.R source file.

    -   Finally, save the data as an RDS file for later use.

-   **1.2_cleaning_growth_data:** This script follows script 1.1 format and workflow to harmonise and clean growth rate data.

-   **1.3_cleaning_respiration_data:** This script follows script 1.1 format and workflow to harmonise and clean respiration rate data.

-   **1.4_cleaning_Q10_data:** This script follows a similar format and workflow to script 1.1-3 to clean historic Q~10~ data.

-   **1.5_assess_taxonomic_coverage:** This script is a work in progress... but is planned to read all datasets in, subset by taxa and perform a taxonomic coverage assessment to resolve zooplankton species bias across rates and zooplankton groups.

-   **2_mapping_the_experiments:** First, this script loads all rate process datasets and historic Q~10~ data. Reads in shapefiles using the *rnaturalearth* package. Converts compiled data lat/lon's to sf for plotting. Finally, generates a series of maps to illustrate the distribution of rate processes and historic temperature dependence data across Earth.

-   **3_table1_historicQ10_data**: This script summarises the compiled historic Q~10~ values for zooplankton and reports them in a table using the *gt* package.

-   **4_Q10_models_Clearance:** This script reads in and filters feeding data by clearance rate and estimates temperature dependence. To do this:

    -   First, we filter the data to exclude zooplankton groups that do not have a suitable number of samples or temperature ranges For example, no less than 20 samples and a temperature range of at least 5 degrees Celcius.

    -   Second, we exclude rate processes that are not biologically reasonable (or are extreme outliers), convert temperature to 1/Kelvin and log the carbon mass-specific rates.

    -   Third, we use bootstrap resampling (n = 9,999) and parellel processing with the package *purrr* to fit Generalised Linear Mixed Models with the Template Model Builder (using the *glmmTMB* package) with zooplankton group as a fixed effect and two random effects (primary author and taxa).

    -   Fourth, we extract the slopes from all models for each zoopGroup and estimate median Q10 coefficients following the Arrhenius equation.

    -   Fifth, we estimate the 95% confidence intervals for Q10 estimates using coefficients from the boostrap models output.

    -   Finally, we save the data for a final visualisation and then generate inverse Arrhenius plots to illustrate the relationship between temperature and carbon mass-specific clearance rates for each zooplankton group available.

-   **5_Q10_models_Ingestion:** This script follows the same workflow as script 4, but for ingestion rate.

-   **6_Q10_models_Growth:** This script follows the same workflow as script 4, but for growth rate.

-   **7_Q10_models_Respiration:** This script follows the same workflow as script 4, but for respiration.

-   **8_figure2_Q10_plots.R:** This script combines all Q10 data for each rate and generates dot plots to illustrate median Q10 (95% confidence intervals) for each zooplankton group with available data.

## Output

A folder to save any outputs such as preliminary figures and tables etc.

## 
