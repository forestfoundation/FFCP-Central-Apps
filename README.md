# FFCP Central Appalachians IFM

This repository contains the matching and calculations code along with data inputs for the VM0045 Family Forest Carbon Program - Central Appalachians IFM project.

## Table of Contents

- [Description](#description)  
- [Project Structure](#project-structure)  
- [Getting Started](#getting-started)  
- [Data](#data)  
- [Usage](#usage)  
- [License](#license)  
- [Acknowledgements](#acknowledgements)  

---

## Description

The code contained in this repository has been validated by Aster Global Inc and Verra under VM0045 v1.2. The repository is referenced in the PD and MR of the Family Forest Carbon Program for the purposes of project validation and verification. The R code and data inputs included herein follow the matching procedure provided in VM0045 Appendix A as well as the calculations and corresponding equations from the methodology. Auxilliary codes are called into the main matching and calculations scripts. 
---

## Project Structure

FFCP-Central-Apps/
├── code/ # Analysis scripts
├── data/
│ ├── input/ # Raw input data
│ └── output/ # Processed outputs, results, figures
├── .gitignore # Files/folders ignored by Git
├── README.md # Project documentation
└── LICENSE # GNU GPL 3.0 license

---

## Data

- Input Data: Raw Datasets
- Output Data: Intermediate and final outputs from the matching and calculations
- Notes:  The project area and plot point shapefiles are not provided with the input data due to confidentiality of landowner property locations. They can be requested separately. The scripts for the first monitoring period are uploaded to the repository for validation. The FIA program has updated the database structure and it is no longer compatible with this version of the code. It should be compatible with the FIA data set downloaded in 2023 for this analysis. Future versions of the script will be updated to accommodate changes to the FIA database.

## Usage

Run the matching script first:
The script requires access to the following:
 
 1. Inventory tree list data (user defined)
 2. Inventory plot list data (user defined)
 3. Stand boundary shapefile (user defined)
 4. Plot point shapefile (user defined)
 5. REF_SPECIES table (from FIA)
 6. SI reference table (from FIA)
 7. SI equation codes (from FIA)
 8. SI location group codes (from FIA)
 9. Ecoregion by state lookup table (customized)
 10. Stocking coefficients (from FIA)

 This script will:
     1) Take initial inventory data from t0 inventories for Cohorts 1 and 2 as inputs.
     2) apply the FFCP matching protocol to those stands, and 
     3) select a donor pool of control (FIA) plots to estimate the baseline carbon stock change
     4) the donor pool plot list will then be input to the FFCP Carbon Accounting Code

---

## License

This project is licensed under the GNU GPL 3.0. 

---

## Acknowledgements

Data Source: USDA Forest Inventory and Analysis (FIA) Program, all participating landowners, American Forest Foundation (AFF), The Nature Conservancy (TNC), Verra, TerraCarbon LLC
Tools and libraries used: rFIA, tidyverse, optmatch, tigris