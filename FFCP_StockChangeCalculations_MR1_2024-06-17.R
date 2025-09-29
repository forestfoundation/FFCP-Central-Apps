######################################################################################################################

# This script calculates emission reductions following the FFCP methodology accounting framework.
# Using the previous and current inventory measurements, stock change is calculated for 
# enrolled project stands and also for baseline controls.

#      Last edited: BR 04-29-2024, BR 02-23-2024, BR 07-19-2023, DTS, 6-29-2022, adapted from EPB, 2-2-21

######################################################################################################################
rm(list = ls())
gc()
library(rFIA)
library(tidyverse)
library(dplyr)
library(tibble)
library(optmatch)
library(RItools)
library(geosphere)
library(readxl)
library(sf)

# Enter cohort, year, and measurement periods

cohort <- "Cohort2021" # <- Enter "Cohort2020" to run for cohort 2020 and "Cohort2021" to run cohort 2021
year <- ifelse(cohort=="Cohort2020", "2020","2021")
tx <- "T0" # Stays the same
ty <- ifelse(cohort=="Cohort2020","T2","T1")

source("NVEL_2024-02-07.R") # Load NVEL estimation function

# Update output.folder name to match the relevant output folder from the Matching script
output.folder <- paste0("OUTPUT_2024-06-27_submitted/")
ifelse(dir.exists(output.folder), NA, dir.create(output.folder))
ifelse(dir.exists(paste0(output.folder,cohort,"/INDIVIDUAL MATCHES")), 
       NA, 
       dir.create(paste0(output.folder,cohort,"/INDIVIDUAL MATCHES"), recursive=T))


####################################################
# Determine nearest ranger districts to each stand
# to pull region, forest number, and district number
# to use the National Volume Estimator Library (NVEL)
# to estimate volumes for standing dead biomass.
####################################################

# Determine nearest ranger district for each stand
  # Import US Ranger districts spatial data (source: https://data.fs.usda.gov/geodata/edw/datasets.php)
ranger_dist <- st_read("INPUT/S_USA.RangerDistrict/S_USA.RangerDistrict.shp")
ranger_dist <- ranger_dist %>% st_transform(crs="EPSG:4326") %>% # transform datum to match plot points
                          st_make_valid()
                          #select(REGION, FORESTNUMB, DISTRICTNU) # select relevant parameters

# Load FFCP stand boundaries
standBdry <- st_transform(st_zm(st_read(paste0("INPUT/", cohort, "_Shapefiles/FFCP_Stands_", cohort, ".shp"))), crs="EPSG:4326")

# Calculate stand centroids
STD_CENT <- st_centroid(st_make_valid(standBdry))

# Calculate ranger district centroids
DISTRICT_CENT <- st_centroid(ranger_dist)


  # assigns the correct UTM zone to each point to account for differnt UTMs
  points_temp <- STD_CENT %>% 
    mutate(utm_zone = floor((st_coordinates(.)[,1] + 180) / 6) + 1, 
            dist_to_ranger = NA, .after=6)
  
  # Create data frame for matrix of stand numbers and region codes
  ranger_code <- data.frame(matrix(ncol=4, nrow=0))
  colnames(ranger_code) <- c("StandNum", "REGION", "FORESTNUMB", "DISTRICTNU")
  
  # loops through stand centroids to determine the min distance to the nearest ranger district
  
  for(i in 1:nrow(points_temp)){
    
    temp <- points_temp[i,]
    temp <- st_transform(temp, crs = paste0("+init=epsg:326",temp$utm_zone)) #project stand centroids to UTM
    
    DISTRICT_CENT_utm <- st_transform(DISTRICT_CENT, crs = paste0("+init=epsg:326",temp$utm_zone)) #project distric centroids to UTM
    temp_dist <- as.vector(st_distance(temp, DISTRICT_CENT_utm)) #calculate distance between stand and ranger centroids
    proj_district <- cbind(DISTRICT_CENT_utm, temp_dist) #add distances to stands to district centroids
    
    proj_district <- proj_district %>% # Find minimum distance to ranger district
      filter(temp_dist == min(temp_dist)) 
      
    
    proj_district <- as.data.frame(proj_district) # Set to dataframe
    proj_district <- dplyr::select(proj_district, REGION, FORESTNUMB, DISTRICTNU)
    temp <- as.data.frame(temp) # Set to dataframe
    
    temp <- cbind(temp[1], proj_district[,1:3]) # combine stand number and ranger codes

    ranger_code <- rbind(ranger_code, temp) # add codes to table
    
  }
  
  # Delete temporary tables
  rm(proj_district,ranger_dist,temp,points_temp,DISTRICT_CENT, DISTRICT_CENT_utm, STD_CENT, standBdry)


####################################################
#         Import treelist and other data
####################################################
  
# States
states <- read.csv(paste0(output.folder, cohort, "_states.csv"))

# Treelist for treatments, times 1 and 2
treelist.treat.t0 <- read_xlsx(paste0("INPUT/FFCP_InventoryData_", tx,"_", cohort,"_TC.xlsx"), sheet = "Tree_Data_Input")
  treelist.treat.t0 <- treelist.treat.t0[which(!is.na(treelist.treat.t0$TREE_INDEX) & treelist.treat.t0$TREEVAL>0),] # remove any blank rows from the dataframe and any TREEVAL records = 0
  treelist.treat.t0 <- rbind(treelist.treat.t0, treelist.treat.t0[treelist.treat.t0$TREEVAL==2,]) # duplicate tree records in treelist for walkthrough trees
  treelist.treat.t0 <- treelist.treat.t0 %>%
    mutate(across(c(TREE_INDEX, TREE, DIA, AZIMUTH, HD), as.numeric),
           across(c(INVYR, TREEVAL, SPCD, STATUSCD, DECAYCLASS), as.integer))
  summary(treelist.treat.t0)
  
treelist.treat.t1 <- read_xlsx(paste0("INPUT/FFCP_InventoryData_", ty, "_", cohort,"_TC.xlsx"), sheet = "Tree_Data_Input")
  treelist.treat.t1 <- treelist.treat.t1[which(!is.na(treelist.treat.t1$TREE_INDEX) & treelist.treat.t1$TREEVAL>0),] # remove any blank rows from the dataframe and any TREEVAL records = 0
  treelist.treat.t1 <- rbind(treelist.treat.t1, treelist.treat.t1[treelist.treat.t1$TREEVAL==2,]) # duplicate tree records in treelist for walkthrough trees
  treelist.treat.t1 <- treelist.treat.t1 %>%
    mutate(across(c(TREE_INDEX, TREE, DIA, AZIMUTH, HD), as.numeric),
           across(c(INVYR, TREEVAL, SPCD, STATUSCD, DECAYCLASS), as.integer))
  summary(treelist.treat.t1)
  

  
# Plotlist for treatments, times 0 and 1
plotlist.treat.t0 <- read_xlsx(paste0("INPUT/FFCP_InventoryData_", tx, "_", cohort,"_TC.xlsx"), sheet = "Plot_Data_Input")
  plotlist.treat.t0 <- plotlist.treat.t0[which(!is.na(plotlist.treat.t0$UNITCD)),]
  plotlist.treat.t0 <- left_join(plotlist.treat.t0, ranger_code, by=c("UNITCD"="StandNum")) # Add ranger codes to plot data
plotlist.treat.t1 <- read_xlsx(paste0("INPUT/FFCP_InventoryData_", ty, "_", cohort,"_TC.xlsx"), sheet = "Plot_Data_Input")
  plotlist.treat.t1 <- plotlist.treat.t1[which(!is.na(plotlist.treat.t1$UNITCD)),]
  plotlist.treat.t1 <- left_join(plotlist.treat.t1, ranger_code, by=c("UNITCD"="StandNum")) # Add ranger codes to plot data
summary(plotlist.treat.t0)
summary(plotlist.treat.t1)
  
# Reformat tree data to match FIA matched tree data
treelist.treat.t0 <- treelist.treat.t0 %>%
  left_join(plotlist.treat.t0, by = "PLOT")
  
# Edit time 1 data
treelist.treat.t0 <- treelist.treat.t0[treelist.treat.t0$SPCD != 0,] # Remove non-tree species from treelist
treelist.treat.t0$BA <- treelist.treat.t0$DIA^2 * 0.005454154 # individual basal area
treelist.treat.t0$TPA_UNADJ <- ifelse(treelist.treat.t0$DIA < 5, 100,(treelist.treat.t0$BAF / treelist.treat.t0$BA)) # Calculate individual tree TPA
treelist.treat.t0 <- treelist.treat.t0 %>%
  mutate(mt = 1,
         Z = 1) %>%
  dplyr::select(mt,
                PLT_CN = PLOT,
                TREE_INDEX,
                UNITCD,
                STATECD,
                INVYR,
                SPCD,
                DIA,
                TOTAL_HT,
                STATUSCD,
                DECAYCD = DECAYCLASS,
                TPA_UNADJ,
                REGION,
                FORESTNUMB,
                DISTRICTNU,
                Z)

# Edit time 2 data
treelist.treat.t1 <- treelist.treat.t1 %>%
  left_join(plotlist.treat.t1, by = "PLOT")

treelist.treat.t1 <- treelist.treat.t1[treelist.treat.t1$SPCD != 0,] # Remove non-tree species from treelist
treelist.treat.t1$BA <- treelist.treat.t1$DIA^2 * 0.005454154 # individual basal area
treelist.treat.t1$TPA_UNADJ <- ifelse(treelist.treat.t1$DIA < 5, 100,(treelist.treat.t1$BAF / treelist.treat.t1$BA)) # Calculate individual tree TPA
treelist.treat.t1 <- treelist.treat.t1 %>%
  mutate(mt = 2,
         Z = 1) %>%
  dplyr::select(mt,
                PLT_CN = PLOT,
                TREE_INDEX,
                UNITCD,
                STATECD,
                INVYR,
                SPCD,
                DIA,
                TOTAL_HT,
                STATUSCD,
                DECAYCD = DECAYCLASS,
                TPA_UNADJ,
                REGION,
                FORESTNUMB,
                DISTRICTNU,
                Z)

# Create combined treelist
treelist.treat <- rbind(treelist.treat.t0, treelist.treat.t1)
treelist.treat[treelist.treat=='NULL'] <- NA # check for null values in data and replace with NA
treelist.treat$TREE_ID <- seq(1, nrow(treelist.treat)) # Add unique tree ID for each tree record

# Upload treatment and matched plots from matching code output
matches <- read.csv(paste0(output.folder, cohort, "_Treatment and Matched Controls.csv"))


#########################################################################
# Create treelist data from matched FIA plots
#########################################################################

# If matched treelist already exists, load matched treelist.
if(file.exists(paste0(output.folder, cohort, "_Match_Treelist.csv"))) { 
  treelist.match <- read.csv(paste0(output.folder, cohort, "_Match_Treelist.csv")) 
  #Otherwise, create matched treelist using FIA data
} else {
    # Re-Load downloaded FIA data if previously downloaded onto hard drive
    data = readFIA(states = states[states!="DC"], 
                   dir = "FIA_data/",
                   tables = c("PLOT", "TREE"),
                   nCores = 14)
    
    ## Create tree lists of selected matches and treatments
    # Select plot data for stock change calculations
    matches.plot.list <- matches %>%
      filter(Z==0) %>% # Select only FIA plots from data table
      dplyr::select(PLT_CN, treat_PLOTID, PLOTID) %>%
      mutate(PLT_CN = as.numeric(as.character(PLT_CN))) %>% # Set PLT_CN to numeric
      left_join(data$PLOT, by = c("PLT_CN"="CN"), keep=T) %>%  # Join FIA plot data
      dplyr::select(PLT_CN, 
                    PREV_PLT_CN1=PREV_PLT_CN, 
                    treat_PLOTID, 
                    MEASYEAR, 
                    STATECD, 
                    PLOTID=PLOTID) %>% # Add PLT_CN from 1 measurement period ago, most recent inventory year, and state code
      left_join(data$PLOT %>% 
                  dplyr::select(CN,
                                PREV_PLT_CN2=PREV_PLT_CN, 
                                PREV_MEASYEAR1=MEASYEAR), 
                by = c("PREV_PLT_CN1"="CN"),
                keep=T) %>% # Add PLT_CN of 2 measurements periods ago and inventory year of previous measurement
      left_join(data$PLOT %>% dplyr::select(CN, 
                                            PREV_PLT_CN3=PREV_PLT_CN, 
                                            PREV_MEASYEAR2=MEASYEAR),
                by = c("PREV_PLT_CN2"="CN")) %>% # Add PLT_CN of 3 measurements ago and inventory year of 2 measurements ago
      dplyr::select(treat_PLOTID, 
                    PLOTID,
                    PLT_CN, 
                    PREV_PLT_CN1, 
                    PREV_PLT_CN2, 
                    PREV_PLT_CN3, 
                    CURR_MEASYEAR=MEASYEAR, 
                    PREV_MEASYEAR1, 
                    PREV_MEASYEAR2)

    
    # Join FIA tree level data to selected plot matches
    treelist.match <- matches.plot.list %>%
      right_join(data$TREE,by = c("PLT_CN")) %>%
      dplyr::select(PLT_CN,
                    PREV_PLT_CN1,
                    PREV_TRE_CN,
                    CN,
                    PLOTID,
                    treat_PLOTID,
                    STATECD,
                    CURR_MEASYEAR,
                    PREV_MEASYEAR1,
                    SPCD,
                    PREVDIA,
                    DIA,
                    STATUSCD,
                    DECAYCD,
                    TPA_UNADJ,
                    VOLCFNET,
                    VOLCFSND,
                    SPGRPCD,
                    REMVCFAL) %>%
      filter(!is.na(treat_PLOTID)) %>%
      left_join(data$TREE %>% dplyr::select(CN, 
                                            PREVDECAYCD=DECAYCD, 
                                            PREVSTATUSCD=STATUSCD, 
                                            PREVTPA_UNADJ=TPA_UNADJ, 
                                            PREVVOLCFNET=VOLCFNET,
                                            PREVDIA_CHECK=DIA), 
                by =c("PREV_TRE_CN"="CN"))
    
    # Write matched treelist to output folder
    write.csv(treelist.match, paste0(output.folder, cohort, "_Match_Treelist.csv"))
    
    # Remove FIA data
    rm(data)
    gc()
}
############################################################################

# Treelist for matches, only plots with 2 measurements
treelist.match$treat_PLOTID <- as.factor(treelist.match$treat_PLOTID)

treelist.match <- treelist.match %>%
  dplyr::select(PLT_CN,
                CN,
                PREV_TRE_CN,
                treat_PLOTID,
                PLOTID,
                STATECD,
                PREV_MEASYEAR1,
                CURR_MEASYEAR,
                SPCD,
                PREVDIA,
                DIA,
                STATUSCD,
                DECAYCD,
                PREVSTATUSCD,
                PREVDECAYCD,
                TPA_UNADJ,
                PREVTPA_UNADJ,
                VOLCFNET,
                PREVVOLCFNET,
                VOLCFSND,
                SPGRPCD,
                REMVCFAL)


# Adjust plot matches
matches$treat_PLOTID <- as.factor(matches$treat_PLOTID)
matches$PLOTID <- as.factor(matches$PLOTID)
matches$PLT_CN <- as.numeric(as.character(matches$PLT_CN))

# Reference table containing
# FIA REF_SPECIES table
# Density Reduction Factors from Harmon et al. 2011
REF_SPECIES <- read.csv(paste0("INPUT/REF_SPECIES.csv"), na.strings = c("", "#N/A"))

# Storage Factors
sf = read.csv(paste0("INPUT/StorageFactors.csv"))
sf$Region = as.factor(sf$Region)
sf$Type = as.factor(sf$Type)

# FIPS codes and region codes
fips = read.csv(paste0("INPUT/StateFIPS.csv"))
fips$State = as.factor(fips$State)
fips$FRA_Region = as.factor(fips$FRA_Region)
fips$FRA_Region_abbr = as.factor(fips$FRA_Region_abbr)

# Stand weights to correct for post-inventory boundary corrections
STAND_WEIGHT = read_xlsx(paste0("INPUT/STAND_WEIGHT.xlsx"), sheet = "STAND_WEIGHT")
acres <- read_xlsx(paste0("INPUT/STAND_WEIGHT.xlsx"), sheet = "Acres")



####################################################
#         Process FIA Data
####################################################

# Add unique tree record ID
treelist.match$TREE_ID <- seq(1, nrow(treelist.match))

## Create a custom TPA value for use in calculating Harvested Wood Products, because some trees have TPA_UNADJ = NA
# This value will be 6.01, unless the tree has a diameter indicating it is a sapling.
# This should be conservative, as correctly measured trees cannot have TPA < 6.01 on subplots
treelist.match$TPA_HWP = ifelse(treelist.match$PREVDIA < 5, 43560/(pi * 6.8^2 * 4), 43560/(pi * 24^2 * 4)) 
treelist.match$TPA_HWP = ifelse(treelist.match$STATUSCD == 3 & is.na(treelist.match$PREVDIA), 43560/(pi * 24^2 * 4), treelist.match$TPA_HWP) # SHOULD THIS BE AN 'AND' or 'OR' STATEMENT?


####################################################
#   Biomass Calculations for FIA Plots (from Jenkins)
####################################################

## Live Tree Biomass
# Append Jenkins coefficients to each treelist entry
treelist.match$JENKINS_TOTAL_B1 = REF_SPECIES[match(treelist.match$SPCD, REF_SPECIES$SPCD), "JENKINS_TOTAL_B1"]
treelist.match$JENKINS_TOTAL_B2 = REF_SPECIES[match(treelist.match$SPCD, REF_SPECIES$SPCD), "JENKINS_TOTAL_B2"]

# Append Jenkins Root to Shoot Ratio coefficients
treelist.match$JENKINS_ROOT_RATIO_B1 = REF_SPECIES[match(treelist.match$SPCD, REF_SPECIES$SPCD), "JENKINS_ROOT_RATIO_B1"]
treelist.match$JENKINS_ROOT_RATIO_B2 = REF_SPECIES[match(treelist.match$SPCD, REF_SPECIES$SPCD), "JENKINS_ROOT_RATIO_B2"]


# Calculate Aboveground Biomass for all live trees - MOST RECENT MEASUREMENT
treelist.match$LAG_t1 = ifelse(treelist.match$STATUSCD == 1, 
                          exp(treelist.match$JENKINS_TOTAL_B1 + (treelist.match$JENKINS_TOTAL_B2*log(treelist.match$DIA * 2.54))) *  # Jenkins eqn
                            (1/1000) * 0.5 * (44/12) * treelist.match$TPA_UNADJ,       # Sapling, kg to Mg, C ratio, C to CO2e, ExpFac
                          0)
# Calculate Belowground Biomass for all live trees - MOST RECENT MEASUREMENT
treelist.match$LBG_t1 = ifelse(treelist.match$STATUSCD == 1, 
                          treelist.match$LAG_t1 * exp(treelist.match$JENKINS_ROOT_RATIO_B1 + (treelist.match$JENKINS_ROOT_RATIO_B2 / (treelist.match$DIA * 2.54))),
                          0)

# Calculate Aboveground Biomass for all live trees - PREVIOUS MEASUREMENT
treelist.match$LAG_t0 = ifelse(treelist.match$PREVSTATUSCD == 1, 
                            exp(treelist.match$JENKINS_TOTAL_B1 + (treelist.match$JENKINS_TOTAL_B2*log(treelist.match$PREVDIA * 2.54))) *  # Jenkins eqn
                              (1/1000) * 0.5 * (44/12) * treelist.match$PREVTPA_UNADJ,       # Sapling, kg to Mg, C ratio, C to CO2e, ExpFac
                            0)

# Calculate Belowground Biomass for all live trees - PREVIOUS MEASUREMENT
treelist.match$LBG_t0 = ifelse(treelist.match$PREVSTATUSCD == 1, 
                            treelist.match$LAG_t0 * exp(treelist.match$JENKINS_ROOT_RATIO_B1 + (treelist.match$JENKINS_ROOT_RATIO_B2 / (treelist.match$PREVDIA * 2.54))),
                            0)


## Deadwood Biomass
# For any dead trees that lack DECAYCD value, conservatively assume decay class 5
treelist.match$DECAYCD = ifelse(is.na(treelist.match$DECAYCD) & treelist.match$STATUSCD == 2, 5, treelist.match$DECAYCD)

# For previous measurement period
treelist.match$PREVDECAYCD = ifelse(is.na(treelist.match$PREVDECAYCD) & treelist.match$PREVSTATUSCD == 2, 5, treelist.match$PREVDECAYCD)


# Add current measurement DRF using Harmon et al decay class factors
treelist.match$DRF <- ifelse(treelist.match$DECAYCD==1, 0.97,
                             ifelse(treelist.match$DECAYCD==2, 0.97,
                                    ifelse(treelist.match$DECAYCD==3, 0.86,
                                           ifelse(treelist.match$DECAYCD==4, 0.53,
                                                  ifelse(treelist.match$DECAYCD==5, 0.53, NA)))))

# Add previous measurement using Harmon et al decay class factors
treelist.match$PREVDRF <- ifelse(treelist.match$PREVDECAYCD==1, 0.97,
                             ifelse(treelist.match$PREVDECAYCD==2, 0.97,
                                    ifelse(treelist.match$PREVDECAYCD==3, 0.86,
                                           ifelse(treelist.match$PREVDECAYCD==4, 0.53,
                                                  ifelse(treelist.match$PREVDECAYCD==5, 0.53, NA)))))

# Append Specific Gravity from REF_SPECIES to treelist
treelist.match$SG = REF_SPECIES[match(treelist.match$SPCD, REF_SPECIES$SPCD), "WOOD_SPGR_GREENVOL_DRYWT"]
treelist.match$StemWood_Ratio_b1 <- REF_SPECIES[match(treelist.match$SPCD, REF_SPECIES$SPCD), "JENKINS_STEM_WOOD_RATIO_B1"]
treelist.match$StemWood_Ratio_b2 <- REF_SPECIES[match(treelist.match$SPCD, REF_SPECIES$SPCD), "JENKINS_STEM_WOOD_RATIO_B2"]
treelist.match$StemWood_Ratio <- exp(treelist.match$StemWood_Ratio_b1 + (treelist.match$StemWood_Ratio_b2 / (2.54 * treelist.match$DIA)))
treelist.match$PREVStemWood_Ratio <- exp(treelist.match$StemWood_Ratio_b1 + (treelist.match$StemWood_Ratio_b2 / (2.54 * treelist.match$PREVDIA)))


# Calculate deadwood biomass using VOLCFNET and CRM for previous measurement
# Adjust VOLCFNET with TPA_UNADJ
treelist.match$VOLCFNET <- treelist.match$VOLCFNET * treelist.match$TPA_HWP
treelist.match$PREVVOLCFNET <- treelist.match$PREVVOLCFNET * treelist.match$PREVTPA_UNADJ
treelist.match_deadvolt0 <- NVEL(treelist.match, 0)
treelist.match <- treelist.match %>% left_join(treelist.match_deadvolt0 %>% dplyr::select(PLT_CN, treat_PLOTID, PREV_TRE_CN, mTCO2e), 
                                                    by=(c("PLT_CN", "treat_PLOTID","PREV_TRE_CN"))) %>%
  rename(DW_t0=mTCO2e)

# Calculate deadwood biomass of the stem portion of each dead tree using sound cubic foot volumes from FIA
# Many smaller trees or those with decaycd = 5 have NA volume
# Biomass calculation is volume * water lbs/ft3 * specific gravity * Density Reduction Factor * lbs to grams * g to kg * kg to Mg * C ratio * CO2e ratio * ExpFac


#treelist.match$DW_t0 = ifelse(treelist.match$PREVSTATUSCD == 2, 
#                              exp(treelist.match$JENKINS_TOTAL_B1 + (treelist.match$JENKINS_TOTAL_B2*log(treelist.match$PREVDIA * 2.54))) *  # Jenkins eqn
#                                (1/1000) * 0.5 * (44/12) * treelist.match$PREVTPA_UNADJ * treelist.match$PREVDRF * treelist.match$PREVStemWood_Ratio, 
#                              0)

# Biomass calculation of dead wood for current measurement period
treelist.match_deadvolt1 <- NVEL(treelist.match, 1) 
treelist.match <- treelist.match %>% left_join(treelist.match_deadvolt1 %>% dplyr::select(PLT_CN, treat_PLOTID, CN, mTCO2e), 
                                                    by=(c("PLT_CN","treat_PLOTID","CN"))) %>%
  rename(DW_t1=mTCO2e)


#treelist.match$DW_t1 = ifelse(treelist.match$STATUSCD == 2, 
#                              exp(treelist.match$JENKINS_TOTAL_B1 + (treelist.match$JENKINS_TOTAL_B2*log(treelist.match$DIA * 2.54))) *  # Jenkins eqn
#                                (1/1000) * 0.5 * (44/12) * treelist.match$TPA_UNADJ * treelist.match$DRF * treelist.match$StemWood_Ratio, 
#                              0)

# Calculate live tree removals (LT parameter)
# Calculate Aboveground Biomass for live tree removals
treelist.match$LT_LAG_bsl = ifelse(treelist.match$STATUSCD == 3 & treelist.match$PREVSTATUSCD == 1, 
                               exp(treelist.match$JENKINS_TOTAL_B1 + (treelist.match$JENKINS_TOTAL_B2*log(treelist.match$PREVDIA * 2.54))) *  # Jenkins eqn
                                 (1/1000) * 0.5 * (44/12) * treelist.match$TPA_HWP,       # Sapling, kg to Mg, C ratio, C to CO2e, ExpFac
                               0)
# Calculate Belowground Biomass for live tree removals
treelist.match$LT_LBG_bsl = ifelse(treelist.match$STATUSCD == 3 & treelist.match$PREVSTATUSCD == 1, 
                               treelist.match$LT_LAG_bsl * exp(treelist.match$JENKINS_ROOT_RATIO_B1 + (treelist.match$JENKINS_ROOT_RATIO_B2 / (treelist.match$PREVDIA * 2.54))),
                               0)
hrvst_removals_check <- filter(treelist.match, STATUSCD==3)

# Calculate total live tree removal (to be used in leakage equation)
treelist.match$LTrem_bsl <- treelist.match$LT_LAG_bsl + treelist.match$LT_LBG_bsl

write.csv(treelist.match, paste0(output.folder, cohort, "_Treelist_Match_Calcs.csv"))

####################################################
#   Biomass Calculations for Treatment Plots (from Jenkins)
####################################################

## Live Tree Biomass
# Append Jenkins coefficients to each treelist entry
treelist.treat$JENKINS_TOTAL_B1 = REF_SPECIES[match(treelist.treat$SPCD, REF_SPECIES$SPCD), "JENKINS_TOTAL_B1"]
treelist.treat$JENKINS_TOTAL_B2 = REF_SPECIES[match(treelist.treat$SPCD, REF_SPECIES$SPCD), "JENKINS_TOTAL_B2"]

# Append Jenkins Sapling Adjustment Factors to trees < 5" dbh
treelist.treat$JENKINS_SAPLING_ADJUSTMENT = REF_SPECIES[match(treelist.treat$SPCD, REF_SPECIES$SPCD), "JENKINS_SAPLING_ADJUSTMENT"]

# Append Jenkins Root to Shoot Ratio coefficients
treelist.treat$JENKINS_ROOT_RATIO_B1 = REF_SPECIES[match(treelist.treat$SPCD, REF_SPECIES$SPCD), "JENKINS_ROOT_RATIO_B1"]
treelist.treat$JENKINS_ROOT_RATIO_B2 = REF_SPECIES[match(treelist.treat$SPCD, REF_SPECIES$SPCD), "JENKINS_ROOT_RATIO_B2"]


# Calculate Aboveground Biomass for all live trees
treelist.treat$LAG = ifelse(treelist.treat$STATUSCD == 1, 
                            exp(treelist.treat$JENKINS_TOTAL_B1 + (treelist.treat$JENKINS_TOTAL_B2*log(treelist.treat$DIA * 2.54))) *  # Jenkins eqn
                              (1/1000) * 0.5 * (44/12) * treelist.treat$TPA_UNADJ,       # Sapling, kg to Mg, C ratio, C to CO2e, ExpFac
                            0)
# Calculate Belowground Biomass for all live trees
treelist.treat$LBG = ifelse(treelist.treat$STATUSCD == 1, 
                            treelist.treat$LAG * exp(treelist.treat$JENKINS_ROOT_RATIO_B1 + (treelist.treat$JENKINS_ROOT_RATIO_B2 / (treelist.treat$DIA * 2.54))),
                            0)

## Deadwood Biomass
# For any dead trees that lack DECAYCD value, conservatively assume decay class 5
#treelist.treat$DECAYCD <- ifelse(is.na(treelist.treat$DECAYCD) & treelist.treat$STATUSCD == 2, 5, treelist.treat$DECAYCD)

# Add Harmon et al decay class factor to treatment treelist
#treelist.treat$DRF <- ifelse(treelist.treat$DECAYCD==1, 0.97,
#                             ifelse(treelist.treat$DECAYCD==2, 0.97,
#                                    ifelse(treelist.treat$DECAYCD==3, 0.86,
#                                           ifelse(treelist.treat$DECAYCD==4, 0.53,
#                                                  ifelse(treelist.treat$DECAYCD==5, 0.53, NA)))))

# Append Specific Gravity from REF_SPECIES to treelist
#treelist.treat$SG <- REF_SPECIES[match(treelist.treat$SPCD, REF_SPECIES$SPCD), "WOOD_SPGR_GREENVOL_DRYWT"]
#treelist.treat$StemWood_Ratio_b1 <- REF_SPECIES[match(treelist.treat$SPCD, REF_SPECIES$SPCD), "JENKINS_STEM_WOOD_RATIO_B1"]
#treelist.treat$StemWood_Ratio_b2 <- REF_SPECIES[match(treelist.treat$SPCD, REF_SPECIES$SPCD), "JENKINS_STEM_WOOD_RATIO_B2"]
#treelist.treat$StemWood_Ratio <- exp(treelist.treat$StemWood_Ratio_b1 + (treelist.treat$StemWood_Ratio_b2 / (2.54 * treelist.treat$DIA)))


# Calculate deadwood biomass
# Calculate volume using Jenkins equations for with project standing dead wood
#treelist.treat$DW <- ifelse(treelist.treat$STATUSCD == 2,
 #                           exp(treelist.treat$JENKINS_TOTAL_B1 + (treelist.treat$JENKINS_TOTAL_B2*log(treelist.treat$DIA * 2.54))) *  # Jenkins eqn
  #                            (1/1000) * 0.5 * (44/12) * treelist.treat$TPA_UNADJ * treelist.treat$DRF * treelist.treat$StemWood_Ratio, 
   #                            0)

# Calculate standing dead wood volume using NVEL and component ratio method

treelist.treat_deadvol <- NVEL(treelist.treat, 3)
# Add calculated dead wood biomass to treelist
treelist.treat <- treelist.treat %>% left_join(treelist.treat_deadvol %>% dplyr::select(TREE_ID, mTCO2e), 
                                                        by=(c("TREE_ID"))) %>%
  rename(DW=mTCO2e)

####################################################
#         Removals Calculations (HWP)
####################################################

## Append storage factors 

# create variables for region and timber type for each tree in treelist
treelist.match$Region = fips[match(treelist.match$STATECD, fips$FIPS), "FRA_Region"]
treelist.match$Type = ifelse(treelist.match$SPGRPCD < 25, "softwood", "hardwood")
treelist.match$Type = as.factor(treelist.match$Type)

# Append SF, based on region and species type
treelist.match$SF = ifelse(treelist.match$STATUSCD == 3, 
                         ifelse(treelist.match$Type == 'softwood', 
                                ifelse(treelist.match$PREVDIA > 9, sf[match(paste(treelist.match$Region, treelist.match$Type), paste(sf$Region, sf$Type)), "SFsaw"], 
                                       sf[match(paste(treelist.match$Region, treelist.match$Type), paste(sf$Region, sf$Type)), "SFpulp"]),
                                ifelse(treelist.match$PREVDIA > 11, sf[match(paste(treelist.match$Region, treelist.match$Type), paste(sf$Region, sf$Type)), "SFsaw"], 
                                       sf[match(paste(treelist.match$Region, treelist.match$Type), paste(sf$Region, sf$Type)), "SFpulp"])), NA)


# Calculate Harvest Wood Products (Equation 9)
treelist.match$HWP = ifelse(treelist.match$STATUSCD == 3 & treelist.match$PREVSTATUSCD == 1,
                          treelist.match$REMVCFAL * 62.4 * treelist.match$SG * 453.5 * (1/1000) * (1/1000) * 0.5 * (44/12) * treelist.match$TPA_HWP * treelist.match$SF, 0)


###################################################
#         Baseline stock change
###################################################

# Create plot summary data for baseline
summ_bsl1 <- treelist.match %>%
  group_by(treat_PLOTID, PLOTID, PREV_MEASYEAR1, CURR_MEASYEAR) %>%
  summarize(LAG_t0 = sum(LAG_t0, na.rm = TRUE),
            LAG_t1 = sum(LAG_t1, na.rm = TRUE),
            LBG_t0 = sum(LBG_t0, na.rm = TRUE),
            LBG_t1 = sum(LBG_t1, na.rm = TRUE),
            DW_t0 = sum(DW_t0, na.rm = TRUE),
            DW_t1 = sum(DW_t1, na.rm = TRUE),
            HWP = sum(HWP, na.rm = TRUE), # Equation 9
            LTrem = sum(LTrem_bsl, na.rm = TRUE)
            ) %>%
  left_join(matches %>% dplyr::select(PLOTID, treat_PLOTID, weight), by = c("PLOTID", "treat_PLOTID")) # add weight to each plot from matching

# Calculate measurement period
summ_bsl1$X_bsl <- summ_bsl1$CURR_MEASYEAR - summ_bsl1$PREV_MEASYEAR1

# Check growth on plots without harvest
nonharvested_plots<- filter(summ_bsl1, LTrem == 0)

# Check % of plots with harvests
perc_hrvst <- nrow(filter(summ_bsl1, LTrem > 0)) / nrow(summ_bsl1) * 100
perc_hrvst
# Equation 3, 4, 5
# Annualized stock change in baseline carbon pools at measurement time
summ_bsl1$ann_LAG_bsl <- (summ_bsl1$LAG_t1 - summ_bsl1$LAG_t0) * (1/summ_bsl1$X_bsl)
summ_bsl1$ann_LBG_bsl <- (summ_bsl1$LBG_t1 - summ_bsl1$LBG_t0) * (1/summ_bsl1$X_bsl)
summ_bsl1$ann_DW_bsl <- (summ_bsl1$DW_t1 - summ_bsl1$DW_t0) * (1/summ_bsl1$X_bsl)
summ_bsl1$ann_HWP_bsl <- summ_bsl1$HWP * (1/summ_bsl1$X_bsl)
summ_bsl1$ann_LTrem_bsl <- summ_bsl1$LTrem * (1/summ_bsl1$X_bsl)

# Equation 6, 7, 8
# Calculate weighted sum of stock change in each carbon pool in the baseline
summ_bsl1 <- summ_bsl1 %>%
  mutate(ann_stk_chg_bsl = ann_LAG_bsl + ann_LBG_bsl + ann_DW_bsl + ann_HWP_bsl,
         wt_CHG_LAG_bsl = ann_LAG_bsl * weight,
         wt_CHG_LBG_bsl = ann_LBG_bsl * weight,
         wt_CHG_DW_bsl = ann_DW_bsl * weight,
         wt_LTrem_bsl = ann_LTrem_bsl * weight
         ) %>%
  left_join(STAND_WEIGHT, by = "treat_PLOTID")

  
write.csv(summ_bsl1, paste0(output.folder, cohort, "_Annualized_Baseline_Stocking.csv"))

# Equation 11 #
summ_bsl <- summ_bsl1 %>%
  group_by(treat_PLOTID) %>%
  summarize(wt_CHG_LAG_bsl = sum(wt_CHG_LAG_bsl),
            wt_CHG_LBG_bsl = sum(wt_CHG_LBG_bsl),
            wt_CHG_DW_bsl = sum(wt_CHG_DW_bsl),
            HWP_bsl = sum(ann_HWP_bsl * weight),
            wt_LTrem_bsl = sum(wt_LTrem_bsl),
            net_stk_chg_bsl = wt_CHG_LAG_bsl + wt_CHG_LBG_bsl + wt_CHG_DW_bsl + HWP_bsl
            )

write.csv(summ_bsl, paste0(output.folder, cohort, "_Composite_Baseline_Stock_Change.csv"))

ggplot() +
  geom_histogram(data=summ_bsl, binwidth=1, color='gray', aes(x=net_stk_chg_bsl)) +
  labs(title = "Baseline Net Stock Change Distribution") +
  xlab("Net Stock Change (tCO2/ac/yr)")

####################################################
#         Project Stock Change
####################################################


# Create plot summary data for actual
summ_wp1 <- treelist.treat %>%
  group_by(PLT_CN, UNITCD, mt, INVYR) %>%
  summarize(LAG = sum(LAG, na.rm=T),
            LBG = sum(LBG, na.rm=T),
            DW = sum(DW, na.rm=T),
            n_trees = n())

# Add both measurements for each stand to one row
# Measure year is added
summ_wp1 <- summ_wp1 %>%
  pivot_wider(id_cols=c(PLT_CN,UNITCD), 
              names_from=mt, 
              values_from =c(LAG, LBG, DW, n_trees)) %>%
  left_join(plotlist.treat.t0 %>% dplyr::select(PLOT,MEASYEAR), by = c("PLT_CN"="PLOT")) %>%
  rename(MEASYEAR_0 = MEASYEAR) %>%
  inner_join(plotlist.treat.t1 %>% dplyr::select(PLOT,MEASYEAR), by = c("PLT_CN"="PLOT")) %>%
  rename(MEASYEAR_1 = MEASYEAR)
# No tally plots appear as NA, replace with 0
summ_wp1[is.na(summ_wp1)] <- 0


# Diameter increment check
#summ_wp1 <- summ_wp1 %>%
#  filter(n_trees_1 == n_trees_2)

# Equation 13, 14, 15 #
summ_wp1 <- summ_wp1 %>%
  mutate(X_wp = MEASYEAR_1 - MEASYEAR_0,
         ann_chg_LAG = (LAG_2 - LAG_1) * (1/X_wp),
         ann_chg_LBG = (LBG_2 - LBG_1) * (1/X_wp),
         ann_chg_DW = (DW_2 - DW_1) * (1/X_wp),
         ann_chg_stk_wp = ann_chg_LAG + ann_chg_LBG + ann_chg_DW)

write.csv(summ_wp1, paste0(output.folder, cohort, "_Annualized_WithProject_Stocking", Sys.Date(), ".csv"))

# Equation 23 #
summ_wp <- summ_wp1 %>%
  group_by(UNITCD) %>%
  summarize(n_plots = n(),
            m = mean(X_wp),
            net_LAG_wp = sum(ann_chg_LAG, na.rm=T) * (1/n_plots),
            net_LBG_wp = sum(ann_chg_LBG, na.rm=T) * (1/n_plots),
            net_DW_wp = sum(ann_chg_DW, na.rm=T) * (1/n_plots),
            net_stk_chg_wp = sum(ann_chg_stk_wp, na.rm=T) * (1/n_plots))

write.csv(summ_wp, paste0(output.folder, cohort, "_Annualized_WithProject_Stock_Change", Sys.Date(), ".csv"))

# Print histogram of net stock change outcomes
ggplot() +
  geom_histogram(data=summ_wp, binwidth=5, color='gray', aes(x=net_stk_chg_wp)) +
  labs(title = "With-Project Net Stock Change Distribution") +
  xlab("Net Stock Change (tCO2/ac/yr)")

###############################################################
#         ERT and CRT Calculations - Equation 26 and 27
###############################################################

# Summary stock change table
# calculate intermediate parameters for equations 30 and 31
summ_stk_chg <- summ_bsl %>%
  left_join(summ_wp, by=c("treat_PLOTID"="UNITCD")) %>% # Combine baseline and with-project summary tables
  left_join(STAND_WEIGHT, by = "treat_PLOTID") %>% # Append stand weights to summary table
  mutate(net_stk_chg = net_stk_chg_wp - net_stk_chg_bsl, # Calculate net stock change
         PE = 0 + 0, # equation 24, project emissions (no fertilizer or biomass burning)
         BE = 0 + 0, # equation 12, baseline emissions (no fertilizer or biomass burning is conservative)
         wt_net_stk_chg_wp = STAND_WEIGHT * net_stk_chg_wp,
         wt_net_stk_chg_bsl = STAND_WEIGHT * net_stk_chg_bsl,
         min_netstkchg_wp = ifelse(wt_net_stk_chg_wp > 0, 0, wt_net_stk_chg_wp),
         min_netstkchg_bsl = ifelse(wt_net_stk_chg_bsl > 0, 0, wt_net_stk_chg_bsl),
         max_netstkchg_wp = ifelse(wt_net_stk_chg_wp < 0, 0, wt_net_stk_chg_wp),
         max_netstkchg_bsl = ifelse(wt_net_stk_chg_bsl < 0, 0, wt_net_stk_chg_bsl),
         net_ER = STAND_WEIGHT * net_stk_chg) # Correct net stock change by weighting using corrected stand areas

write.csv(summ_stk_chg, paste0(output.folder, cohort, "_Summary_Stock_Change_Table.csv"))

# Equation 30
mean_ERt <- summ_stk_chg %>%
  reframe(mean_ERt = ifelse(sum(wt_net_stk_chg_wp * m) > 0, 1, 0) 
          * (1/sum(STAND_WEIGHT)) * sum(-PE - min_netstkchg_bsl + min_netstkchg_wp) 
          + (1 - ifelse(sum(wt_net_stk_chg_wp * m) > 0, 1, 0)) 
          * (1/sum(STAND_WEIGHT)) 
          * sum(PE 
                - BE 
                - min_netstkchg_bsl 
                + min_netstkchg_wp 
                + max_netstkchg_wp 
                - max_netstkchg_bsl))
mean_ERt <- as.numeric(mean_ERt)

# Non-permanence risk percentage calculated using the VCS NPRT Tool v4.3
NPR = 0.16 # to be updated if NPRT score changes

# Equation 31
mean_CRt <- summ_stk_chg %>%
  reframe(mean_CRt = ifelse(sum(net_stk_chg_wp * m) > 0, 1, 0) 
          * (1/sum(STAND_WEIGHT)) 
          * sum(max_netstkchg_wp - max_netstkchg_bsl))
mean_CRt <- as.numeric(mean_CRt)

# Equation 33 - Buffer Credits
Bu_CRt <- summ_stk_chg %>%
  reframe(Bu_CRt = ifelse(sum(wt_net_stk_chg_wp) > 0, 1, 0)
          * ifelse(sum(wt_net_stk_chg_wp) > 0, 1, 0) 
          * (1/sum(STAND_WEIGHT))
          * sum(max_netstkchg_wp - max_netstkchg_bsl)
          * NPR)
Bu_CRt <- as.numeric(Bu_CRt)

# Equation 34 
Bu_ERt <- summ_stk_chg %>%
  reframe(Bu_CRt = ifelse(sum(wt_net_stk_chg_wp) > 0, 1, 0)
          * ifelse(sum(wt_net_stk_chg_wp) > 0, 1, 0) 
          * (1/sum(STAND_WEIGHT))
          * sum(min_netstkchg_wp - min_netstkchg_bsl)
          * NPR
          + (1 - ifelse(sum(wt_net_stk_chg_wp) > 0, 1, 0))
          * sum(min_netstkchg_wp
          - min_netstkchg_bsl
          + max_netstkchg_wp
          -max_netstkchg_bsl)
          * NPR)
Bu_ERt <- as.numeric(Bu_ERt)

###############################################################
#         Leakage - Equation 25
###############################################################

# Calculate total acreage
acres <- acres %>%
  filter(Cohort==year) %>%
  group_by(Strata)

acres$Acres <- as.numeric(acres$Acres)

# No harvest took place in wp scenario during 1st monitoring period, 
# therefore, LTremoved_wp = 0
# The project activity does not involve permanent reduction in timber supply,
# therefore LFt = 0.10
LK_t_avg <- summ_stk_chg %>%
  mutate(net_LTrem_bsl = (0 - wt_LTrem_bsl),
         net_LTrem = (0 - wt_LTrem_bsl) * 0.10)

LK_t <- LK_t_avg %>%
  group_by(Strata) %>%
  summarize(sum_net_LTrem = sum(net_LTrem),
            n_stands = n()) %>%
  left_join(acres, by = "Strata") %>%
  mutate(LK_t = Acres * sum_net_LTrem * (1/n_stands)) %>%
  mutate(LK_t = ifelse(LK_t > 0, 0, LK_t)) %>%
  dplyr::select(-n_stands, -Acres)

# Equation 28
LK_ERt <- LK_t$LK_t * (mean_ERt/(mean_ERt+mean_CRt))
# Equation 29
LK_CRt <- LK_t$LK_t * (mean_CRt/(mean_ERt+mean_CRt))

write.csv(LK_t_avg, paste0(output.folder,cohort,"_LeakageTable.csv"))

###############################################################
#         Uncertainty calculations - Equation 32
###############################################################

# Calculate variance in treatment plots
summ_var_wp <- summ_stk_chg %>%
        group_by(Strata) %>%
        summarize(sum_net_stk_chg_wp_sq = sum(net_stk_chg_wp^2),
                  sum_net_stk_chg_wp = sum(net_stk_chg_wp),
                  n_stands = n(), # Calculate the number of stands per stratum
                  var_wp = (sum_net_stk_chg_wp_sq 
                            - sum_net_stk_chg_wp^2 
                            * (1/n_stands)) 
                            * (1/(n_stands - 1)) 
                            * (1/n_stands), # Calculate variance per Shiver and Borders eqs 8.25 and 8.26
                  tot_ER_wp = sum(wt_net_stk_chg_wp)/sum(STAND_WEIGHT),
                  tot_ER_bsl = sum(wt_net_stk_chg_bsl)/sum(STAND_WEIGHT),
                  ER_ac_unwtd = mean(net_stk_chg),
                  ERt_ac = sum(net_ER)/sum(STAND_WEIGHT)) %>% # Sum net emission reductions
  mutate(mean_ERt = mean_ERt,
         mean_CRt = mean_CRt)

# Calculate variance in control plots
variance_bsl <- summ_bsl1 %>%
  group_by(Strata) %>%
  distinct(PLOTID, .keep_all = T) %>% # Remove all duplicate plots
  summarize(sum_ann_stk_chg_bsl_sq = sum(ann_stk_chg_bsl^2),
            summ_ann_stk_chg_bsl = sum(ann_stk_chg_bsl),
            n_FIA_plots = nrow(.),
            var_bsl = (sum_ann_stk_chg_bsl_sq - summ_ann_stk_chg_bsl^2 * (1/n_FIA_plots)) * (1/(n_FIA_plots - 1)) * (1/n_FIA_plots)) 
              # Calculate variance on individual baseline plots per Shiver and Borders eqs 8.25 and 8.26

# Add stratum level variance to control plot data
summ_bsl2 <- summ_bsl1 %>%
  left_join(variance_bsl, by = "Strata") %>%
  mutate(wt_var_bsl = weight^2 * var_bsl) # Calculate weighted variance by plot

# Sum weighted variance by stratum
summ_var_bsl <- summ_bsl2 %>%
  group_by(Strata) %>%
  summarize(sum_wt_var_bsl = sum(wt_var_bsl,na.rm=TRUE)) # Sum weighted variance

# Join control and with-project variance estimates
summ_var <- summ_var_wp %>% 
  left_join(summ_var_bsl, by = "Strata") %>%
  left_join(acres, by = "Strata")

# Calculate 95% confidence interval as % of mean ER by stratum
# Calculate total uncertainty (propagate across strata)
summ_var <- summ_var %>%
  mutate(SE = sqrt(var_wp),
         CI = qt(0.975, n_stands - 1) * ((1 / n_stands) * var_wp + (1/(n_stands^2)) * sum_wt_var_bsl)^(1/2),
         UNC = CI *(1/(mean_ERt+mean_CRt)))


################################################################
#         Net emission reductions and removals - Equation 35
################################################################
  
# Output final summary table
summ <- summ_var %>%
  left_join(LK_t, by = "Strata") %>%
  mutate(ERt = Acres * mean_ERt,
         CRt = Acres * mean_CRt,
         Bu_CRt = Acres * Bu_CRt,
         Bu_ERt = Acres * Bu_ERt,
         LK_ERt = LK_ERt,
         LK_CRt = LK_CRt,
         net_ERt = (ERt  + LK_ERt),
         net_CRt = (CRt + LK_CRt),
         UNC = UNC * 100,
         VCU_CRt = CRt - Bu_CRt,
         VCU_ERt = ERt - Bu_ERt) %>%
  dplyr::select(Strata,
                "# of Stands" = n_stands,
                "Area (At) (Acres)" = Acres,
                "Annualized Net Stock Change WP (tCO2e/ac)"=tot_ER_wp,
                "Annualized Net Stock Change BSL (tCO2e/ac)"=tot_ER_bsl,
                "Annualized Mean Emission Reductions (ERt) (tCO2e/ac)" = mean_ERt,
                "Annualized Mean Emission Removals (CRt) (tCO2e/ac)" = mean_CRt,
                "Gross Emission Reductions (ERt) (tCO2e)"=ERt,
                "Gross Emission Removals (CRt) (tCO2e)"=CRt,
                "Variance of WP Stock Change"= var_wp,
                "Weighted Variance of BSL Stock Change"=sum_wt_var_bsl,
                "Standard Error"=CI,
                "Leakage (LKt) (tCO2e)"=LK_t,
                "Leakge (CRt) (tCO2e)" = LK_CRt,
                "Leakage (ERt) (tCO2e)" = LK_ERt,
                "Uncertainty (UNCt) (%)"=UNC,
                "Buffer Units (CRt)" = Bu_CRt,
                "Buffer Units (ERt)" = Bu_ERt,
                "Verified Carbon Units (CRt)" = VCU_CRt,
                "Verified Carbon Units (ERt)" = VCU_ERt,
                "Net Emission Reductions (ERt) (tCO2e)"=net_ERt)


write.csv(summ, paste0(output.folder, cohort, "_ER_Summary.csv"))

###############################################################
# END CODE HERE
##############################################################

