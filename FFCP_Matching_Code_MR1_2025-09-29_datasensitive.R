#***********************************************************************
############## Introduction ###################

# This script runs the matching protocol from Appendix A of the validated Verra Methodology
# for Improved Forest Management Using Dynamic Matched Baselines from National Forest Inventories (VM00045).
#
# The code is composed of 3 main sections. First, it calculates the required covariates
# from data inputs specific to the project area. Second, it prepares the
# project data and national inventory data to ready it for the matching algorithm.
# Third, it runs the matching algorithm using the 'optmatch' package and outputs
# a set of control plots to be used for tracking stock change in the baseline.
#
# The matched plots data and corresponding tree data are archived and used to run 
# stock change calculations in a next section of code, currently saved in "FFCP_StockChangeCalculations_*".
#
# This code was built for and tested with the FFCP Central Appalachians Carbon Project.
#
# **Important** - data tables need to be updated in line with the region in which data is collected.
#
# The script requires access to the following:
# 
# 1. Inventory tree list data (user defined)
# 2. Inventory plot list data (user defined)
# 3. Stand boundary shapefile (user defined)
# 4. Plot point shapefile (user defined)
# 5. REF_SPECIES table (from FIA)
# 6. SI reference table (from FIA)
# 7. SI equation codes (from FIA)
# 8. SI location group codes (from FIA)
# 9. Ecoregion by state lookup table (customized)
# 10. Stocking coefficients (from FIA)
#
# This script will:
#     1) Take initial inventory data from t0 inventories for Cohorts 1 and 2 as inputs.
#     2) apply the FFCP matching protocol to those stands, and 
#     3) select a donor pool of control (FIA) plots to estimate the baseline carbon stock change
#     4) the donor pool plot list will then be input to the FFCP Carbon Accounting Code

#      Last edited: BDSR 2025-09-29
#
# Note: This version of the code has been edited for data privacy reasons and will 
# run using intermediate datasets that do not include personally identifiable information.

#***********************************************************************


#**********************************************
## Libraries, User Entered Values, Calling Data        ----
#**********************************************

#### Libraries ####
  rm(list = setdiff(ls(), "data")) #use to remove all objects except FIA data
  gc() # reset memory
  library(rFIA) # automated access to FIA servers
  library(tidyverse) # data organization
  library(optmatch) #matching package
  library(RItools)
  library(readxl) # sub-package of tidyverse
  library(geosphere) # geospatial package
  library(raster) # geospatial package for rasters
  library(terra) # used for visualizing maps
  library(elevatr) #accessing elevation data
  library(sf) # spatial package
  library(lwgeom)
  library(dplyr)
  library(tigris) # remotely loading roads data from US census
  library(rstac)
  library(gdalcubes)


#### User Specified Values  ####
  cohort <- "Cohort2020" # Specifies the cohort which will be used for labeling input and output files and selecting the appropriate field data
  k <- 10 # The number of control plots to match to each treatment plot
  ownership <- "Private" # enter "Private" or "Public" or "All"
  std_origin <- "Natural" # enter "Natural" or "Plantation"

  # Create new output folder
  output.folder <- paste0("OUTPUT_", Sys.Date(), "/")
  ifelse(dir.exists(output.folder), NA, dir.create(output.folder))
  ifelse(dir.exists(paste0(output.folder,cohort,"/INDIVIDUAL MATCHES")), 
         NA, 
         dir.create(paste0(output.folder,cohort,"/INDIVIDUAL MATCHES"), recursive=T))

#### Import Data ####
  # 
  # ### Measured field data ###
  # 
  # # Download treelist and plotlist inventory field data from treatment plots
  # # A stand level table will be created from the tree and plot level data
  # treelist <- read_xlsx(paste0("INPUT/FFCP_InventoryData_T0_", cohort,"_TC.xlsx"), sheet = "Tree_Data_Input") #tree data
  #   treelist <- treelist[which(!is.na(treelist$TREE_INDEX) & treelist$TREEVAL>0),] # remove any blank rows from the dataframe and any TREEVAL records = 0
  #   treelist <- rbind(treelist, treelist[treelist$TREEVAL==2,]) # duplicate tree records in treelist for walkthrough trees
  # plotlist <- read_xlsx(paste0("INPUT/FFCP_InventoryData_T0_", cohort,"_TC.xlsx"), sheet = "Plot_Data_Input") #plot data
  #   plotlist <- plotlist[which(!is.na(plotlist$UNITCD)),]
  # 
  # # Create vector of all states included in cohort
  # STATECD <- unique(plotlist$STATECD)
  # 
  # # Pull list of unique species codes of site index trees from inventory
  # speciescomp <- as.data.frame(unique(treelist$SPCD[treelist$SITREE>0]))
  # colnames(speciescomp) <- "SPCD"
  # 
  # # Quick check treelist
  # summary(treelist)
  # 
  # # Quick check plotlist
  # summary(plotlist)
  # 
  # 
  ### Reference data sets ###
  
  # Download Reference Species Table from FIADB (access here: https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html)
  # **Check for periodic updates to REF_SPECIES table**
  REF_SPECIES <- read.csv(paste0("INPUT/REF_SPECIES.csv"), na.strings = c("", "#N/A"))
  
  # Download Stocking coefficient table for calculating FORTYPGRP
  stk_coef <- read.csv(paste0("INPUT/STOCKING_COEFFICIENTS.csv"), na.strings = c("", "#N/A"))
  
  # Download Chapman-Richards coefficients to REF_SPECIES table (coefficients derived specifically for site trees present in inventory)
  # The coefficients were consolidated based on the site trees in the inventory.
  # New cohorts will require re-assessment of site tree species and the appropriate 
  # stand index equations and coefficients will need to be collated.
  SI_LOCGRP <- read.csv(paste0("INPUT/SI_LOCGRP.csv"))
  # Select the SI location group code based on states in which plots were sampled
  # Location code is used to determine the appropriate SI coefficients to use.
  # The Northeastern states all use region code S24 and apply the Westfall et al coefficients
  # Where coefficients for a species is not included in Westfall et al, the Scott and Voorhis 
  # coefficients will be applied.
  LOCATIONCD <- SI_LOCGRP$LOC_GRP[which(SI_LOCGRP$STATECD %in% plotlist$STATECD)]
  # Reference table with species specific coefficients from Westfall et al (only appropriate for NE US)
  SI_REFTABLE_WESTFALL <- read.csv(paste0("INPUT/EQ_COEFF_WESTFALL.csv"))
  # Reference table with species specific coefficients from Scott and Voorhis
  SI_REFTABLE_SV <- read.csv(paste0("INPUT/SI_REFTABLE_SV.csv"))
  # Reference table of CONFIG_ID (based on location), SPECIES_NUM (SPCD), and SI_EQ (species specific SI equation code)
  NIMS_SPECIES_CONFIG <- read.csv(paste0("INPUT/NIMS_SPECIES_CONFIG.csv"))
  NIMS_SPECIES_CONFIG <- filter(NIMS_SPECIES_CONFIG, CONFIG_ID%in%LOCATIONCD)
  # Evaluate species coeffs available in Westfall and replace with S&V if unavailable
  SI_EQCODE <- speciescomp %>%
    left_join(NIMS_SPECIES_CONFIG, by = c("SPCD"="SPECIES_NUM"))
  # Reference table of site class bins based on location and species
  NIMS_REF_RANGE_SITE_CLASS <- read.csv(paste0("INPUT/NIMS_REF_RANGE_SITE_CLASS.csv"))
  NIMS_REF_RANGE_SITE_CLASS <- filter(NIMS_REF_RANGE_SITE_CLASS, LOC_GRP%in%LOCATIONCD)
    # Subset table to region
  # Join SI tables to create SI lookup table for Westfall coeffs and for Scott and Voorhis
  SI_LOOKUP_WESTFALL <- SI_REFTABLE_WESTFALL %>%
    inner_join(NIMS_REF_RANGE_SITE_CLASS, by = "SPCD")
  
  SI_LOOKUP_SV <- SI_REFTABLE_SV %>%
    inner_join(NIMS_REF_RANGE_SITE_CLASS, by = "SPCD")

  # Download list of states with intersection eco-subsections for cross-referencing states to gather FIA data.
  # Eco-subsection will be overlaid on a US state map to determine which states to select FIA data from.
  ecolist <- read.csv(paste0("INPUT/EcoregionByState_Table_2021_05_03.csv"), 
                    header=TRUE, 
                    sep = ',')

  
  # ### Spatial Data sets ###
  # 
  # # Download plot points and stand boundary shapefiles of the cohort being matched
  # plotPts <- st_read(paste0("INPUT/", cohort, "_Shapefiles/FFCP_Plots_", cohort, ".shp")) # plot point shapefiles
  # standBdry <- st_transform(st_zm(st_read(paste0("INPUT/", cohort, "_Shapefiles/FFCP_Stands_", cohort, ".shp"))), crs="EPSG:4326")
  # standBdry <- standBdry[standBdry$Selected == "Yes",] # subset to only properties selected for sampling, Note: this only works if
  #                                                     # a selection attribute is added to the data
  # # create temporary table with unit codes and state codes
  # plotlist_state <- plotlist %>%
  #   group_by(UNITCD) %>%
  #   summarize(STATECD = min(STATECD))
  # 
  # # add state code to stand boundary shapefile
  # standBdry <- standBdry %>%
  #   left_join(plotlist_state %>% dplyr::select(UNITCD, STATECD), by = c("StandNum"="UNITCD"))
  # 
  # # Calculate stand centroids
  # sf_use_s2(FALSE)
  # STD_CENT <- st_centroid(standBdry)
  # STD_CENT <- filter(STD_CENT, !is.na(STATECD))
  # 
  # # US states for map viewing
  # us_states <- states(cb=TRUE)
  # us_states <- st_transform(us_states, st_crs(standBdry))
  # 
  # # Display properties
  # ggplot() +
  #   geom_sf(data=us_states[us_states$STATEFP %in% STATECD,]) +
  #   geom_sf(data=standBdry) +
  #   geom_sf(data=STD_CENT)

  # Download Digital Elevation Model data for measuring elevation of each plot
  # Uses the 'elevatr' package to download elevation data from the AWS server
  # Projection: "EPSG:4326" = WGS84
  # zoom = 13, resolution ~ 13.5m at 45 deg lat
  # It is important to re-download the DEM if the property boundary shapefiles are updated, otherwise, the saved DEM can be used.
  #for(i in nrow(1:standBdry)) {
   #DEM <- get_elev_raster(standBdry, # Download DEM tiles for all stands
    #                   z=12, 
     #                  prj="EPSG:4326",
      #                 clip="bbox",
       #                src="aws")
#    
#   #writeRaster(DEM, paste0("INPUT/", cohort, "_stand_DEM_bbox"), overwrite=T) # Save downloaded DEM raster to avoid redownloading and reprocessing with each run
#   DEM <- raster(paste0("INPUT/", cohort, "_stand_DEM_bbox.grd")) # Read already downloaded and processed DEM into R script
# #}
#   
#   # Download ecosubsection shapefile of the US
#   # This shapefile will be matched against US states to determine which states to download FIA data from
#   ecosub <- st_read(paste0("INPUT/S_USA.EcomapSubsections/S_USA.EcomapSubsections.shp"))
#   ecosub <- st_transform(ecosub, st_crs(standBdry)) #transform projection from NAD83 to WGS84
#   
#   # Identify states in the ecological subsection
#   sf_use_s2(FALSE) # turns off spherical topology to allow st_join to run
#   ecosub.plot <- st_join(standBdry, ecosub, largest=TRUE) # spatial join of stand boundaries to ecosubsections
#   ecosub.plot <- as.data.frame(ecosub.plot) # creates vector of unique map units
#   plotlist$ECOSUBCD <- ecosub.plot[match(plotlist$UNITCD, ecosub.plot$StandNum), "MAP_UNIT_S"] # add ECOSUBCD to plotlist
#   plotlist$ECOSEC <- substr(plotlist$ECOSUBCD, 1, nchar(as.character(plotlist$ECOSUBCD)) - 1) # create eco section code
#   plotlist$ECOPROV <- substr(plotlist$ECOSUBCD, 1, nchar(as.character(plotlist$ECOSUBCD)) - 2) # create eco province code
#   states <- as.data.frame(unique(ecolist[which(ecolist$ecological.province %in% plotlist$ECOPROV), "State"])) # matches selected ecoprovinces to states to download correct FIA data
#   colnames(states) <- "states"
#   write_csv(states, paste0(output.folder, cohort, "_states.csv"))
states <- read.csv(paste0("INPUT/",cohort,"_states.csv"))
  # 
  # # Calculate slope at the stand level
  # # Uses 'raster' package to analyze DEMs
  # # Slope is measured in degrees for the DEM and converted to percent later
  # # Method analyzes slope based on the 8 surrounding raster pixels ('Queen method')
  # slopeMap <- terrain(DEM,
  #                     opt='slope',
  #                     unit='degrees',
  #                     neighbors=8)

  gc()
  ### Download control plot data ###
  
  # Download FIA data from DataMart, if not yet downloaded, to hard drive
  # Downloads the plot, condition, and tree tables for selected states
  # Runs parallel processing to speed download times, may need to be adjusted depending on machine
  # Load downloaded FIA data if previously downloaded onto hard drive
  data = readFIA(states = states[states!="DC"], 
                dir = "FIA_Data/",
                tables = c("PLOT", "COND", "TREE", "SITETREE"),
                inMemory = TRUE,
                nCores = 14)
  #data = data.full

  ## Define functions
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  
  source("fortypgrp_calc.R")
  
  # Remove tables that are no longer needed
  rm(ecolist, SI_LOCGRP, ecosub, NIMS_REF_RANGE_SITE_CLASS)
  gc()
  
# #**********************************************
# ## Produce covariates for treatment data     ----
# #**********************************************
# 
#   ### Spatial data analysis ###
#   # Input slope from plot list
#   #slope <- plotlist %>%
#    # dplyr::select(UNITCD, Slope) %>%
#     #group_by(UNITCD) %>%
#     #summarize(Slope = mean(Slope, na.rm = T)*100)
#   
#   # Run zonal statistics (calculate avg slope)
#   # Uses 'raster' package to calculate mean slope for each stand
#   slope_zone <- raster::extract(slopeMap, 
#                                plotPts,
#                                method="simple",
#                                fun=mean,
#                                exact=T,
#                                df=T,
#                                sp=T,
#                                na.rm=T)
#   slope <- slope_zone@data # subset shapefile data to data.frame
#   slope$slopePer <- tan(slope$slope*(pi/180))*100 # convert slope from degrees to percent slope
#   slope <- slope %>%
#     group_by(StandNum) %>%
#     summarize(slopePer=mean(slopePer, na.rm=T))
# 
#   # create dataframe of centroid lat/lon
#   STD_LAT_LONG <- STD_CENT %>%
#     mutate(LON = unlist(map(STD_CENT$geometry,1)),
#            LAT = unlist(map(STD_CENT$geometry,2)))
# 
#   points <- data.frame()
# #### Road Distance ####
# # Loops through each state that a property is in to estimate distance to nearest road from property centroids
# for(m in 1:length(STATECD)) {
#   
#   j = STATECD[m]
#   
#   # assigns the correct UTM zone to each point to account for differnt UTMs
#   points_temp <- STD_CENT[STD_CENT$STATECD==j,] %>% 
#     mutate(utm_zone = floor((st_coordinates(.)[,1] + 180) / 6) + 1, 
#            RDDIST = NA, .after=6)
# 
# 
#   # grabs list of county FIPS codes in the project states
#   proj_counties <- tigris::counties(state = j)
#   
#   # loops through stand centroids to determine the min distance to the nearest road
#   
#   for(i in 1:nrow(points_temp)){
#     
#     temp <- points_temp[i,]
#     temp <- st_transform(temp, crs = paste0("+init=epsg:326",temp$utm_zone))
#     
#     temp_buff <- st_buffer(temp, dist = 10000)
#     proj_counties <- st_transform(proj_counties, crs = paste0("+init=epsg:326",temp$utm_zone))
#     
#     buff_county_int <- st_intersection(temp_buff, proj_counties)$GEOID
#     
#     roads <- tigris::roads(county = str_sub(buff_county_int, start = 3), state = j)
#     roads <- st_transform(roads, crs = paste0("+init=epsg:326",temp$utm_zone))
#     
#     
#     points_temp$RDDIST[i] <- min(st_distance(temp, roads))
#     
#   }
# 
#   points <- rbind(points, points_temp)
# }
#   
#   
#    #Assign code for distance to road by distance ranges in feet
#    #Source (FIADB)
#   RDDISTCD <- points %>%
#     mutate(RDDIST = RDDIST * 3.28084,
#       RDDISTCD = case_when(
#       RDDIST <= 100 ~ 1,
#       RDDIST > 100 & RDDIST <= 300 ~ 2,
#       RDDIST > 300 & RDDIST <= 500 ~ 3,
#       RDDIST > 500 & RDDIST <= 1000 ~ 4,
#       RDDIST > 1000 & RDDIST <= 2640 ~ 5,
#       RDDIST > 2640 & RDDIST <= 5280 ~ 6,
#       RDDIST > 5280 & RDDIST <= 15840 ~ 7,
#       RDDIST > 15840 & RDDIST <= 26400 ~ 8,
#       RDDIST > 26400 ~ 9
#     ), .after=10)
#   
#   
#   
#   # Create stand buffer to remove FIA plots within 1.6km of a stand
#   # prep stand boundary buffers for evaluating distance and location of FIA plots relative to stands
#   
#   standBdry_buff <- data.frame() # create empty dataframe for stand boundary buffers
#   
#   # for loop that creates 1.6 km buffers around each stand boundary
#   # will be used in the matching algorithm to evaluate if an FIA plot falls within the project area
#   for(i in 1:nrow(RDDISTCD)){
#     
#     temp <- RDDISTCD[i,]
#     temp <- st_transform(temp, crs = paste0("+init=epsg:326",temp$utm_zone)) # project in UTM so we can buffer in meters
#     
#     buff <- st_buffer(temp, dist = 1600) # add 1600 m buffer to stand boundary
#     buff <- st_transform(buff, crs = "EPSG:4269") # transform to NAD83 to be compatible FIA plot dist analysis
#     
#     standBdry_buff <- rbind(standBdry_buff, buff)
#     
#   }
#   
# #### Add parameters to treelist ####
#   plotlist <- plotlist %>%
#       group_by(UNITCD) %>%
#       mutate(STAND_N = n())
#     
#   # Join plotlist table and REF_SPECIES table to treelist table
# 
#   treelist <- treelist %>%
#     left_join(by = "PLOT", plotlist) %>%
#     left_join(by = "SPCD", REF_SPECIES) %>%
#     left_join(by = "STOCKING_SPGRPCD", stk_coef) %>%
#     mutate(b1 = ifelse(is.na(b1), 0, b1),
#            b0 = ifelse(is.na(b0), 0, b0))
# 
# treelist_calcs <- function(treelist) {
#   treelist <- treelist[treelist$SPCD != 0,] # Remove non-tree species from treelist
#   treelist$BA <- treelist$DIA^2 * 0.005454154 # individual basal area
#   treelist$TPA_UNADJ <- ifelse(treelist$DIA < 5, 100,(treelist$BAF / treelist$BA)) # Calculate individual tree TPA
#   treelist$BA_AC <- treelist$BA * treelist$TPA_UNADJ # Calculate per acre basal area (same as BAF for VRP)               
#   treelist$RD <- 2.47 * (0.00015 + (0.00218 * treelist$WOOD_SPGR_GREENVOL_DRYWT)) * ((treelist$DIA/10)^1.6) * (1/treelist$STAND_N) # Individual Ducey-Knapp relative density - English units 
#   treelist$RD_AC <- treelist$RD * treelist$TPA_UNADJ # per acre Ducey-Knapp relative density
#   #generate stocking value. Set stockability proportion to 1. Arner et al Step 1
#   treelist$STOCKING_ALT <- ifelse(treelist$STATUSCD == 1 & treelist$DIA >= 5, 
#                                   treelist$b0 * treelist$DIA ^ (treelist$b1) * treelist$TPA_UNADJ/1,
#                                   0) 
#   
#   return(treelist)
# }
# 
# treelist <- treelist_calcs(treelist)
# 
# ###### FOREST TYPE GROUP ####### 
# 
# # Create table with total stocking values of each species group
#   # This process is based on the FIA process for assigning forest type group codes by species composition
#   # ADD SOURCE LOCATION FOR FIA STOCKING CALCULATIONS
# 
# fortypgrp_list <- fortypcode(treelist)
# 
# # Output a copy of the forest type table
# write.csv(fortypgrp_list, paste0(output.folder, "forest_type_groupcode_table", cohort, ".csv"))
# 
# #### Create stand list table ####
# 
#   # Site index estimation
# 
#     # Pull species and age of site tree on each plot
#     #Add method for SI calculation to treelist
#     treelist <- mutate(treelist, 
#                        SI_METHOD = case_when(SITREE>0 & SPCD %in% SI_LOOKUP_WESTFALL$SPCD ~ "WESTFALL", 
#                                              SITREE>0 ~ "SV",
#                                              ))
#     treelist$SITREE <- treelist$SITREE + 5 # Add 5 years to tree core for all eastern species
#     treelist$TOTAL_HT <- as.numeric(treelist$TOTAL_HT)
#     SiteIndex_WESTFALL <- treelist %>%
#       filter(SITREE > 0 & SI_METHOD=="WESTFALL") %>%
#       dplyr::select(PLOT, SPCD, TOTAL_HT, SITREE, SI_METHOD) %>%
#       left_join(SI_LOOKUP_WESTFALL, by = "SPCD")
#     # Prep site index using Scott and Voorhis method for species not included in Westfall
#     SiteIndex_SV <- treelist %>%
#       filter(SITREE > 0 & SI_METHOD=="SV") %>%
#       left_join(SI_EQCODE, by = ("SPCD")) %>%
#       dplyr::select(PLOT, SPCD, SI_EQ, TOTAL_HT, SITREE, SI_METHOD) %>%
#       left_join(SI_LOOKUP_SV, by = "SI_EQ")
# 
#     # Run SI equation for each site tree with parameters: ht, coefficients, age
#     # Note: FIA changed the SI equation from Scott and Voorhis to Westfall et al in 2022
#     # Site class is calculated using the same method for donor and treatment plots
#     # Equation source: Westfall et al. 2017. Site Index Models for Tree Species in the Northeastern United States. Forest Science, Vol 63:3. 283-290.
#     # Ht at age 50
# 
#     SiteIndex_WESTFALL <- SiteIndex_WESTFALL %>%
#                                 mutate(RES = TOTAL_HT - (4.5 + b0 * (1 - exp(-b1 * SITREE))^b2),
#                                        Z = (1 - exp(-b1 * SITREE))^b2 + (1 - exp(-b1 * SITREE))^b2 * b3 * log((1 - exp(-b1 * SITREE))) * b0,
#                                        R = exp(L * SITREE),
#                                        X = D * Z * ((Z * D * Z + R)^(-1)) * RES,
#                                        SI = if_else(TOTAL_HT < MAX_HT &
#                                                     TOTAL_HT > MIN_HT &
#                                                     SITREE > MIN_AGE &
#                                                     SITREE < MAX_AGE, 
#                                             4.5 + (b0 + X) * (1 - exp(-b1 * 50))^(b2 + b3 * X), NA))
# 
#     # For species that do not have a Westfall coeffient, Scott and Voorhis is used
#     SiteIndex_SV <- SiteIndex_SV %>%
#                             mutate(SI = if_else(SITREE > MIN_AGE &
#                                                 SITREE < MAX_AGE, 
#                                       TOTAL_HT / (b1 * (1.0 - exp(-b2 * SITREE))^b3), NA))
#     
#     
#     # Bin site index into site class code based on species and min SI values (source: NIMS_REF_RANGE_SITECLASS, LOC_GRP = S24)
#     SiteIndex_WESTFALL <- SiteIndex_WESTFALL %>%
#       mutate(SITECLCD = case_when(
#         SI < SITECL5_LOW | is.na(SITECL5_LOW) & SI < SITECL6_LOW ~ 6,
#         SI < SITECL4_LOW & SI >= SITECL5_LOW | is.na(SITECL4_LOW) & SI >= SITECL5_LOW ~ 5,
#         SI < SITECL3_LOW & SI >= SITECL4_LOW | is.na(SITECL3_LOW) & SI >= SITECL4_LOW ~ 4,
#         SI < SITECL2_LOW & SI >= SITECL3_LOW | is.na(SITECL2_LOW) & SI >= SITECL3_LOW ~ 3,
#         SI < SITECL1_LOW & SI >= SITECL2_LOW | is.na(SITECL1_LOW) & SI >= SITECL2_LOW ~ 2,
#         SI >= SITECL1_LOW & is.numeric(SITECL1_LOW) ~ 1))
#     
#     SiteIndex_SV <- SiteIndex_SV %>%
#       mutate(SITECLCD = case_when(
#         SI < SITECL5_LOW | is.na(SITECL5_LOW) & SI < SITECL6_LOW ~ 6,
#         SI < SITECL4_LOW & SI >= SITECL5_LOW | is.na(SITECL4_LOW) & SI >= SITECL5_LOW ~ 5,
#         SI < SITECL3_LOW & SI >= SITECL4_LOW | is.na(SITECL3_LOW) & SI >= SITECL4_LOW ~ 4,
#         SI < SITECL2_LOW & SI >= SITECL3_LOW | is.na(SITECL2_LOW) & SI >= SITECL3_LOW ~ 3,
#         SI < SITECL1_LOW & SI >= SITECL2_LOW | is.na(SITECL1_LOW) & SI >= SITECL2_LOW ~ 2,
#         SI >= SITECL1_LOW & is.numeric(SITECL1_LOW) ~ 1))
# 
#     # Combine the two site index tables into one
#     SiteIndex <- SiteIndex_WESTFALL %>%
#       dplyr::select(PLOT, SI, SITECLCD) %>%
#       bind_rows(SiteIndex_SV %>% dplyr::select(PLOT, SI, SITECLCD))
#     
#     # Add Site Index and Site Class Code to plotlist
#     SiteIndex <- SiteIndex %>% dplyr::select(c("PLOT","SI","SITECLCD"))
#     plotlist <- plotlist %>%
#       dplyr::select(-SITECLCD) %>%
#       left_join(SiteIndex, by="PLOT")
# 
# #### Summarize and organize stand list data ####    
#     
#   # Summarize tree list data to create stand level covariates
#   # Included: 
#     # Quadratic Mean Diameter (QMD),
#     # Relative density of commercial species (RD.comm)
#     # Relative density of saplings (RD.sap)
#     # Stand Age (STDAGE)
#   standlist1 <- treelist%>%
#     group_by(UNITCD) %>%
#     summarize(QMD = sqrt( (sum(BA_AC[DIA >= 5 & STATUSCD == 1]) / sum(TPA_UNADJ[DIA >= 5 & STATUSCD == 1])) / 0.005454), # Calculate QMD 
#               RD.comm = sum(RD_AC[DIA >= 5 & STATUSCD == 1 & E_SPGRPCD != 23 & E_SPGRPCD != 43 & E_SPGRPCD != 48], na.rm = TRUE), # calculate relative density for commercial trees
#               RD.sap = sum(RD_AC[DIA < 5 & STATUSCD == 1 & E_SPGRPCD != 23 & E_SPGRPCD != 43 & E_SPGRPCD != 48], na.rm = TRUE), # calculate relative density for commercial saplings
#               STDAGE = mean(SITREE, na.rm = TRUE))
#     
#   # Summarize plot list data to create stand level covariates
#   # Included:
#     # State Code (STATECD)
#     # County Code (COUNTYCD)
#     # Iventory year (MEASYEAR)
#     # Basal Area Factor (BAF)
#     # Site class code (SITECLCD)
#   standlist2 <- plotlist %>%
#     group_by(UNITCD) %>%
#     summarize(STATECD = mean(STATECD),
#               COUNTYCD = mean(COUNTYCD),
#               MEASYEAR = mean(MEASYEAR),
#               BAF = mean(BAF),
#               SITECLCD = mean(SITECLCD, na.rm=T),
#               ECOSUBCD = min(ECOSUBCD),
#               ECOSEC = min(ECOSEC),
#               ECOPROV = min(ECOPROV))
#   
#   # Consolidate stand list data and add remaining parameters
#   standlist <- standlist1 %>%
#     left_join(standlist2, by = "UNITCD") %>%
#     left_join(slope, by = c("UNITCD"="StandNum")) %>%
#     #left_join(elevTable, by = "UNITCD") %>% # add mean elevation to stand list
#     left_join(fortypgrp_list, by = "UNITCD") %>% # add forest type group codes
#     left_join(RDDISTCD, by = c("UNITCD"="StandNum")) %>% # add road distance code
#     left_join(STD_LAT_LONG, by = c("UNITCD" ="StandNum")) %>% # add stand centroid lat long
#     mutate(PLT_CN = paste0(UNITCD,"_",MEASYEAR), # create unique identifier for plot and time
#            OWNGRPCD = ifelse(ownership=="Private", 40, 30), # add ownership code (40 = Private)
#            STDORGCD = 0, # Stand origin code for naturally regenerated stands = 0
#            Z = 1) %>% # add code for treatment plots
#     dplyr::select(PLT_CN,
#                   UNITCD, 
#                   STATECD, 
#                   COUNTYCD,
#                   STDAGE,
#                   SITECLCD,
#                   RD.sap,
#                   RD.comm,
#                   #ELEV,
#                   SLOPE = slopePer,
#                   QMD,
#                   RDDISTCD,
#                   MEASYEAR, 
#                   OWNGRPCD, 
#                   FORTYPGRP,
#                   STDORGCD,
#                   ECOSUBCD,
#                   ECOSEC,
#                   ECOPROV,
#                   LAT,
#                   LON,
#                   Z) %>%
#     add_column(PLOTID = NA)

# Read the stand list into the script which contains all covariates used for matching
standlist <-  read.csv(paste0("INPUT/",cohort,"_standlist.csv"))

  
#**********************************************
## Process FIA Data           ----
#**********************************************

#### Set forest group code ####
    # Forest group code is determined using forest type code
      data$COND$FORTYPGRP = as.factor(ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 106, 100,
                                                  ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 130, 120,
                                                         ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 143, 140,
                                                                ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 170, 160,
                                                                       ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 180, 170,
                                                                              ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 400, 380,
                                                                                     ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 500, 400,
                                                                                            ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 600, 500,
                                                                                                   ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 700, 600,
                                                                                                          ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 800, 700,
                                                                                                                 ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 900, 800,
                                                                                                                        ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 906, 900,
                                                                                                                               ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 913, 910,
                                                                                                                                      ifelse(as.numeric(as.character(data$COND$FORTYPCD)) < 999, 990, 
                                                                                                                                             999)))))))))))))))

#### Calculate key individual tree variables ####
    data$TREE$ba = data$TREE$DIA^2 * 0.005454 # Calculate individual tree basal
    data$TREE$BA_AC = data$TREE$ba * data$TREE$TPA_UNADJ # Calculate per acre basal area                      
    data$TREE$volbfnet.ac = data$TREE$VOLBFNET * data$TREE$TPA_UNADJ # Calculate per acre net cubic volume
    data$TREE$SPEC_GRAV = REF_SPECIES[match(data$TREE$SPCD, REF_SPECIES$SPCD), "WOOD_SPGR_GREENVOL_DRYWT"] # Specific gravity, used to calculate RD
    data$TREE$RD = 2.47 * (0.00015 + (0.00218 * data$TREE$SPEC_GRAV)) * ((data$TREE$DIA/10)^1.6) # Calc individual Ducey-Knapp RD - English units  
    data$TREE$RD_AC = data$TREE$RD * data$TREE$TPA_UNADJ # Calc per acre Ducey-Knapp RD
    
    # create stocking variables for trees in FIA data
    data$TREE$SFTWD_HRDWD = REF_SPECIES[match(data$TREE$SPCD, REF_SPECIES$SPCD), "SFTWD_HRDWD"] # Add softwood/hardwood code for forest typing
    data$TREE$STOCKING_SPGRPCD = REF_SPECIES[match(data$TREE$SPCD, REF_SPECIES$SPCD), "STOCKING_SPGRPCD"] # Add species group code from REF_SPECIES table to match stocking coefficients to trees
    data$TREE$FOREST_TYPE_SPGRPCD = REF_SPECIES[match(data$TREE$SPCD, REF_SPECIES$SPCD), "FOREST_TYPE_SPGRPCD"] # Add forest type species group from REF_SPECIES for forest typing
    data$TREE$b1 = stk_coef[match(data$TREE$STOCKING_SPGRPCD, stk_coef$STOCKING_SPGRPCD), "b1"] # add stocking coefficient b1
    data$TREE$b1 = ifelse(is.na(data$TREE$b1), 0, data$TREE$b1) # Set stocking coefficient to zero if there is no coefficient
    data$TREE$b0 = stk_coef[match(data$TREE$STOCKING_SPGRPCD, stk_coef$STOCKING_SPGRPCD), "b0"] # add stocking coefficient b0
    data$TREE$b0 = ifelse(is.na(data$TREE$b0), 0, data$TREE$b0)
    data$TREE$STOCKING_ALT <- ifelse(data$TREE$STATUSCD == 1 & data$TREE$DIA >= 5 & !is.na(data$TREE$b0), 
                                     data$TREE$b0 * data$TREE$DIA ^ (data$TREE$b1) * data$TREE$TPA_UNADJ/1,
                                     0)   #generate stocking value. Set stockability proportion to 1. Arner et al Step 1

    ###### Add forest type group assignment to FIA data ####### 
    
    # Create table with total stocking values of each species group
    # This process is based on the FIA process for assigning forest type group codes by species composition
    
    FIA_fortypgrp_list <- data$TREE %>%
      group_by(PLT_CN) %>%
      summarise(STOCKING_ALT_SFTWD = sum(STOCKING_ALT * (SFTWD_HRDWD == "S"), na.rm=T),
                STOCKING_ALT_HRDWD = sum(STOCKING_ALT * (SFTWD_HRDWD == "H"), na.rm=T),
                STOCKING_ALT_ERCEDAR = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD == 64), na.rm=T),
                STOCKING_ALT_OAKPIN = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(41, 42, 44:54, 64)), na.rm=T),
                STOCKING_ALT_EXHWDS = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(144, 145, 146, 148)), na.rm=T),
                STOCKING_ALT_EXSFTWD = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(70,71,72,73)), na.rm=T),
                STOCKING_ALT_NROBLWA = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(85, 108)), na.rm=T),
                STOCKING_ALT_RMWASHETC = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(95, 103, 94, 102, 121)), na.rm=T),
                STOCKING_ALT_MBB = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(66, 96, 98, 107, 124, 95, 103, 94, 102, 121)), na.rm=T),
                STOCKING_ALT_ASPBRCH = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(99, 117, 119, 154, 224)), na.rm=T),
                STOCKING_ALT_ELMASHCW = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(91, 97, 100, 104, 115, 116, 118, 123, 129, 135, 137, 208, 222)), na.rm=T),
                STOCKING_ALT_ALDRMAPL = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(130,131)), na.rm=T),
                STOCKING_ALT_LOBSHRTP = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(44, 45, 47, 49:52, 54)), na.rm=T),
                STOCKING_ALT_LNGLFSLH = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(46,48)), na.rm=T),
                STOCKING_ALT_OKGMCYP = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(59, 61, 87, 90, 111, 112, 114, 127, 128, 143, 151, 229)), na.rm=T),
                STOCKING_ALT_ASHBCHBLCH = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(103, 105, 102, 121)), na.rm=T),
                STOCKING_ALT_OAKHCK = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(81:86, 88, 89, 92, 93, 101, 108, 110, 120, 122, 124, 202, 205, 206, 207, 211, 125, 201, 203, 204, 109, 209, 94, 106, 113, 95, 105)), na.rm=T),
                STOCKING_ALT_ESPRFIR = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(16, 17, 55, 56, 58, 60, 65)), na.rm=T),
                STOCKING_ALT_RWJPIN = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD %in% c(41, 42, 53, 66)), na.rm=T),
                STOCKING_ALT_TOTAL = sum(STOCKING_ALT, na.rm=T))
    

    # Assign to groups by species composition
    
    
    FIA_fortypgrp_list$STOCKING_ALT_OAKHCK= ifelse(FIA_fortypgrp_list$STOCKING_ALT_ASHBCHBLCH < (0.5 * FIA_fortypgrp_list$STOCKING_ALT_TOTAL) & FIA_fortypgrp_list$STOCKING_ALT_OAKHCK > (0.05 * FIA_fortypgrp_list$STOCKING_ALT_TOTAL),
                                               FIA_fortypgrp_list$STOCKING_ALT_OAKHCK + FIA_fortypgrp_list$STOCKING_ALT_ASHBCHBLCH,
                                               FIA_fortypgrp_list$STOCKING_ALT_OAKHCK)
    
    FIA_fortypgrp_list$STOCKING_ALT_MBB = ifelse(FIA_fortypgrp_list$STOCKING_ALT_NROBLWA < (0.5 * FIA_fortypgrp_list$STOCKING_ALT_TOTAL) & FIA_fortypgrp_list$STOCKING_ALT_MBB > (0.05 * FIA_fortypgrp_list$STOCKING_ALT_TOTAL),  
                                             FIA_fortypgrp_list$STOCKING_ALT_MBB + FIA_fortypgrp_list$STOCKING_ALT_NROBLWA,  
                                             FIA_fortypgrp_list$STOCKING_ALT_MBB)
    
    FIA_fortypgrp_list <- FIA_fortypgrp_list %>%
      rowwise() %>%
      mutate(STOCKING_MAX = ifelse(
        STOCKING_ALT_HRDWD > STOCKING_ALT_SFTWD,
        max(STOCKING_ALT_EXHWDS, STOCKING_ALT_MBB, STOCKING_ALT_ASPBRCH, STOCKING_ALT_ELMASHCW, STOCKING_ALT_ALDRMAPL, STOCKING_ALT_OKGMCYP, STOCKING_ALT_OAKHCK), 
        max(STOCKING_ALT_EXSFTWD, STOCKING_ALT_LOBSHRTP, STOCKING_ALT_LNGLFSLH, STOCKING_ALT_ESPRFIR, STOCKING_ALT_RWJPIN, STOCKING_ALT_ERCEDAR))) %>%
      ungroup()
    
    # Assign forest type groups based on max stocking values
    FIA_fortypgrp_list$FORTYPGRP <- ifelse(FIA_fortypgrp_list$STOCKING_ALT_HRDWD > FIA_fortypgrp_list$STOCKING_ALT_SFTWD & FIA_fortypgrp_list$STOCKING_ALT_OAKPIN >= (0.25 * FIA_fortypgrp_list$STOCKING_ALT_TOTAL), 400,
                                           ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_EXHWDS, 990, 
                                                  ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_MBB, 800, 
                                                         ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_ASPBRCH, 900,
                                                                ifelse(FIA_fortypgrp_list$STOCKING_MAX ==  FIA_fortypgrp_list$STOCKING_ALT_ELMASHCW, 700,
                                                                       ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_ALDRMAPL, 910,
                                                                              ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_OKGMCYP, 600,
                                                                                     ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_OAKHCK, 500,
                                                                                            ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_EXSFTWD, 380,
                                                                                                   ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_LOBSHRTP,160,
                                                                                                          ifelse(FIA_fortypgrp_list$STOCKING_MAX ==  FIA_fortypgrp_list$STOCKING_ALT_LNGLFSLH, 140,
                                                                                                                 ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_ESPRFIR, 120,
                                                                                                                        ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_ERCEDAR, 170,
                                                                                                                               ifelse(FIA_fortypgrp_list$STOCKING_MAX == FIA_fortypgrp_list$STOCKING_ALT_RWJPIN, 100,
                                                                                                                        999))))))))))))))
    FIA_fortypgrp_list$FORTYPGRP <- as.numeric(FIA_fortypgrp_list$FORTYPGRP)

    
    # Join forest type codes to FIA data
    data$COND$FORTYPGRP_ALT <- FIA_fortypgrp_list[match(data$COND$PLT_CN, FIA_fortypgrp_list$PLT_CN), "FORTYPGRP"]
    
    total_stocking <- data$TREE %>%
      group_by(PLT_CN) %>%
      summarize(STOCKING = sum(STOCKING, na.rm=T),
                STOCKING_ALT = sum(STOCKING_ALT))
    
#### Organize FIA data ####
    
  ## Create columns with key data from previous measurement for each tree
  # This will facilitate comparison of current measurement with previous.

## COND TABLE
  ## Create previous data for COND table
    # link COND table to PREV and SUBS PLT_CN
    data$COND$PREV_PLT_CN = data$PLOT[match(data$COND$PLT_CN, data$PLOT$CN), "PREV_PLT_CN"]
    data$COND$SUBS_PLT_CN = data$PLOT[match(data$COND$PLT_CN, data$PLOT$PREV_PLT_CN), "CN"]

  # Append previous and subsequent data to COND table
    data$COND$PLOT_STATUS_CD = data$PLOT[match(data$COND$PLT_CN, data$PLOT$CN), "PLOT_STATUS_CD"]
    data$COND$plot_status_cd_prev = data$COND[match(data$COND$PREV_PLT_CN, data$COND$PLT_CN), "PLOT_STATUS_CD"]
    data$COND$plot_status_cd_subs = data$COND[match(data$COND$PLT_CN, data$COND$PREV_PLT_CN), "PLOT_STATUS_CD"]

  # Append MEASYEAR to COND table
    data$COND$MEASYEAR = data$PLOT[match(data$COND$PLT_CN, data$PLOT$CN), "MEASYEAR"]
    
## TREE TABLE
  ## Create previous data for TREE table
    data$TREE$PREV_PLT_CN = data$PLOT[match(data$TREE$PLT_CN, data$PLOT$CN), "PREV_PLT_CN"]
    data$TREE$SUBS_PLT_CN = data$PLOT[match(data$TREE$PLT_CN, data$PLOT$PREV_PLT_CN), "CN"]

  # Append columns with previous and subsequent data
  # Inventory Year
    data$TREE$MEASYEAR = data$PLOT[match(data$TREE$PLT_CN, data$PLOT$CN), "MEASYEAR"]
    data$TREE$MEASYEAR.prev = data$TREE[match(data$TREE$PREV_PLT_CN, data$TREE$PLT_CN), "MEASYEAR"]
    data$TREE$MEASYEAR.subs = data$PLOT[match(data$TREE$PLT_CN, data$PLOT$PREV_PLT_CN), "MEASYEAR"]
    
  # Append previous individual tree measurements
    data$TREE$BA_AC.prev = data$TREE[match(data$TREE$PREV_TRE_CN, data$TREE$CN), "BA_AC"]
    data$TREE$RD_AC.prev = data$TREE[match(data$TREE$PREV_TRE_CN, data$TREE$CN), "RD_AC"]
    data$TREE$diam.prev = data$TREE[match(data$TREE$PREV_TRE_CN, data$TREE$CN), "DIA"]
    data$TREE$TPA.prev = data$TREE[match(data$TREE$PREV_TRE_CN, data$TREE$CN), "TPA_UNADJ"]
    data$TREE$statuscd.prev = data$TREE[match(data$TREE$PREV_TRE_CN, data$TREE$CN), "STATUSCD"]
    data$TREE$treeclcd.prev = data$TREE[match(data$TREE$PREV_TRE_CN, data$TREE$CN), "TREECLCD"]
    data$TREE$volbfnet.ac.prev = data$TREE[match(data$TREE$PREV_TRE_CN, data$TREE$CN), "volbfnet.ac"]


## Add ID variable which is unique for all physical plot locations, but not across time
  # Concatenates state code, county code, and plot ID to create unique code
    data$PLOT$PLOTID = as.factor(paste0(data$PLOT$STATECD, ".", data$PLOT$COUNTYCD, ".", data$PLOT$PLOT))
    data$COND$PLOTID = as.factor(paste0(data$COND$STATECD, ".", data$COND$COUNTYCD, ".", data$COND$PLOT))
    data$TREE$PLOTID = as.factor(paste0(data$TREE$STATECD, ".", data$TREE$COUNTYCD, ".", data$TREE$PLOT))


## Adjust ECOSUBCD variables
  # Create eco section and province
    data$PLOT$ECOSEC = as.factor(substr(data$PLOT$ECOSUBCD, 1, nchar(as.character(data$PLOT$ECOSUBCD)) - 1))
    data$PLOT$ECOPROV = as.factor(substr(data$PLOT$ECOSUBCD, 1, nchar(as.character(data$PLOT$ECOSUBCD)) - 2))

  # Remove preceding or trailing blank spaces
    data$PLOT$ECOSEC = str_trim(data$PLOT$ECOSEC)
    data$PLOT$ECOPROV = str_trim(data$PLOT$ECOPROV)
    
  # Append to COND and TREE table
    data$COND$ECOSEC = data$PLOT[match(data$COND$PLT_CN, data$PLOT$CN), "ECOSEC"]
    data$COND$ECOPROV = data$PLOT[match(data$COND$PLT_CN, data$PLOT$CN), "ECOPROV"]
    data$TREE$ECOSEC = data$PLOT[match(data$TREE$PLT_CN, data$PLOT$CN), "ECOSEC"]
    data$TREE$ECOPROV = data$PLOT[match(data$TREE$PLT_CN, data$PLOT$CN), "ECOPROV"]


#**********************************************
## Subset Data           ----
#**********************************************

# Backup full data
    data.full = data
#    data = data.full
    #data.tree.save <- as.data.frame(data.full$TREE)
    #data.cond.save <- as.data.frame(data.full$COND)
    #data.cond.save <- as.data.frame(data.full$PLOT)
    #write.csv(data.tree.save, paste0(output.folder, "/", cohort, "_FIA_tree_dataprep.csv"))
    #write.csv(data.cond.save, paste0(output.folder, "/", cohort, "_FIA_cond_dataprep.csv"))
    #write.csv(data.plot.save, paste0(output.folder, "/", cohort, "_FIA_plot_dataprep.csv"))
    

    ## Assess which year all states have remeasured plots
    table(data$COND$INVYR, data$COND$STATECD)
    
##### Subset data to include only plots within the defined donor pool, removing all:

#     1. Periodic inventory points;
#     2. Single Condition plots;
#     3. Plots without accessible, sampled forestland;
#     4. Publicly owned land;
#     5. Land reserved from timber harvesting (as identified in FIA dataset);
#     6. Plots that do not meet program eligibility requirements. Here, this is defined for the 
#         purpose of the trial as being within the specified ecological section.
#     7. Plots with prior and subsequent measurements, to enable C accting.

## 1. Remove periodic inventory plots (this is unnecessary, but makes subsequent code run faster, by eliminating a large portion of data)
    # Append KINDCD to cond table
    data$COND$KINDCD <- data$PLOT[match(data$COND$PLT_CN, data$PLOT$CN), "KINDCD"]
    # Subset by KINDCD
    data$COND <- data$COND[which(data$COND$KINDCD %in% c(1,2,3)), ]

## 2. Subset to remove plots that span multiple conditions
    # Subset the Condition table to remove plots with multiple condition status codes.
    # The condprop_unadj variable shows the proportion of a plot that is in a given condition.
    # Thus if condprop_unadj = 1, the plot is entirely within a single condition.
    # Remove plots from cond table that have multiple condition status codes (split cover types/owners/landuses)
    data$COND <- data$COND[which(data$COND$CONDPROP_UNADJ==1), ]
    
## 3. Remove trees from plots with PLOT_STATUS_CD indicating non-sampled plots (PLOT_STATUS_CD == 3)
    # Remove natively  from PLOT table
    data$COND$PLOT_STATUS_CD <- data$PLOT[match(data$COND$PLT_CN, data$PLOT$CN), "PLOT_STATUS_CD"]
    data$COND <- data$COND[which(data$COND$PLOT_STATUS_CD == 1), ]

## 4. Subset to match treatment stands
    # Remove natively  from COND table
    data$COND <- data$COND[which(data$COND$OWNGRPCD %in% standlist$OWNGRPCD), ]

## 5. Subset to remove plots outside the specified program eligibility requirements
    # Remove natively  from COND table
    data$COND <- data$COND[which(data$COND$FORTYPGRP %in% standlist$FORTYPGRP), ]

## 6. Subset to include only plots with previous and subsequent measurements
    # NOTE: this is important, as we need pre-harvest data to perform matching and post-harvest to allow rebound of C within monitoring period
    # subset to include only condition classes that had sampled forest area 
    #   both current and previous measurement interval
    data$COND <- data$COND[which(data$COND$PLOT_STATUS_CD == 1 &
                                data$COND$plot_status_cd_prev == 1), ]
    
## 7. Remove plot measurements that are older than 7 years or less than 5 years from the inventory date
    data$COND <- data$COND[which(data$COND$MEASYEAR >= (min(standlist$MEASYEAR) - 7)),]
    
## 8. Select only most recent measurements of plots.   
    data$COND <- data$COND %>%
            group_by(PLOTID) %>%
            filter(MEASYEAR == max(MEASYEAR)) %>%
      ungroup()
      
## Apply each of the previous subsets to the TREE and PLOT tables
    data$PLOT <- data$PLOT[which(data$PLOT$CN %in% unique(data$COND$PLT_CN)), ]
    data$TREE <- data$TREE[which(data$TREE$PLT_CN %in% unique(data$COND$PLT_CN)), ]
    
## Remove noncommercial species from assessment pool
    data$TREE <- data$TREE[which(data$TREE$SPGRPCD != 43 | data$TREE$SPGRPCD != 28 | data$TREE$SPGRPCD != 48), ]
    
    plot.pool <- unique(data$COND[, "PLT_CN"]) # Create list of PLT_CN for donor pool plots
    plot.pool <- pull(plot.pool) # Convert plot.pool to vector of numeric values
    
    # Remove donor pool plots that are within 1.6 km of any stand boundaries
    fia_lat_lon <- subset(data$PLOT[which(data$PLOT$CN %in% plot.pool),], select = c("CN", "LAT", "LON")) # subset FIA plot data based on selected plots and pull lat longs
    fia_lat_lon <- st_as_sf(fia_lat_lon, coords=c("LON", "LAT"), crs="EPSG:4269") # convert FIA plots to sf object
      print(paste0(length(fia_lat_lon$CN), " donor plots selected before buffering"))
    fia_buffered <- st_difference(standBdry_buff, fia_lat_lon) # select all plots that do not intersect with buffered stand boundaries
    fia_buffered_filt <-  unique(fia_buffered$CN)# subset dataframe of remaining plots
    fia_lat_lon <- filter(fia_lat_lon, CN %in% fia_buffered_filt)
      print(paste0(length(fia_lat_lon$CN), " donor plots selected after buffering"))
    
## Use final iteration of COND table to subset TREE table
    tree.pool <- data$TREE[which(data$TREE$PLT_CN %in% fia_buffered_filt), ]
    
    # subset plots where T2 measurement is more than 7 years after initial measurement  ######
    tree.pool <- tree.pool[which((tree.pool$MEASYEAR.prev - tree.pool$MEASYEAR) < 8),] 
    
    
## Assess which year all states have remeasured plots
    table(data$COND$MEASYEAR, data$COND$STATECD)
    
    
    # plot eligible donor pool plots (fuzzed FIA plots)
    #png(filename = paste0(output.folder, "Figures/",cohort,"EligibleDonorPlots", Sys.Date(),".png"),
        #width = 700, height = 400)
    ggplot() + 
      geom_sf(data=us_states[us_states$STUSPS %in% states,]) + 
      geom_sf(data=fia_lat_lon, 
              aes(fill='black'), 
              show.legend = "point") +
      geom_sf(data=st_as_sf(STD_CENT), 
              aes(fill='red'), 
              show.legend = "point") +
      scale_fill_manual(values=c("A"='black', "B"='red'), 
                        labels=c("Donor Plots", "Selected Stands"))  
    dev.off()
#**********************************************
## Summarize Plot Data           ----
#**********************************************

    
## Create plot summaries for assessment pool plots
    control <- tree.pool %>%
      group_by(PLT_CN, PREV_PLT_CN, SUBS_PLT_CN, PLOTID, STATECD, MEASYEAR, 
               ECOPROV, ECOSEC) %>%
      summarize(RD.comm = sum(RD_AC[DIA >= 5 & STATUSCD == 1], na.rm = TRUE),
                RD.sap = sum(RD_AC[DIA < 5 & STATUSCD == 1], na.rm = TRUE),
                BA_AC = sum(BA_AC[STATUSCD == 1 & DIA >= 5], na.rm = TRUE),
                TPA = sum(TPA_UNADJ[STATUSCD == 1 & DIA >= 5], na.rm = TRUE),
                QMD = sqrt((BA_AC/TPA)/0.005454),
                MEASYEAR.prev = mean(MEASYEAR.prev, na.rm = TRUE))
    control <- as.data.frame(control)
    
    # Calculate Site Class independently to standardize site class codes between treatment plots and FIA plots
    # Remove plots with multiple condition classes from SITETREE table and filter out plots with estimated SI (Method = 3 or 4)
    SITETREE <- data$SITETREE %>%
      filter(!duplicated(PLT_CN),
             METHOD==1|METHOD==2)
    
    # Join site tree table to selected plots
    data_SI_WESTFALL <- SITETREE %>%
      filter(PLT_CN %in% control$PLT_CN) %>%
      left_join(SI_LOOKUP_WESTFALL, by = "SPCD") %>%
      dplyr::select(CN, 
                    PLT_CN, 
                    SPCD, 
                    DIA, 
                    HT, 
                    AGEDIA, 
                    SI_SP_GRP, 
                    MIN_HT, 
                    MAX_HT, 
                    MIN_AGE, 
                    MAX_AGE, 
                    b0, 
                    b1, 
                    b2, 
                    b3,
                    L,
                    D,
                    SITECL1_LOW, 
                    SITECL2_LOW, 
                    SITECL3_LOW, 
                    SITECL4_LOW, 
                    SITECL5_LOW,
                    SITECL6_LOW,
                    SITREE)
    
    data_SI_SV <- data_SI_WESTFALL %>%
      filter(is.na(SI_SP_GRP)) %>%
      dplyr::select(CN, PLT_CN, SPCD, DIA, HT, AGEDIA, SITREE) %>%
      inner_join(SI_EQCODE, by = "SPCD") %>%
      left_join(SI_LOOKUP_SV, by = "SI_EQ")
    
    # Run SI equation for each site tree with parameters: ht, coefficients, age
    data_SI_WESTFALL <- data_SI_WESTFALL %>%
      filter(!is.na(SI_SP_GRP)) %>%
      mutate(RES = HT - (4.5 + b0 * (1 - exp(-b1 * AGEDIA))^b2),
             Z = (1 - exp(-b1 * AGEDIA))^b2 + (1 - exp(-b1 * AGEDIA))^b2 * b3 * log((1 - exp(-b1 * AGEDIA))) * b0,
             R = exp(L * AGEDIA),
             X = D * Z * ((Z * D * Z + R)^(-1)) * RES,
             SI = if_else(HT < MAX_HT &
                            HT > MIN_HT &
                            AGEDIA > MIN_AGE &
                            AGEDIA < MAX_AGE, 
                          4.5 + (b0 + X) * (1 - exp(-b1 * 50))^(b2 + b3 * X), NA))
    
    
    # Equation source: Scott, C.T.; Voorhis, N.G.  1986.  Northeastern Forest Survey site index equations and site productivity classes.  North. J. Appl. For. 3:144-148
    data_SI_SV <- data_SI_SV %>%
                    mutate(SI = if_else(AGEDIA > MIN_AGE &
                                        AGEDIA < MAX_AGE, 
                                HT / (b1 * (1.0 - exp(-b2 * AGEDIA))^b3), NA))
    
    data_SI <- data_SI_WESTFALL %>%
      dplyr::select(PLT_CN, DIA, HT, SITECL1_LOW, 
                    SITECL2_LOW, 
                    SITECL3_LOW, 
                    SITECL4_LOW, 
                    SITECL5_LOW,
                    SITECL6_LOW, 
                    SITREE, 
                    SI) %>%
      bind_rows(data_SI_SV %>% dplyr::select(PLT_CN, DIA, HT, 
                                             SITECL1_LOW, 
                                             SITECL2_LOW, 
                                             SITECL3_LOW, 
                                             SITECL4_LOW, 
                                             SITECL5_LOW,
                                             SITECL6_LOW, 
                                             SITREE, 
                                             SI))
    
    
    # Compare FIA generated SI to manually computed SI
    ggplot(data=data_SI) +
      geom_point(aes(SITREE, SI))
    arrange(data_SI, desc(SI))
    
    # Bin site index into site class code based on species and min SI values (source: NIMS_REF_RANGE_SITECLASS, LOC_GRP = S24)
    data_SI_bins <- data_SI %>%
      mutate(SITECLCD = case_when(
        SI < SITECL5_LOW | is.na(SITECL5_LOW) & SI < SITECL6_LOW ~ 6,
        SI < SITECL4_LOW & SI >= SITECL5_LOW | is.na(SITECL4_LOW) & SI >= SITECL5_LOW ~ 5,
        SI < SITECL3_LOW & SI >= SITECL4_LOW | is.na(SITECL3_LOW) & SI >= SITECL4_LOW ~ 4,
        SI < SITECL2_LOW & SI >= SITECL3_LOW | is.na(SITECL2_LOW) & SI >= SITECL3_LOW ~ 3,
        SI < SITECL1_LOW & SI >= SITECL2_LOW | is.na(SITECL1_LOW) & SI >= SITECL2_LOW ~ 2,
        SI >= SITECL1_LOW & is.numeric(SITECL1_LOW) ~ 1))
    
    # Add Site Index and Site Class Code to donor pool list
    control_SI <- data_SI_bins %>% 
      dplyr::select(c("PLT_CN","SI","SITECLCD"))
    
    control <- inner_join(control, control_SI %>% dplyr::select(PLT_CN, SITECLCD), by = "PLT_CN")
      
    
    ## Add covariates
    control$UNITCD = NA
    control$PREV_PLT_CN = data.full$COND[match(control$PLT_CN, data.full$PLOT$PLT_CN), "PREV_PLT_CN"]
    control$STATECD = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "STATECD"]
    control$COUNTYCD = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "COUNTYCD"]
    control$STDAGE = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "STDAGE"]
    #control$SITECLCD = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "SITECLCD"]
    control$FORTYPGRP = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "FORTYPGRP"]
    control$STDORGCD = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "STDORGCD"]
    #control$ELEV = data.full$PLOT[match(control$PLT_CN, data.full$PLOT$CN), "ELEV"]
    control$SLOPE = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "SLOPE"]
    control$RDDISTCD = data.full$PLOT[match(control$PLT_CN, data.full$PLOT$CN), "RDDISTCD"]
    control$LAT = data.full$PLOT[match(control$PLT_CN, data.full$PLOT$CN), "LAT"]
    control$LON = data.full$PLOT[match(control$PLT_CN, data.full$PLOT$CN), "LON"]
    control$CONDPROP_UNADJ = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "CONDPROP_UNADJ"]
    control$OWNGRPCD = data.full$COND[match(control$PLT_CN, data.full$COND$PLT_CN), "OWNGRPCD"]
    # Assign treatment variables
    control$Z = 0

  
  # Subset to exclude any plots with split conditions or "ineligible" forest type
  # consider also excluding "ineligible" stocking levels, i.e. <2000 bdft/ac
  # This is necessary, despite above subsetting for CONDPROP_UNADJ, because there are plots which
  # were unsplit at the harvest measurement, but were split at previous, pre-harvest measurement.
  control = control[which(control$CONDPROP_UNADJ == 1 &
                        control$FORTYPGRP %in% standlist$FORTYPGRP), ]
  
  print(paste0("There are ", nrow(control), " selected in the donor pool for matching"))



## Reshape control dataset

control = control %>%
  ungroup() %>%
  dplyr::select(PLT_CN,
                PLOTID,
                MEASYEAR,
                ECOSEC,
                ECOPROV,
                OWNGRPCD,
                FORTYPGRP,
                STDORGCD,
                RDDISTCD,
                SLOPE,
                QMD,
                RD.comm,
                RD.sap,
                STDAGE,
                SITECLCD,
                #ELEV,
                LAT,
                LON,
                Z)

## Write summaries to file
write.csv(control, paste0(output.folder, cohort, "_DonorPool_.csv"), row.names = FALSE)

## Reshape treatment dataset to mirror controls so they can be compiled
standlist = standlist %>%
  dplyr::select(PLT_CN,
                PLOTID = UNITCD,
                MEASYEAR,
                ECOSEC,
                ECOPROV,
                OWNGRPCD,
                FORTYPGRP,
                STDORGCD,
                RDDISTCD,
                SLOPE,
                QMD,
                RD.comm,
                RD.sap,
                STDAGE,
                SITECLCD,
                #ELEV,
                LAT,
                LON,
                Z)




## Make sure all treatments and controls have a value in each variable
control = control[complete.cases(control),] 
standlist = standlist[complete.cases(standlist),]


# Upload standlist to hard drive for records
write.csv(standlist, paste0(output.folder, cohort, "_standlist.csv"))

rm(tree.pool)


# Convert to dataframe and write to file
    write.csv(control, paste0(output.folder, cohort, "_Control_plot summary_at enrollment.csv"), row.names = FALSE)
    #control = read.csv(paste0(output.folder, cohort, "_Control_plot summary_at enrollment.csv"), header = TRUE, sep = ",")

    #***********************************************************************#


#**********************************************
## Prep matching dataframes    ----
#**********************************************
  
## Create dataframes to receive output
  # Define dataframe treat.comp.controls to accept treatments and their respective composite controls
      treat.comp.controls <- data.frame(PLT_CN = as.factor(character()), 
                           PREV_PLT_CN = as.factor(character()),
                           STATECD = as.factor(character()),
                           COUNTYCD = as.factor(character()),
                           PLOTID = as.factor(character()),
                           MEASYEAR = as.integer(character()), 
                           ECOSEC = as.factor(character()),
                           ECOPROV = as.factor(character()),
                           OWNGRPCD = as.factor(character()),
                           FORTYPGRP = as.factor(character()),
                           RDDISTCD = as.factor(character()),
                           SLOPE = as.numeric(character()),
                           STDORGCD = as.factor(character()),
                           QMD = as.numeric(character()),
                           RD.comm = as.numeric(character()),
                           RD.sap = as.numeric(character()),
                           STDAGE = as.numeric(character()),
                           SITECLCD = as.numeric(character()),
                           #ELEV = as.numeric(character()),
                           LAT = as.numeric(character()),
                           LON = as.numeric(character()),
                           Z = as.numeric(character()),
                           Mahalanobis.Distance = as.numeric(character()),
                           inv.Mahalanobis.Distance = as.numeric(character()),
                           weight = as.numeric(character()),
                           stringsAsFactors=FALSE) 

  # Define dataframe matches to accept treatment units and all 10 of their best matches
    matches <- treat.comp.controls
      
#**********************************************
## Matching                     ----
#**********************************************

# Matching will be performed individually for each plot. That is, 
# each treatment will be isolated, matched and archived.
# A for loop is used to iterate through the various treatments in a project.

# set starting value for treatment
    treatment <- 1
    control.plot <- 2

# backup control data, so it can be reset at each iteration
    control.backup <- control
    
# create data.frame for listing exact matching requirements
    donor.pool.check <- as.data.frame(matrix(ncol=4, nrow=0))
    colnames(donor.pool.check) <- c("PLOTID", "GROUPING","GROUPCODE", "PLOTCOUNT")
    

## For loop 
    # 1. defines a unique donor pool for each treatment plot,
    # 2. creates physical distance variable from treatment to each potential donor
    # 3. performs matching
    # 4. generates treelist
    # 5. adds treelist, matches, etc. to various output/archive dataframes

for (treatment in 1:nrow(standlist)) {
  ## Define initial dataset as the treatment to be matched and the full list of control plots
      # Backup control data
      control <- control.backup
      
      # Row bind the treatment unit to the donor pool
      data.matching <- rbind(standlist[treatment, ], control)
      
      # Remove any plots that are missing covariates necessary for matching
      data.matching <- data.matching[complete.cases(data.matching), ]
      
      # First define the parameters specific to the treatment being matched
      FORTYPGRP = data.matching$FORTYPGRP[1]
      ECOSEC = data.matching$ECOSEC[1]
      ECOPROV = data.matching$ECOPROV[1]
      STDORGCD = data.matching$STDORGCD[1]
      
      # Match to selected grouping level
      # Define final donor pool as an recursive if statement that verifies there are at least 50 available donor plots for matching

      d = 
        if (nrow(data.matching[which(data.matching$ECOSEC == ECOSEC & 
                                     data.matching$FORTYPGRP == FORTYPGRP), ]) 
            > 100) {
          c(standlist$PLT_CN[treatment],"ECOSEC", 1, nrow(data.matching[which(data.matching$ECOSEC == ECOSEC & 
                                                                        data.matching$FORTYPGRP == FORTYPGRP), ]))
        } else if (nrow(data.matching[which(data.matching$ECOPROV == ECOPROV &
                                            data.matching$FORTYPGRP == FORTYPGRP), ])
                   > 100) {
          c(standlist$PLT_CN[treatment],"ECOPROV", 2, nrow(data.matching[which(data.matching$ECOPROV == ECOPROV &
                                                                         data.matching$FORTYPGRP == FORTYPGRP), ]))
        } else if (nrow(data.matching[which(data.matching$FORTYPGRP == FORTYPGRP), ])
                   > 100) {
          c(standlist$PLT_CN[treatment],"FORTYPGRP", 3, nrow(data.matching[which(data.matching$FORTYPGRP == FORTYPGRP), ]))
        } else print(paste0("Insufficient plots available in donor pool (n = ", (nrow(data.matching)), ") for plot ", 
                            standlist$PLT_CN[treatment ] ,", please adjust covariates and try again.")) 
      
      donor.pool.check[nrow(donor.pool.check)+1,] <- d
      
      data.matching = 
        if (donor.pool.check$GROUPING[treatment] == "ECOSEC") {
          data.matching[which(data.matching$ECOSEC == ECOSEC &
                                data.matching$FORTYPGRP == FORTYPGRP),]
        } else if (donor.pool.check$GROUPING[treatment] == "ECOPROV") {
          data.matching[which(data.matching$ECOPROV == ECOPROV &
                                data.matching$FORTYPGRP == FORTYPGRP),]
        } else if (donor.pool.check$GROUPING[treatment] == "FORTYPGRP") {
          data.matching[which(data.matching$FORTYPGRP == FORTYPGRP),]
        } else NA
      
      
      # Match to Stand Origin Code
      data.matching = data.matching[which(data.matching$STDORGCD==STDORGCD),]
      
  ## Calculate physical distance
      # Create blank vector for distances between each control plot and the stand centroid
      data.matching$dist = as.numeric(NA)
      data.matching[which(data.matching$Z == 1), "dist"] = 0
      
      ## A second, nested for loop is used to calculate distance from treatment to 
      # each control plot in the donor pool using distGeo from the geosphere package. 
      for (control.plot in 2:nrow(data.matching)) {
        
        #*** Note:  all lon/lat values should be in NAD83 GCS to reflect FIA standards. Arguments a and f below refer to
        #           radius and flattening values for GRS80 ellipsoid used in NAD83.
        dist = distGeo(p1 = c(as.character(data.matching[1, "LON"]), as.character(data.matching[1, "LAT"])), 
                       p2 = c(as.character(data.matching[control.plot, "LON"]), as.character(data.matching[control.plot, "LAT"])),
                       a = 6378137, f = 1/298.257222101)
    
        
        # Default distance is in meters. Convert to miles.
        dist.mi = dist/1609.34
        
        # Store resulting distance in data.matching dataframe
        data.matching[control.plot, "dist"] = dist.mi
        
        # Advance counter for control plot
        control.plot = control.plot + 1
      } 

  
  ## Calculate Mahalanobis Distances for each treatment plot
  
      # Make the FIA PLT_CN the row name for each row
      # This keeps the PLT_CN available through matching protocol
      rownames(data.matching) = data.matching$PLT_CN
      
      # Use the optmatch package in R to calculate Mahalanobis distances for each pairwise comparison of treat x control plot
      # this does not include ANY other columns besides the MDs
      m.dists = as.data.frame(t(match_on(Z ~ dist + SLOPE + RD.sap + QMD + RDDISTCD + STDAGE + SITECLCD + RD.comm, 
                                         data = data.matching, method = "mahalanobis")))
      
      
      ## Create column with FIA plot CNs from row.names
      m.dists$PLT_CN = row.names(m.dists)
    
      ## Order m.dists by Mahalanobis Distance
      colnames(m.dists)[1] = "Mahalanobis.Distance"
      m.dists = m.dists[order(m.dists[,"Mahalanobis.Distance"]), ]
  
  
      
  ## Create new dataframe with the treatment and its donor plots
      # Create plot list from the plot identifier for the treatment plot, and the top k matches (by Mahalanobis distance)
      plot.list = as.vector(c(data.matching$PLT_CN[1], m.dists$PLT_CN[1:k]))
    
      # Create a dataframe with all data corresponding to that plot list
      match = data.matching[which(data.matching$PLT_CN %in% plot.list), ]
  
  ## Calculate constituent plot weights
      # Append MDs to appropriate plots in match dataframe
      match$Mahalanobis.Distance = m.dists$Mahalanobis.Distance[match(match$PLT_CN, m.dists$PLT_CN)]
      
      # calculate inverse MDs
      match$inv.Mahalanobis.Distance = (1/match$Mahalanobis.Distance) * 100
      
      # calculate weight
      match$weight = match$inv.Mahalanobis.Distance / sum(match$inv.Mahalanobis.Distance, na.rm = TRUE)
      
  ## Create composite control 
      # fill basic info from treatment plot, setting all non-identical values to NA
      # calculate weighted mean values for numeric covariates
      # revert to dataframe
      match=as.data.frame(match)
      
      ## assign a few variables that have issues in the assigned row above. not sure why
      # Create new level for PLT_CN factor, called Composite Control
      match$PLT_CN = as.factor(match$PLT_CN)
      levels(match$PLT_CN) <- c(levels(match$PLT_CN), "Composite Control")
    
      donors = match[2:nrow(match),]
      match[nrow(match)+1, ] = c("PLT_CN" = "Composite Control", 
                      "PLOTID" = NA,
                      "MEASYEAR" = round(weighted.mean(donors$MEASYEAR, donors$weight),0),
                      "ECOSEC" = NA,
                      "ECOPROV" = NA,
                      "OWNGRPCD" = match[1, "OWNGRPCD"],
                      "FORTYPGRP" = NA,
                      "STDORGCD" = weighted.mean(donors$STDORGCD, donors$weight),
                      "RDDISTCD" = weighted.mean(donors$RDDISTCD, donors$weight),
                      "SLOPE" = weighted.mean(donors$SLOPE, donors$weight),
                      "QMD" = weighted.mean(donors$QMD, donors$weight),
                      "RD.comm" = weighted.mean(donors$RD.comm, donors$weight),
                      "RD.sap" = weighted.mean(donors$RD.sap, donors$weight),
                      "STDAGE" = weighted.mean(donors$STDAGE, donors$weight),
                      "SITECLCD" = weighted.mean(donors$SITECLCD, donors$weight),
                      #"ELEV" = weighted.mean(donors$ELEV, donors$weight),
                      "LAT" = weighted.mean(donors$LAT, donors$weight),
                      "LON" = weighted.mean(donors$LON, donors$weight),
                      "Z" = NA,
                      "dist" = weighted.mean(donors$dist, donors$weight),
                      "Mahalanobis.Distance" = NA,
                      "inv.Mahalanobis.Distance" = NA,
                      "weight" = NA
      )
      
    
      
      # Compute ECOSEC or ECOPROV or FORTYPGRP
      # from the data.matching dataframe. Importantly, this will help reflect whether the
      # data used for matching was from the same ecological section, province, etc. 
      match[nrow(match),"ECOSEC"] = if(length(getmode(data.matching$ECOSEC)) < 2) {getmode(data.matching$ECOSEC)} else{NA}
      match[nrow(match),"ECOPROV"] = if(length(getmode(data.matching$ECOPROV)) < 2) {getmode(data.matching$ECOPROV)} else{NA}
      match[nrow(match),"FORTYPGRP"] = mean(as.numeric(as.character(data.matching$FORTYPGRP)))
      
  
      ## create column denoting which treatment plot each control plot is associated with
      match$treat_PLT_CN = match[1, "PLT_CN"]
      match$treat_PLOTID = match[1, "PLOTID"]
  
      
  ## Store data and archive
      ## store treatment plot and composite control in treat.comp.controls dataframe
      treat.comp.controls = rbind(treat.comp.controls, match[c(1, nrow(match)), ])
      
      ## store treatment plot and k best matches in matches dataframe
      matches = rbind(matches, match[c(1:nrow(match)-1), ])
      ## Archive individual treatment-best matches dataframe
      write.csv(match, paste0(output.folder,cohort, "/INDIVIDUAL MATCHES/", cohort, "_match_", match[1, "PLT_CN"],".csv"), row.names = FALSE)
      
  
  ## Advance counter for treatment
  treatment = treatment + 1
}

    
# Convert numeric columns from character to numeric
matches[,8:23] <- as.data.frame(sapply(matches[,8:23], as.numeric))
treat.comp.controls[,8:23] <- as.data.frame(sapply(treat.comp.controls[,8:23], as.numeric))

# Create shapefile object of matched plot points for mapping
match_pts <- st_as_sf(matches, coords=c("LON", "LAT"), crs = "EPSG:4269")

# Create plot of treatments and control plots overlaid on map of US States
#png(filename = paste0(output.folder, "Figures/",cohort,"_MatchedPlots_", Sys.Date(),".png"),
 #   width = 700, height = 400)
ggplot() +
  geom_sf(data=us_states[us_states$STATEFP %in% c(STATECD),]) +
  geom_sf(data=match_pts, aes(color=as.factor(match_pts$Z))) +
  labs(title="Treatment and Matched Control Plots", color="Plots") +
  scale_color_manual(labels=c("Control", "Treatment"), values =c("dark blue", "red"))
dev.off()
    
# Remove temporary objects
rm(data.matching, donors, m.dists, match, control.backup, control.plot, dist, dist.mi, 
   ECOPROV, ECOSEC, FORTYPGRP, plot.list, treatment)


## Create dataframe containing only composite control plots, to compare to treatment units
comp.controls = treat.comp.controls[which(treat.comp.controls$PLT_CN == "Composite Control"), ]
write.csv(comp.controls, paste0(output.folder, cohort, "_Composite Controls.csv"), row.names = FALSE)


## Archive results and data used in matching
write.csv(treat.comp.controls, paste0(output.folder, cohort, "_Treatment and Composite Controls.csv"), row.names = FALSE)
write.csv(matches, paste0(output.folder, cohort, "_Treatment and Matched Controls.csv"), row.names = FALSE)
write.csv(data$TREE, paste0(output.folder, cohort, "_Donor Pool_Tree Table.csv"), row.names = FALSE)
write.csv(data$COND, paste0(output.folder, cohort, "_Donor Pool_Condition Table.csv"), row.names = FALSE)
write.csv(data$PLOT, paste0(output.folder, cohort, "_Donor Pool_Plot Table.csv"), row.names = FALSE)
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#############************END CODE HERE*****************################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





