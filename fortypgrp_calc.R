###############################
# Function to assign forest type coding based on 
# FIA forest type coding decision framework
# for the Northeastern US
# 
# The function requires treelist with the following parameters:
# 1. A stocking estimate
# 2. A plot ID named 'UNITCD'
# 3. FOREST_TYPE_SPGRPCD from FIA 
# 4. Softwood/Hardwood designation
#
# Edited by Ben Rifkin 5-23-2023
###############################

# Create table with total stocking values of each species group
# This process is based on the FIA process for assigning forest type group codes by species composition

fortypcode <- function(treelist) {
  
  fortypgrp <- treelist %>%
    group_by(UNITCD) %>%
    summarise(STOCKING_ALT_SFTWD = sum(STOCKING_ALT * (SFTWD_HRDWD == "S"), na.rm=T),
              STOCKING_ALT_HRDWD = sum(STOCKING_ALT * (SFTWD_HRDWD == "H"), na.rm=T),
              STOCKING_ALT_ERCEDAR = sum(STOCKING_ALT * (FOREST_TYPE_SPGRPCD == 64), na.rm=T),
              STOCKING_ALT_OAKPIN = sum(STOCKING_ALT *(FOREST_TYPE_SPGRPCD %in% c(41, 42, 44:54, 64)), na.rm=T),
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
              STOCKING_ALT_TOTAL = sum(STOCKING_ALT, na.rm=T)
    )
  
  # Assign to groups by species composition
  fortypgrp$STOCKING_ALT_OAKHCK <- ifelse(fortypgrp$STOCKING_ALT_ASHBCHBLCH < (0.5 * fortypgrp$STOCKING_ALT_TOTAL) & fortypgrp$STOCKING_ALT_OAKHCK > (0.05 * fortypgrp$STOCKING_ALT_TOTAL),
                                          fortypgrp$STOCKING_ALT_OAKHCK + fortypgrp$STOCKING_ALT_ASHBCHBLCH,
                                          fortypgrp$STOCKING_ALT_OAKHCK)
  
  fortypgrp$STOCKING_ALT_MBB <- ifelse(fortypgrp$STOCKING_ALT_NROBLWA < (0.5 * fortypgrp$STOCKING_ALT_TOTAL) & fortypgrp$STOCKING_ALT_MBB > (0.05 * fortypgrp$STOCKING_ALT_TOTAL),  
                                       fortypgrp$STOCKING_ALT_MBB + fortypgrp$STOCKING_ALT_NROBLWA,  
                                       fortypgrp$STOCKING_ALT_MBB)
  
  fortypgrp <- fortypgrp %>%
    rowwise() %>%
    mutate(STOCKING_MAX = ifelse(
      STOCKING_ALT_HRDWD > STOCKING_ALT_SFTWD,
      max(STOCKING_ALT_EXHWDS, STOCKING_ALT_MBB, STOCKING_ALT_ASPBRCH, STOCKING_ALT_ELMASHCW, STOCKING_ALT_ALDRMAPL, STOCKING_ALT_OKGMCYP, STOCKING_ALT_OAKHCK), 
      max(STOCKING_ALT_EXSFTWD, STOCKING_ALT_LOBSHRTP, STOCKING_ALT_LNGLFSLH, STOCKING_ALT_ESPRFIR, STOCKING_ALT_RWJPIN, STOCKING_ALT_ERCEDAR))) %>%
    ungroup()
  
  fortypgrp$FORTYPGRP <- ifelse(fortypgrp$STOCKING_ALT_HRDWD > fortypgrp$STOCKING_ALT_SFTWD & fortypgrp$STOCKING_ALT_OAKPIN >= (0.25 * fortypgrp$STOCKING_ALT_TOTAL), 400,
                                ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_EXHWDS, 990, 
                                       ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_MBB, 800, 
                                              ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_ASPBRCH, 900,
                                                     ifelse(fortypgrp$STOCKING_MAX ==  fortypgrp$STOCKING_ALT_ELMASHCW, 700,
                                                            ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_ALDRMAPL, 910,
                                                                   ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_OKGMCYP, 600,
                                                                          ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_OAKHCK, 500,
                                                                                 ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_EXSFTWD, 380,
                                                                                        ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_LOBSHRTP,160,
                                                                                               ifelse(fortypgrp$STOCKING_MAX ==  fortypgrp$STOCKING_ALT_LNGLFSLH, 140,
                                                                                                      ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_ESPRFIR, 120,
                                                                                                             ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_ERCEDAR, 170,
                                                                                                                    ifelse(fortypgrp$STOCKING_MAX == fortypgrp$STOCKING_ALT_RWJPIN, 100,
                                                                                                                           999))))))))))))))
  
  print(fortypgrp$FORTYPGRP)
  
  return(fortypgrp)
}