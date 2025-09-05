#NVEL Script
##Aaron Holley - 06 November 2023
## Ben Rifkin EDITED - 14 February 2024

NVEL <- function(treelist, t){

  nm <- deparse(substitute(treelist)) # Pull name of object input
  
if(t==3){

# Prep tree list
  
  # Subset dead trees
  treelist <- treelist[which(treelist$STATUSCD==2),]
  # Subset relevant columns
  treelist <- subset(treelist, select=c("mt","TREE_ID","TREE_INDEX", "PLT_CN","SPCD","DIA","TOTAL_HT","REGION","FORESTNUMB","DISTRICTNU"))
  #Set parameters in treelist table
  colnames(treelist)[which(colnames(treelist)=="REGION")] <- "regn"
  treelist$regn <- as.numeric(treelist$regn)
  colnames(treelist)[which(colnames(treelist)=="FORESTNUMB")] <- "forst"
  colnames(treelist)[which(colnames(treelist)=="DISTRICTNU")] <- "dist"
  colnames(treelist)[which(colnames(treelist)=="DIA")] <- "dbhob"
  colnames(treelist)[which(colnames(treelist)=="TOTAL_HT")] <- "httot"
  treelist$httot <- as.numeric(treelist$httot)
  colnames(treelist)[which(colnames(treelist)=="SPCD")] <- "spec"
  treelist$spec <- as.numeric(treelist$spec)
} else {
  # Prep tree list
  if(t==1) { # choose parameters for current measurement
  # Subset dead trees
  treelist <- treelist[which(treelist$STATUSCD==2),]
  
  # Subset relevant columns
  treelist <- subset(treelist, select=c("PLT_CN","CN","TREE_ID","treat_PLOTID","SPCD","DIA", "VOLCFNET"))
  #Set parameters in treelist table
  colnames(treelist)[which(colnames(treelist)=="DIA")] <- "dbhob"
  colnames(treelist)[which(colnames(treelist)=="SPCD")] <- "spec"
  treelist$spec <- as.numeric(treelist$spec)
  colnames(treelist)[which(colnames(treelist)=="VOLCFNET")] <- "TotCuFtVol_FIA"
  } else if(t==0) { #choose parameters for previous measurement
    # Subset dead trees
    treelist <- treelist[which(treelist$PREVSTATUSCD==2),]
    
    # Subset relevant columns
    treelist <- subset(treelist, select=c("PLT_CN","PREV_TRE_CN","treat_PLOTID","SPCD","PREVDIA", "PREVVOLCFNET"))
    #Set parameters in treelist table
    colnames(treelist)[which(colnames(treelist)=="PREVDIA")] <- "dbhob"
    colnames(treelist)[which(colnames(treelist)=="SPCD")] <- "spec"
    treelist$spec <- as.numeric(treelist$spec)
    colnames(treelist)[which(colnames(treelist)=="PREVVOLCFNET")] <- "TotCuFtVol_FIA"
  }
}

if(nm=="treelist.treat") {
  #Calculate Stem Volume using NVEL
  #Calling the NVEL .dll. Must be set to location of latest NVEL .dll on local machine
  dyn.load("INPUT/VolLibDll20230818/VolLibDll20230818/vollib-64bits/vollib.dll")
  
  #Check .dll version number. Make sure this is latest downloaded from USFS.
  .Fortran("vernum_r",vernum=integer(1))
  
  #Assigning variables
  
  #Calculate CuFt. Vol.
  i=1
  while (i <=nrow(treelist)) {
    regn<-treelist$regn[i]
    forst<-treelist$forst[i]
    dist<-treelist$dist[i]
    spec<-treelist$spec[i]
    voleq<-""
    errflg<-0
    dbhob<-treelist$dbhob[i]
    httot<-treelist$httot[i]
    mtopp<-0
    mtops<-0
    ht1prd<-0
    ht2prd<-0
    upsht1<-0
    upsd1<-0
    stump<-1
    fclass<-0
    dbtbh<-0
    btr<-0
    vol<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    fiavol<-0
    voltype<-""
    
    #Obtain VOLEQ for a species/site combination
    VOLEQ<-as.character(.Fortran("getvoleq_r",as.integer(regn),as.character(forst),as.character(dist),as.integer(spec),as.character(voleq),as.integer(errflg))[5])
    treelist$VOLEQ[i]<-VOLEQ
    
    #Implement FIA Volume Equation
    FIAVOL_Calced<-.Fortran("vollib_r",as.character(VOLEQ),as.integer(regn),as.character(forst),as.character(dist),as.integer(spec),as.double(dbhob),as.double(httot),as.double(mtopp),as.double(mtops),as.double(ht1prd),as.double(ht2prd),as.double(upsht1),as.double(upsd1),as.double(stump),as.integer(fclass),as.double(dbtbh),as.double(btr),as.double(vol),as.integer(errflg))[18]
  
    TotCuFtVol_FIA<-FIAVOL_Calced[[1]]
    treelist$TotCuFtVol_FIA[i]<-TotCuFtVol_FIA[1]
    
    i=i+1
  }
  
}
  
#CRM
#Following ARB/FIA method https://ww2.arb.ca.gov/sites/default/files/cap-and-trade/protocols/usforest/2014/biomass_estimation_component_ratio_method.pdf
  #Import REF_SPECIES and merge with treelist
  REF_SPECIES<-read.csv("INPUT/REF_SPECIES.csv")
  treelist<-merge(treelist,REF_SPECIES, by.x="spec", by.y="SPCD")
  
  #Specific Gravity
  treelist$WdSpGr<-treelist$WOOD_SPGR_GREENVOL_DRYWT*62.43
  treelist$BkSpGr<-treelist$BARK_SPGR_GREENVOL_DRYWT*62.43
  
  #Jenkins
    #Total AGB
    treelist$AGB_lbs_Jenkins_Total<-exp(treelist$JENKINS_TOTAL_B1+treelist$JENKINS_TOTAL_B2*log(treelist$dbhob*2.54)) *2.2046
    #stem_ratio
    treelist$stem_ratio<-exp(treelist$JENKINS_STEM_WOOD_RATIO_B1 + treelist$JENKINS_STEM_WOOD_RATIO_B2 / (treelist$dbhob*2.54))
    #bark_ratio
    treelist$bark_ratio<-exp(treelist$JENKINS_STEM_BARK_RATIO_B1 + treelist$JENKINS_STEM_BARK_RATIO_B2 / (treelist$dbhob*2.54) )
    #foliage_ratio
    treelist$foliage_ratio<-exp(treelist$JENKINS_FOLIAGE_RATIO_B1 + treelist$JENKINS_FOLIAGE_RATIO_B2 / (treelist$dbhob*2.54) ) 
    #root_ratio
    treelist$root_ratio<- exp(treelist$JENKINS_ROOT_RATIO_B1 + treelist$JENKINS_ROOT_RATIO_B2 / (treelist$dbhob*2.54) ) 
    #stem_biomass_jenkins
    treelist$stem_biomass_jenkins_lbs<-treelist$AGB_lbs_Jenkins_Total * treelist$stem_ratio 
    #bark_biomass_jenkins
    treelist$bark_biomass_jenkins_lbs<-treelist$AGB_lbs_Jenkins_Total * treelist$bark_ratio 
    #bole_biomass_jenkins
    treelist$bole_biomass_jenkins_lbs<-treelist$stem_biomass_jenkins_lbs + treelist$bark_biomass_jenkins_lbs 
    #foliage_biomass_jenkins
    treelist$foliage_biomass_jenkins_lbs<-treelist$AGB_lbs_Jenkins_Total * treelist$foliage_ratio 
    #root_biomass_jenkins
    treelist$root_biomass_jenkins_lbs<-treelist$AGB_lbs_Jenkins_Total * treelist$root_ratio 
    #stump_biomass
      #Inside Bark
      treelist$stump_dib0<-(treelist$dbhob*treelist$RAILE_STUMP_DIB_B1)+(treelist$dbhob*treelist$RAILE_STUMP_DIB_B2*(4.5-0)/(1+0)) #DIB at 0'
      treelist$stump_dib1<-(treelist$dbhob*treelist$RAILE_STUMP_DIB_B1)+(treelist$dbhob*treelist$RAILE_STUMP_DIB_B2*(4.5-1)/(1+1)) #DIB at 1'
      treelist$stump_dib<-(treelist$stump_dib0+treelist$stump_dib1)/2 #Average DIB
      treelist$stump_IBVol_CuFt<-(pi*(treelist$stump_dib/2)^2)/144 #Average Area IB. Since it is for a 1' stump, the sqft = cuft
      #Outside Bark
      treelist$stump_dob0<-treelist$dbhob+(treelist$dbhob*treelist$RAILE_STUMP_DOB_B1*(4.5-0)/(0+1)) #DOB at 0'
      treelist$stump_dob1<-treelist$dbhob+(treelist$dbhob*treelist$RAILE_STUMP_DOB_B1*(4.5-1)/(1+1)) #DOB at 1'
      treelist$stump_dob<-(treelist$stump_dob0+treelist$stump_dob1)/2 #Average DOB
      treelist$stump_OBVol_CuFt<-(pi*(treelist$stump_dob/2)^2)/144 #Average Area OB. Since it is for a 1' stump, the sqft = cuft
      
      #Stump Bio
      treelist$stump_biomass_bark_lbs<-(treelist$stump_OBVol_CuFt-treelist$stump_IBVol_CuFt)*treelist$BkSpGr
      treelist$stump_biomass_wood_lbs<-treelist$stump_IBVol_CuFt*treelist$WdSpGr
      treelist$stump_biomass_total_lbs<-treelist$stump_biomass_bark_lbs+treelist$stump_biomass_wood_lbs
      
      
    #top_biomass_jenkins
    treelist$top_biomass_jenkins_lbs<-treelist$AGB_lbs_Jenkins_Total-treelist$stem_biomass_jenkins_lbs-treelist$bark_biomass_jenkins_lbs-treelist$foliage_biomass_jenkins_lbs-treelist$stump_biomass_total_lbs
    
  #DryBio
    #DRYBIOBOLE using NVEL
    treelist$drybio_bole_lbs<-(treelist$TotCuFtVol_FIA*(treelist$BARK_VOL_PCT/100.0) * (treelist$BARK_SPGR_GREENVOL_DRYWT * 62.4) ) + (treelist$TotCuFtVol_FIA * (treelist$WOOD_SPGR_GREENVOL_DRYWT * 62.4) ) 
    
    #Adjustment Factor
    treelist$AdjFact<-treelist$drybio_bole_lbs/treelist$bole_biomass_jenkins_lbs
  
    #DryBio Top
    treelist$drybio_top_lbs<-treelist$top_biomass_jenkins_lbs*treelist$AdjFact
  
    #DryBio Stump
    treelist$drybio_stump_lbs<-treelist$stump_biomass_total_lbs*treelist$AdjFact
  
    #DryBio BG
    treelist$drybio_bg_lbs<-treelist$root_biomass_jenkins_lbs*treelist$AdjFact
    
#AG mTCO2e
    #Parms tables of VM0045 state that "standing dead wood is restricted here to aboveground stem (bole) biomass"
    treelist$mTCO2e<-treelist$drybio_bole_lbs/2204.6*0.5*(44/12)
    
    treelist_deadvol <- treelist
    
    return(treelist_deadvol)

}
