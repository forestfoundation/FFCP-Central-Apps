#############################
# 
#     Goodness of Fit 
#
# This code is meant to evaluate the goodness
# of fit (GOF) of the matching control points
# selected for the dynamic baseline as prescribed
# in VM0045 under the VCS Standard. The code evaluates
# the data outputs from the "FFCP_Matching_Code".

# This code has been updated to follow v1.2 of the methodology (2025-10-16)
#############################

## Goodness of Fit for matches is determined at the project scale, 
## as the Standardized Mean Difference between treatment and control's covariates.
## SDM = Absolute Value(mean with project - mean baseline) / standard deviation with project, for each covariate
##**NOTE: This will not work for physical distance, as all distances for treatment plots are necessarily 0.

# For the first monitoring period of the FFCP Central Apps project, the goodness of 
# fit is evaluated for the two cohorts (2020 and 2021) in aggregate.

# load libraries
rm(list = (ls())) #use to remove all objects except FIA data
gc() # reset memory
library(tidyverse) # data organization
library(readxl) # sub-package of tidyverse
library(cobalt) # evaluates matching balance

# Set output directory to pull output data from
output.folder <- "OUTPUT_2024-06-19/"

# Load and combine cohorts for GOF analysis
comp.controls_1 = read.csv(paste0(output.folder, "Cohort2020_Composite Controls.csv"))
standlist_1 = read.csv(paste0(output.folder, "Cohort2020_standlist.csv"))
standlist_1 = standlist_1[,-1]
comp.controls_2 = read.csv(paste0(output.folder, "Cohort2021_Composite Controls.csv"))
standlist_2 = read.csv(paste0(output.folder, "Cohort2021_standlist.csv"))
standlist_2 = standlist_2[,-1]
control_1 = read.csv(paste0(output.folder, "Cohort2020_Control_plot summary_at enrollment.csv"), header = TRUE, sep = ",")
control_2 = read.csv(paste0(output.folder, "Cohort2021_Control_plot summary_at enrollment.csv"), header = TRUE, sep = ",")
matches_1 = read.csv(paste0(output.folder, "Cohort2020_Treatment and Matched Controls.csv"))
matches_2 = read.csv(paste0(output.folder, "Cohort2021_Treatment and Matched Controls.csv"))
comp.controls = rbind(comp.controls_1,comp.controls_2)
standlist = rbind(standlist_1, standlist_2)
control = rbind(control_1, control_2)
matches = rbind(matches_1, matches_2)

# Equation A2
## Calculate SDMs as part of a vector
SDMs = c(abs(mean(standlist$SLOPE) - mean(as.numeric(comp.controls$SLOPE))) / sqrt((sd(standlist$SLOPE)^2+sd(comp.controls$SLOPE)^2)/2),
         abs(mean(standlist$QMD) - mean(as.numeric(comp.controls$QMD))) / sqrt((sd(standlist$QMD)^2+sd(comp.controls$QMD)^2)/2),
         abs(mean(standlist$RDDISTCD) - mean(as.numeric(comp.controls$RDDISTCD))) / sqrt((sd(standlist$RDDISTCD)^2+sd(comp.controls$RDDISTCD)^2)/2),
         abs(mean(standlist$RD.comm) - mean(as.numeric(comp.controls$RD.comm))) / sqrt((sd(standlist$RD.comm)^2+sd(comp.controls$RD.comm)^2)/2),
         abs(mean(standlist$RD.sap) - mean(as.numeric(comp.controls$RD.sap))) / sqrt((sd(standlist$RD.sap)^2+sd(comp.controls$RD.sap)^2)/2),
         abs(mean(standlist$STDAGE) - mean(as.numeric(comp.controls$STDAGE))) / sqrt((var(standlist$STDAGE)+var(comp.controls$STDAGE))/2),
         abs(mean(standlist$SITECLCD) - mean(as.numeric(comp.controls$SITECLCD))) / sqrt((sd(standlist$SITECLCD)^2+sd(comp.controls$SITECLCD)^2)/2))

## Create Goodness of fit dataframe, to view and output
GOF = data.frame(Covariate = c("SLOPE", "QMD", "RDDISTCD", "RD.comm", "RD.sap", "STDAGE", "SITECLCD"),
                 SDM = SDMs)

GOF





# Evaluate covariate balance using Cobalt package
# Create table of all donor plots with treatments and m.dists for evaluations
all.plots <- rbind(standlist, control)
all.plots <- all.plots %>%
  left_join(matches %>% dplyr::select("PLT_CN", "Mahalanobis.Distance", "weight"), by = "PLT_CN")

# Input 1 for weight of all treatment plots
all.plots$weight <- ifelse(all.plots$Z==1 ,1, all.plots$weight)

# Input 0 for weight of unselected donor plots
all.plots$weight <- ifelse(is.na(all.plots$weight), 0, all.plots$weight)

# Input 0 for mahalanobis distance to treatment plots
all.plots$Mahalanobis.Distance <- ifelse(all.plots$Z==1, 0, all.plots$Mahalanobis.Distance)

# Input 10 for m.dist for unselected donor plots
all.plots$Mahalanobis.Distance <- ifelse(is.na(all.plots$Mahalanobis.Distance), 10, all.plots$Mahalanobis.Distance)

# Create data frame of covariates
covs <- all.plots[,9:15]

# Run summary balanace stats in "cobalt" package
balance <-bal.tab(covs, 
                  treat=all.plots$Z, 
                  weights=all.plots$weight, 
                  stats="mean.diffs", 
                  disp=c("means","sds"),
                  s.d.denom="treated",
                  thresholds=c(m=0.25),
                  un=TRUE)
# Create covariate balance plot showing adjusted vs unadjusted variables
#png(filename = paste0(output.folder, "Figures/","Cohort2020_2021","_GoodnessOfFit_", Sys.Date(),".png"),
 #   width = 700, height = 400)
love.plot(balance,
          abs=TRUE)
dev.off()

write.csv(balance[1], paste0(output.folder, "combined_cohort2020&cohort2021_", "Covariate_balance.csv"))

#Combine cohort 2020 and 2021 for GOF analysis
#treat.comp.controls_1 = read.csv(paste0(output.folder, "Cohort2020_Treatment and Composite Controls.csv"))
#treat.comp.controls = rbind(treat.comp.controls,treat.comp.controls_1)
# Plot stand covariates against composite covariates to assess goodness of fit
png(filename = paste0(output.folder, "Figures/","Cohort2020_2021","_MatchedCovariateScatterPlots_", Sys.Date(),".png"),
    width = 800, height = 500)
par(mfrow=c(2, 4))
plot(treat.comp.controls$SLOPE[which(treat.comp.controls$Z==1)], treat.comp.controls$SLOPE[which(is.na(treat.comp.controls$Z))],
     main="SLOPE",
     xlab="Treatment",
     ylab="Control")
abline(0,1)
plot(treat.comp.controls$QMD[which(treat.comp.controls$Z==1)], treat.comp.controls$QMD[which(is.na(treat.comp.controls$Z))],
     main="QMD",
     xlab="Treatment",
     ylab="Control")
abline(0,1)
plot(treat.comp.controls$SITECLCD[which(treat.comp.controls$Z==1)], treat.comp.controls$SITECLCD[which(is.na(treat.comp.controls$Z))],
     main="SITE CLASS CODE",
     xlab="Treatment",
     ylab="Control")
abline(0,1)
plot(treat.comp.controls$RDDISTCD[which(treat.comp.controls$Z==1)], treat.comp.controls$RDDISTCD[which(is.na(treat.comp.controls$Z))],
     main="ROAD DISTANCE",
     xlab="Treatment",
     ylab="Control")
abline(0,1)
plot(treat.comp.controls$STDAGE[which(treat.comp.controls$Z==1)], treat.comp.controls$STDAGE[which(is.na(treat.comp.controls$Z))],
     main="STAND AGE",
     xlab="Treatment",
     ylab="Control")
abline(0,1)
plot(treat.comp.controls$RD.comm[which(treat.comp.controls$Z==1)], treat.comp.controls$RD.comm[which(is.na(treat.comp.controls$Z))],
     main="RD COMMERCIAL",
     xlab="Treatment",
     ylab="Control")
abline(0,1)
plot(treat.comp.controls$RD.sap[which(treat.comp.controls$Z==1)], treat.comp.controls$RD.sap[which(is.na(treat.comp.controls$Z))],
     main="RD SAPLING",
     xlab="Treatment",
     ylab="Control")
abline(0,1)
dev.off()


## Generate warning if goodness of fit metrics indicate a poor quality match
if(max(SDMs) >= 0.25) {warning("Goodness of Fit statistics indicate the data did not achieve a valid match. Of the Standardized Difference of Means calculated, ", length(SDMs[SDMs >= 0.25]), " were greater than the 0.25 threshold required by the Verra Methodology for Improved Forest Management. Redefine the treat.comp.controls and matches datasets in the Prep Matching Data section of the code, then repeat the donor pool selection and matching steps above with progressively smaller k values until a valid overall match is achieved.")} else{message("The project achieved a valid match. Congratulations!")}

## Output goodness of fit statistics
write.csv(GOF, paste0(output.folder, "_Goodness of Fit Statistics.csv"), row.names = FALSE)

