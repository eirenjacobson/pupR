#This script provides a complete pup production
#estimation of the 2012 Harp and Hood survey
#carried out in the West Ice.

library(mgcv)
library(ggplot2)

#Source functions needed
source("R/functions.r")
#source("ADMButils.R")

#Load the data and clean it up
returnList = loadAndPrepareData(survey = "WestIce2012",
                                fname = "WestIce2012.csv",
                                population = c("harp","hood"))


#Correct pup counts from photos for readers error
correctedList = correctReading(dataList = returnList)

#Estimate number of pups - without reader correction
EstHarp = kingsley(returnList$xycord$x,
                   returnList$data$HarpFinalCounts,
                   returnList$data$Area,
                   returnList$data$Transect)

EstHooded = kingsley(returnList$xycord$x,
                     returnList$data$HoodedFinalCounts,
                     returnList$data$Area,
                     returnList$data$Transect)

#Estimate number of pups - with reader correction - no area correction
EstHarpCorr = kingsley(correctedList$xycord$x,
                       correctedList$data$CountsHarpCorr,
                       correctedList$data$Area,
                       correctedList$data$Transect)

EstHoodedCorr = kingsley(correctedList$xycord$x,
                         correctedList$data$CountsHoodedCorr,
                         correctedList$data$Area,
                         correctedList$data$Transect)

#Total uncertainty estimate of pup estimation + uncertainty of reader errors
VarHarpCorr = EstHarpCorr$SE^2+correctedList$VmeasHarp
VarHoodedCorr = EstHoodedCorr$SE^2+correctedList$VmeasHooded

#Correction factors for pups not born and pups that has left the ice prior photographing
#These estimates are obtained using the BirthDistR.R code
#Look at function for this
QHarp = 0.99
SEQHarp = 0.00014
QHooded = 0.82
SEQHooded = 0.0134

estimates = list()
estimates$EstHarpFinal = EstHarp$N/QHarp
estimates$EstHarpFinalCorr = EstHarpCorr$N/QHarp
estimates$EstHoodedFinal = EstHooded$N/QHooded
estimates$EstHoodedFinalCorr = EstHoodedCorr$N/QHooded

estimates$VarTotalHarp = (1/QHarp)^2*EstHarp$SE^2+(EstHarp$N/QHarp^2)*SEQHarp^2
estimates$VarTotalHooded = (1/QHooded)^2*EstHooded$SE^2+(EstHooded$N/QHooded^2)*SEQHooded^2
estimates$VarTotalHarpCorr = (1/QHarp)^2*VarHarpCorr+(EstHarpCorr$N/QHarp^2)*SEQHarp^2
estimates$VarTotalHoodedCorr = (1/QHooded)^2*VarHoodedCorr+(EstHoodedCorr$N/QHooded^2)*SEQHooded^2

#Print pup production estimates based on the Kingsley method to screen
printKingsley2screen(estimates)

############################################################
#GAM based pup production estimates
############################################################
#GAM estimates without readers correction
EstGAM = GAMestimate(harpcounts = returnList$data$HarpFinalCounts,
                     hoodedcounts = returnList$data$HoodedFinalCounts,
                     area = returnList$data$Area,
                     xycord = returnList$xycord,
                     transect = returnList$data$Transect,
                     distr = "negbin")

#GAM estimates with readers correction
EstGAMCorr = GAMestimate(harpcounts = correctedList$data$CountsHarpCorr,
                         hoodedcounts = correctedList$data$CountsHoodedCorr,
                         area = correctedList$data$Area,
                         xycord = correctedList$xycord,
                         transect = correctedList$data$Transect,
                         distr = "negbin")

#Total uncertainty estimate of pup estimation + uncertainty of reader errors
VarHarpGAMCorr = EstGAMCorr$SEHarpGAM^2+correctedList$VmeasHarp
VarHoodedGAMCorr = EstGAMCorr$SEHoodedGAM^2+correctedList$VmeasHooded

estimatesGAM = list()
estimatesGAM$EstHarpGAMFinal = EstGAM$NHarpGAM/QHarp
estimatesGAM$EstHarpGAMFinalCorr = EstGAMCorr$NHarpGAM/QHarp
estimatesGAM$EstHoodedGAMFinal = EstGAM$NHoodedGAM/QHooded
estimatesGAM$EstHoodedGAMFinalCorr = EstGAMCorr$NHoodedGAM/QHooded

estimatesGAM$VarTotalHarpGAM = (1/QHarp)^2*EstGAM$VarHarpGAM+(EstGAM$NHarpGAM/QHarp^2)*SEQHarp^2
estimatesGAM$VarTotalHarpGAMCorr = (1/QHarp)^2*VarHarpGAMCorr+(EstGAMCorr$NHarpGAM/QHarp^2)*SEQHarp^2
estimatesGAM$VarTotalHoodedGAM = (1/QHooded)^2*EstGAM$VarHoodedGAM+(EstGAM$NHoodedGAM/QHooded^2)*SEQHooded^2
estimatesGAM$VarTotalHoodedGAMCorr = (1/QHooded)^2*VarHoodedGAMCorr+(EstGAMCorr$NHoodedGAM/QHooded^2)*SEQHooded^2

#Print GAM based pup production estimates on screen
printGAM2screen(estimatesGAM)

