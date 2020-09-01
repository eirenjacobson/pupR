# This script show an example of a complete pup production
# estimation using data from the 2012 West Ice survey.


# Loading the pupR package
library(pupR)


########################
# Loading the data
########################

# If using the demo data set - loading the demo data set
data("pupsDemoData")

# Break up the elements of the list in the demo data
dataList = demoData$dataList
harpStages = demoData$HarpStages
hoodedStages = demoData$HoodedStages

# If using the complete data set - loading the data
dataList = loadAndPrepareData(survey = "WestIce2012",
                              fname = "WestIce2012.csv",
                              population = c("harp","hood"))


# Correct for readers errors
correctedList = correctReading(dataList = dataList)


###################################
# Abundance estimation
###################################

# The Kingsley method
#*********************
#
# The harp seal population
#----------------------------
#
# Using the uncorrected pup counts
EstHarp = kingsley(dataList$xycord$x,
                   dataList$data$HarpFinalCounts,
                   dataList$data$Area,
                   dataList$data$Transect)

# Using readers corrected pup counts
EstHarpCorr = kingsley(correctedList$xycord$x,
                       correctedList$data$CountsHarpCorr,
                       correctedList$data$Area,
                       correctedList$data$Transect)


# The hooded seal population
#------------------------------
#
# Using the uncorrected pup counts
EstHooded = kingsley(dataList$xycord$x,
                     dataList$data$HoodedFinalCounts,
                     dataList$data$Area,
                     dataList$data$Transect)


# Using readers corrected pup counts
EstHoodedCorr = kingsley(correctedList$xycord$x,
                         correctedList$data$CountsHoodedCorr,
                         correctedList$data$Area,
                         correctedList$data$Transect)


# Calculating the total uncertainty estimate of the abundance estimate
VarHarpCorr = EstHarpCorr$SE^2+correctedList$VmeasHarp
VarHoodedCorr = EstHoodedCorr$SE^2+correctedList$VmeasHooded

# The GAM method
#*****************
#
# Both the harp and the hooded populationS
# Using the pup counts not corrected for readers errors
EstGAM = GAMestimate(harpcounts =
                       dataList$data$HarpFinalCounts,
                     hoodedcounts =
                       dataList$data$HoodedFinalCounts,
                     area = dataList$data$Area,
                     xycord = dataList$xycord,
                     transect = dataList$data$Transect,
                     distr = "negbin")

# Using readers corrected pup counts
EstGAMCorr = GAMestimate(harpcounts =
                           correctedList$data$CountsHarpCorr,
                         hoodedcounts =correctedList$data$CountsHoodedCorr,
                         area = correctedList$data$Area,
                         xycord = correctedList$xycord,
                         transect = correctedList$data$Transect,
                         distr = "negbin")

# Calculating the total uncertainty
VarHarpGAMCorr = EstGAMCorr$SEHarpGAM^2+correctedList$VmeasHarp
VarHoodedGAMCorr = EstGAMCorr$SEHoodedGAM^2+correctedList$VmeasHooded


# Correct for seals not born yet or seals which has left the ice
#*****************************************************************

#Estimating the birth distribution for the harp seals
bdistHarp = birthDist(data = harpStages,
                      harpLengthStages = c(2.4,4.42,11.39),
                      harpKappa = 12.4,
                      population = "harp",
                      datePhoto = 28)


#Estimating the birth distribution for the hooded seals
bdistHooded = birthDist(data = hoodedStages,
                        hoodedLengthStages = c(2,1,4),
                        hoodedKappa = 8.6,
                        population = "hooded",
                        datePhoto = 28)

# Extracting the correction factor and the uncertainty
QHarp = bdistHarp$PropIce            # Correction factor for harp seals
SEQHarp = bdistHarp$PropIceSD        # Uncertainty of the harp seal correction factor
QHooded = bdistHooded$PropIce        # Correction factor for harp seals
SEQHooded = bdistHooded$PropIceSD    # Uncertainty of the harp seal correction factor

# Correcting for missing seals
estimates = list()


# Correcting the Kingsley method estimates
#------------------------------------------

# Harp seal population
# Estimated pup abundance based on raw pup counts
estimates$EstHarpFinal = EstHarp$N/QHarp

# Estimated pup abundance based on readers corrected pup counts
estimates$EstHarpFinalCorr = EstHarpCorr$N/QHarp

# Hooded seal population
# Estimated pup abundance based on raw pup counts
estimates$EstHoodedFinal = EstHooded$N/QHooded

# Estimated pup abundance based on readers corrected pup counts
estimates$EstHoodedFinalCorr = EstHoodedCorr$N/QHooded

# Calculating the total uncertainty
estimates$VarTotalHarp = (1/QHarp)^2*EstHarp$SE^2+(EstHarp$N/QHarp^2)*SEQHarp^2
estimates$VarTotalHooded = (1/QHooded)^2*EstHooded$SE^2+(EstHooded$N/QHooded^2)*SEQHooded^2
estimates$VarTotalHarpCorr = (1/QHarp)^2*VarHarpCorr+(EstHarpCorr$N/QHarp^2)*SEQHarp^2
estimates$VarTotalHoodedCorr = (1/QHooded)^2*VarHoodedCorr+(EstHoodedCorr$N/QHooded^2)*SEQHooded^2


# Correcting the GAM estimates
#-------------------------------
estimatesGAM = list()
estimatesGAM$EstHarpGAMFinal = EstGAM$NHarpGAM/QHarp
estimatesGAM$EstHarpGAMFinalCorr = EstGAMCorr$NHarpGAM/QHarp
estimatesGAM$EstHoodedGAMFinal = EstGAM$NHoodedGAM/QHooded
estimatesGAM$EstHoodedGAMFinalCorr = EstGAMCorr$NHoodedGAM/QHooded

estimatesGAM$VarTotalHarpGAM = (1/QHarp)^2*EstGAM$VarHarpGAM+(EstGAM$NHarpGAM/QHarp^2)*SEQHarp^2
estimatesGAM$VarTotalHarpGAMCorr = (1/QHarp)^2*VarHarpGAMCorr+(EstGAMCorr$NHarpGAM/QHarp^2)*SEQHarp^2
estimatesGAM$VarTotalHoodedGAM = (1/QHooded)^2*EstGAM$VarHoodedGAM+(EstGAM$NHoodedGAM/QHooded^2)*SEQHooded^2
estimatesGAM$VarTotalHoodedGAMCorr = (1/QHooded)^2*VarHoodedGAMCorr+(EstGAMCorr$NHoodedGAM/QHooded^2)*SEQHooded^2


# Print the results to screen
#******************************
#
# Pup production estimates based on the Kingsley method:
printKingsley2screen(estimates)

# Pup production estimates based on the GAM method
printGAM2screen(estimatesGAM)

