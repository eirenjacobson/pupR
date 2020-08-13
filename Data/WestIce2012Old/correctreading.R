
#1: both readers together, 2: separate correction for each reader, 3: both readers together with intercept
correction_type = 2 
#Only for option 1, runs MC simulations to get empirical estimate of reading error
run_mc_sim = FALSE

if (correction_type == 1){

#Both readers pooled
indcorrHarp = which(1-is.na(data$WPTrue)==1)			#Find which photos have been checked by both readers
lincorrHarp = lm(data$WPTrue[indcorrHarp]~data$HarpFinalCounts[indcorrHarp] - 1)
plot(data$HarpFinalCounts[indcorrHarp],data$WPTrue[indcorrHarp])
indcorrHooded = which(1-is.na(data$BBTrue)==1)			#Find which photos have been checked by both readers
lincorrHooded = lm(data$BBTrue[indcorrHooded]~data$HoodedFinalCounts[indcorrHooded] - 1)
plot(data$HoodedFinalCounts[indcorrHooded],data$BBTrue[indcorrHooded])

var_bHarp = coef(summary(lincorrHarp))[1,2]; var_bHarp = var_bHarp^2
var_bHooded = coef(summary(lincorrHooded))[1,2]; var_bHooded = var_bHooded^2
var_uHarp = sd(lincorrHarp$residuals); var_uHarp = var_uHarp^2
var_uHooded = sd(lincorrHooded$residuals); var_uHooded = var_uHooded^2

data$CountsHarpCorr <- lincorrHarp$coefficients*data$HarpFinalCounts
data$CountsHoodedCorr <- lincorrHooded$coefficients*data$HoodedFinalCounts

#Get uncertainty associatet with reader errors
VmeasHarp = vmeasest(xycord$x,data$HarpFinalCounts,data$Area,data$Transect,Spacing,gap,var_bHarp,var_uHarp)
VmeasHooded = vmeasest(xycord$x,data$HoodedFinalCounts,data$Area,data$Transect,Spacing,gap,var_bHooded,var_uHooded)

if (run_mc_sim == TRUE){
#Compare this uncertainty with simulation based uncertainty.
Msim = 10000
EstimatesHarp = array(NA,Msim)
EstimatesHooded = array(NA,Msim)
for (i in 1:Msim){
	data$simcountsHarp = rnorm(1,lincorrHarp$coefficients,sqrt(var_bHarp))*data$HarpFinalCounts+1*rnorm(1,0,sqrt(var_uHarp))
	EstHarpSim = kingsley(xycord$x,data$simcountsHarp,data$Area,data$Transect,Spacing,gap)
    EstimatesHarp[i] = EstHarpSim$N
    data$simcountsHooded = rnorm(1,lincorrHooded$coefficients,var_bHooded)*data$HoodedFinalCounts+1*rnorm(1,0,var_uHooded)
	EstHoodedSim = kingsley(xycord$x,data$simcountsHooded,data$Area,data$Transect,Spacing,gap)
    EstimatesHooded[i] = EstHoodedSim$N
          }
     }
}


if (correction_type == 2){

##############################
#Separate fit for each reader

#Reader 1
indcorrHarpLL = which((1-is.na(data$WPTrue)==1) & data$Reader == "LL")			#Find which photos have been checked by LL
lincorrHarpLL = lm(WPTrue[indcorrHarpLL]~HarpFinalCounts[indcorrHarpLL] - 1,data = data)
plot(data$HarpFinalCounts[indcorrHarpLL],data$WPTrue[indcorrHarpLL])
abline(a = 0, b = lincorrHarpLL$coefficients,lwd = 3,col = "red")

indcorrHoodedLL = which(1-is.na(data$BBTrue)==1 & data$Reader == "LL")			#Find which photos have been checked by LL
lincorrHoodedLL = lm(data$BBTrue[indcorrHoodedLL]~data$HoodedFinalCounts[indcorrHoodedLL] - 1)
plot(data$HoodedFinalCounts[indcorrHoodedLL],data$BBTrue[indcorrHoodedLL])
abline(a = 0, b = lincorrHoodedLL$coefficients,lwd = 3,col = "red")


var_bHarpLL = coef(summary(lincorrHarpLL))[1,2]; var_bHarpLL = var_bHarpLL^2
var_bHoodedLL = coef(summary(lincorrHoodedLL))[1,2]; var_bHoodedLL = var_bHoodedLL^2
var_uHarpLL = sd(lincorrHarpLL$residuals); var_uHarpLL = var_uHarpLL^2
var_uHoodedLL = sd(lincorrHoodedLL$residuals); var_uHoodedLL = var_uHoodedLL^2

#Korriger bare de som Lotta har lest
indLL = which(data$Reader == "LL")
data$CountsHarpCorr[indLL] <- lincorrHarpLL$coefficients*data$HarpFinalCounts[indLL]
data$CountsHoodedCorr[indLL] <- lincorrHoodedLL$coefficients*data$HoodedFinalCounts[indLL]

#Get uncertainty associatet with reader errors
VmeasHarpLL = vmeasest(xycord$x[indLL],data$HarpFinalCounts[indLL],data$Area[indLL],data$Transect[indLL],Spacing,gap,var_bHarpLL,var_uHarpLL)
VmeasHoodedLL = vmeasest(xycord$x[indLL],data$HoodedFinalCounts[indLL],data$Area[indLL],data$Transect[indLL],Spacing,gap,var_bHoodedLL,var_uHoodedLL)



#Reader 2
indcorrHarpMP = which(1-is.na(data$WPTrue)==1 & data$Reader == "MP")			#Find which photos have been checked by MP
lincorrHarpMP = lm(data$WPTrue[indcorrHarpMP]~data$HarpFinalCounts[indcorrHarpMP] - 1)
plot(data$HarpFinalCounts[indcorrHarpMP],data$WPTrue[indcorrHarpMP])
abline(a = 0, b = lincorrHarpMP$coefficients,lwd = 3,col = "red")

indcorrHoodedMP = which(1-is.na(data$BBTrue)==1 & data$Reader == "MP")			#Find which photos have been checked by MP
lincorrHoodedMP = lm(data$BBTrue[indcorrHoodedMP]~data$HoodedFinalCounts[indcorrHoodedMP] - 1)
#plot(data$HoodedFinalCounts[indcorrHoodedMP],data$BBTrue[indcorrHoodedMP])
#abline(a = 0, b = lincorrHoodedMP$coefficients,lwd = 3,col = "red")

var_bHarpMP = coef(summary(lincorrHarpMP))[1,2]; var_bHarpMP = var_bHarpMP^2
var_bHoodedMP = coef(summary(lincorrHoodedMP))[1,2]; var_bHoodedMP = var_bHoodedMP^2
var_uHarpMP = sd(lincorrHarpMP$residuals); var_uHarpMP = var_uHarpMP^2
var_uHoodedMP = sd(lincorrHoodedMP$residuals); var_uHoodedMP = var_uHoodedMP^2

#Korriger bare de som Michael har lest
indMP = which(data$Reader == "MP")
data$CountsHarpCorr[indMP] <- lincorrHarpMP$coefficients*data$HarpFinalCounts[indMP]
data$CountsHoodedCorr [indMP]<- lincorrHoodedMP$coefficients*data$HoodedFinalCounts[indMP]

#Get uncertainty associatet with reader errors
VmeasHarpMP = vmeasest(xycord$x[indMP],data$HarpFinalCounts[indMP],data$Area[indMP],data$Transect[indMP],Spacing,gap,var_bHarpMP,var_uHarpMP)
VmeasHoodedMP = vmeasest(xycord$x[indMP],data$HoodedFinalCounts[indMP],data$Area[indMP],data$Transect[indMP],Spacing,gap,var_bHoodedMP,var_uHoodedMP)

VmeasHarp = VmeasHarpLL + VmeasHarpMP

VmeasHooded = VmeasHoodedLL + VmeasHoodedMP
}

if (correction_type == 3){

#Both readers pooled
indcorrHarp = which(1-is.na(data$WPTrue)==1)			#Find which photos have been checked by both readers
lincorrHarp = lm(data$WPTrue[indcorrHarp]~data$HarpFinalCounts[indcorrHarp])
plot(data$HarpFinalCounts[indcorrHarp],data$WPTrue[indcorrHarp])
indcorrHooded = which(1-is.na(data$BBTrue)==1)			#Find which photos have been checked by both readers
lincorrHooded = lm(data$BBTrue[indcorrHooded]~data$HoodedFinalCounts[indcorrHooded])
plot(data$HoodedFinalCounts[indcorrHooded],data$BBTrue[indcorrHooded])

var_aHarp = coef(summary(lincorrHarp))[1,2]; var_aHarp = var_aHarp^2
var_aHooded = coef(summary(lincorrHooded))[1,2]; var_aHooded = var_aHooded^2
var_bHarp = coef(summary(lincorrHarp))[2,2]; var_bHarp = var_bHarp^2
var_bHooded = coef(summary(lincorrHooded))[2,2]; var_bHooded = var_bHooded^2
var_uHarp = sd(lincorrHarp$residuals); var_uHarp = var_uHarp^2
var_uHooded = sd(lincorrHooded$residuals); var_uHooded = var_uHooded^2

data$CountsHarpCorr <- lincorrHarp$coefficients*data$HarpFinalCounts
data$CountsHoodedCorr <- lincorrHooded$coefficients*data$HoodedFinalCounts

#Get uncertainty associatet with reader errors
VmeasHarp = vmeasest(xycord$x,data$HarpFinalCounts,data$Area,data$Transect,Spacing,gap,var_bHarp,var_uHarp,var_aHarp,cov_ab = 0,intcpt = TRUE)
VmeasHooded = vmeasest(xycord$x,data$HoodedFinalCounts,data$Area,data$Transect,Spacing,gap,var_bHooded,var_uHooded,var_aHooded,cov_ab = 0,intcpt = TRUE)

}