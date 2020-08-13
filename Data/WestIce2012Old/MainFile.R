library(mgcv)

#Source functions needed
source("functions.r")
source("ADMButils.R")



Spacing = 3 		#Distance between transects in NM, needed for Kingsley method
gap = 1					#Maximum distance in NM between allowed, for removing holes where camera has been stopped due to
						    #flying over open water


#Complete readings with readers error
data <- read.table("WestIce2012_MichaelAndLotta_final.csv",header = TRUE,sep = ",")

#Remove some photos marked in the comment column
indRemove = which(data$Comments == "Remove")
data <- data[-indRemove,]

dimframe <- dim(data)
Nphotos <- dimframe[1]


#Sjekk bilder som er lest to ganger og bruk tall fra 2. lesning
ind2countHarp = which(data$WP2Reading != "NA")
ind2countHooded = which(data$BB2Reading != "NA")

data$HarpFinalCounts <- data$WP1Reading
data$HarpFinalCounts[ind2countHarp] <- data$WP2Reading[ind2countHarp]

#Set empty photos to zero
indNAHarp = which(is.na(data$HarpFinalCounts))
data$HarpFinalCounts[indNAHarp] = 0

#Same as above for hooded seals, i.e., use 2nd count and set empty photos to zero
data$HoodedFinalCounts <- data$BB1Reading
data$HoodedFinalCounts[ind2countHooded] <- data$BB2Reading[ind2countHooded]
indNAHooded = which(is.na(data$HoodedFinalCounts))
data$HoodedFinalCounts[indNAHooded] = 0

# #For selection of photos to use for readers error
# indLL = which((data$WP1Reading >0) & (data$BB1Reading > 0) & (data$Reader == "LL"))
# indMP = which((data$WP1Reading >0) & (data$BB1Reading > 0) & (data$Reader == "MP"))

# selectionLL = indMP[sample.int(length(indMP),50)]
# selectionLL = sort(selectionLL)

# plot(data$HarpFinalCounts[selectionLL])
# plot(data$HoodedFinalCounts[selectionLL])

# selectionMP = indLL[sample.int(length(indLL),50)]
# selectionMP = sort(selectionMP)

# X11("")
# par(mfrow = c(1,2))
# plot(data$HarpFinalCounts[selectionMP])
# plot(data$HoodedFinalCounts[selectionMP])


# #data$PhotoNo[selection]


#Find the number of zero photos
indzeroharp = which(data$HarpFinalCounts == 0)
propzeroharp = length(indzeroharp)/length(data$HarpFinalCounts)
indzerohooded = which(data$HoodedFinalCounts == 0)
propzerohooded = length(indzerohooded)/length(data$HoodedFinalCounts)
cat(paste("Proportion of zero counts for harp seals:",propzeroharp),"\n")
cat(paste("Proportion of zero counts for hooded seals:",propzerohooded),"\n")

#Convert to cartesian coordinates
xycord = lb2xy(data$FrameLon,data$FrameLat,mean(data$FrameLon),mean(data$FrameLat))

#Find the number of photos in each transect
uTrans = sort(unique(data$Transect))
nphotos = list()
ncountsHarp = list()
totalHarps = 0
ncountsHooded = list()
totalHooded = 0
meanypos = array(NA,27)
for (i in 1:length(uTrans)){
	nphotos[[i]] = length(which(data$Transect == uTrans[i]))
	ncountsHarp[[i]] = sum(data$HarpFinalCounts[data$Transect == uTrans[i]])
	totalHarps = totalHarps + ncountsHarp[[i]]
	ncountsHooded[[i]] = sum(data$HoodedFinalCounts[data$Transect == uTrans[i]])
	totalHooded = totalHooded + ncountsHooded[[i]]
	meanypos[i] = mean(xycord$y[data$Transect == uTrans[i]])
}
mean(abs(diff(meanypos)))

# #Mean height and range of flying
# indhigh = which(data$Altitude > 320)
# indlow = which(data$Altitude < 320)
# mean(data$Altitude[indhigh])
# min(data$Altitude[indhigh])
# max(data$Altitude[indhigh])

# mean(data$Altitude[indlow])
# min(data$Altitude[indlow])
# max(data$Altitude[indlow])

#Find area covered of each photo
data$Area <- (0.06786*data$Altitude/0.1005)*(0.10386*data$Altitude/0.1005)

# mean(data$Area[indhigh])
# min(data$Area[indhigh])
# max(data$Area[indhigh])

# mean(data$Area[indlow])
# min(data$Area[indlow])
# max(data$Area[indlow])

#Find length covered of each photo (in Nautical Miles)
data$length <- (0.06786*data$Altitude/0.1005) / 1852
# mean(data$length[indhigh])*1852
# min(data$length[indhigh])*1852
# max(data$length[indhigh])*1852

# mean(data$length[indlow])*1852
# min(data$length[indlow])*1852
# max(data$length[indlow])*1852

#Find transect width of each photo (in m)
data$trwidth <- (0.10386*data$Altitude/0.1005)
data$trwidthNM <- (0.10386*data$Altitude/0.1005) / 1852
# mean(data$trwidth[indhigh])
# min(data$trwidth[indhigh])
# max(data$trwidth[indhigh])

# mean(data$trwidth[indlow])
# min(data$trwidth[indlow])
# max(data$trwidth[indlow])


# Adjust area for overaping photos  ---- CHECK TRANSECT 15 WHERE TRANSECT IS IN BOTH DIRECTIONS (done)
# warnings can be ignored
source("adjustareaforoverlap.R")

#Correcting for readers errors
source("correctreading.R")

#Estimate number of pups - without reader correction
EstHarp = kingsley(xycord$x,data$HarpFinalCounts,data$Area,data$Transect,Spacing,gap)
EstHooded = kingsley(xycord$x,data$HoodedFinalCounts,data$Area,data$Transect,Spacing,gap)

#Estimate number of pups - with reader correction - no area correction
EstHarpCorr = kingsley(xycord$x,data$CountsHarpCorr,data$Area,data$Transect,Spacing,gap)
EstHoodedCorr = kingsley(xycord$x,data$CountsHoodedCorr,data$Area,data$Transect,Spacing,gap)

#Estimate number of pups - with reader correction - with area correction
EstHarpCorr = kingsley(xycord$x,data$CountsHarpCorr,data$AreaAdj,data$Transect,Spacing,gap)
EstHoodedCorr = kingsley(xycord$x,data$CountsHoodedCorr,data$AreaAdj,data$Transect,Spacing,gap)


#Uncertainty of pup estimation + uncertainty of reader errors
VarHarpCorr = EstHarpCorr$SE^2+VmeasHarp
VarHoodedCorr = EstHoodedCorr$SE^2+VmeasHooded

#Correct for pups not born and pups that has left the ice prior photographing
#These estimates are obtained using the BirthDistR.R code
QHarp = 0.99
SEQHarp = 0.00014
QHooded = 0.82
SEQHooded = 0.0134

EstHarpFinal = EstHarp$N/QHarp
EstHarpFinalCorr = EstHarpCorr$N/QHarp
EstHoodedFinal = EstHooded$N/QHooded
EstHoodedFinalCorr = EstHoodedCorr$N/QHooded

VarTotalHarp = (1/QHarp)^2*EstHarp$SE^2+(EstHarp$N/QHarp^2)*SEQHarp^2
VarTotalHooded = (1/QHooded)^2*EstHarp$SE^2+(EstHooded$N/QHooded^2)*SEQHooded^2
VarTotalHarpCorr = (1/QHarp)^2*VarHarpCorr+(EstHarpCorr$N/QHarp^2)*SEQHarp^2
VarTotalHoodedCorr = (1/QHooded)^2*VarHoodedCorr+(EstHoodedCorr$N/QHooded^2)*SEQHooded^2

#Print the results to screen
source("PrintToScreenKingsley.R")

############################################################
#Plot various graphics
############################################################
	
#Plot Transects - not working as the columns StartLonT and StopLatT ar not there
uTrans = unique(data$Transect)
X11("",9,7)
par(mar=c(6,5,4,5),bg = "white")
plot(data$StartLonT,data$StopLatT,type = "n",,xlab = "Longitude",ylab = "Latitude",cex.lab = 1.5,cex.main = 1.5,bty = "l")
for (i in 1:length(uTrans)){
	ind = which(data$Transect == uTrans)
	lines(c(data$StartLonT[ind[1]],data$StopLonT[ind[1]]),c(data$StartLatT[ind[1]],data$StopLatT[ind[1]]),col = "blue",lwd = 3)	
}
	
#Plot counts
X11("",9,7)
par(mar=c(6,5,4,5),bg = "white")
symbols(data$FrameLon,data$FrameLat, circles=sqrt(data$HarpFinalCounts/3.141593), inches=0.35, fg="white", bg="red", xlab="Longitude", ylab="Latitude",main = "Spatial distribution of Harp seal pup counts",cex.lab = 1.5,cex.main = 1.5,bty = "l")
X11("",9,7)
par(mar=c(6,5,4,5),bg = "white")
symbols(data$FrameLon,data$FrameLat, circles=sqrt(data$HoodedFinalCounts/3.141593), inches=0.35, fg="white", bg="red", xlab="Longitude", ylab="Latitude",main = "Spatial distribution of Hooded seal pup counts",cex.lab = 1.5,cex.main = 1.5,bty = "l")

#Harps only
X11("",9,7)
par(mar=c(6,5,4,5),bg = "white")
symbols(data$FrameLon,data$FrameLat, circles=sqrt(data$HarpFinalCounts/3.141593), inches=0.2, fg="white", bg="grey10", xlab="Longitude", ylab="Latitude",cex.lab = 1.5,cex.main = 1.5,bty = "l")
#text(-19,72.15,"(A)",cex = 1.5)
#symbols(c(-16,-16,-16),c(71.0,70.9,70.8),circles = sqrt(seq(1,max(data$HoodedFinalCounts),length = 4)/3.141593), inches=0.2, fg="white", bg="red",,add = TRUE)
symbols(c(-16-0.5,-16-0.5,-16-0.5,-16-0.5),c(71.0+0.1,70.85+0.1,70.7+0.1,70.55+0.1),circles = sqrt(c(150,100,50,1)/3.141593), inches=0.2, fg="white", bg="grey10",,add = TRUE)
text(c(-16-0.2,-16-0.2,-16-0.2,-16-0.2),c(71.0+0.1,70.85+0.1,70.7+0.1,70.55+0.1),c("150","100","50","1"))



#harps and hoods
X11("",13,7)
par(mfrow = c(1,2))
par(mar=c(6,5,4,5),bg = "white")
symbols(data$FrameLon,data$FrameLat, circles=sqrt(data$HarpFinalCounts/3.141593), inches=0.2, fg="white", bg="red", xlab="Longitude", ylab="Latitude",main = "Harp seal pups",cex.lab = 1.5,cex.main = 1.5,bty = "l")
text(-19,72.15,"(A)",cex = 1.5)
#symbols(c(-16,-16,-16),c(71.0,70.9,70.8),circles = sqrt(seq(1,max(data$HoodedFinalCounts),length = 4)/3.141593), inches=0.2, fg="white", bg="red",,add = TRUE)
symbols(c(-16-0.5,-16-0.5,-16-0.5,-16-0.5),c(71.0+0.1,70.85+0.1,70.7+0.1,70.55+0.1),circles = sqrt(c(150,100,50,1)/3.141593), inches=0.2, fg="white", bg="red",,add = TRUE)
text(c(-16-0.2,-16-0.2,-16-0.2,-16-0.2),c(71.0+0.1,70.85+0.1,70.7+0.1,70.55+0.1),c("150","100","50","1"))


symbols(data$FrameLon,data$FrameLat, circles=sqrt(data$HoodedFinalCounts/3.141593), inches=0.2, fg="white", bg="red", xlab="Longitude", ylab="Latitude",main = "Hooded seal pups",cex.lab = 1.5,cex.main = 1.5,bty = "l")
text(-19,72.15,"(B)",cex = 1.5)
symbols(c(-16-0.5,-16-0.5,-16-0.5,-16-0.5),c(71.0+0.1,70.85+0.1,70.7+0.1,70.55+0.1),circles = sqrt(c(12,8,4,1)/3.141593), inches=0.2, fg="white", bg="red",,add = TRUE)
text(c(-16-0.2,-16-0.2,-16-0.2,-16-0.2),c(71.0+0.1,70.85+0.1,70.7+0.1,70.55+0.1),c("12","8","4","1"))



#Plot photos taken - I
X11("",9,7)
par(mar=c(6,5,4,5),bg = "white")
recsize = matrix(NA,length(data$Area),2)
recsize[,1] = data$length
recsize[,2] = data$trwidthNM
symbols(xycord$x,xycord$y, rectangles=recsize, inches=FALSE, fg="white", bg="black", xlab="Relative position in NM (East-West)", ylab="Relative position in NM (North-South)",main = "Area covered by photos",cex.lab = 1.5,cex.main = 1.5,bty = "l",add = TRUE)

#Plot photos taken - II
X11("",9,7)
par(mar=c(6,5,4,5),bg = "white")
plot(xycord$x,xycord$y,type = "n",xlim = c(-10,-5),ylim = c(-17,-10))
for (i in 1:500){
	lines(c((xycord$x[i]-data$length[i]/2),(xycord$x[i]+data$length[i]/2)),c((xycord$y[i]-data$trwidthNM[i]/2),(xycord$y[i]-data$trwidthNM[i]/2)),lwd = 1,col = "red")
	lines(c((xycord$x[i]+data$length[i]/2),(xycord$x[i]+data$length[i]/2)),c((xycord$y[i]-data$trwidthNM[i]/2),(xycord$y[i]+data$trwidthNM[i]/2)),lwd = 1,col = "red")
	lines(c((xycord$x[i]+data$length[i]/2),(xycord$x[i]-data$length[i]/2)),c((xycord$y[i]+data$trwidthNM[i]/2),(xycord$y[i]+data$trwidthNM[i]/2)),lwd = 1,col = "red")
	lines(c((xycord$x[i]-data$length[i]/2),(xycord$x[i]-data$length[i]/2)),c((xycord$y[i]+data$trwidthNM[i]/2),(xycord$y[i]-data$trwidthNM[i]/2)),lwd = 1,col = "red")
}


#Find degree of overlap
indphotos = list()
indbb = list()
indwp = list()
indoverlap = list()
indoverlapmean = array(NA,14)
indoverlapmin = array(NA,14)
indoverlapmax = array(NA,14)
overlapprop = array(NA,14)
count = 1
for (Tr in 1:14){
	indTr = which(data$Transect == Tr)
	XcordLower = xycord$x-data$length/2
	XcordUpper = xycord$x+data$length/2

	XcordDiff = array(NA,length(indTr)-1)
	for (i in 2:length(XcordDiff)){
		if (Tr%%2 == 0) {XcordDiff[i] = XcordLower[indTr[i+1]]-XcordUpper[indTr[i]]
		} else {XcordDiff[i] = XcordUpper[indTr[i+1]]-XcordLower[indTr[i]]}
		
	}

	XcordDiffM = XcordDiff*1852
	if (Tr%%2 == 0) {ind = which(XcordDiffM<0)} else {ind = which(XcordDiffM>0)}
	indphotos[[Tr]] = indTr[ind]
	indbb[[Tr]] = data$BB1Reading[indTr[ind]]
	indwp[[Tr]] = data$WP1Reading[indTr[ind]]
	indoverlap[[Tr]] = XcordDiffM[ind]
	indoverlapmean[count] = abs(mean(XcordDiffM[ind]))
	indoverlapmin[count] = abs(min(XcordDiffM[ind]))
	indoverlapmax[count] = abs(max(XcordDiffM[ind]))
	overlapprop[count] = length(ind)/length(indTr)*100
	count = count + 1
}

indzero = which(overlapprop == 0)
indoverlapmean[indzero] = 0
indoverlapmin[indzero] = 0
indoverlapmax[indzero] = 0

overlapsum = data.frame(Transect = 1:14,Percentage = round(100*overlapprop)/100,Min = round(100*indoverlapmin)/100,Max = round(100*indoverlapmax)/100,Mean = round(100*indoverlapmean)/100)



############################################################
#GAM analysis
############################################################
#Proportion zero observations
cat(paste("Proportion zero counts harp seal pups: ",length(which(data$HarpFinalCounts == 0))/length(data$HarpFinalCounts)))
cat(paste("Proportion zero counts hooded seal pups: ",length(which(data$HoodedFinalCounts == 0))/length(data$HoodedFinalCounts)))
data$log.area <- log(data$AreaAdj)
data$x <- xycord$x
data$y <- xycord$y

#Use these for GAM estimates from uncorrected counts
#HarpGam <- gam(HarpFinalCounts ~ s(xycord$x,xycord$y)+offset(log.area),family=negbin(c(0.1,10)),method=gam.method(gam="outer"),gamma = 1.4)
HarpGam <- gam(HarpFinalCounts ~ s(x,y)+offset(log.area),data = data,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
HarpGam$family$getTheta() #Check estimate
HarpGamP <- gam(HarpFinalCounts ~ s(x,y)+offset(log.area),data = data,family="poisson",gamma = 1.4)


#HoodedGam <- gam(HoodedFinalCounts ~ s(xycord$x,xycord$y)+offset(log.area),family=negbin(c(0.1,10)),method=gam.method(gam="outer"),gamma = 1.4)
HoodedGam <- gam(HoodedFinalCounts ~ s(x,y)+offset(log.area),data = data,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
HoodedGam$family$getTheta() #Check estimate
HoodedGamP <- gam(HoodedFinalCounts ~ s(x,y)+offset(log.area),data = data,family="poisson",gamma = 1.4)


#Use these for GAM estimates from corrected counts
#HarpGamCorr <- gam(data$CountsHarpCorr ~ s(xycord$x,xycord$y)+offset(log.area),family=negbin(c(0.1,10)),method=gam.method(gam="outer"),gamma = 1.4)
HarpGamCorr <- gam(CountsHarpCorr ~ s(x,y)+offset(log.area),data = data,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
HarpGamCorr$family$getTheta() #Check estimate

#HoodedGamCorr <- gam(data$CountsHoodedCorr ~ s(xycord$x,xycord$y)+offset(log.area),family=negbin(c(0.1,10)),method=gam.method(gam="outer"),gamma = 1.4)
HoodedGamCorr <- gam(CountsHoodedCorr ~ s(x,y)+offset(log.area),data = data,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
HoodedGamCorr$family$getTheta() #Check estimate


#Finne riktig areal for patchen
meshvalue = ceiling(max(c(abs(min(xycord$x)),abs(max(xycord$x)),abs(min(xycord$y)),abs(max(xycord$y)))));
ds = 0.1 																		#Check this!!!!
meshvector = seq(-meshvalue,meshvalue,ds)
XI = outer(meshvector*0,meshvector,FUN="+")
YI = outer(meshvector,0*meshvector,FUN="+")
dimX = dim(XI)
inmatrix = matrix(0,dimX[1],dimX[2])
empty = inmatrix 

dxnm = 1852			#Unit for nautical miles

tr_label = sort(unique(data$Transect));
for (k in 1:length(tr_label)) {
	tr = tr_label[k]
	postr = which(data$Transect == tr);
	x_tr = sort(xycord$x[postr])
	y_tr = sort(xycord$y[postr])
	
	xv_left = x_tr[1]-1/2*mean(sqrt(data$Area))/dxnm
	xv_right = x_tr[length(x_tr)]+1/2*mean(sqrt(data$Area))/dxnm

	yv_left = y_tr[1]-1/2*Spacing*1.0
	yv_right = y_tr[length(y_tr)]+1/2*Spacing*1.0
	
	indx <- which(YI >= xv_left & YI <= xv_right & XI >= yv_left & XI <= yv_right)
	
	inmatrix[indx] = 1
	
	}

#plot this for checking the area selected
#image(meshvector,meshvector,inmatrix)	

area.index = which(inmatrix > 0)
area.latlon = which(inmatrix > 0,arr.ind=TRUE)
lon = meshvector[area.latlon[,1]]
lat = meshvector[area.latlon[,2]]


#plot(lon,lat)			#Use for check area
#lines(xycord$x,xycord$y)

#Create new data frame
sealdatap <- data.frame(x=lon,y=lat,log.area=0*lon)

#Predict the GAm model on the new surface - uncorrected counts
zzHarp <- array(NA,length(meshvector)^2)
zzHarp[area.index] <- predict(HarpGam,sealdatap)
image(meshvector,meshvector,matrix(zzHarp,dimX[1],dimX[2]),xlim=c(-25,25),ylim = c(-25,25),xlab = "Relative position (nm)",ylab = "Relative position (nm)")
contour(meshvector,meshvector,matrix(zzHarp,dimX[1],dimX[1]),add=TRUE)

zzHooded <- array(NA,length(meshvector)^2)
zzHooded[area.index] <- predict(HoodedGam,sealdatap)
image(meshvector,meshvector,matrix(zzHooded,dimX[1],dimX[2]),xlim=c(-25,25),ylim = c(-25,25),xlab = "Relative position (nm)",ylab = "Relative position (nm)")
contour(meshvector,meshvector,matrix(zzHooded,dimX[1],dimX[1]),add=TRUE)

#Both together
X11("",13,7)
par(mfrow = c(1,2))
par(mar=c(6,5,4,5),bg = "white")
image(meshvector,meshvector,matrix(zzHarp,dimX[1],dimX[2]),xlim=c(-25,25),ylim = c(-25,25),xlab = "Relative position (nm)",ylab = "Relative position (nm)",main = "Harp seals",cex.lab = 1.5,cex.main = 1.5,bty = "l")
contour(meshvector,meshvector,matrix(zzHarp,dimX[1],dimX[1]),add=TRUE)
text(-22,22,"(A)",cex = 1.5)
image(meshvector,meshvector,matrix(zzHooded,dimX[1],dimX[2]),xlim=c(-25,25),ylim = c(-25,25),xlab = "Relative position (nm)",ylab = "Relative position (nm)",main = "Hooded seals",cex.lab = 1.5,cex.main = 1.5,bty = "l")
contour(meshvector,meshvector,matrix(zzHooded,dimX[1],dimX[1]),add=TRUE)
text(-22,22,"(B)",cex = 1.5)




#Predict the GAm model on the new surface - corrected counts
zzHarpCorr <- array(NA,length(meshvector)^2)
zzHarpCorr[area.index] <- predict(HarpGamCorr,sealdatap)
image(meshvector,meshvector,matrix(zzHarpCorr,dimX[1],dimX[1]),xlim=c(-25,25),ylim = c(-25,25),xlab = "Relative position (nm)",ylab = "Relative position (nm)")
contour(meshvector,meshvector,matrix(zzHarpCorr,dimX[1],dimX[1]),add=TRUE)

zzHoodedCorr <- array(NA,length(meshvector)^2)
zzHoodedCorr[area.index] <- predict(HoodedGamCorr,sealdatap)
image(meshvector,meshvector,matrix(zzHoodedCorr,dimX[1],dimX[1]),xlim=c(-25,25),ylim = c(-25,25),xlab = "Relative position (nm)",ylab = "Relative position (nm)")
contour(meshvector,meshvector,matrix(zzHoodedCorr,dimX[1],dimX[1]),add=TRUE)



#Number of pups
NHarpGAM = sum(exp(zzHarp),na.rm=TRUE)*(ds*dxnm)^2
NHoodedGAM = sum(exp(zzHooded),na.rm=TRUE)*(ds*dxnm)^2
NHarpGAMCorr = sum(exp(zzHarpCorr),na.rm=TRUE)*(ds*dxnm)^2
NHoodedGAMCorr = sum(exp(zzHoodedCorr),na.rm=TRUE)*(ds*dxnm)^2

XpHarp <- predict(HarpGam,newdata=sealdatap,type="lpmatrix")
XpHooded <- predict(HoodedGam,newdata=sealdatap,type="lpmatrix")
XpHarpCorr <- predict(HarpGamCorr,newdata=sealdatap,type="lpmatrix")
XpHoodedCorr <- predict(HoodedGamCorr,newdata=sealdatap,type="lpmatrix")

dimXp = dim(XpHarp)
splinesHarp = array(0,dimXp[1])
splinesHooded = array(0,dimXp[1])
splinesHarpCorr = array(0,dimXp[1])
splinesHoodedCorr = array(0,dimXp[1])

for(ell in 1:dimXp[1]){
	splinesHarp[ell]=XpHarp[ell,]%*%coef(HarpGam)
	splinesHooded[ell]=XpHooded[ell,]%*%coef(HoodedGam)
	splinesHarpCorr[ell]=XpHarpCorr[ell,]%*%coef(HarpGamCorr)
	splinesHoodedCorr[ell]=XpHoodedCorr[ell,]%*%coef(HoodedGamCorr)

	}

VarHarpGAM = exp(t(splinesHarp))%*%XpHarp%*%HarpGam$Vp%*%t(XpHarp)%*%exp(splinesHarp)*(ds*dxnm)^4
SEHarpGAM = sqrt(VarHarpGAM)
CVHarpGAM = SEHarpGAM/NHarpGAM*100
VarHoodedGAM = exp(t(splinesHooded))%*%XpHooded%*%HoodedGam$Vp%*%t(XpHooded)%*%exp(splinesHooded)*(ds*dxnm)^4
SEHoodedGAM = sqrt(VarHoodedGAM)
CVHoodedGAM = SEHoodedGAM/NHoodedGAM*100

VarHarpGAMCorr = exp(t(splinesHarpCorr))%*%XpHarpCorr%*%HarpGamCorr$Vp%*%t(XpHarpCorr)%*%exp(splinesHarpCorr)*(ds*dxnm)^4
SEHarpGAMCorr = sqrt(VarHarpGAMCorr)
CVHarpGAMCorr = SEHarpGAMCorr/NHarpGAMCorr*100
VarHoodedGAMCorr = exp(t(splinesHoodedCorr))%*%XpHoodedCorr%*%HoodedGamCorr$Vp%*%t(XpHoodedCorr)%*%exp(splinesHoodedCorr)*(ds*dxnm)^4
SEHoodedGAMCorr = sqrt(VarHoodedGAMCorr)
CVHoodedGAMCorr = SEHoodedGAMCorr/NHoodedGAMCorr*100


NHarpGAMFinal = NHarpGAM/QHarp
NHarpGAMCorr = NHarpGAMCorr/QHarp
NHoodedGAMFinal = NHoodedGAM/QHooded
NHoodedGAMCorr = NHoodedGAMCorr/QHooded


VarTotalHarpGAM = (1/QHarp)^2*VarHarpGAM+(NHarpGAM/QHarp^2)*SEQHarp^2
VarTotalHoodedGAM = (1/QHooded)^2*VarHoodedGAM+(NHoodedGAM/QHooded^2)*SEQHooded^2

VarTotalHarpGAMCorr = (1/QHarp)^2*VarHarpGAMCorr+(NHarpGAMCorr/QHarp^2)*SEQHarp^2
VarTotalHoodedGAMCorr = (1/QHooded)^2*VarHoodedGAMCorr+(NHoodedGAMCorr/QHooded^2)*SEQHooded^2

TotalAreaCovered = sum(inmatrix)*(ds*dxnm)^2/(1000*1000)

#Print results to screen
cat("------------------------------------------------------------------------------------------------")
cat("Harp seals \n")
cat(paste("Total estimates of Harp seal pups without correction: ",NHarpGAMFinal," , SE = ",sqrt(VarTotalHarpGAM),", CV = ",100*sqrt(VarTotalHarpGAM)/NHarpGAMFinal,"%\n"))
cat(paste("Total estimates of Harp seal pups with correction: ",NHarpGAMCorr," , SE = ",sqrt(VarTotalHarpGAMCorr),", CV = ",100*sqrt(VarTotalHarpGAMCorr)/NHarpGAMCorr,"%\n"))
cat(paste("  * Estimate without reader correction:",NHarpGAM," SE = ",SEHarpGAM,", CV = ",CVHarpGAM,"% \n"))
cat(paste("  * Estimate with reader correction:   ",NHarpGAMCorr," SE = ",SEHarpGAMCorr,", CV = ",CVHarpGAMCorr,"% \n"))
cat("------------------------------------------------------------------------------------------------")
cat("Hooded seals \n")
cat(paste("Total estimates of Hooded seal pups without correction: ",NHoodedGAMFinal," , SE = ",sqrt(VarTotalHoodedGAM),", CV = ",100*sqrt(VarTotalHoodedGAM)/NHoodedGAMFinal,"%\n"))
cat(paste("Total estimates of Hooded seal pups with correction: ",NHoodedGAMCorr," , SE = ",sqrt(VarTotalHoodedGAMCorr),", CV = ",100*sqrt(VarTotalHoodedGAMCorr)/NHoodedGAMCorr,"%\n"))
cat(paste("  * Estimate without reader correction:",NHoodedGAM," SE = ",SEHoodedGAM,", CV = ",CVHoodedGAM,"% \n"))
cat(paste("  * Estimate with reader correction:   ",NHoodedGAMCorr," SE = ",SEHoodedGAMCorr,", CV = ",CVHoodedGAMCorr,"% \n"))
cat("------------------------------------------------------------------------------------------------")



