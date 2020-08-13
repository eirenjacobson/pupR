##########################################
indphotos = list()
indbb = list()
indwp = list()
indoverlap = list()
indoverlapmean = array(NA,14)
indoverlapmin = array(NA,14)
indoverlapmax = array(NA,14)
overlapprop = array(NA,14)
count = 1

data$AreaAdj = data$Area
XcordLower = xycord$x-data$length/2
XcordUpper = xycord$x+data$length/2
YcordLower = xycord$y-data$trwidthNM/2
YcordUpper = xycord$y+data$trwidthNM/2

utrans = sort(unique(data$Transect))

#They had to redo transect 24
for (Tr in 1:23){
	indTr = which(data$Transect == Tr)
	
	XcordDiff = array(NA,length(indTr)-1)
	YcordDiff = array(NA,length(indTr)-1)

	for (i in 2:length(XcordDiff)){
		if (Tr%%2 == 0) {XcordDiff[i] = XcordLower[indTr[i+1]]-XcordUpper[indTr[i]]
		} else {XcordDiff[i] = XcordUpper[indTr[i+1]]-XcordLower[indTr[i]]}
		
	}

	XcordDiffM = XcordDiff*1852
	if (Tr%%2 == 0) {ind = which(XcordDiffM<0)} else {ind = which(XcordDiffM>0)}
	indphotos[[Tr]] = indTr[ind]-1
	indbb[[Tr]] = data$BB1Reading[indTr[ind]-1]
	indwp[[Tr]] = data$WP1Reading[indTr[ind]-1]
	indoverlap[[Tr]] = XcordDiffM[ind]
	indoverlapmean[count] = abs(mean(XcordDiffM[ind]))
	indoverlapmin[count] = abs(min(XcordDiffM[ind]))
	indoverlapmax[count] = abs(max(XcordDiffM[ind]))
	overlapprop[count] = length(ind)/length(indTr)*100
	count = count + 1
	
	PhotosWithOverlap = indTr[ind]-1
	
	#Check if there is any overlap and if on Tr 15 do nothing since they fly both directions there
	if (length(PhotosWithOverlap)>0 & Tr != 15){
	DeltaX = abs(XcordDiffM[ind])
	DeltaY = array(NA,length(DeltaX))
	
	for (j in 1:(length(PhotosWithOverlap)-1)){
		#Check if the next photo in line comes right after
		if ((PhotosWithOverlap[j+1]-PhotosWithOverlap[j])==1){
			#Direction of flight irrelevant
			#Find out if next photo is higher or lower
			if (xycord$y[PhotosWithOverlap[j]+1]<xycord$y[PhotosWithOverlap[j]]){
				DeltaY[j] = 1852*(YcordUpper[PhotosWithOverlap[j]+1]-YcordLower[PhotosWithOverlap[j]])
				a = 0			
			} else {
				DeltaY[j] = 1852*(YcordUpper[PhotosWithOverlap[j]]-YcordLower[PhotosWithOverlap[j]+1])		
				a = 1	
			}
						
			data$AreaAdj[PhotosWithOverlap[j]] = data$AreaAdj[PhotosWithOverlap[j]] - DeltaX[j]*DeltaY[j]
		} 
	}
	#Do what's needed for last photo
	j = length(PhotosWithOverlap)
	#Find out if next photo is higher or lower
	if (xycord$y[PhotosWithOverlap[j]+1]<xycord$y[PhotosWithOverlap[j]]){
			DeltaY[j] = 1852*(YcordUpper[PhotosWithOverlap[j]+1]-YcordLower[PhotosWithOverlap[j]])			
		} else {
			DeltaY[j] = 1852*(YcordUpper[PhotosWithOverlap[j]]-YcordLower[PhotosWithOverlap[j]+1])			
		}
						
	data$AreaAdj[PhotosWithOverlap[j]] = data$AreaAdj[PhotosWithOverlap[j]] - DeltaX[j]*DeltaY[j]
	}
}

for (Tr in 24:27){
	indTr = which(data$Transect == Tr)
	
	XcordDiff = array(NA,length(indTr)-1)
	YcordDiff = array(NA,length(indTr)-1)

	for (i in 2:length(XcordDiff)){
		if (Tr%%2 != 0) {XcordDiff[i] = XcordLower[indTr[i+1]]-XcordUpper[indTr[i]]
		} else {XcordDiff[i] = XcordUpper[indTr[i+1]]-XcordLower[indTr[i]]}
		
	}

	XcordDiffM = XcordDiff*1852
	if (Tr%%2 != 0) {ind = which(XcordDiffM<0)} else {ind = which(XcordDiffM>0)}
	indphotos[[Tr]] = indTr[ind]-1
	indbb[[Tr]] = data$BB1Reading[indTr[ind]-1]
	indwp[[Tr]] = data$WP1Reading[indTr[ind]-1]
	indoverlap[[Tr]] = XcordDiffM[ind]
	indoverlapmean[count] = abs(mean(XcordDiffM[ind]))
	indoverlapmin[count] = abs(min(XcordDiffM[ind]))
	indoverlapmax[count] = abs(max(XcordDiffM[ind]))
	overlapprop[count] = length(ind)/length(indTr)*100
	count = count + 1
	
	PhotosWithOverlap = indTr[ind]-1
	#Check if there is any overlap 
	if (length(PhotosWithOverlap)>0){
	DeltaX = abs(XcordDiffM[ind])
	DeltaY = array(NA,length(DeltaX))
	
	for (j in 1:(length(PhotosWithOverlap)-1)){
		#Check if the next photo in line comes right after
		if ((PhotosWithOverlap[j+1]-PhotosWithOverlap[j])==1){
			#Direction of flight irrelevant
			#Find out if next photo is higher or lower
			if (xycord$y[PhotosWithOverlap[j]+1]<xycord$y[PhotosWithOverlap[j]]){
				DeltaY[j] = 1852*(YcordUpper[PhotosWithOverlap[j]+1]-YcordLower[PhotosWithOverlap[j]])
				a = 0			
			} else {
				DeltaY[j] = 1852*(YcordUpper[PhotosWithOverlap[j]]-YcordLower[PhotosWithOverlap[j]+1])		
				a = 1	
			}
						
			data$AreaAdj[PhotosWithOverlap[j]] = data$AreaAdj[PhotosWithOverlap[j]] - DeltaX[j]*DeltaY[j]
		} 
	}
	#Do what's needed for last photo
	j = length(PhotosWithOverlap)
	#Find out if next photo is higher or lower
	if (xycord$y[PhotosWithOverlap[j]+1]<xycord$y[PhotosWithOverlap[j]]){
			DeltaY[j] = 1852*(YcordUpper[PhotosWithOverlap[j]+1]-YcordLower[PhotosWithOverlap[j]])			
		} else {
			DeltaY[j] = 1852*(YcordUpper[PhotosWithOverlap[j]]-YcordLower[PhotosWithOverlap[j]+1])			
		}
						
	data$AreaAdj[PhotosWithOverlap[j]] = data$AreaAdj[PhotosWithOverlap[j]] - DeltaX[j]*DeltaY[j]
	}
}



data$AreaDiff = data$Area-data$AreaAdj
ind = which(data$AreaDiff > 0)
min(data$AreaDiff[ind])
max(data$AreaDiff[ind])
sum(data$AreaDiff[ind])

nphotos = 0
for (i in 1:27){nphotos = nphotos + length(indphotos[[i]])}
