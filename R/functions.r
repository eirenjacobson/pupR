#' Convert latitude and longitude to Cartesian koordinates
#'
#' @param survey Specify which survey to analyze
#' @param fname File name of data file
#' @param population Specify which population to be analysed, harp, hood or both.
#' @return
#' @keywords
#' @export
#' @examples
#' lb3xy()

loadAndPrepareData <- function(survey = "WestIce2012",
                               fname = "WestIce2012_MichaelAndLotta_final.csv",
                               population = c("harp","hood")
                               )
{

  data <- read.table(paste0("data/",survey,"/",fname),header = TRUE,sep = ",")

  #cat("\n Removing faulty photos....")
  #Remove some photos marked in the comment column
  indRemove = which(data$Comments == "Remove")
  data <- data[-indRemove,]

  dimframe <- dim(data)

  Nphotos <- dimframe[1]

  xycord = lb2xy(data$FrameLon,data$FrameLat,mean(data$FrameLon),mean(data$FrameLat))

  uTrans = sort(unique(data$Transect))
  nphotos = list()
  meanypos = array(NA,27)


  if("harp" %in% population){
    #Sjekk bilder som er lest to ganger og bruk tall fra 2. lesning
    ind2countHarp = which(data$WP2Reading != "NA")
    ind2countHooded = which(data$BB2Reading != "NA")

    data$HarpFinalCounts <- data$WP1Reading
    data$HarpFinalCounts[ind2countHarp] <- data$WP2Reading[ind2countHarp]

    #Set empty photos to zero
    indNAHarp = which(is.na(data$HarpFinalCounts))
    data$HarpFinalCounts[indNAHarp] = 0

    #Find the number of zero photos
    indzeroharp = which(data$HarpFinalCounts == 0)
    propzeroharp = length(indzeroharp)/length(data$HarpFinalCounts)
    cat(paste("Proportion of zero counts for harp seals:",round(propzeroharp,digits = 2)),"\n")

    #Find the number of photos in each transect
    ncountsHarp = list()
    totalHarps = 0
    for (i in 1:length(uTrans)){
      nphotos[[i]] = length(which(data$Transect == uTrans[i]))
      ncountsHarp[[i]] = sum(data$HarpFinalCounts[data$Transect == uTrans[i]])
      totalHarps = totalHarps + ncountsHarp[[i]]
      meanypos[i] = mean(xycord$y[data$Transect == uTrans[i]])
    }

  }

  mean(abs(diff(meanypos)))


  if("hood" %in% population){
    #Sjekk bilder som er lest to ganger og bruk tall fra 2. lesning
    ind2countHooded = which(data$BB2Reading != "NA")

    #Set empty photos to zero
    data$HoodedFinalCounts <- data$BB1Reading
    data$HoodedFinalCounts[ind2countHooded] <- data$BB2Reading[ind2countHooded]
    indNAHooded = which(is.na(data$HoodedFinalCounts))
    data$HoodedFinalCounts[indNAHooded] = 0

    #Find the number of zero photos
    indzerohooded = which(data$HoodedFinalCounts == 0)
    propzerohooded = length(indzerohooded)/length(data$HoodedFinalCounts)
    cat(paste("Proportion of zero counts for hooded seals:",round(propzerohooded,digits = 2)),"\n")

    ncountsHooded = list()
    totalHooded = 0

    for (i in 1:length(uTrans)){
      nphotos[[i]] = length(which(data$Transect == uTrans[i]))
      ncountsHooded[[i]] = sum(data$HoodedFinalCounts[data$Transect == uTrans[i]])
      totalHooded = totalHooded + ncountsHooded[[i]]
      meanypos[i] = mean(xycord$y[data$Transect == uTrans[i]])
    }

  }

  #Find area covered of each photo
  data$Area <- (0.06786*data$Altitude/0.1005)*(0.10386*data$Altitude/0.1005)

  #Find length covered of each photo (in Nautical Miles)
  data$length <- (0.06786*data$Altitude/0.1005) / 1852

  #Find transect width of each photo (in m)
  data$trwidth <- (0.10386*data$Altitude/0.1005)
  data$trwidthNM <- (0.10386*data$Altitude/0.1005) / 1852

  returnList = list()
  returnList$data = data
  returnList$xycord = xycord

  return(returnList)
}

#' Convert latitude and longitude to Cartesian koordinates
#'
#' @param lon Longitude
#' @param lat Latitude
#' @param lon0 Mean longitude
#' @param lat0 Mean latitude
#' @return
#' @keywords
#' @export
#' @examples
#' lb2xy()

lb2xy <- function(lon,lat,lon0,lat0)
{
    n = length(lon);
    l0 = lon0/180*pi;
    b0 = lat0/180*pi;
    R=6360/1.852;
    l = lon/180*pi;
    b = lat/180*pi;

    x=array(0,length(l));
    y=array(0,length(b));
    cb0=cos(b0);
    sb0=sin(b0);
    A = rbind(c(0,1,0),c(-sb0,0,cb0),c(cb0,0,sb0))
    cl=cos(l-l0);
    sl=sin(l-l0);
    cb=cos(b);
    sb=sin(b);
    xg=rbind(cb*cl,cb*sl,sb)
    xp=A%*%xg;
    x=xp[1,]*R;
    y=xp[2,]*R;
    NM = data.frame(x=x,y=y)
    return(NM)
}

#' Estimates the pup production using the Kingsley method
#'
#' @param x X Coordinates in Cartesian system
#' @param count Vector containing the pup counts for each photo
#' @param area Vector containing area of each photo
#' @param transect Vector containing the transect each photo comes from
#' @param Spacing Transect width in NM, default 3
#' @param gap Maximum distance in NM between allowed, for removing holes where camera has been stopped due to flying over open water, default = 1
#' @return
#' @keywords
#' @export
#' @examples
#' kingsley()

kingsley <- function(x,count,area,transect,Spacing = 3,gap = 1)
{

    N = length(x);
    tr_label = sort(unique(transect));
    N_tr      = length(tr_label);

    disttr = array(0,N_tr)
    areatr = array(0,N_tr)
    count_ref = array(0,N_tr)

    for(i in 1:N_tr){
        tr = tr_label[i];
        ptr = which(transect == tr);
        xtr = x[ptr];
        distance = sqrt(diff(sort(xtr))^2);                                	   #Find the distance between each photo
        ind1nm_vec= which(distance > gap);                                   	   #Find area where water occurs
        #disttr[i] = 1852*abs(xtr[1]-xtr[length(xtr)]) - sum(distance[ind1nm_vec]);  #Transect length minus length of gaps
         disttr[i] = 1852*abs(min(xtr)-max(xtr)) - 0*sum(distance[ind1nm_vec]);  #Transect length minus length of gaps


        areatr[i] = sum(area[ptr]);                                        	   #Sampled area of transect
        count_ref[i] = sum(count[ptr])/areatr[i];                              #Estimate number of pups along transect

    }


    N_new = (1852*Spacing)*sum(disttr*count_ref);
    V_new = (1852*Spacing*(1852*Spacing-sum(area)/sum(disttr)))/(2*(N_tr-1))*sum(disttr^2)*sum(diff(count_ref)^2);
    V_new = (1852*Spacing*N_tr*(1852*Spacing-sum(area)/sum(disttr)))/(2*(N_tr-1))*sum(diff(disttr*count_ref)^2)
    SE_new = sqrt(V_new);
    CV_new = 100*SE_new/N_new;
    estimates <- data.frame(N = round(N_new),SE = SE_new,CV = CV_new)
    return(estimates)
	}


#' Plot GAM fit do data
#'
#' @param seal.mod
#' @param x Vector containing x coordinates in Cartesian system
#' @param y Vector containing y coordinates in Cartesian system
#' @param transect Transect each photo beongs to
#' @param area Area of each photo
#' @param spacing Distance between transects
#' @param patch ?
#' @return
#' @keywords
#' @export
#' @examples
#' plotGAM()


plotGAM <- function(seal.mod,x,y,transect,area,spacing,patch){

	meshvalue = ceiling( max(c(abs(min(x)),abs(max(x)),abs(min(y)),abs(max(y)))));
	dxnm = 1852
	ds = 0.1
	meshvector = seq(-meshvalue,meshvalue,ds)
	XI = outer(meshvector*0,meshvector,FUN="+")
	YI = outer(meshvector,0*meshvector,FUN="+")
	dimX = dim(XI)
	inmatrix = matrix(0,dimX[1],dimX[2])
	empty = inmatrix



	tr_label = unique(transect);
	for (k in 1:length(tr_label)) {
		tr = tr_label[k]
		postr = which(transect == tr);
		x_tr = sort(x[postr])
		y_tr = sort(y[postr])

		xv_left = x_tr[1]-1/2*mean(sqrt(area))/dxnm
		xv_right = x_tr[length(x_tr)]+1/2*mean(sqrt(area))/dxnm

		yv_left = y_tr[1]-1/2*Spacing*1.0
		yv_right = y_tr[length(y_tr)]+1/2*Spacing*1.0

		indx <- which(YI >= xv_left & YI <= xv_right & XI >= yv_left & XI <= yv_right)

		inmatrix[indx] = 1

		}

#image(meshvector,meshvector,inmatrix)

	area.index = which(inmatrix > 0)
	area.latlon = which(inmatrix > 0,arr.ind=TRUE)
	lon = meshvector[area.latlon[,1]]
	lat = meshvector[area.latlon[,2]]

#plot(lon,lat)
#lines(x,y)

	sealdatap <- data.frame(x=lon,y=lat,log.area=0*lon)

	zz <- array(NA,length(meshvector)^2)
	zz[area.index] <- predict(seal.mod,sealdatap)

	#old.par.settings <- par(cex.lab=1.5)


	image(meshvector,meshvector,matrix(zz,dimX[1],dimX[1]),main= paste("Estimated GAM surface for Patch ",patch),xlab = "Relative position (nm)",ylab = "Relative position (nm)")
	contour(meshvector,meshvector,matrix(zz,dimX[1],dimX[1]),add=TRUE)
	#par(old.par.settings)
}

#' Estimate uncertainty from measurement errors
#'
#' @param x X Coordinates in Cartesian system
#' @param count Vector containing the pup counts for each photo
#' @param area Vector containing area of each photo
#' @param transect Vector containing the transect each photo comes from
#' @param Spacing Distance between transect
#' @param gap Distance threshold for length if flown over open water
#' @param var_b
#' @param var_u
#' @param var_a
#' @param cov_ab
#' @param intcpt
#' @return
#' @keywords
#' @export
#' @examples
#' vmeasest()

vmeasest <- function(x,
                     count,
                     area,
                     transect,
                     Spacing,
                     gap,
                     var_b,
                     var_u,
                     var_a,
                     cov_ab = NA,
                     intcpt = FALSE)
	{
		N = length(x)
		tr_label = sort(unique(transect))
		N_tr = length(tr_label)

	    distr = array(0,N_tr)
        Fj = array(0,N_tr)
        term1 = array(0,N_tr)
        term2 = array(0,N_tr)
        term3 = array(0,N_tr)
        term4 = array(0,N_tr)

		for(i in 1:N_tr){
        	tr = tr_label[i];
        	ptr = which(transect == tr);
        	xtr = x[ptr];
        	distance = sqrt(diff(sort(xtr))^2);                                	   #Find the distance between each photo
        	ind1nm_vec= which(distance > gap);                                   	   #Find area where water occurs
        	#disttr[i] = abs(xtr[1]-xtr[length(xtr)]) - sum(distance[ind1nm_vec]);  #Transect length minus length of gaps
            distr[i] = 1*abs(min(xtr)-max(xtr)) - 0*sum(distance[ind1nm_vec])
        	K_tr = length(ptr)
        	Fj[i] = (1852*distr[i])/sum(area[ptr])
        	term1[i] = Fj[i]*K_tr
        	term2[i] = Fj[i]*sum(count[ptr])
        	term3[i] = Fj[i]^2*K_tr
        	term4[i] = sum(var_b*count[ptr]^2)

    	}
		Spacing = 1852*Spacing;

		if (intcpt == TRUE) {
			V_meas = Spacing^2*(sum(term1)^2*var_a+2*cov_ab*sum(term1)*sum(term2)+sum(term2)^2*var_b+sum(term3)*var_u)
			#cat("\n Intercept used\n")
			} else {
			  V_meas = Spacing^2*(sum(term2)^2*var_b+sum(term3)*var_u)
				#cat("\n No intercept used\n")
			  }

		#V_meas = T*sum(Fj*term4)
		return(V_meas)

}

#' Estimate uncertainty from readers errors
#'
#' @param data The data set with photo counts
#' @param type Type of correction. Use 1 for separate correction for each reader (assumes two readers) and 2 for the same correction for both readers
#' @param population Specify which population or populations to estimate the uncertainty from
#' @param readers Specify the ID of each reader. Assumes 2 readers
#' @param transect Vector containing the transect each photo comes from
#' @param Spacing Transect width
#' @param gap Distance threshold for length if flown over open water
#' @param var_b
#' @param var_u
#' @param var_a
#' @param cov_ab
#' @param intcpt
#' @return
#' @keywords
#' @export
#' @examples
#' correctReading()

correctReading <- function(dataList,
                           type = 1,
                           population = c("harp","hood"),
                           readers = c("LL","MP"),
                           plotFig = TRUE,
                           Spacing = 3,
                           gap = 1,
                           grDev = TRUE)

  {

  data = dataList$data
  xycord = dataList$xycord
  if (type == 1){

    ##############################
    #Separate fit for each reader

    #Reader 1
    if("harp" %in% population){
      indcorrHarpLL = which((1-is.na(data$WPTrue)==1) & data$Reader == readers[1])			#Find which photos have been checked by LL
      lincorrHarpLL = lm(WPTrue[indcorrHarpLL]~HarpFinalCounts[indcorrHarpLL] - 1,data = data)
      if(plotFig){
        #windows(width = 11,height = 7)
        if(grDev) graphDev(width = width,height = height)
        plot(data$HarpFinalCounts[indcorrHarpLL],data$WPTrue[indcorrHarpLL],xlab = "Harp seal pup counts",ylab = "True number of harp seal pups",main = "Linear fit of uncorrected counts to true number of pups",sub = "Reader 1")
        abline(a = 0, b = lincorrHarpLL$coefficients,lwd = 3,col = "red")
      }
    }

    if("hood" %in% population){
      indcorrHoodedLL = which(1-is.na(data$BBTrue)==1 & data$Reader == readers[1])			#Find which photos have been checked by LL
      lincorrHoodedLL = lm(data$BBTrue[indcorrHoodedLL]~data$HoodedFinalCounts[indcorrHoodedLL] - 1)
      if(plotFig){
        #windows(width = 11,height = 7)
        if(grDev) graphDev(width = width,height = height)

        plot(data$HoodedFinalCounts[indcorrHoodedLL],data$BBTrue[indcorrHoodedLL],xlab = "Hooded seal pup counts",ylab = "True number of hooded pups",main = "Linear fit of uncorrected counts to true number of pups",sub = "Reader 1")
        abline(a = 0, b = lincorrHoodedLL$coefficients,lwd = 3,col = "red")
      }
    }

    if("harp" %in% population){
      var_bHarpLL = coef(summary(lincorrHarpLL))[1,2]; var_bHarpLL = var_bHarpLL^2
      var_uHarpLL = sd(lincorrHarpLL$residuals); var_uHarpLL = var_uHarpLL^2
    }

    if("hood" %in% population){
      var_bHoodedLL = coef(summary(lincorrHoodedLL))[1,2]; var_bHoodedLL = var_bHoodedLL^2
      var_uHoodedLL = sd(lincorrHoodedLL$residuals); var_uHoodedLL = var_uHoodedLL^2

    }

    #Only correct those reader 1 has read
    indLL = which(data$Reader == readers[1])
    if("harp" %in% population) data$CountsHarpCorr[indLL] <- lincorrHarpLL$coefficients*data$HarpFinalCounts[indLL]
    if("hood" %in% population) data$CountsHoodedCorr[indLL] <- lincorrHoodedLL$coefficients*data$HoodedFinalCounts[indLL]

    #Get uncertainty associatet with reader errors
    if("harp" %in% population) VmeasHarpLL = vmeasest(xycord$x[indLL],data$HarpFinalCounts[indLL],data$Area[indLL],data$Transect[indLL],Spacing,gap,var_bHarpLL,var_uHarpLL)
    if("hood" %in% population) VmeasHoodedLL = vmeasest(xycord$x[indLL],data$HoodedFinalCounts[indLL],data$Area[indLL],data$Transect[indLL],Spacing,gap,var_bHoodedLL,var_uHoodedLL)



    #Reader 2
    if("harp" %in% population){

      indcorrHarpMP = which(1-is.na(data$WPTrue)==1 & data$Reader == readers[2])			#Find which photos have been checked by MP
      lincorrHarpMP = lm(data$WPTrue[indcorrHarpMP]~data$HarpFinalCounts[indcorrHarpMP] - 1)
      if(plotFig){
        #windows(width = 11,height = 7)
        if(grDev) graphDev(width = width,height = height)
        plot(data$HarpFinalCounts[indcorrHarpMP],data$WPTrue[indcorrHarpMP],xlab = "Harp seal pup counts",ylab = "True number of harp seal pups",main = "Linear fit of uncorrected counts to true number of pups",sub = "Reader 2")
        abline(a = 0, b = lincorrHarpMP$coefficients,lwd = 3,col = "red")
      }
    }

    if("hood" %in% population){
      indcorrHoodedMP = which(1-is.na(data$BBTrue)==1 & data$Reader == readers[2])			#Find which photos have been checked by MP
      lincorrHoodedMP = lm(data$BBTrue[indcorrHoodedMP]~data$HoodedFinalCounts[indcorrHoodedMP] - 1)
      #plot(data$HoodedFinalCounts[indcorrHoodedMP],data$BBTrue[indcorrHoodedMP])
      #abline(a = 0, b = lincorrHoodedMP$coefficients,lwd = 3,col = "red")
    }


    if("harp" %in% population){
      var_bHarpMP = coef(summary(lincorrHarpMP))[1,2]; var_bHarpMP = var_bHarpMP^2
      var_uHarpMP = sd(lincorrHarpMP$residuals); var_uHarpMP = var_uHarpMP^2

      #Only correct photos ready by reader 2
      indMP = which(data$Reader == readers[2])
      data$CountsHarpCorr[indMP] <- lincorrHarpMP$coefficients*data$HarpFinalCounts[indMP]

      #Get uncertainty associatet with reader errors
      VmeasHarpMP = vmeasest(xycord$x[indMP],data$HarpFinalCounts[indMP],data$Area[indMP],data$Transect[indMP],Spacing,gap,var_bHarpMP,var_uHarpMP)

    }

    if("hood" %in% population){
      var_bHoodedMP = coef(summary(lincorrHoodedMP))[1,2]; var_bHoodedMP = var_bHoodedMP^2
      var_uHoodedMP = sd(lincorrHoodedMP$residuals); var_uHoodedMP = var_uHoodedMP^2

      #Only correct photos ready by reader 2
      indMP = which(data$Reader == readers[2])
      data$CountsHoodedCorr [indMP]<- lincorrHoodedMP$coefficients*data$HoodedFinalCounts[indMP]

      #Get uncertainty associatet with reader errors
      VmeasHoodedMP = vmeasest(xycord$x[indMP],data$HoodedFinalCounts[indMP],data$Area[indMP],data$Transect[indMP],Spacing,gap,var_bHoodedMP,var_uHoodedMP)


    }

    correctedData = list()
    correctedData$data = data

    if("harp" %in% population){
      VmeasHarp = VmeasHarpLL + VmeasHarpMP
      correctedData$VmeasHarp = VmeasHarp
    }

    if("hood" %in% population){
      VmeasHooded = VmeasHoodedLL + VmeasHoodedMP
      correctedData$VmeasHooded = VmeasHooded
    }

    correctedData$xycord = xycord

    return(correctedData)
  }

  if (type == 2){

    #Both readers pooled
    if("harp" %in% population){
      indcorrHarp = which(1-is.na(data$WPTrue)==1)			#Find which photos have been checked by both readers
      lincorrHarp = lm(data$WPTrue[indcorrHarp]~data$HarpFinalCounts[indcorrHarp])
      if(plotFig){
        windows(width = 11,height = 7)
        plot(data$HarpFinalCounts[indcorrHarp],data$WPTrue[indcorrHarp],xlab = "Harp seal pup counts",ylab = "True number of harp seal pups",main = "Linear fit of uncorrected counts to true number of harp seal pups")
        abline(a = 0, b = lincorrHarp$coefficients,lwd = 3,col = "red")
      }

      var_aHarp = coef(summary(lincorrHarp))[1,2]; var_aHarp = var_aHarp^2
      var_bHarp = coef(summary(lincorrHarp))[2,2]; var_bHarp = var_bHarp^2
      var_uHarp = sd(lincorrHarp$residuals); var_uHarp = var_uHarp^2

      data$CountsHarpCorr <- lincorrHarp$coefficients*data$HarpFinalCounts

    }

    if("hood" %in% population){
      indcorrHooded = which(1-is.na(data$BBTrue)==1)			#Find which photos have been checked by both readers
      lincorrHooded = lm(data$BBTrue[indcorrHooded]~data$HoodedFinalCounts[indcorrHooded])
      if(plotFig){
        windows(width = 11,height = 7)
        plot(data$HoodedFinalCounts[indcorrHooded],data$BBTrue[indcorrHooded],xlab = "Hooded seal pup counts",ylab = "True number of hooded seal pups",main = "Linear fit of uncorrected counts to true number of hooded seal pups")
        abline(a = 0, b = lincorrHooded$coefficients,lwd = 3,col = "red")
      }

      var_aHooded = coef(summary(lincorrHooded))[1,2]; var_aHooded = var_aHooded^2
      var_bHooded = coef(summary(lincorrHooded))[2,2]; var_bHooded = var_bHooded^2
      var_uHooded = sd(lincorrHooded$residuals); var_uHooded = var_uHooded^2

      data$CountsHarpCorr <- lincorrHarp$coefficients*data$HarpFinalCounts
      data$CountsHoodedCorr <- lincorrHooded$coefficients*data$HoodedFinalCounts

    }

    #Get uncertainty associatet with reader errors
    correctedData = list()
    correctedData$newdata = data

    if("harp" %in% population){
      VmeasHarp = vmeasest(xycord$x,data$HarpFinalCounts,data$Area,data$Transect,Spacing,gap,var_bHarp,var_uHarp,var_aHarp,cov_ab = 0,intcpt = TRUE)
      correctedData$VmeasHarp = VmeasHarp
    }

    if("hood" %in% population){
      VmeasHooded = vmeasest(xycord$x,data$HoodedFinalCounts,data$Area,data$Transect,Spacing,gap,var_bHooded,var_uHooded,var_aHooded,cov_ab = 0,intcpt = TRUE)
      correctedData$VmeasHooded = VmeasHooded
    }

  }

}


#' Print pup production estimates using the Kingsley method to screen
#'
#' @param data The data set with photo counts
#' @param type Type of correction. Use 1 for separate correction for each reader (assumes two readers) and 2 for the same correction for both readers
#' @param population Specify which population or populations to estimate the uncertainty from
#' @param readers Specify the ID of each reader. Assumes 2 readers
#' @param transect Vector containing the transect each photo comes from
#' @param Spacing Transect width
#' @param gap Distance threshold for length if flown over open water
#' @param var_b
#' @param var_u
#' @param var_a
#' @param cov_ab
#' @param intcpt
#' @return
#' @keywords
#' @export
#' @examples
#' print2screen()

printKingsley2screen <- function(estimates,
                                 population = c("harp","hood"))
{

  #Print results to screen
  if("harp" %in% population){
    #Print results to screen
    cat("---------------------------------------------------------------------\n")
    cat("Estimates of harp seal pup production using the Kingsley method \n")
    cat(paste("Estimate without readers error correction: ",round(estimates$EstHarpFinal,digits = 0)," , SE = ",round(sqrt(estimates$VarTotalHarp),digits = 0),", CV = ",round(100*sqrt(estimates$VarTotalHarp)/estimates$EstHarpFinal,digits = 1),"%\n"))
    cat(paste("Estimate with readers error correction: ",round(estimates$EstHarpFinalCorr,digits = 0)," , SE = ",round(sqrt(estimates$VarTotalHarpCorr),digits = 0),", CV = ",round(100*sqrt(estimates$VarTotalHarpCorr)/estimates$EstHarpFinalCorr,digits = 1),"%\n"))
    cat("---------------------------------------------------------------------\n")
  }
  if("hood" %in% population){
    cat("---------------------------------------------------------------------\n")
    cat("Estimates of hooded seal pup production using the Kingsley method \n")
    cat(paste("Estimate without readers error correction: ",round(estimates$EstHoodedFinal,digits = 0)," , SE = ",round(sqrt(estimates$VarTotalHooded),digits = 0),", CV = ",round(100*sqrt(estimates$VarTotalHooded)/estimates$EstHarpFinal,digits = 1),"%\n"))
    cat(paste("Estimate with readers error correction: ",round(estimates$EstHoodedFinalCorr,digits = 0)," , SE = ",round(sqrt(estimates$VarTotalHoodedCorr),digits = 0),", CV = ",round(100*sqrt(estimates$VarTotalHoodedCorr)/estimates$EstHarpFinalCorr,digits = 1),"%\n"))
    cat("---------------------------------------------------------------------\n")
  }
}

#' Print pup production estimates using the GAM method to screen
#'
#' @param data The data set with photo counts
#' @param type Type of correction. Use 1 for separate correction for each reader (assumes two readers) and 2 for the same correction for both readers
#' @param population Specify which population or populations to estimate the uncertainty from
#' @param readers Specify the ID of each reader. Assumes 2 readers
#' @param transect Vector containing the transect each photo comes from
#' @param Spacing Transect width
#' @param gap Distance threshold for length if flown over open water
#' @param var_b
#' @param var_u
#' @param var_a
#' @param cov_ab
#' @param intcpt
#' @return
#' @keywords
#' @export
#' @examples
#' print2screen()

printGAM2screen <- function(estimates,
                            population = c("harp","hood"))
{

  #Print results to screen
  if("harp" %in% population){
    #Print results to screen
    cat("-------------------------------------------------------------")
    cat("Estimates of harp seal pup production using the GAM method \n")
    cat(paste("Estimate without readers error correction: ",round(estimates$EstHarpGAMFinal,digits = 0)," , SE = ",round(sqrt(estimates$VarTotalHarpGAM),digits = 0),", CV = ",round(100*sqrt(estimates$VarTotalHarpGAM)/estimates$EstHarpFinalGAM,digits = 1),"%\n"))
    cat(paste("Estimate with readers error correction: ",round(estimates$EstHarpGAMFinalCorr,digits = 0)," , SE = ",round(sqrt(estimates$VarTotalHarpGAMCorr),digits = 0),", CV = ",round(100*sqrt(estimates$VarTotalHarpGAMCorr)/estimates$EstHarpFinalGAMCorr,digits = 1),"%\n"))
    cat("-----------------------------------------------------------\n")
  }
  if("hood" %in% population){
    cat("-----------------------------------------------------------\n")
    cat("Estimates of hooded seal pup production using the GAM method \n")
    cat(paste("Estimate without readers error correction: ",round(estimates$EstHoodedGAMFinal,digits = 0)," , SE = ",round(sqrt(estimates$VarTotalHoodedGAM),digits = 0),", CV = ",round(100*sqrt(estimates$VarTotalHoodedGAM)/estimates$EstHoodedGAMFinal,digits = 1),"%\n"))
    cat(paste("Estimate with readers error correction: ",round(estimates$EstHoodedGAMFinalCorr,digits = 0)," , SE = ",round(sqrt(estimates$VarTotalHoodedGAMCorr),digits = 0),", CV = ",round(100*sqrt(estimates$VarTotalHoodedGAMCorr)/estimates$EstHoodedGAMFinalCorr,digits = 1),"%\n"))
    cat("-----------------------------------------------------------")
  }
}


#' Estimate the pup production using GAMs
#'
#' @param data The data set with photo counts
#' @param type Type of correction. Use 1 for separate correction for each reader (assumes two readers) and 2 for the same correction for both readers
#' @param population Specify which population or populations to estimate the uncertainty from
#' @param readers Specify the ID of each reader. Assumes 2 readers
#' @param transect Vector containing the transect each photo comes from
#' @param Spacing Transect width
#' @param gap Distance threshold for length if flown over open water
#' @param var_b
#' @param var_u
#' @param var_a
#' @param cov_ab
#' @param intcpt
#' @return
#' @keywords
#' @export
#' @examples
#' GAMestimate()

GAMestimate <- function(harpcounts = NA,
                        hoodedcounts = NA,
                        area,
                        xycord,
                        transect,
                        population = c("harp","hood"),
                        distr = "negbin",
                        Spacing = 3,
                        ds = 0.1,
                        plotMesh = TRUE,
                        grDev = TRUE)
{

  #data = dataList$data
  #xycord = dataList$xycord

  #Proportion zero observations
  if("harp" %in% population) cat(paste("Proportion zero counts harp seal pups: ",round(length(which(harpcounts == 0))/length(harpcounts),digits = 2)),"\n")
  if("hood" %in% population) cat(paste("Proportion zero counts hooded seal pups: ",round(length(which(hoodedcounts == 0))/length(hoodedcounts),digits = 2)),"\n")

  data = data.frame(log.area = log(area),
                    x = xycord$x,
                    y = xycord$y)
  if("harp" %in% population) data$harpcounts = harpcounts
  if("hood" %in% population) data$hoodedcounts = hoodedcounts
  #data$log.area <- log(area)
  #data$x <- xycord$x
  #data$y <- xycord$y

  #Use these for GAM estimates from uncorrected counts
  #HarpGam <- gam(HarpFinalCounts ~ s(xycord$x,xycord$y)+offset(log.area),family=negbin(c(0.1,10)),method=gam.method(gam="outer"),gamma = 1.4)
  if("harp" %in% population){
    cat("\n Fitting GAM to harp pup counts...\n")
    if("negbin" %in% distr) HarpGam <- bam(harpcounts ~ s(x,y)+offset(log.area),data = data,family=nb(theta = NULL,link = "log"),gamma = 1.4)

      #HarpGam <- gam(HarpFinalCounts ~ s(x,y)+offset(log.area),data = data,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)
    #HarpGam$family$getTheta() #Check estimate
    if("poisson" %in% distr) HarpGam <- gam(harpcounts ~ s(x,y)+offset(log.area),data = data,family="poisson",gamma = 1.4)
  }

  #HoodedGam <- gam(HoodedFinalCounts ~ s(xycord$x,xycord$y)+offset(log.area),family=negbin(c(0.1,10)),method=gam.method(gam="outer"),gamma = 1.4)
  if("hood" %in% population){
    cat("\n Fitting GAM to hooded pup counts...\n")

    if("negbin" %in% distr) HoodedGam <- bam(hoodedcounts ~ s(x,y)+offset(log.area),data = data,family=nb(theta = NULL,link = "log"),gamma = 1.4)

      #HoodedGam <- gam(HoodedFinalCounts ~ s(x,y)+offset(log.area),data = data,family=negbin(c(0.1,10)),optimizer = "perf",gamma = 1.4)

    #HoodedGam$family$getTheta() #Check estimate
    if("poisson" %in% distr) HoodedGam <- gam(hoodedcounts ~ s(x,y)+offset(log.area),data = data,family="poisson",gamma = 1.4)
  }

  cat("\n Creating area mesh grid...\n")
  meshList = areaMesh(xycord = xycord)
  meshvector = meshList$meshvector
  inmatrix = meshList$inmatrix
  XI = meshList$XI
  YI = meshList$YI
  dimX = dim(XI)

  dxnm = 1852			#Unit for nautical miles

  tr_label = sort(unique(transect));
  for (k in 1:length(tr_label)) {
    tr = tr_label[k]
    postr = which(transect == tr);
    x_tr = sort(xycord$x[postr])
    y_tr = sort(xycord$y[postr])

    xv_left = x_tr[1]-1/2*mean(sqrt(area))/dxnm
    xv_right = x_tr[length(x_tr)]+1/2*mean(sqrt(area))/dxnm

    yv_left = y_tr[1]-1/2*Spacing*1.0
    yv_right = y_tr[length(y_tr)]+1/2*Spacing*1.0

    indx <- which(YI >= xv_left & YI <= xv_right & XI >= yv_left & XI <= yv_right)

    inmatrix[indx] = 1

  }

  #plot this for checking the area selected
  if(plotMesh){
    #windows(height = 7,width = 9)
    if(grDev) graphDev(width = width,height = height)
    image(meshvector,meshvector,inmatrix,main = "Area to estimate density surface of")
  }

  area.index = which(inmatrix > 0)
  area.latlon = which(inmatrix > 0,arr.ind=TRUE)
  lon = meshvector[area.latlon[,1]]
  lat = meshvector[area.latlon[,2]]

  #Create new data frame
  sealdatap <- data.frame(x=lon,y=lat,log.area=0*lon)

  Nmesh <- length(meshvector)

  GAMestimates = list()

  #Predict the GAm model on the new surface - uncorrected counts
  if("harp" %in% population){
    cat("\n Estimating number of harp seal pups...\n")
    zzHarp <- array(NA,length(meshvector)^2)
    zzHarp[area.index] <- predict(HarpGam,sealdatap)
    NHarpGAM = sum(exp(zzHarp),na.rm=TRUE)*(ds*dxnm)^2

    #median=qnbinom(p=0.5,size = shapeNbHarps,mu = zzHarp)
    GAMestimates$NHarpGAM = NHarpGAM


    if(plotMesh){
      #windows(height = 7,width = 9)
      harpDF = data.frame(x = rep(meshvector,Nmesh),
                          y = rep(meshvector,each=Nmesh),
                          mean = (zzHarp))
      indNA = is.na(harpDF$mean)
      harpDF = harpDF[-which(indNA),]
      #harpDF$mean[is.na(harpDF$mean)] = 255
      if(grDev) graphDev(width = width,height = height)
      #image(meshvector,meshvector,matrix(zzHarp,dimX[1],dimX[2]),xlim=c(-30,30),ylim = c(-45,40),xlab = "Relative position (nm)",ylab = "Relative position (nm)",main = "Spacial distribution of harp seal pups")
      #contour(meshvector,meshvector,matrix(zzHarp,dimX[1],dimX[1]),add=TRUE)
      ff <- ggplot(harpDF,aes(x=x,y=y)) + geom_raster(aes(fill=mean)) + ggtitle("Spacial distribution of harp seal pups") + xlab("Relative position (nm)") +  ylab("Relative position (nm)")
      ff <- ff +  scale_fill_gradientn(colours=heat.colors(256,rev=T))
      ff <- ff + theme_bw() + theme(plot.title = element_text(size = 18,
                                                              face = "bold",
                                                              vjust = 5),axis.title.y = element_text(vjust = 5), axis.title.x = element_text(vjust = -5), axis.text=element_text(size=14),
                                    axis.title=element_text(size=16),
                                    legend.text = element_text(size=12),
                                    legend.title = element_text(size=16) ,plot.margin=unit(c(1,1,1.5,1.2),"cm"),panel.border = element_blank())
      print(ff)
    }

    XpHarp <- predict(HarpGam,newdata=sealdatap,type="lpmatrix")

    dimXp = dim(XpHarp)
    splinesHarp = array(0,dimXp[1])

    for(ell in 1:dimXp[1]){
      splinesHarp[ell]=XpHarp[ell,]%*%coef(HarpGam)
    }

    VarHarpGAM = exp(t(splinesHarp))%*%XpHarp%*%HarpGam$Vp%*%t(XpHarp)%*%exp(splinesHarp)*(ds*dxnm)^4
    SEHarpGAM = sqrt(VarHarpGAM)
    CVHarpGAM = SEHarpGAM/NHarpGAM*100

    GAMestimates$VarHarpGAM = VarHarpGAM
    GAMestimates$SEHarpGAM = SEHarpGAM
    GAMestimates$CVHarpGAM = CVHarpGAM

  }

  if("hood" %in% population){
    cat("\n Estimating number of hooded seal pups...\n")

    zzHooded <- array(NA,length(meshvector)^2)
    zzHooded[area.index] <- predict(HoodedGam,sealdatap)
    NHoodedGAM = sum(exp(zzHooded),na.rm=TRUE)*(ds*dxnm)^2

    GAMestimates$NHoodedGAM = NHoodedGAM

    if(plotMesh){
      hoodedDF = data.frame(x = rep(meshvector,Nmesh),
                          y = rep(meshvector,each=Nmesh),
                          mean = zzHooded)
      indNA = is.na(hoodedDF$mean)
      hoodedDF = hoodedDF[-which(indNA),]
      #windows(height = 7,width = 9)
      if(grDev) graphDev(width = width,height = height)
      #image(meshvector,meshvector,matrix(zzHooded,dimX[1],dimX[2]),xlim=c(-30,30),ylim = c(-45,40),xlab = "Relative position (nm)",ylab = "Relative position (nm)",main = "Spatial distribution of hooded seal pups")
      #contour(meshvector,meshvector,matrix(zzHooded,dimX[1],dimX[1]),add=TRUE)
      ff <- ggplot(hoodedDF,aes(x=x,y=y)) + geom_raster(aes(fill=mean)) + ggtitle("Spacial distribution of hooded seal pups") + xlab("Relative position (nm)") +  ylab("Relative position (nm)")
      ff <- ff +  scale_fill_gradientn(colours=heat.colors(256,rev=T))
      ff <- ff + theme_bw() + theme(plot.title = element_text(size = 18,
                                                              face = "bold",
                                                              vjust = 5),axis.title.y = element_text(vjust = 5), axis.title.x = element_text(vjust = -5), axis.text=element_text(size=14),
                                    axis.title=element_text(size=16),
                                    legend.text = element_text(size=12),
                                    legend.title = element_text(size=16) ,plot.margin=unit(c(1,1,1.5,1.2),"cm"),panel.border = element_blank())
      print(ff)
    }

    XpHooded <- predict(HoodedGam,newdata=sealdatap,type="lpmatrix")

    dimXp = dim(XpHooded)
    splinesHooded = array(0,dimXp[1])

    for(ell in 1:dimXp[1]){
      splinesHooded[ell]=XpHooded[ell,]%*%coef(HoodedGam)
    }

    VarHoodedGAM = exp(t(splinesHooded))%*%XpHooded%*%HoodedGam$Vp%*%t(XpHooded)%*%exp(splinesHooded)*(ds*dxnm)^4
    SEHoodedGAM = sqrt(VarHoodedGAM)
    CVHoodedGAM = SEHoodedGAM/NHoodedGAM*100

    GAMestimates$VarHoodedGAM = VarHoodedGAM
    GAMestimates$SEHoodedGAM = SEHoodedGAM
    GAMestimates$CVHoodedGAM = CVHoodedGAM

  }

  TotalAreaCovered = sum(inmatrix)*(ds*dxnm)^2/(1000*1000)

  GAMestimates$TotalAreaCovered = TotalAreaCovered
  return(GAMestimates)
}


#' Find the area of the Patch and create a mesh grid
#'
#' @param data The data set with photo counts
#' @param type Type of correction. Use 1 for separate correction for each reader (assumes two readers) and 2 for the same correction for both readers
#' @param population Specify which population or populations to estimate the uncertainty from
#' @param readers Specify the ID of each reader. Assumes 2 readers
#' @param transect Vector containing the transect each photo comes from
#' @param Spacing Transect width
#' @param gap Distance threshold for length if flown over open water
#' @param var_b
#' @param var_u
#' @param var_a
#' @param cov_ab
#' @param intcpt
#' @return
#' @keywords
#' @export
#' @examples
#' areaMesh()

areaMesh <- function(xycord,
                     ds = 0.1)
{

  #Find the correct area to create the mesh
  meshvalue = ceiling(max(c(abs(min(xycord$x)),abs(max(xycord$x)),abs(min(xycord$y)),abs(max(xycord$y))))) 																		#Check this!!!!
  meshvector = seq(-meshvalue,meshvalue,ds)
  XI = outer(meshvector*0,meshvector,FUN="+")
  YI = outer(meshvector,0*meshvector,FUN="+")
  dimX = dim(XI)
  inmatrix = matrix(0,dimX[1],dimX[2])
  empty = inmatrix

  meshReturn = list()
  meshReturn$meshvector = meshvector
  meshReturn$inmatrix = inmatrix
  meshReturn$XI = XI
  meshReturn$YI = YI
  return(meshReturn)

}


# Prephare graphical device for a given OS
# @param width Width of the graphical window
# @param height Height of the grapchical window
# @return The an OS dependent graphical window
graphDev = function(width = 7,height = 5) {

  system = Sys.info()
  if(system['sysname']=="Windows"){
    windows(width = 7,height = 5)
  }

  if(system['sysname']=="Linux"){
    X11(width = 7,height = 5)
  }

  if(system['sysname']=="Darwin"){
    quartz("",width = 7,height = 5)
  }
}


#' Estimate the temporal distribution of births
#'
#' Estimating the temporal distribution of births to correct for pups not born yet or have left the ice
#' @param data The data set with photo counts
#' @param type Type of correction. Use 1 for separate correction for each reader (assumes two readers) and 2 for the same correction for both readers
#' @param population Specify which population or populations to estimate the uncertainty from
#' @param readers Specify the ID of each reader. Assumes 2 readers
#' @param transect Vector containing the transect each photo comes from
#' @param Spacing Transect width
#' @param gap Distance threshold for length if flown over open water
#' @param var_b
#' @param var_u
#' @param var_a
#' @param cov_ab
#' @param intcpt
#' @return
#' @keywords
#' @export
#' @examples
#' birtDist()

birthDist <- function(harpfname = "HarpStages2012.txt",
                      harpLengthStages = c(2.4,4.42,11.39),
                      harpKappa = 12.4,
                      hoodedfname = "HoodStages2012.txt",
                      hoodedLengthStages = c(2,1,4),
                      hoodedKappa = 8.6,
                      survey = "WestIce2012",
                      datePhoto = 28,
                        transect,
                        population = c("harp","hood"),
                        distr = "negbin",
                        Spacing = 3,
                        ds = 0.1,
                        plotMesh = TRUE,
                        grDev = TRUE)
{


  if("harp" %in% population){

    filename = paste0("/data/",survey,"/",harpfname)

    data = read.table(filename,sep = "",header = TRUE)
    days = data$Date
    ndays = length(days)


    length_stage1 = harpLengthStages[1]
    length_stage2 = harpLengthStages[2]
    length_stage3 = harpLengthStages[3]

    kappa = harpKappa


    rho1 = length_stage1 / kappa;
    rho2 = length_stage2 / kappa;
    rho3 = length_stage3 / kappa;


    staging =  as.matrix(data[1:ndays,2:8])

    #Combine newborn/yellow and fat/grey stages
    staging[,2] = staging[,1] + staging[,2]
    staging[,4] = staging[,4] + staging[,5]

    #Remove the newborn, grey, and ragged stages
    #Only three stages are used in this model
    #Can be expanded
    staging = staging[,c(-1,-5,-6,-7)]

    xmin = 15
    xmax = 45



  }

  if("hood" %in% population){

    filename = paste0("data/",survey,"/",harpfname,"/")
    data = read.table(filename,sep = "",header = TRUE)
    days = data$Date
    ndays = length(days)



    #rho1 = 0.18
    #rho2 = 0.24
    #rho3 = 0.46

    kappa = hoodedKappa


    # These give a good fit and is used in Øigård et al 2013
    rho1 = hoodedLengthStages[1]/kappa
    rho2 = hoodedLengthStages[2]/kappa
    rho3 = hoodedLengthStages[3]/kappa

    staging =  as.matrix(data[1:ndays,2:5])

    #Combine newborn and thin stages
    staging[,2] = staging[,1] + staging[,2]
    staging = staging[,-1]

    xmin = 15
    xmax = 35

  }

  DimStaging = dim(staging)
  nstages = DimStaging[2]
  stageprop = staging/rowSums(staging)

  spacing = 1/24
  spacing = 0.05
  tau = seq(0,10,by = spacing)

  phi1 = dgamma(tau[-1],kappa,1/rho1)
  phi2 = dgamma(tau[-1],kappa,1/rho2)
  phi3 = dgamma(tau[-1],kappa,1/rho3)

  tau = seq(spacing,30,by = spacing)
  tm1 = seq(spacing,40,by = spacing)
  tm2 = seq(tau[2]+tm1[1],tau[length(tau)]+tm1[length(tm1)],length = (length(tm1)+length(phi1)-1))
  tm3 = seq(tau[2]+tm2[1],tau[length(tau)]+tm2[length(tm2)],length = (length(tm2)+length(phi2)-1))

  cphi1 = 1-pgamma(tau,kappa,1/rho1)
  cphi2 = 1-pgamma(tau,kappa,1/rho2)
  cphi3 = 1-pgamma(tau,kappa,1/rho3)

  tn1 = seq(tm1[1]+tau[1],tm1[length(tm1)]+tau[length(tau)], length = (length(tm1)+length(cphi1)-1))
  tn2 = seq(tm2[1]+tau[1],tm2[length(tm2)]+tau[length(tau)], length = (length(tm2)+length(cphi2)-1))
  tn3 = seq(tm3[1]+tau[1],tm3[length(tm3)]+tau[length(tau)], length = (length(tm3)+length(cphi3)-1))

  t_min = min(c(tn1[1],tn2[1],tn3[1]))-1
  t_max = max(c(tn1[length(tn1)],tn2[length(tn2)],tn3[length(tn3)]))+1
  t_tot = seq(t_min,t_max,by = spacing)
  Nt_tot = length(t_tot)
  #Not used any more
  #priors = matrix(NA,2,2)
  #priors[1,] = c(21,5); priors[2,] = c(1,1)

  datoind = round(mean(which((t_tot>datePhoto) & (t_tot<(datePhoto+1)))))

  data = list()
  data$ndays = ndays
  data$nstages = nstages
  data$spacing = spacing
  data$days = days
  data$staging = staging
  data$stageprop = stageprop
  #data$nphi = length(phi1)
  data$phi1=phi1
  data$phi2=phi2
  #data$ncphi=length(cphi1)
  data$cphi1=cphi1
  data$cphi2=cphi2
  data$cphi3=cphi3
  #data$ntm1=length(tm1)
  data$tm1 = tm1
  #data$ntn1=length(tn1)
  data$tn1 = tn1
  #data$ntn2=length(tn2)
  data$tn2 = tn2
  #data$ntn3 = length(tn3)
  data$tn3=tn3
  #data$nttot = length(t_tot)
  data$datoind = datoind
  data$ttot = t_tot

  parameters = list()
  parameters$pmub = 0
  parameters$logsigmab = log(1)


  #############################
  # TMB part
  #############################
  library(TMB)

  compile("birthDist.cpp","-O1 -g",DLLFLAGS="")
  dyn.load(dynlib("birthDist"))

  #load("birthidstData.RDat")

  obj <- MakeADFun(data,parameters,DLL="birthDist",checkParameterOrder = FALSE)

  obj$fn()
  obj$gr()
  system.time(opt <- nlminb(obj$par,obj$fn,obj$gr,control = list(eval.max = 1e6,maxit = 1e6)))

  rep<-sdreport(obj, getJointPrecision=TRUE)
  rep.matrix <- summary(rep)
  rep.rnames = rownames(rep.matrix)

  Report = obj$report()
  nn1 = Report$Nout[,1]
  nn2 = Report$Nout[,2]
  nn3 = Report$Nout[,3]
  nn = Report$Nout[,4]

  #Build data frame
  dfDist = data.frame(x = rep(t_tot,3), dist = NA,group = NA)
  dfDist$dist[1:Nt_tot] = nn1/nn
  dfDist$dist[(Nt_tot+1):(2*Nt_tot)] = nn2/nn
  dfDist$dist[((2*Nt_tot+1)):(3*Nt_tot)] = nn3/nn

  dfCum = data.frame(x = rep(t_tot,3),propIce = NA, group = NA)
  dfCum$propIce[1:Nt_tot] = nn1
  dfCum$propIce[(Nt_tot+1):(2*Nt_tot)] = nn1 + nn2
  dfCum$propIce[((2*Nt_tot+1)):(3*Nt_tot)] = nn1 + nn2 + nn3

  #dfPropCI = data.frame(x = t_tot,LL = (rep.matrix[3,1]-1.96*rep.matrix[3,2]), UL = (rep.matrix[3,1]+1.96*rep.matrix[3,2]))

  if("harp" %in% population){
    dfDist$group[1:Nt_tot] = "Newborn/Yellow"
    dfDist$group[(Nt_tot+1):(2*Nt_tot)] = "Thin white"
    dfDist$group[((2*Nt_tot+1)):(3*Nt_tot)] = "Fat white/Greycoat"

    dfCum$group[1:Nt_tot] = "Newborn/Yellow"
    dfCum$group[(Nt_tot+1):(2*Nt_tot)] = "Thin white"
    dfCum$group[((2*Nt_tot+1)):(3*Nt_tot)] = "Fat white/Greycoat"
  }

  if("hood" %in% population){
    dfDist$group[1:Nt_tot] = "Newborn and Thin"
    dfDist$group[(Nt_tot+1):(2*Nt_tot)] = "Fat"
    dfDist$group[((2*Nt_tot+1)):(3*Nt_tot)] = "Solitary"

    dfCum$group[1:Nt_tot] = "Newborn and Thin"
    dfCum$group[(Nt_tot+1):(2*Nt_tot)] = "Fat"
    dfCum$group[((2*Nt_tot+1)):(3*Nt_tot)] = "Solitary"

  }

  indNa = which(is.na(dfDist$dist))
  dfDist$dist[indNa] = 0

  dfObs = data.frame(x = days)
  dfObs$stage1 = staging[,1]/rowSums(staging)
  dfObs$stage2 = staging[,2]/rowSums(staging)
  dfObs$stage3 = staging[,3]/rowSums(staging)


  # windows(height = 7,width = 9)
  # par(mar=c(6,5,4,5),bg = "white")
  # plot(t_tot,nn1/nn,type = "l",col = colvcol[1],lwd = 4,xlim = c(xmin,xmax),ylim = c(0,1),xlab = "Days since 1. March 2012",ylab = "Proportion",cex.lab = 1.5,cex.main = 1.5,bty = "l")
  # lines(t_tot,nn2/nn,col = colvcol[2],lwd = 4)
  # lines(t_tot,nn3/nn,col= colvcol[3],lwd = 4)
  # lines(days,staging[,1]/rowSums(staging),type = "p",col = colvcol[1],pch = 15,lwd = 5)
  # lines(days,staging[,2]/rowSums(staging),type = "p",col = colvcol[2],pch = 19,lwd = 5)
  # lines(days,staging[,3]/rowSums(staging),type = "p",col = colvcol[3],pch = 17,lwd = 5)

  library(ggplot2)

  theCols <- RColorBrewer::brewer.pal(3, 'Dark2')

  windows(height = 7,width = 9)
  pl <- ggplot() +
    geom_line(data = dfDist ,aes(x=x,y=dist,group = group,color = group),
              size = 1.3,
              linetype = 1) +
    coord_cartesian(xlim = c(xmin,xmax),ylim=c(0,1))
    #xlim(xmin,xmax) + ylim(0,1)

  pl <- pl + geom_point(data = dfObs, aes(x = days, y = stage1),size = 2,color = theCols[2],shape=15) +
    geom_point(data = dfObs, aes(x = days, y = stage2),size = 2,col = theCols[3],shape=16) +
    geom_point(data = dfObs, aes(x = days, y = stage3),size = 2,col = theCols[1],shape=17) +
    geom_vline(xintercept = datePhoto, linetype="dashed",
               color = "lightgrey", size=0.8)


  pl <- pl + theme_classic() +
    theme(text = element_text(size=20), plot.margin = unit(c(1,2,1,1), "cm"), axis.text.y = element_text(angle = 90,margin = margin(t = 0, r = 20, b = 0, l = 30),vjust = 0.5),axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),vjust = 1),
                   legend.title = element_blank(),
                   legend.position = "top") +
    xlab("Days since 1. March") + ylab("Proportion") +
    scale_fill_manual(values = c(theCols[1], theCols[2], theCols[3])) +
    scale_colour_manual(values = c(theCols[1], theCols[2], theCols[3]))

  pl


  windows(height = 7,width = 9)
  pl <- ggplot() +
    geom_line(data = dfCum ,aes(x=x,y=propIce,color = group),
              size = 1.3,
              linetype = 1,
              ) +
   # geom_ribbon(data=dfPropCI,aes(x = x, ymin=LL,ymax=UL),
    #            alpha=0.3,
    #            fill = "lightgrey") +
    geom_vline(xintercept = datePhoto, linetype="dashed",
               color = "lightgrey", size=0.8) +
    coord_cartesian(xlim = c(xmin,xmax),ylim=c(0,1))

    # geom_line(data = dfCum ,aes(x=x,y=nn12),
    #           size = 1.3,
    #           linetype = 1,
    #           color = "blue") +
    # geom_line(data = dfCum ,aes(x=x,y=nn123),
    #           size = 1.3,
    #           linetype = 1,
    #           color = "green") +


  pl <- pl + theme_classic() +
    theme(text = element_text(size=20), plot.margin = unit(c(1,2,1,1), "cm"), axis.text.y = element_text(angle = 90,margin = margin(t = 0, r = 20, b = 0, l = 30),vjust = 0.5),axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),vjust = 1),
          legend.title = element_blank(),
          legend.position = "top") +
    xlab("Days since 1. March") + ylab("Proportion seals on ice")

pl


returnList = list()
returnList$mub = rep.matrix[4,1]
returnList$mubSD = rep.matrix[4,2]
returnList$sigma = rep.matrix[5,1]
returnList$sigmaSD = rep.matrix[5,2]
returnList$PropIce = rep.matrix[3,1]
returnList$PropIceSD = rep.matrix[3,2]
returnList$xaxis = t_tot
returnList$nn1 = nn1
returnList$nn2 = nn2
returnList$nn3 = nn3
returnList$nn = nn

return(returnList)

}
