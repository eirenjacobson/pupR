
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
#' lb3xy()

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
#' @param Twdt Transect width
#' @param gap Distance threshold for length if flown over open water
#' @return
#' @keywords
#' @export
#' @examples
#' kingsley()

kingsley <- function(x,count,area,transect,Twdt,gap)
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


    N_new = (1852*Twdt)*sum(disttr*count_ref);
    V_new = (1852*Twdt*(1852*Twdt-sum(area)/sum(disttr)))/(2*(N_tr-1))*sum(disttr^2)*sum(diff(count_ref)^2);
    V_new = (1852*Twdt*N_tr*(1852*Twdt-sum(area)/sum(disttr)))/(2*(N_tr-1))*sum(diff(disttr*count_ref)^2)
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

		yv_left = y_tr[1]-1/2*spacing*1.0
		yv_right = y_tr[length(y_tr)]+1/2*spacing*1.0

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
#' @param Twdt Transect width
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
                     Twdt,
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
		Twdt = 1852*Twdt;

		if (intcpt == TRUE) {
			V_meas = Twdt^2*(sum(term1)^2*var_a+2*cov_ab*sum(term1)*sum(term2)+sum(term2)^2*var_b+sum(term3)*var_u)
			cat("Intercept used")} else {V_meas = Twdt^2*(sum(term2)^2*var_b+sum(term3)*var_u)
				cat("No intercept used")}

		#V_meas = T*sum(Fj*term4)
		return(V_meas)

	}


