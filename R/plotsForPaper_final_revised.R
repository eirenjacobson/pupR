### Plots for the paper


#### Plotting seal data ####

library(SDMTools)
require(ggplot2)
library(reshape2)
library(tidyr)
library(scales)
library(xtable)
library(rgdal)
library(spatstat)
library(INLA)
library(fields)
library(sp)
library(boot)


rm(list=ls())

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


paperFigureFolder <- "/disk/home/jullum/Prosjekter_lokalt/Fritanek_Point/CoxProcessesINLA_new/Latex/2ndRevisionJRSSC/figures/"


#seals=read.csv(file="/nr/project/stat/PointProcess/Data/WestIce2012.csv")
seals <- readRDS("/nr/project/stat/PointProcess/Data/Seals/OriginalSealsKm.rds")
load("/nr/project/stat/PointProcess/Data/Seals/domain.count.inside.RData")


#> latidude
#> longitude
#> x		transformer til cartesiske koordinater
#> y		samme som ovenfor
#> harps		antall harp pups/grønlandsselunger i et bilde
#> hoods 		antall hooded pups/klappmyssunger i et bilde
#> area		areal dekt av et bilde (i m^2)
#> length		lengden av et bilde (i m)
#> lengthNm	lengden av et bilde (i nautiske mil)
#> width		bredden av et bilde (i m)
#> widthNm	bredden av et bilde i (nm)

library(SealCountINLA)

DefToGlob(INLAPPSealsPoissonReg)
sealPhotoDataFile = "/nr/project/stat/PointProcess/Data/Seals/OriginalSealsKm.rds"#"M:\\PointProcess/Data/Seals/OriginalSealsKm.rds"
## Loading the seal photo data
if (tools::file_ext(sealPhotoDataFile)=="rds"){
  seals <- readRDS(file=sealPhotoDataFile) # RDS-file
} else {
  load(file=sealPhotoDataFile) # RData, with a seal object.
}

# Creating a list where all important data to be used later are stored explanatory names
dataList <- list()
dataList$org <- list()
dataList$org$noPhoto <- dim(seals)[1]   # Total number of photo taken (with or without seals)
dataList$org$coordPhotoX <- seals$xkm   # X coordinate of center of each photo, given in kilometers
dataList$org$coordPhotoY <- seals$ykm   # Y coordinate of center of each photo,  given in kilometers
dataList$org$photoWidth <- seals$lengthkm # The width (X-direction) of each of the photos, given in kilometers
dataList$org$photoHeight <- seals$widthkm # The height (Y-direction) of each of the photos, given in kilometers

rectangleCentersX = dataList$org$coordPhotoX
rectangleCentersY = dataList$org$coordPhotoY
rectangleWidth = dataList$org$photoWidth
rectangleHeight = dataList$org$photoHeight

sealTransectDataFile = "/nr/project/stat/PointProcess/Data/Seals/OigardTablesTransformed.rds"#"M:\\PointProcess/Data/Seals/OigardTablesTransformed.rds"
transectData <- readRDS(sealTransectDataFile)

dataList$transect <- list()
dataList$transect$noTransects <- nrow(transectData)
dataList$transect$transectStartCoordX <- transectData$x.start
dataList$transect$transectStartCoordY <- transectData$y.start
dataList$transect$transectEndCoordX <- transectData$x.end
dataList$transect$transectEndCoordY <- transectData$y.end

theseTransInCountDomain <- 1:dataList$transect$noTransects

countingDomain <- CreateCountDomainPolygon(transectStartCoordX = dataList$transect$transectStartCoordX,
                                           transectStartCoordY = dataList$transect$transectStartCoordY,
                                           transectEndCoordX = dataList$transect$transectEndCoordX,
                                           transectEndCoordY = dataList$transect$transectEndCoordY,
                                           coordPhotoX = dataList$org$coordPhotoX,
                                           coordPhotoY = dataList$org$coordPhotoY,
                                           photoWidth = dataList$org$photoWidth,
                                           photoHeight = dataList$org$photoHeight,
                                           transectYSpan = 1.5*1.852,
                                           theseTransInCountDomain=theseTransInCountDomain)
countingDomain$LineType <- "Whelping region"


seals$harps.pos=seals$harps
seals$harps.pos[seals$harps==0]=NA
seals$hooded.pos=seals$hooded
seals$hooded.pos[seals$hooded==0]=NA

seals2=gather(seals,value="no",key="type",harps.pos,hooded.pos)

seals2$type = factor(seals2$type,
                     levels = c("harps.pos","hooded.pos"),
                     labels = c("Harps","Hooded"))

### First plot: Harps and hooded seals with different colors,
gg_color_hue <- function(n) {hcl(h=seq(15, 375, length=n+1), l=65, c=100)[1:n]}


g1=ggplot(seals2,aes(x=xkm,y=ykm,group=1))
g1 <- g1+
  geom_point(size=0.3,shape=15,color="grey")+
  geom_point(aes(size=no,color=type),alpha=I(0.2))+
  scale_size_continuous(name="Seal count",range=c(1,8))+
  scale_color_discrete(name="Seal type")
g1 <- g1 + theme_bw()
g1 <- g1 + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16))
g1 <- g1 + labs(x = "x (km)",
                y = "y (km)")
g1

g2 <- g1 + geom_path(data = countingDomain, mapping=aes(x=x,y=y,linetype=LineType),color="black",size=1)
g2 <- g2 + geom_path(data = data.frame(x=0,y=0,LineType="Transect"), mapping=aes(x=x,y=y,linetype=LineType),color="grey")
g2 <- g2+ scale_linetype_manual("Lines",values=c(1,1),guide = guide_legend(override.aes = list(colour=c("grey","black"))))
g2

pdf(file=file.path(paperFigureFolder,"seal_data.pdf"),width=7,height=6)
g2
dev.off()


#### ##########


### Testing Tor Arnes GAM-script


seals <- readRDS("/nr/project/stat/PointProcess/Data/Seals/OriginalSealsKm.rds")

sealTransectDataFile = "/nr/project/stat/PointProcess/Data/Seals/OigardTablesTransformed.rds"#"M:\\PointProcess/Data/Seals/OigardTablesTransformed.rds"
transectData <- readRDS(sealTransectDataFile)

trans=rep(NA,nrow(seals))
for(i in 1:nrow(seals)){
  trans[i]=which.min((seals$ykm[i]-(transectData$y.start+transectData$y.end)/2)^2)
}



#### Functions

#Transfers longitude-latitude to x and y in nautical miles.
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


##### KINGSLEY  ##########
Spacing = 3 			#Distance between transects in NM
gap = 1

(KingsleyHarps = kingsley(seals$x,seals$harps,seals$area,trans,Spacing,gap))
(KinglseyHooded = kingsley(seals$x,seals$hooded,seals$area,trans,Spacing,gap))


######### NEW #####

#### Script for fitting the Bayesian homogenous Poisson model

library(SealCountINLA)

#### First just getting the data on the right form
sealPhotoDataFile = "/nr/project/stat/PointProcess/Data/Seals/OriginalSealsKm.rds"
sealTransectDataFile = "/nr/project/stat/PointProcess/Data/Seals/OigardTablesTransformed.rds"

seals <- readRDS(file=sealPhotoDataFile) # RDS-file

# Creating a list where all important data to be used later are stored explanatory names
dataList <- list()
dataList$org <- list()
dataList$org$noPhoto <- dim(seals)[1]   # Total number of photo taken (with or without seals)
dataList$org$coordPhotoX <- seals$xkm   # X coordinate of center of each photo, given in kilometers
dataList$org$coordPhotoY <- seals$ykm   # Y coordinate of center of each photo,  given in kilometers
dataList$org$photoWidth <- seals$lengthkm # The width (X-direction) of each of the photos, given in kilometers
dataList$org$photoHeight <- seals$widthkm # The height (Y-direction) of each of the photos, given in kilometers

## Loading the transect seal data
transectData <- readRDS(sealTransectDataFile)

dataList$transect <- list()
dataList$transect$noTransects <- nrow(transectData)
dataList$transect$transectStartCoordX <- transectData$x.start
dataList$transect$transectStartCoordY <- transectData$y.start
dataList$transect$transectEndCoordX <- transectData$x.end
dataList$transect$transectEndCoordY <- transectData$y.end


#### Prepare running for the full counting domain ####
countingDomain <- CreateCountDomainPolygon(transectStartCoordX = dataList$transect$transectStartCoordX,
                                           transectStartCoordY = dataList$transect$transectStartCoordY,
                                           transectEndCoordX = dataList$transect$transectEndCoordX,
                                           transectEndCoordY = dataList$transect$transectEndCoordY,
                                           coordPhotoX = dataList$org$coordPhotoX,
                                           coordPhotoY = dataList$org$coordPhotoY,
                                           photoWidth = dataList$org$photoWidth,
                                           photoHeight = dataList$org$photoHeight,
                                           transectYSpan = 1.5*1.852,
                                           theseTransInCountDomain=1:27)
areaCountDomain <- rgeos::gArea(coo2sp(countingDomain)) # Should be close to 4569 which Tor-Arne uses in his papers


#### The data #####

omegaSize <- areaCountDomain
observationsHooded <- seals$hooded
observationsHarps <- seals$harps
ASize <- dataList$org$photoHeight*dataList$org$photoWidth
####################


#### Specifying the priors ####
# Using the Kingsley estimator to specify the mean of the prior
meanPriorHooded <- 10#10928/omegaSize)
meanPriorHarps <- 10#85968/omegaSize)

# Approx 10 times Kingsley estimated variance -- might use an even more vague prior
variancePriorHooded <- 10 #10*10^6/omegaSize^2
variancePriorHarps <- 10#10*10^8/omegaSize^2

(aPriorHooded <- (meanPriorHooded^2)/variancePriorHooded)
(aPriorHarps <- (meanPriorHarps^2)/variancePriorHarps)

(bPriorHooded <- (meanPriorHooded)/variancePriorHooded)
(bPriorHarps <- (meanPriorHarps)/variancePriorHarps)

#########

(aPosteriorHooded <- aPriorHooded + sum(observationsHooded))
(aPosteriorHarps <- aPriorHarps + sum(observationsHarps))

(bPosteriorHooded <- bPriorHooded + sum(ASize))
(bPosteriorHarps <- bPriorHarps + sum(ASize))

meanPosteriorHooded <- aPosteriorHooded/bPosteriorHooded
meanPosteriorHarps <- aPosteriorHarps/bPosteriorHarps

varPosteriorHooded <- aPosteriorHooded/bPosteriorHooded^2
varPosteriorHarps <- aPosteriorHarps/bPosteriorHarps^2

##### This gives posterior predictice distribution for the complete area equal to

shapeNbHooded <- aPosteriorHooded # Equal to size in dbinom
shapeNbHarps <- aPosteriorHarps # Equal to size in dbinom

muNbHooded <- omegaSize*aPosteriorHooded/bPosteriorHooded
muNbHarps <- omegaSize*aPosteriorHarps/bPosteriorHarps

#pHooded <- muNbHooded/(shapeNbHooded+muNbHooded)
#pHooded*(shapeNbHarps-1)/(1-pHooded)

intensityHooded <- aPosteriorHooded/bPosteriorHooded
intensityHarps <- aPosteriorHarps/bPosteriorHarps

set.seed(123)
evalFull <- 1:(2*10^5)
posteriorDist <- dnbinom(x=evalFull,size = shapeNbHarps,mu = muNbHarps)
samp <- rnbinom(n = 10^7,size = shapeNbHarps,mu = muNbHarps)

HomoPoisHarps <- list()

HomoPoisHarps$postDens <- list(x = evalFull,
                               y = posteriorDist)

HomoPoisHarps$abundance.est = c(mean=muNbHarps,
                                median=qnbinom(p=0.5,size = shapeNbHarps,mu = muNbHarps),
                                mode=evalFull[which.max(posteriorDist)],
                                IQR = diff(qnbinom(p=c(0.25,0.75),size = shapeNbHarps,mu = muNbHarps)),
                                CI = qnbinom(p=c(0.025,0.975),size = shapeNbHarps,mu = muNbHarps))

set.seed(123)
evalFull <- 1:(2*10^4)
posteriorDist <- dnbinom(x=evalFull,size = shapeNbHooded,mu = muNbHooded)
samp <- rnbinom(n = 10^7,size = shapeNbHooded,mu = muNbHooded)

HomoPoisHooded <- list()

HomoPoisHooded$postDens <- list(x = evalFull,
                               y = posteriorDist)

# HomoPoisHooded$abundance.est = c(mean=muNbHooded,
#                                 median=median(samp),
#                                 mode=Mode(samp),
#                                 IQR = diff(quantile(samp,c(0.25,0.75))),
#                                 CI = quantile(samp,c(0.025,0.975)))

HomoPoisHooded$abundance.est = c(mean=muNbHooded,
                                 median=qnbinom(p=0.5,size = shapeNbHooded,mu = muNbHooded),
                                 mode=evalFull[which.max(posteriorDist)],
                                 IQR = diff(qnbinom(p=c(0.25,0.75),size = shapeNbHooded,mu = muNbHooded)),
                                 CI = qnbinom(p=c(0.025,0.975),size = shapeNbHooded,mu = muNbHooded))


######### NEW ENDING #####




#################


resBaseFolder <- "/home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results"




#### First some plots showing the total posterior predicitve distributions for the various methods

res.list <- list()

RdataREGHoodedPath <- file.path(resBaseFolder,"finalHooded/REG/runType=PoissonReg seal=hooded cov=linearAndSpatial extra.iid=FALSE samp=10000 comment=SpatialAsFixed_TransectAsCountingDomain_140/output.RData")
load(RdataREGHoodedPath)

res.list$REGHooded <- list()
res.list$REGHooded$x <- gridList$gridvalX
res.list$REGHooded$y <- gridList$gridvalY
res.list$REGHooded$mean <- finalResList$mean.field.domain.samp
res.list$REGHooded$sd <- finalResList$sd.field.domain.samp

res.list$REGHooded$params <- c(finalResList$intercept.mean,finalResList$covariate.mean,finalResList$spatialX.mean,finalResList$spatialY.mean,finalResList$spatialXY.mean,finalResList$mean.range.param)
names(res.list$REGHooded$params) = c("intercept","covariate","spatialX","spatialY","spatialXY","range")
res.list$REGHooded$abundance.est = c(mean=finalResList$posteriorMean,
                                     median=finalResList$posteriorMedian,
                                     mode=finalResList$posteriorMode,
                                     IQR = finalResList$posteriorIQR,
                                     CI = finalResList$posteriorCI)

res.list$REGHooded$postDens <- list(x = finalResList$posteriorEvalPoints,
                                    y = finalResList$posteriorDist)



RdataGAMHoodedPath <- file.path(resBaseFolder,"finalHooded/GAM/runType=GAM seal=hooded cov=linearAndSpatial comment=Linear_199/output.RData")
load(RdataGAMHoodedPath)

res.list$GAMHooded <- list()
res.list$GAMHooded$x <- gridList$gridvalX
res.list$GAMHooded$y <- gridList$gridvalY
res.list$GAMHooded$mean <- finalResList$mean.field.domain.samp
res.list$GAMHooded$sd <- finalResList$sd.field.domain.samp

res.list$GAMHooded$fixed.effects <- c(finalResList$intercept.mean,finalResList$covariate.mean,finalResList$spatialX.mean,finalResList$spatialY.mean,finalResList$spatialXY.mean)
names(res.list$GAMHooded$fixed.effects) = c("intercept","covariate","spatialX","spatialY","spatialXY")
res.list$GAMHooded$abundance.est = c(mean=finalResList$posteriorMean,
                                     median=finalResList$posteriorMedian,
                                     mode=finalResList$posteriorMode,
                                     IQR = finalResList$posteriorIQR,
                                     CI = finalResList$posteriorCI)

res.list$GAMHooded$postDens <- list(x = finalResList$posteriorEvalPoints,
                                    y = finalResList$posteriorDist)

RdataGAMPoissonHoodedPath <- file.path(resBaseFolder,"finalHooded/GAMPoisson/runType=GAM seal=hooded cov=linearAndSpatial comment=Linear_223/output.RData")
load(RdataGAMPoissonHoodedPath)

res.list$GAMPoissonHooded <- list()
res.list$GAMPoissonHooded$x <- gridList$gridvalX
res.list$GAMPoissonHooded$y <- gridList$gridvalY
res.list$GAMPoissonHooded$mean <- finalResList$mean.field.domain.samp
res.list$GAMPoissonHooded$sd <- finalResList$sd.field.domain.samp
res.list$GAMPoissonHooded$fixed.effects <- c(finalResList$intercept.mean,finalResList$covariate.mean,finalResList$spatialX.mean,finalResList$spatialY.mean,finalResList$spatialXY.mean)
names(res.list$GAMPoissonHooded$fixed.effects) = c("intercept","covariate","spatialX","spatialY","spatialXY")
res.list$GAMPoissonHooded$abundance.est = c(mean=finalResList$posteriorMean,
                                     median=finalResList$posteriorMedian,
                                     mode=finalResList$posteriorMode,
                                     IQR = finalResList$posteriorIQR,
                                     CI = finalResList$posteriorCI)

res.list$GAMPoissonHooded$postDens <- list(x = finalResList$posteriorEvalPoints,
                                           y = finalResList$posteriorDist)

res.list$HomoPoisHooded <- HomoPoisHooded

#### Adding Kingsley to the res.list as well

res.list$KingsleyHooded$abundance.est <- c(mean=KinglseyHooded$N,
                                     median=KinglseyHooded$N,
                                     mode=KinglseyHooded$N,
                                     IQR = 2*qnorm(0.75)*KinglseyHooded$SE,
                                     CI = c(KinglseyHooded$N-qnorm(0.975)*KinglseyHooded$SE,KinglseyHooded$N+qnorm(0.975)*KinglseyHooded$SE))



#### For simplicity we sample from the posteriors to better match the histogram function of ggplot

####
break1 <- 8000
break2 <- 16000
nobreaks <- 6
logbreaks <-round(10^((log10(break2)-log10(break1))*seq(0,1,length.out=nobreaks)+log10(break1)),-2)

set.seed(123)
dens.df <- data.frame(REG = sample(x = res.list$REGHooded$postDens$x,size = 10^6,prob = res.list$REGHooded$postDens$y,replace = T),
                      GAM = sample(x = res.list$GAMHooded$postDens$x,size = 10^6,prob = res.list$GAMHooded$postDens$y,replace = T),
                      GAMPoisson = sample(x = res.list$GAMPoissonHooded$postDens$x,size = 10^6,prob = res.list$GAMPoissonHooded$postDens$y,replace = T),
                      HomoPois = sample(x = res.list$HomoPoisHooded$postDens$x,size = 10^6,prob = res.list$HomoPoisHooded$postDens$y,replace = T))
dens.df.melt <- melt(dens.df,variable.name = "Model")
pp = ggplot(dens.df.melt,aes(x=value,fill=Model)) + geom_histogram(aes(y=..count../sum(..count..)),alpha=0.3,bins=150,position="identity") +
  coord_cartesian(xlim = c(break1, break2)) +
  xlab("Hooded pup count") +
  ylab("Probability")
pp = pp + geom_errorbarh(data=data.frame(KinglseyHooded,Model="GAM",Mod2="Kingsley"),mapping=aes(x=N,y=-0.0025,xmax=N+1.96*SE,xmin=N-1.96*SE,color=Mod2),
                         height = .0025,size=1)
pp = pp + geom_point(data.frame(KinglseyHooded,Model="GAM",Mod2="Kingsley"),mapping=aes(x=N,y=-0.0025,color=Mod2),size=2.5)
#pp = pp + geom_point(data.frame(KinglseyHooded,Model="GAM"),mapping=aes(x=N,y=-0.005),size=2.5)


pp <- pp + theme_bw()
pp <- pp + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16))
pp = pp +
  scale_color_manual(name="",values=gg_color_hue(4)[4],guide=guide_legend(order=2)) +
  scale_fill_discrete(guide=guide_legend(order=1,override.aes = list(color=NA,alpha=0.3)),
                      labels=c("LGCP","GAM NB","GAM Po", "Hom Po"))

pp
#pdf(file=file.path(paperFigureFolder,"postPred_hooded3.pdf"),width=12,height=6)
#pp
#dev.off()
pdf(file=file.path(paperFigureFolder,"postPred_hooded4.pdf"),width=12,height=6)
pp
dev.off()


pplog <- pp + scale_x_log10(breaks = logbreaks,
                            labels = logbreaks)
pplog

#pdf(file=file.path(paperFigureFolder,"postPred_hooded_log3.pdf"),width=12,height=6)
#pplog
#dev.off()

#### Mean and sd of latent field

no.x <- length(res.list$REGHooded$x)
no.y <- length(res.list$REGHooded$y)


latentfield.df <- data.frame(x=rep(res.list$REGHooded$x,no.y),
                             y=rep(res.list$REGHooded$y,each=no.x),
                             mean=c(res.list$REGHooded$mean),
                             sd = c(res.list$REGHooded$sd))
latentfield.df <- latentfield.df[complete.cases(latentfield.df),]

ff <- ggplot(latentfield.df,aes(x=x,y=y)) + geom_raster(aes(fill=mean)) +   xlab("x (km)") +  ylab("y (km)")
ff <- ff + theme_bw()
ff <- ff + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16))
ff <- ff +  scale_fill_gradientn(colours=heat.colors(256,rev=T))
ff
pdf(file=file.path(paperFigureFolder,"mean_latent_field_Hooded.pdf"),width=6,height=6)
ff
dev.off()

ff <- ggplot(latentfield.df,aes(x=x,y=y)) + geom_raster(aes(fill=sd)) +   xlab("x (km)") +  ylab("y (km)")
ff <- ff + theme_bw()
ff <- ff + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16))
ff <- ff +  scale_fill_gradientn(colours=heat.colors(256,rev=T))
ff
pdf(file=file.path(paperFigureFolder,"sd_latent_field_Hooded.pdf"),width=6,height=6)
ff
dev.off()





#### Konverterer resultattabeller

#####

abundanceDf.Hooded <- as.data.frame(rbind(res.list$REGHooded$abundance.est,res.list$GAMHooded$abundance.est,
                                          res.list$GAMPoissonHooded$abundance.est,res.list$HomoPoisHooded$abundance.est,
                                          res.list$KingsleyHooded$abundance.est))
rownames(abundanceDf.Hooded) <- c("LGCP","GAM NB","GAM Po", "Hom Po", "Kingsley")
colnames(abundanceDf.Hooded) <- c("mean","median","mode","IQR","0.025-quantile","0.975-quantile")


res0 <- xtable(abundanceDf.Hooded,digits = 0,
              align=c("l","c","c","c","c","c","c"),
              caption = "Summary tables for the posterior predictive distributions for the total area counts of hooded seals.",
              label = "tab:summary_hooded")

print.xtable(res0,caption.placement = "top")

# % latex table generated in R 3.6.1 by xtable 1.8-4 package
# % Mon Nov 11 10:39:40 2019
# \begin{table}[ht]
# \centering
# \caption{Summary tables for the posterior predictive distributions for the total area counts of hooded seals.}
# \label{tab:summary_hooded}
# \begin{tabular}{lcccccc}
# \hline
# & mean & median & mode & IQR & 0.025-quantile & 0.975-quantile \\
# \hline
# LGCP & 11649 & 11503 & 11472 & 1699 & 9472 & 14741 \\
# GAM NB & 11178 & 11157 & 11093 & 807 & 10075 & 12395 \\
# GAM Po & 11296 & 11292 & 11093 & 572 & 10467 & 12147 \\
# Hom Po & 11494 & 11489 & 11479 & 571 & 10678 & 12338 \\
# Kingsley & 10928 & 10928 & 10928 & 1969 & 8067 & 13789 \\
# \hline
# \end{tabular}
# \end{table}
#


print(res.list$REGHooded$params)
#> print(res.list$REGHooded$params)
#intercept   covariate    spatialX    spatialY   spatialXY       range
#-1.37152575  9.07356067  0.07472763 -0.04865848 -0.01543162  3.63154806

### HERE !

resultList <- readRDS(file = paste0(resBaseFolder,"/resultListnew_with_BayesHomoPois.rds"))


## Hooded, Photo level

tab <- resultList$finalHooded$CVFold$CRPS.photo
theseMod <- c("REGSpatialCVFold10","GAMSpatialCVFold10","GAMPoissonSpatialCVFold10","BayesHomoPois")
useNames <- c("LGCP","GAM NB","GAM Po","Hom Po")
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.photo.hooded <- data.frame(a = paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$lower,2),nsmall=2),", ",format(round(subtab$upper,2),nsmall=2),")",sep=""))
rownames(df.photo.hooded) <- useNames
colnames(df.photo.hooded) <- "CRPS"


tab <- -resultList$finalHooded$CVFold$logScore.photo # Reverse the sign here
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.photo.hooded$b <- paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$upper,2),nsmall=2),", ",format(round(subtab$lower,2),nsmall=2),")",sep="")
colnames(df.photo.hooded)[2] <- "logScore"


tab <- resultList$finalHooded$LeaveOut$CRPS.photo
theseMod <- c("REGSpatialLeaveOutAsCV","GAMSpatialLeaveOutAsCV","GAMPoissonSpatialLeaveOutAsCV","BayesHomoPois") # Needs to be in same order
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.photo.hooded$c <-  paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$lower,2),nsmall=2),", ",format(round(subtab$upper,2),nsmall=2),")",sep="")
colnames(df.photo.hooded)[3] <- "CRPS"# (5%, 95%) "


tab <- -resultList$finalHooded$LeaveOut$logScore.photo # Reverse the sign here
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.photo.hooded$d <- paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$upper,2),nsmall=2),", ",format(round(subtab$lower,2),nsmall=2),")",sep="")
colnames(df.photo.hooded)[4] <- "logScore"# (5%, 95%) "


df.photo.hooded
res <- xtable(df.photo.hooded,
              align=c("l","c","c","c","c"),
              caption = "Validation results on photo level. Values in parenthesis are shows the lower and upper bounds of 90% bootstrapped confidence intervals for the performance measure. Cells shown in italics are best (smallest) per column. Those which are significantly smaller than the others (defined as having non-overlapping confidence intervals) are also bolded.",
              label = "tab:Val_hooded_photo")

print.xtable(res,caption.placement = "top")

# % latex table generated in R 3.6.1 by xtable 1.8-4 package
# % Mon Nov 11 09:58:47 2019
# \begin{table}[ht]
# \centering
# \caption{Validation results on photo level. Values in parenthesis are shows the lower and upper bounds of 90% bootstrapped confidence intervals for the performance measure. Cells shown in italics are best (smallest) per column. Those which are significantly smaller than the others (defined as having non-overlapping confidence intervals) are also bolded.}
# \label{tab:Val_hooded_photo}
# \begin{tabular}{lcccc}
# \hline
# & CRPS & logScore & CRPS & logScore \\
# \hline
# LGCP & 0.18 (0.16, 0.19) & 0.47 (0.44, 0.49) & 0.22 (0.20, 0.25) & 0.54 (0.51, 0.57) \\
# GAM NB & 0.21 (0.19, 0.23) & 0.51 (0.47, 0.53) & 0.22 (0.20, 0.24) & 0.53 (0.50, 0.56) \\
# GAM Po & 0.22 (0.20, 0.24) & 0.54 (0.51, 0.58) & 0.24 (0.22, 0.26) & 0.58 (0.54, 0.62) \\
# Hom Po & 0.26 (0.24, 0.28) & 0.77 (0.71, 0.84) & 0.26 (0.24, 0.29) & 0.78 (0.72, 0.85) \\
# \hline
# \end{tabular}
# \end{table}


# Manually add this row above the rownames: 	\multicolumn{1}{c}{} &	\multicolumn{2}{c}{Random 10-fold CV} &	\multicolumn{2}{c}{Leave-out full transect}\\ %%%% Manually added


### NOw on transect level

tab <- resultList$finalHooded$CVFold$CRPS.trans
theseMod <- c("REGSpatialCVFold10","GAMSpatialCVFold10","GAMPoissonSpatialCVFold10","BayesHomoPois")
useNames <- c("LGCP","GAM NB","GAM Po","Hom Po")
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.trans.hooded <- data.frame(a = paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$lower,2),nsmall=2),", ",format(round(subtab$upper,2),nsmall=2),")",sep=""))
rownames(df.trans.hooded) <- useNames
colnames(df.trans.hooded) <- "CRPS"


tab <- -resultList$finalHooded$CVFold$logScore.trans # Reverse the sign here
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.trans.hooded$b <- paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$upper,2),nsmall=2),", ",format(round(subtab$lower,2),nsmall=2),")",sep="")
colnames(df.trans.hooded)[2] <- "logScore"


tab <- resultList$finalHooded$LeaveOut$CRPS.trans
theseMod <- c("REGSpatialLeaveOutAsCV","GAMSpatialLeaveOutAsCV","GAMPoissonSpatialLeaveOutAsCV","BayesHomoPois") # Needs to be in same order
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.trans.hooded$c <-  paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$lower,2),nsmall=2),", ",format(round(subtab$upper,2),nsmall=2),")",sep="")
colnames(df.trans.hooded)[3] <- "CRPS "# (5% , 95%) "


tab <- -resultList$finalHooded$LeaveOut$logScore.trans # Reverse the sign here
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.trans.hooded$d <- paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$upper,2),nsmall=2),", ",format(round(subtab$lower,2),nsmall=2),")",sep="")
colnames(df.trans.hooded)[4] <- "logScore"# (5%, 95%) "


df.trans.hooded
res <- xtable(df.trans.hooded,
              align=c("l","c","c","c","c"),
              caption = "Validation results on transect level Values in parenthesis are shows the lower and upper bounds of 90% bootstrapped confidence intervals for the performance measure. Cells shown in italics are best (smallest) per column. Those which are significantly smaller than the others (defined as having non-overlapping confidence intervals) are also bolded.",
              label = "tab:Val_hooded_transect")

print.xtable(res,caption.placement = "top")

# % latex table generated in R 3.6.1 by xtable 1.8-4 package
# % Mon Nov 11 10:00:45 2019
# \begin{table}[ht]
# \centering
# \caption{Validation results on transect level Values in parenthesis are shows the lower and upper bounds of 90% bootstrapped confidence intervals for the performance measure. Cells shown in italics are best (smallest) per column. Those which are significantly smaller than the others (defined as having non-overlapping confidence intervals) are also bolded.}
# \label{tab:Val_hooded_transect}
# \begin{tabular}{lcccc}
# \hline
# & CRPS & logScore & CRPS & logScore \\
# \hline
# LGCP & 5.43 (4.04,  6.99) & 3.68 (3.51, 3.86) &  9.91 (5.99, 14.80) &  3.67 (3.26,  4.09) \\
# GAM NB & 5.93 (4.95,  7.00) & 3.79 (3.68, 3.91) &  9.37 (5.66, 13.63) &  3.68 (3.11,  4.27) \\
# GAM Po & 5.90 (4.49,  7.42) & 3.72 (3.50, 3.96) & 10.14 (5.86, 15.09) &  4.14 (3.32,  5.01) \\
# Hom Po & 4.89 (2.21, 11.06) & 3.58 (3.14, 4.65) & 20.70 (1.29, 55.68) & 12.81 (2.46, 37.28) \\
# \hline
# \end{tabular}
# \end{table}


####################### HARPS ################




RdataREGHarpsPath <- file.path(resBaseFolder,"newfinalHarps/REG/runType=PoissonReg seal=harps cov=linearAndSpatial extra.iid=FALSE samp=10000 comment=SpatialAsFixedTransectAsCountingDomain_347/output.RData")
load(RdataREGHarpsPath)

res.list$REGHarps <- list()
res.list$REGHarps$x <- gridList$gridvalX
res.list$REGHarps$y <- gridList$gridvalY
res.list$REGHarps$mean <- finalResList$mean.field.domain.samp
res.list$REGHarps$sd <- finalResList$sd.field.domain.samp

res.list$REGHarps$params <- c(finalResList$intercept.mean,finalResList$covariate.mean,finalResList$spatialX.mean,finalResList$spatialY.mean,finalResList$spatialXY.mean,finalResList$mean.range.param)
names(res.list$REGHarps$params) = c("intercept","covariate","spatialX","spatialY","spatialXY","range")
res.list$REGHarps$abundance.est = c(mean=finalResList$posteriorMean,
                                     median=finalResList$posteriorMedian,
                                     mode=finalResList$posteriorMode,
                                     IQR = finalResList$posteriorIQR,
                                     CI = finalResList$posteriorCI)

res.list$REGHarps$postDens <- list(x = finalResList$posteriorEvalPoints,
                                    y = finalResList$posteriorDist)



RdataGAMHarpsPath <- file.path(resBaseFolder,"newfinalHarps/GAM/runType=GAM seal=harps cov=linearAndSpatial comment=LinearAndSpatial_920/output.RData")
load(RdataGAMHarpsPath)

res.list$GAMHarps <- list()
res.list$GAMHarps$x <- gridList$gridvalX
res.list$GAMHarps$y <- gridList$gridvalY
res.list$GAMHarps$mean <- finalResList$mean.field.domain.samp
res.list$GAMHarps$sd <- finalResList$sd.field.domain.samp

res.list$GAMHarps$fixed.effects <- c(finalResList$intercept.mean,finalResList$covariate.mean,finalResList$spatialX.mean,finalResList$spatialY.mean,finalResList$spatialXY.mean)
names(res.list$GAMHarps$fixed.effects) = c("intercept","covariate","spatialX","spatialY","spatialXY")
res.list$GAMHarps$abundance.est = c(mean=finalResList$posteriorMean,
                                     median=finalResList$posteriorMedian,
                                     mode=finalResList$posteriorMode,
                                     IQR = finalResList$posteriorIQR,
                                     CI = finalResList$posteriorCI)

res.list$GAMHarps$postDens <- list(x = finalResList$posteriorEvalPoints,
                                    y = finalResList$posteriorDist)

RdataGAMPoissonHarpsPath <- file.path(resBaseFolder,"finalHarps/GAMPoisson/runType=GAM seal=harps cov=linearAndSpatial comment=LinearAndSpatial_81/output.RData")
load(RdataGAMPoissonHarpsPath)

res.list$GAMPoissonHarps <- list()
res.list$GAMPoissonHarps$x <- gridList$gridvalX
res.list$GAMPoissonHarps$y <- gridList$gridvalY
res.list$GAMPoissonHarps$mean <- finalResList$mean.field.domain.samp
res.list$GAMPoissonHarps$sd <- finalResList$sd.field.domain.samp
res.list$GAMPoissonHarps$fixed.effects <- c(finalResList$intercept.mean,finalResList$covariate.mean,finalResList$spatialX.mean,finalResList$spatialY.mean,finalResList$spatialXY.mean)
names(res.list$GAMPoissonHarps$fixed.effects) = c("intercept","covariate","spatialX","spatialY","spatialXY")
res.list$GAMPoissonHarps$abundance.est = c(mean=finalResList$posteriorMean,
                                            median=finalResList$posteriorMedian,
                                            mode=finalResList$posteriorMode,
                                            IQR = finalResList$posteriorIQR,
                                            CI = finalResList$posteriorCI)

res.list$GAMPoissonHarps$postDens <- list(x = finalResList$posteriorEvalPoints,
                                           y = finalResList$posteriorDist)

res.list$HomoPoisHarps <- HomoPoisHarps
#> print(res.list$REGHarps$params)
#intercept    covariate     spatialX     spatialY    spatialXY        range
#-2.767656740 14.701701293  0.026821323 -0.011515962 -0.003776375  2.895524643

#### Adding Kingsley to the res.list as well

res.list$KingsleyHarps$abundance.est <- c(mean=KingsleyHarps$N,
                                          median=KingsleyHarps$N,
                                          mode=KingsleyHarps$N,
                                          IQR = 2*qnorm(0.75)*KingsleyHarps$SE,
                                          CI = c(KingsleyHarps$N-qnorm(0.975)*KingsleyHarps$SE,KingsleyHarps$N+qnorm(0.975)*KingsleyHarps$SE))


#### For simplicity we sample from the posteriors to better match the density function of ggplot

####
break1 <- 60000
break2 <- 350000
nobreaks <- 6
logbreaks <-round(10^((log10(break2)-log10(break1))*seq(0,1,length.out=nobreaks)+log10(break1)),-4)

set.seed(123)

dens.df <- data.frame(REG = sample(x = res.list$REGHarps$postDens$x,size = 10^6,prob = res.list$REGHarps$postDens$y,replace = T),
                      GAM = sample(x = res.list$GAMHarps$postDens$x,size = 10^6,prob = res.list$GAMHarps$postDens$y,replace = T),
                      GAMPoisson = sample(x = res.list$GAMPoissonHarps$postDens$x,size = 10^6,prob = res.list$GAMPoissonHarps$postDens$y,replace = T),
                      HomoPois = sample(x = res.list$HomoPoisHarps$postDens$x,size = 10^6,prob = res.list$HomoPoisHarps$postDens$y,replace = T))
dens.df.melt <- melt(dens.df,variable.name = "Model")
pp = ggplot(dens.df.melt,aes(x=value,fill=Model)) + geom_histogram(aes(y=..count../sum(..count..)),alpha=0.3,bins=150,position="identity") +
  coord_cartesian(xlim = c(break1, break2)) +
  xlab("Harp pup count") +
  ylab("Probability")
pp = pp + geom_errorbarh(data=data.frame(KingsleyHarps,Model="GAM",Mod2="Kingsley"),mapping=aes(x=N,y=-0.0075,xmax=N+1.96*SE,xmin=N-1.96*SE,color=Mod2),
                         height = .0075,size=1)
pp = pp + geom_point(data.frame(KingsleyHarps,Model="GAM",Mod2="Kingsley"),mapping=aes(x=N,y=-0.0075,color=Mod2),size=2.5)

pp <- pp + theme_bw()
pp <- pp + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16))
pp = pp +
  scale_color_manual(name="",values=gg_color_hue(4)[4],guide=guide_legend(order=2)) +
  scale_fill_discrete(guide=guide_legend(order=1,override.aes = list(color=NA,alpha=0.3)),
                      labels=c("LGCP","GAM NB","GAM Po","Hom Po"))


pp

#pdf(file=file.path(paperFigureFolder,"postPred_Harps4.pdf"),width=12,height=6)
#pp
#dev.off()


pplog <- pp + scale_x_log10(breaks = logbreaks,
                            labels = logbreaks)
pplog

#pdf(file=file.path(paperFigureFolder,"postPred_Harps_log4.pdf"),width=12,height=6)
#pplog
#dev.off()

pdf(file=file.path(paperFigureFolder,"postPred_Harps_log5.pdf"),width=12,height=6)
pplog
dev.off()



#### Mean and sd of latent field

no.x <- length(res.list$REGHarps$x)
no.y <- length(res.list$REGHarps$y)


latentfield.df <- data.frame(x=rep(res.list$REGHarps$x,no.y),
                             y=rep(res.list$REGHarps$y,each=no.x),
                             mean=c(res.list$REGHarps$mean),
                             sd = c(res.list$REGHarps$sd))
latentfield.df <- latentfield.df[complete.cases(latentfield.df),]


#### Including posterior mean plot for all methods

LGCP_field <- latentfield.df


no.x <- length(res.list$GAMHarps$x)
no.y <- length(res.list$GAMHarps$y)


GAM_NB_field <- data.frame(x=rep(res.list$GAMHarps$x,no.y),
                           y=rep(res.list$GAMHarps$y,each=no.x),
                           mean=c(res.list$GAMHarps$mean),
                           sd = c(res.list$GAMHarps$sd))
GAM_NB_field <- GAM_NB_field[complete.cases(GAM_NB_field),]

GAM_Po_field <- data.frame(x=rep(res.list$GAMPoissonHarps$x,no.y),
                           y=rep(res.list$GAMPoissonHarps$y,each=no.x),
                           mean=c(res.list$GAMPoissonHarps$mean),
                           sd = c(res.list$GAMPoissonHarps$sd))
GAM_Po_field <- GAM_Po_field[complete.cases(GAM_Po_field),]

Hom_Po_field <- GAM_Po_field
Hom_Po_field$mean <- log(intensityHarps)
Hom_Po_field$sd <- NA

library(data.table)
mean_field_dt <- as.data.table(rbind(cbind(LGCP_field,type="LGCP"),
                       cbind(GAM_NB_field,type="GAM NB"),
                       cbind(GAM_Po_field,type="GAM Po"),
                       cbind(Hom_Po_field,type="Hom Po")))



limits <- range(mean_field_dt$mean)

mean_field_dt <- mean_field_dt[type!="LGCP",]

ff <- ggplot(mean_field_dt,aes(x=x,y=y)) + geom_raster(aes(fill=mean)) +   xlab("x (km)") +  ylab("y (km)")
ff <- ff + theme_bw()
ff <- ff + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16),
                 strip.text = element_text(size=16))
ff <- ff +  scale_fill_gradientn(colours=heat.colors(256,rev=T),
                                 limits = c(-7.5,9), #"=c(-8,9),
                                 breaks = c(-8,-4,0,4,8))
#
#                                 values = rescale(mean_field_df$mean))


ff <- ff + facet_wrap(~type)

pdf(file=file.path(paperFigureFolder,"mean_field_other_methods_Harps.pdf"),width=14,height=6)
ff
dev.off()



### AND just here we create the original mean+sd latent field Harps plot
ff <- ggplot(latentfield.df,aes(x=x,y=y)) + geom_raster(aes(fill=mean))+   xlab("x (km)") +  ylab("y (km)")
ff <- ff + theme_bw()
ff <- ff + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16))
ff <- ff +  scale_fill_gradientn(colours=heat.colors(256,rev=T),
                                 limits = c(-7.5,9), #"=c(-8,9),
                                 breaks = c(-4,0,4,8))
ff
pdf(file=file.path(paperFigureFolder,"mean_latent_field_Harps.pdf"),width=6,height=6)
ff
dev.off()

ff <- ggplot(latentfield.df,aes(x=x,y=y)) + geom_raster(aes(fill=sd)) +   xlab("x (km)") +  ylab("y (km)")
ff <- ff + theme_bw()
ff <- ff + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16))
ff <- ff +  scale_fill_gradientn(colours=heat.colors(256,rev=T))
ff
pdf(file=file.path(paperFigureFolder,"sd_latent_field_Harps.pdf"),width=6,height=6)
ff
dev.off()




#### Konverterer resultattabeller

#####

abundanceDf.Harps <- as.data.frame(rbind(res.list$REGHarps$abundance.est,res.list$GAMHarps$abundance.est,
                                         res.list$GAMPoissonHarps$abundance.est,res.list$HomoPoisHarps$abundance.est,
                                         res.list$KingsleyHarps$abundance.est))
rownames(abundanceDf.Harps) <- c("LGCP","GAM NB","GAM Po", "Hom Po","Kingsley")
colnames(abundanceDf.Harps) <- c("mean","median","mode","IQR","0.025-quantile","0.975-quantile")


res0 <- xtable(abundanceDf.Harps,digits = 0,
               align=c("l","c","c","c","c","c","c"),
               caption = "Summary tables for the posterior predictive distributions for the total area counts of harp seals.",
               label = "tab:summary_Harps")

print.xtable(res0,caption.placement = "top")

# % latex table generated in R 3.6.1 by xtable 1.8-4 package
# % Mon Nov 11 10:51:14 2019
# \begin{table}[ht]
# \centering
# \caption{Summary tables for the posterior predictive distributions for the total area counts of harp seals.}
# \label{tab:summary_Harps}
# \begin{tabular}{lcccccc}
# \hline
# & mean & median & mode & IQR & 0.025-quantile & 0.975-quantile \\
# \hline
# LGCP & 147919 & 127965 & 110996 & 72347 & 69267 & 357185 \\
# GAM NB & 98617 & 98035 & 91876 & 12895 & 81023 & 119349 \\
# GAM Po & 84852 & 84852 & 84910 & 1536 & 82681 & 87094 \\
# Hom Po & 88272 & 88267 & 88257 & 1583 & 85986 & 90587 \\
# Kingsley & 85968 & 85968 & 85968 & 15926 & 62829 & 109107 \\
# \hline
# \end{tabular}
# \end{table}


print(res.list$REGHarps$params)

#intercept    covariate     spatialX     spatialY    spatialXY        range
#-2.767656740 14.701701293  0.026821323 -0.011515962 -0.003776375  2.895524643






####### Validation tables also for Harps


tab <- resultList$finalHarps$CVFold$CRPS.photo
theseMod <- c("REGSpatialCVFold10","GAMSpatialCVFold10","GAMPoissonSpatialCVFold10","BayesHomoPois")
useNames <- c("LGCP","GAM NB","GAM Po","Hom Po")
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.photo.harps <- data.frame(a = paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$lower,2),nsmall=2),", ",format(round(subtab$upper,2),nsmall=2),")",sep=""))
rownames(df.photo.harps) <- useNames
colnames(df.photo.harps) <- "CRPS"


tab <- -resultList$finalHarps$CVFold$logScore.photo # Reverse the sign here
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.photo.harps$b <- paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$upper,2),nsmall=2),", ",format(round(subtab$lower,2),nsmall=2),")",sep="")
colnames(df.photo.harps)[2] <- "logScore"


tab <- resultList$finalHarps$LeaveOut$CRPS.photo
theseMod <- c("REGSpatialLeaveOutAsCV","GAMSpatialLeaveOutAsCV","GAMPoissonSpatialLeaveOutAsCV","BayesHomoPois") # Needs to be in same order
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.photo.harps$c <-  paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$lower,2),nsmall=2),", ",format(round(subtab$upper,2),nsmall=2),")",sep="")
colnames(df.photo.harps)[3] <- "CRPS"# (5%, 95%) "


tab <- -resultList$finalHarps$LeaveOut$logScore.photo # Reverse the sign here
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.photo.harps$d <- paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$upper,2),nsmall=2),", ",format(round(subtab$lower,2),nsmall=2),")",sep="")
colnames(df.photo.harps)[4] <- "logScore"# (5%, 95%) "


df.photo.harps
res <- xtable(df.photo.harps,
              align=c("l","c","c","c","c"),
              caption = "Validation results on photo level. Values in parenthesis are shows the lower and upper bounds of 90% bootstrapped confidence intervals for the performance measure. Cells shown in italics are best (smallest) per column. Those which are significantly smaller than the others (defined as having non-overlapping confidence intervals) are also bolded.",
              label = "tab:Val_harps_photo")

print.xtable(res,caption.placement = "top")

# % latex table generated in R 3.6.1 by xtable 1.8-4 package
# % Mon Nov 11 10:05:16 2019
# \begin{table}[ht]
# \centering
# \caption{Validation results on photo level. Values in parenthesis are shows the lower and upper bounds of 90% bootstrapped confidence intervals for the performance measure. Cells shown in italics are best (smallest) per column. Those which are significantly smaller than the others (defined as having non-overlapping confidence intervals) are also bolded.}
# \label{tab:Val_harps_photo}
# \begin{tabular}{lcccc}
# \hline
# & CRPS & logScore & CRPS & logScore \\
# \hline
# LGCP & 1.14 (1.01, 1.27) & 0.95 (0.91, 1.00) & 1.96 (1.72, 2.20) & 1.28 (1.22, 1.33) \\
# GAM NB & 1.78 (1.58, 2.00) & 1.17 (1.11, 1.22) & 1.90 (1.67, 2.13) & 1.27 (1.21, 1.33) \\
# GAM Po & 2.32 (2.10, 2.55) & 2.09 (2.00, 2.17) & 2.46 (2.22, 2.71) & 2.17 (2.08, 2.26) \\
# Hom Po & 2.64 (2.40, 2.90) & 3.47 (3.28, 3.67) & 2.66 (2.42, 2.92) & 3.49 (3.30, 3.69) \\
# \hline
# \end{tabular}
# \end{table}


# Manually add this row above the rownames: 	\multicolumn{1}{c}{} &	\multicolumn{2}{c}{Random 10-fold CV} &	\multicolumn{2}{c}{Leave-out full transect}\\ %%%% Manually added


### NOw on transect level

tab <- resultList$finalHarps$CVFold$CRPS.trans
theseMod <- c("REGSpatialCVFold10","GAMSpatialCVFold10","GAMPoissonSpatialCVFold10","BayesHomoPois")
useNames <- c("LGCP","GAM NB","GAM Po","Hom Po")
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.trans.harps <- data.frame(a = paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$lower,2),nsmall=2),", ",format(round(subtab$upper,2),nsmall=2),")",sep=""))
rownames(df.trans.harps) <- useNames
colnames(df.trans.harps) <- "CRPS"


tab <- -resultList$finalHarps$CVFold$logScore.trans # Reverse the sign here
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.trans.harps$b <- paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$upper,2),nsmall=2),", ",format(round(subtab$lower,2),nsmall=2),")",sep="")
colnames(df.trans.harps)[2] <- "logScore"


tab <- resultList$finalHarps$LeaveOut$CRPS.trans
theseMod <- c("REGSpatialLeaveOutAsCV","GAMSpatialLeaveOutAsCV","GAMPoissonSpatialLeaveOutAsCV","BayesHomoPois") # Needs to be in same order
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.trans.harps$c <-  paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$lower,2),nsmall=2),", ",format(round(subtab$upper,2),nsmall=2),")",sep="")
colnames(df.trans.harps)[3] <- "CRPS"# (5%, 95%) "


tab <- -resultList$finalHarps$LeaveOut$logScore.trans # Reverse the sign here
theseRows <- match(theseMod,rownames(tab))
subtab <- tab[theseRows,]

df.trans.harps$d <- paste(format(round(subtab$mean,2),nsmall=2)," (",format(round(subtab$upper,2),nsmall=2),", ",format(round(subtab$lower,2),nsmall=2),")",sep="")
colnames(df.trans.harps)[4] <- "logScore"# (5%, 95%) "


df.trans.harps
res <- xtable(df.trans.harps,
              align=c("l","c","c","c","c"),
              caption = "Validation results on transect level Values in parenthesis are shows the lower and upper bounds of 90% bootstrapped confidence intervals for the performance measure. Cells shown in italics are best (smallest) per column. Those which are significantly smaller than the others (defined as having non-overlapping confidence intervals) are also bolded.",
              label = "tab:Val_harps_transect")

print.xtable(res,caption.placement = "top")

# % latex table generated in R 3.6.1 by xtable 1.8-4 package
# % Mon Nov 11 10:13:51 2019
# \begin{table}[ht]
# \centering
# \caption{Validation results on transect level Values in parenthesis are shows the lower and upper bounds of 90% bootstrapped confidence intervals for the performance measure. Cells shown in italics are best (smallest) per column. Those which are significantly smaller than the others (defined as having non-overlapping confidence intervals) are also bolded.}
# \label{tab:Val_harps_transect}
# \begin{tabular}{lcccc}
# \hline
# & CRPS & logScore & CRPS & logScore \\
# \hline
# LGCP &  95.98 (51.20, 148.78) &  6.57 (5.93,   7.33) & 152.95 (111.93, 198.59) &  6.49 (5.98,   6.94) \\
# GAM NB &  88.33 (70.94, 106.87) &  6.45 (6.31,   6.62) &  96.70 ( 61.05, 139.91) &  7.00 (6.24,   7.76) \\
# GAM Po &  60.93 (42.57,  79.87) &  8.03 (6.83,   9.22) &  55.96 ( 39.62,  74.51) &  8.42 (7.34,   9.39) \\
# Hom Po & 102.62 ( 8.87, 296.09) & 62.92 (4.31, 300.60) & 149.38 (  7.70, 454.53) & 82.63 (4.08, 251.00) \\
# \hline
# \end{tabular}
# \end{table}

###HERE

#### Mesh ####


countingDomainExpanded <- CreateCountDomainPolygon(transectStartCoordX = dataList$transect$transectStartCoordX,
                                                   transectStartCoordY = dataList$transect$transectStartCoordY,
                                                   transectEndCoordX = dataList$transect$transectEndCoordX,
                                                   transectEndCoordY = dataList$transect$transectEndCoordY,
                                                   coordPhotoX = dataList$org$coordPhotoX,
                                                   coordPhotoY = dataList$org$coordPhotoY,
                                                   photoWidth = dataList$org$photoWidth,
                                                   photoHeight = dataList$org$photoHeight,
                                                   transectYSpan = 2*1.852,
                                                   transectXSpan = 0.5*1.852,
                                                   theseTransInCountDomain=theseTransInCountDomain)

### NEW: We use the mesh loaded when extracting the REGHarps results above instead as the one below here is not what we have used in practice!!! ###

# domain.outer  <- inla.nonconvex.hull(points=cbind(rectangleCentersX,rectangleCentersY),
#                                      convex=convHullVar.convex+0.03,
#                                      concave=convHullVar.concave+0.03,
#                                      resolution=convHullVar.resolution)
#
# domain.outer.final  <- inla.nonconvex.hull(points=cbind(rectangleCentersX,rectangleCentersY),
#                                            convex=convHullVar.convex,
#                                            concave=convHullVar.concave,
#                                            resolution=convHullVar.resolution)
# domain.inner <- INLA::inla.mesh.segment(loc = as.matrix(countingDomainExpanded))
#
# obligMeshLoc <- as.data.frame(cbind(x=rectangleCentersX,y=rectangleCentersY))
#
# mesh <- inla.mesh.2d(loc=obligMeshLoc,
#                      boundary = list(domain.inner,domain.outer),
#                      max.edge=c(1.5,10), # Try c(1.5,10)
#                      offset=meshVar.offset,
#                      cutoff=0.17) # Try 0.17

#plot(mesh)
#axis(side=1)
#axis(side=2)

# Formatting data


countingDomain$LineType <- "Whelping region"

meshdata <- data.frame(a=rbind(mesh$loc[mesh$graph$tv[,1],c(1,2)],mesh$loc[mesh$graph$tv[,2],c(1,2)],mesh$loc[mesh$graph$tv[,1],c(1,2)]),
                       b=rbind(mesh$loc[mesh$graph$tv[,2],c(1,2)],mesh$loc[mesh$graph$tv[,3],c(1,2)],mesh$loc[mesh$graph$tv[,3],c(1,2)]),
                       LineType = "Mesh")
boundaryData <- data.frame(mesh$loc[mesh$segm$bnd$idx,c(1,2)],
                           LineType="Modeling boundary")


gg = ggplot() +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16)) +
  xlab("x (km)") +
  ylab("y (km)")
gg <- gg + geom_path(data = countingDomain, aes_string(x="x",y="y",color="LineType"))
gg <- gg + geom_path(data = boundaryData, aes_string(x="X1",y="X2",color="LineType"))
gg_nomesh <- gg + scale_color_manual(name="Line type",
                                     breaks = c("Modeling boundary","Whelping region"),
                                     values = c(gg_color_hue(3)[2],"red"))
gg <- gg + geom_segment(data = meshdata,
                        aes_string(x="a.1",y="a.2",xend="b.1",yend="b.2",color="LineType"),size=0.1)
gg <- gg + scale_color_manual(name="Line type",
                              breaks = c("Mesh","Modeling boundary","Whelping region"),
                              values = c("black",gg_color_hue(3)[3],"red"))
#gg_nomesh
#gg

#pdf(file=file.path(paperFigureFolder,"mesh_new.pdf"),width=7,height=6)
#gg
#dev.off()

### Trying to a dd a subplot as well to show finer detail

xmin = 32
xmax = 39
ymin = -9
ymax = 2

ggsub <- gg + coord_cartesian(ylim=c(ymin, ymax),xlim=c(xmin,xmax)) +
  theme(legend.position="none",
        axis.text = element_text(size=8),
        axis.title = element_blank()) +
  theme(plot.background = element_rect(color = gg_color_hue(3)[2],size=1))

ggsub_nomesh <- gg_nomesh + coord_cartesian(ylim=c(ymin, ymax),xlim=c(xmin,xmax)) +
  theme(legend.position="none",
        axis.text = element_text(size=8),
        axis.title = element_blank()) +
  theme(plot.background = element_rect(color = gg_color_hue(3)[2],size=1))


ggfinal0 <- gg + geom_rect(data=meshdata,xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, alpha=0.006,fill = gg_color_hue(3)[2]) +
  geom_segment(data=meshdata,x=xmin,y=ymin,xend=12,yend=-40.5,size=0.3,col=gg_color_hue(3)[2]) +
  geom_segment(data=meshdata,x=xmax,y=ymin,xend=69,yend=-40.5,size=0.3,col=gg_color_hue(3)[2])
layers <- ggfinal0$layers
#ggfinal0$layers[[1]] <- layers$layers[[4]]
#ggfinal0$layers[[2]] <- layers$layers[[1]]
#ggfinal0$layers[[3]] <- layers$layers[[2]]
#ggfinal0$layers[[4]] <- layers$layers[[3]]


ggfinal <- ggfinal0 + annotation_custom(ggplotGrob(ggsub),xmin=12,xmax=69,ymin=-108,ymax=-40)
ggfinal

pdf(file=file.path(paperFigureFolder,"mesh_new2.pdf"),width=7,height=6)
ggfinal
dev.off()
ggsave(filename = file.path(paperFigureFolder,"mesh_new2.pdf"),
         width = 7,
         height = 6)
 ggsave(filename = file.path(paperFigureFolder,"mesh_new2.png"),
        width = 7,
        height = 6)
 ggsave(filename = file.path(paperFigureFolder,"mesh_new2v2.png"),
        width = 7,
        height = 6,
        dpi = 600)
 ggsave(filename = file.path(paperFigureFolder,"mesh_new2v3.png"),
        width = 7,
        height = 6,
        dpi = 1200)

#### Extending the ggfinal plot with photo extents ###

 library(data.table)

 photo_polygons <- data.table(x = c(dataList$org$coordPhotoX-0.5*dataList$org$photoWidth,
                                    dataList$org$coordPhotoX+0.5*dataList$org$photoWidth,
                                    dataList$org$coordPhotoX+0.5*dataList$org$photoWidth,
                                    dataList$org$coordPhotoX-0.5*dataList$org$photoWidth),
                              y = c(dataList$org$coordPhotoY-0.5*dataList$org$photoHeight,
                                    dataList$org$coordPhotoY-0.5*dataList$org$photoHeight,
                                    dataList$org$coordPhotoY+0.5*dataList$org$photoHeight,
                                    dataList$org$coordPhotoY+0.5*dataList$org$photoHeight),
                              id = rep(1:length(dataList$org$coordPhotoX),4))
 setkey(photo_polygons,"id")
 #photo_polygons[,value:=id]

  id.inside.subfig <- which(dataList$org$coordPhotoX > (xmin-1) & dataList$org$coordPhotoX < (xmax+1) &
                     dataList$org$coordPhotoY > (ymin-1) & dataList$org$coordPhotoY < (ymax+1))

 photo_polygons <-  photo_polygons[id %in% id.inside.subfig,]
 photo_polygons[,col:=(id %%2)+1]
 photo_polygons[,id:=as.factor(id)]


 # set.seed(12)
 # nopic <- nrow(photo_polygons)/4
 # newid <- data.table(id=unique(photo_polygons$id),newid=sample(1:nopic,nopic,replace=F))
 # photo_polygons <- merge(photo_polygons,newid,by="id")
 # photo_polygons[,id:=newid]
 #
 #
 # photo_polygons[,id:=as.factor(id)]


# ggplot(photo_polygons,aes(x=x,y=y)) + geom_polygon(aes(x=x,y=y,group=id),alpha=0.5)




 photos <- list()
 photos$pic.x.left <- dataList$org$coordPhotoX-0.5*dataList$org$photoWidth
 photos$pic.x.right <- dataList$org$coordPhotoX+0.5*dataList$org$photoWidth
 photos$pic.y.bottom <- dataList$org$coordPhotoY-0.5*dataList$org$photoHeight
 photos$pic.y.top <- dataList$org$coordPhotoY+0.5*dataList$org$photoHeight
 orgPhotos <- as.data.frame(photos)
 modPhotos <- orgPhotos

 modObservationDomain <- as.data.frame(MergePolygons(modPhotos))


 voronoiTess <- CreateVoronoiTessellation(locationsCoordX = mesh$loc[,1],
                                          locationsCoordY = mesh$loc[,2],
                                          observationDomain = modObservationDomain,
                                          parallelize.numCores = 10) ### Might consider replacing countingDomian with a somewhat larger modelling domain here instead

 # insidePhotos.vec <- c()
 # for (i in 1:mesh$n){
 #   insidePhotos.vec <- c(insidePhotos.vec,which(as.logical(point.in.polygon(dataList$mod$coordPhotoX,dataList$mod$coordPhotoY,
 #                                                     voronoiTess$tiles[[i]]$x,voronoiTess$tiles[[i]]$y))))
 # }
 #
 # table(sapply(insidePhotos.list,length))

 voronoiTess_polygons <- rbindlist(lapply(voronoiTess$tiles,FUN = function(x)data.table(x=x$x,y=x$y)),idcol = "id")
 voronoiTess_polygons[,value:=id]

 voronoiTess_polygons[,id:=as.factor(id)]

# ggplot(voronoiTess_polygons) + geom_polygon(aes(x=x,y=y,group=id,fill=value),alpha=0.5)

 ggsub1 <-  ggsub +geom_polygon(data = photo_polygons,aes(x=x,y=y,group=id,fill=col),alpha=0.5) +
   theme_bw() +
   theme(legend.position = "none") +
   xlab("x (km)") +
   ylab("y (km)") + ggtitle("Mesh and photos")

ggsub2 <- ggsub_nomesh + geom_polygon(data=voronoiTess_polygons,aes(x=x,y=y,group=id),color="purple",alpha=0,size=0.2) +
   geom_polygon(data = photo_polygons,aes(x=x,y=y,group=id,fill=col),alpha=0.5) +
  theme_bw() +
theme(legend.position = "none") +
  xlab("x (km)") + ggtitle("Voronoi tesselations and photos")
#  ylab("y (km)")


#ggsub1
#ggsub2
library(gridExtra)


pdf(file=file.path(paperFigureFolder,"mesh_details.pdf"),width=8,height=5)
grid.arrange(ggsub1,ggsub2,ncol=2)
dev.off()


 #
# #
#  png(filename=file.path(paperFigureFolder,"mesh_new2.png"),
#      width=480*7/6,height=480*6/6)
#  ggfinal
#  dev.off()
#  png(filename=file.path(paperFigureFolder,"mesh_new2v2.png"),
#      width=480*7/6,height=480*6/6,pointsize = 12)
#  ggfinal
#  dev.off()


 # jpeg(filename=file.path(paperFigureFolder,"mesh_new2.jpeg"),
#     width=7*100,height=6*100)
# ggfinal
# dev.off()
# bmp(filename=file.path(paperFigureFolder,"mesh_new2.bmp"),
#      width=7*100,height=6*100)
# ggfinal
# dev.off()
# tiff(filename=file.path(paperFigureFolder,"mesh_new2.tiff"),
#     width=7*100,height=6*100)
# ggfinal
# dev.off()


#
#
# ,b=mesh$loc[mesh$graph$tv[,2],c(1,2)]),
#                   data.frame(a=mesh$loc[mesh$graph$tv[,2],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
#                   data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]))
#
#
#
# gg = gg + geom_segment(data = ),
#                        aes_string(x="a.1",y="a.2",xend="b.1",yend="b.2"),size=0.1)
#
# gg = gg + geom_segment(data = data.frame(a=mesh$loc[mesh$graph$tv[,2],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)]),
#                        aes_string(x="a.1",y="a.2",xend="b.1",yend="b.2"),size=0.1)
#
# gg = gg + geom_segment(data = data.frame(a=mesh$loc[mesh$graph$tv[,1],c(1,2)],b=mesh$loc[mesh$graph$tv[,3],c(1,2)],LineType="Mesh"),
#                        aes_string(x="a.1",y="a.2",xend="b.1",yend="b.2",color="LineType"),size=0.1)
#
# gg <- gg + geom_path(data = countingDomain2, aes_string(x="x",y="y",color="LineType"))
#
#
# points(dataList$mod$coordPhotoX[dataList$mod$noObsPerPhoto==0],dataList$mod$coordPhotoY[dataList$mod$noObsPerPhoto==0],col="red",cex=0.5,pch=16)
# points(dataList$mod$coordPhotoX[dataList$mod$noObsPerPhoto>0],dataList$mod$coordPhotoY[dataList$mod$noObsPerPhoto>0],col="green",cex=0.5,pch=16)
#
#
#
#
# pdf(file=file.path(savingFolder,"mesh.pdf"),width=6,height=6)
# par(mar=c(0,0,0,0))
# plot(mesh, asp=1, main='')
# legend("bottomright",c("y=0","y>0"),col=c("red","green"),pch=16)
# legend("topleft",c("Model boundary","counting domain"),col=c("blue",6),lty=c(1),lwd=3)
#

# predPhotoMeshPointList <- MeshPointsInPredPhotos(mesh = mesh,
#                                                  predPhotos = orgPhotos)
#
# noMeshPoints = rep(NA,length(predPhotoMeshPointList))
# for (i in 1:length(predPhotoMeshPointList)){
#   noMeshPoints[i] <- predPhotoMeshPointList[[i]]$noMeshPoints
# }
# print(paste("ALL prediction photos matching exactly 1 mesh point? ",all.equal(noMeshPoints,rep(1,length(predPhotoMeshPointList))),sep=""))
#
# thisMeshPoint <- rep(NA,length(predPhotoMeshPointList))
# for (i in 1:length(predPhotoMeshPointList)){
#   thisMeshPoint[i] <- which(predPhotoMeshPointList[[i]]$logical)
# }
#




#### HERE !!!!

####################### Satellite imagery ######

### COPIED FROM PREVIOUS FOLDER !!!! #

satelliteFile = "/nr/project/stat/PointProcess/Data/Seals/Satellite/cov_grid_band1.rds"#"M:\\PointProcess/Data/Seals/Satellite/cov_grid_band1.rds"
satelliteFile = "/nr/project/stat/PointProcess/Data/Seals/Satellite/cov_grid_band1_5000.rds"#"M:\\PointProcess/Data/Seals/Satellite/cov_grid_band1.rds"

 covGrid <- readRDS(satelliteFile)
#plot(cov)

  rangeX <- range(countingDomain$x)
  rangeY <- range(countingDomain$y)
#} else {
  rangeX <- range(mesh$loc[,1]) # range of x-coord of grid
  rangeY <- range(mesh$loc[,2]) # range of y-coord of grid
#}


nxy <- round(c(diff(rangeX),diff(rangeY))/grid.pixelsize) # The number of points of the grid in x- and y-direction
projgrid <- inla.mesh.projector(mesh,xlim=rangeX,ylim=rangeY,dims=nxy) # Defines the projection
gridvalX <- projgrid$x
gridvalY <- projgrid$y


indGridX <- covGrid$xcol
indGridY <- covGrid$yrow

alloldGridX <- rep(indGridX,times=covGrid$dim[2])
alloldGridY <- rep(indGridY,each=covGrid$dim[1])

covGrid$v2 <- covGrid$v

covGrid$v2[covGrid$v==mean(covGrid$v)] <- NA

covNewGrid <- fields::as.image(Z=c(t(covGrid$v2)), x = cbind(x=alloldGridX,y=alloldGridY), grid = list(x=gridvalX,y=gridvalY),na.rm=T)
covNewGridval <- covNewGrid$z
covNewGridval[covNewGridval<0.01] <- NA # Added!
covNewGrid$z <- covNewGridval # Added

#### Mean and sd of latent field

no.x <- length(covNewGrid$x)
no.y <- length(covNewGrid$y)


cov.df <- data.frame(x=rep(covNewGrid$x,no.y),
                     y=rep(covNewGrid$y,each=no.x),
                     z=c(covNewGrid$z))

cov.df <- cov.df[complete.cases(cov.df),]

bb=pnt.in.poly(cbind(cov.df$x,cov.df$y),cbind(boundaryData$X1,boundaryData$X2))
#aa=point.in.polygon(cov.df$x,cov.df$y,boundaryData$X1,boundaryData$X2)
#cov.df[as.logical(aa),]
cov.df.new <- cov.df[which(bb$pip==1),]


ff <- ggplot(cov.df.new,aes(x=x,y=y)) + geom_raster(aes(fill=z))+   xlab("x (km)") +  ylab("y (km)")
ff <- ff + theme_bw()
ff <- ff + theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=16),
                 legend.text = element_text(size=12),
                 legend.title = element_text(size=16))
#ff <- ff +  scale_fill_gradientn(colours=tim.colors(128),name="Density")
ff <- ff +  scale_fill_gradientn(colours=two.colors(n=256,start="darkblue",middle="white",end="red",alpha=1.0),name="Density")

#ff <- ff + scale_fill_discrete()
ff
ff2 <- ff + geom_path(data = countingDomain, mapping=aes(x=x,y=y,linetype=LineType),color="black",size=1)
#ff2 <- ff2 + geom_path(data = boundaryData, aes(x=X1,y=X2,linetype=LineType),color=gg_color_hue(4)[3],size=1.5)
ff2 <- ff2 + guides(fill = guide_colorbar(order = 1)) +
  scale_linetype_manual("",values=c(1))#,guide = guide_legend(override.aes = list(color=c(gg_color_hue(5)[5])),
                                                                #   order = 0))
ff2

pdf(file=file.path(paperFigureFolder,"satellite_data.pdf"),width=7,height=6)
ff2
dev.off()




# g1=ggplot(seals2,aes(x=xkm,y=ykm,group=1))
# g1 <- g1+
#   geom_point(size=0.3,shape=15,color="grey")+
#   geom_point(aes(size=no,color=type),alpha=I(0.2))+
#   scale_size_continuous(name="Seal count",range=c(1,8))+
#   scale_color_discrete(name="Seal type")
# g1 <- g1 + theme_bw()
# g1 <- g1 + theme(axis.text=element_text(size=14),
#                  axis.title=element_text(size=16),
#                  legend.text = element_text(size=12),
#                  legend.title = element_text(size=16))
# g1 <- g1 + labs(x = "x (km)",
#                 y = "y (km)")
# g1
#
# g2 <- g1 + geom_path(data = countingDomain, mapping=aes(x=x,y=y,linetype=LineType),color=gg_color_hue(4)[4])
# g2 <- g2 + geom_path(data = data.frame(x=0,y=0,LineType="Transect"), mapping=aes(x=x,y=y,linetype=LineType),color="grey")
# g2 <- g2+ scale_linetype_manual("Other",values=c(1,1),guide = guide_legend(override.aes = list(colour=c(gg_color_hue(4)[4],"grey"))))
# g2
#
# pdf(file=file.path(paperFigureFolder,"seal_data.pdf"),width=7,height=6)
# g2
# dev.off()



#dfdf=data.frame(x=res.list$REGHooded$x,y=res.list$REGHooded$y,z=res.list$REGHooded$z)#
#
#ggplot(dfdf) + geom_raster(aes(x,y,z))



# pg <- ggplot_build(pp)
#
#
# scale_y_continuous(breaks = c(0,0.00025,0.00050,0.00075,0.001),
#                    labels = c(0,0.00025,0.00050,0.00075,0.001)*)
#
#
# dens.df <- data.frame(x = 1:50000,
#                       y.REGHooded = 0,
#                       y.GAMHooded = 0,
#                       y.GAMPoissonHooded = 0)
#
# for (i in 1:length(res.list$REGHooded$postDens$x)){
#   xx <- res.list$REGHooded$postDens$x[i]
#   this.x <- which(dens.df$x==xx)
#   dens.df$y.REGHooded[this.x]= res.list$REGHooded$postDens$y[i]
# }
#
# for (i in 1:length(res.list$GAMHooded$postDens$x)){
#   xx <- res.list$GAMHooded$postDens$x[i]
#   this.x <- which(dens.df$x==xx)
#   dens.df$y.GAMHooded[this.x] = res.list$GAMHooded$postDens$y[i]
# }
#
# for (i in 1:length(res.list$GAMPoissonHooded$postDens$x)){
#   xx <- res.list$GAMPoissonHooded$postDens$x[i]
#   this.x <- which(dens.df$x==xx)
#   dens.df$y.GAMPoissonHooded[this.x] = res.list$GAMPoissonHooded$postDens$y[i]
# }
#
# ggplot(dens.df,aes(x)) +
#   geom_line(aes(y=y.REGHooded,colour="REG")) +
#   geom_line(aes(y=y.GAMHooded,colour="GAM")) +
#   geom_line(aes(y=y.GAMPoissonHooded,colour="GAMPoisson")) +
#   scale_x_continuous(limits = c(8000, 16000))
#
#
#
#
# stacked <- with(dens.df,
#                 data.frame(value = c(y.REGHooded, y.GAMHooded,y.GAMPoissonHooded),
#                            variable = factor(rep(c("REG","GAM","GAMPoisson"),
#                                                  each = NROW(dens.df))),
#                            x = rep(x, 3)))
#
#
# require(ggplot2)
# p <- ggplot(stacked, aes(x, value, colour = variable))
# p + geom_line()
#
#
# library(ggplot2)
#


############################# Individual transect figures #######


### Bootstrapping function

library(boot)


finalHoodedFolders <- c("REGSpatialLeaveOutAsCV","GAMSpatialLeaveOutAsCV") # Needs to be in same order

resultsBaseFolderVec1 <- paste("/home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/",finalHoodedFolders,"/",sep="")

finalHarpsFolders <- c("REGSpatialLeaveOutAsCV","GAMSpatialLeaveOutAsCV") # Needs to be in same order

resultsBaseFolderVec2 <- paste("/home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/",finalHarpsFolders,"/",sep="")

resultsBaseFolderVec <- c(resultsBaseFolderVec1,resultsBaseFolderVec2)

save.resultList <- TRUE

transDfList <- list()
covarageList50 <- list()
covarageList90 <- list()

for (m in 1:length(resultsBaseFolderVec)){

  resultsBaseFolder <- resultsBaseFolderVec[m]

  #### Preparing re-ordering of CV-photos, only used for those cases ####
  sealPhotoDataFile = "/nr/project/stat/PointProcess/Data/Seals/OriginalSealsKm.rds"
  seals <- readRDS(file=sealPhotoDataFile) # RDS-file

  n <- nrow(seals)
  set.seed(123)
  randomOrder <- sample(1:n,n)
  splittedRandomOrder=split(randomOrder,ceiling((1:n)/n*10))

  thisOrder <- NULL
  for (i in 1:10){
    thesePhotos <- (1:n %in% splittedRandomOrder[[i]])
    thisOrder <- c(thisOrder,which(thesePhotos))
  }
  usedOrder <- order(thisOrder)

  split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))
  sealFolder <- split_path(resultsBaseFolder)[2]
  lastFolder <- split_path(resultsBaseFolder)[1]

  ## Number of observed seals of correct type for each photo
  if (sealFolder=="finalHooded") sealCount <- seals$hooded
  if (sealFolder=="finalHarps")  sealCount <- seals$harps


  sealTransectDataFile = "/nr/project/stat/PointProcess/Data/Seals/OigardTablesTransformed.rds"
  transectData <- readRDS(sealTransectDataFile)

  # Assigning each photo to a single transect
  transectYmid <- (transectData$y.start + transectData$y.end)/2

  photoinTransectVec <- rep(NA,n)
  for (i in 1:n){
    photoinTransectVec[i] <-which.min(abs(seals$ykm[i]-transectYmid))
  }




  ################## Here the actual function starts ####

  folders <- list.dirs(resultsBaseFolder,recursive = F)
  folders0 <- list.dirs(resultsBaseFolder,recursive = F,full.names=F)

  noTrans <- length(folders)

  transVec <- rep(NA,length(noTrans))
  for (i in 1:noTrans){
    transVec[i] <- as.numeric(getstr(folders0[i]))
  }

  transVec <- transVec[!is.na(transVec)]

  transExists <- rep(NA,noTrans)
  comparisonList <- list()
  for (i in 1:length(transVec)){
    thisfile <- paste(folders[i],"/photo_pred_comparisons_transect_",transVec[i],".rds",sep="")

    if (file.exists(thisfile)){
      comparisonList[[i]] <- readRDS(file = thisfile)
      transExists[i] <- T
    } else {
      transExists[i] <- F
    }
  }
  transVec <- transVec[transExists]



  ### In the meanwhile...

  quant.finder <- function(vec,q,evalvec){
    vec.q = vec-q
    evalvec[vec.q>0][1]
  }




  transDf <- NULL
  infiniteLogScoreTrans <- NULL

  for (i in order(transVec)){

    countsTrans <- comparisonList[[i]]$trans.list$posteriorevalFullTransect
    posteriorDistTransectNew <- comparisonList[[i]]$trans.list$posteriorDistTransect
    FhatTrans <- comparisonList[[i]]$trans.list$FhatTrans

    trueCountTrans <- sum(comparisonList[[i]][[5]]$photoCounts)

    initial <- FhatTrans - (trueCountTrans <= countsTrans)*1

    CRPSTransL <- sum(initial[initial>0]^2) # Underestimating
    CRPSTransU <- sum(initial[initial<=0]^2) # Over estimation

    CRPSTrans <- sum(initial^2)

    meanTrans <- c(t(countsTrans)%*%posteriorDistTransectNew)
    medTrans <- quant.finder(vec=FhatTrans,q=0.5,evalvec = countsTrans)

    quant005Trans <- quant.finder(vec=FhatTrans,q=0.05,evalvec = countsTrans)
    quant025Trans <- quant.finder(vec=FhatTrans,q=0.25,evalvec = countsTrans)
    quant075Trans <- quant.finder(vec=FhatTrans,q=0.75,evalvec = countsTrans)
    quant095Trans <- quant.finder(vec=FhatTrans,q=0.95,evalvec = countsTrans)


    thisEval <- which(countsTrans==trueCountTrans)
    if (length(thisEval)==1){
      if (thisEval==1){
        logScoreTrans <- log(FhatTrans[1])
      } else {
        logScoreTrans <- log(FhatTrans[thisEval]-FhatTrans[thisEval-1])
      }
    } else {
      logScoreTrans <- -Inf
    }

    infiniteLogScore <- is.infinite(logScoreTrans)
    if (infiniteLogScore){
      allf <- diff(c(0,FhatTrans))
      logScoreTrans <- log(min(allf[allf>0]))
    }
    infiniteLogScoreTrans <- c(infiniteLogScoreTrans,infiniteLogScore)

    transDf <- rbind(transDf,
                     c(trueCountTrans,
                       quant005Trans,
                       quant025Trans,
                       medTrans,
                       quant075Trans,
                       quant095Trans,
                       meanTrans,
                       CRPSTrans,
                       CRPSTransL,
                       CRPSTransU,
                       logScoreTrans))

  }


  transDf <- as.data.frame(transDf)
  names(transDf) <- c("TrueCount","Pred_q_0.05","Pred_q_0.25","Pred_q_0.5","Pred_q_0.75","Pred_q_0.95","Pred_mean","CRPS","CRPS_Underest","CRPS_Overest","logScore")

  s.mean <- function(x,i){mean(x[i])}

  (PropTrueCountIn050CI <- mean(transDf$Pred_q_0.25 <= transDf$TrueCount & transDf$Pred_q_0.75 >= transDf$TrueCount))
  set.seed(123)
  (PropTrueCountIn050CI.transBoot <- boot(transDf$Pred_q_0.25 <= transDf$TrueCount & transDf$Pred_q_0.75 >= transDf$TrueCount,s.mean,R=10^4,stype="i"))
  (PropTrueCountIn050CI.uncer <- quantile(PropTrueCountIn050CI.transBoot$t,c(0.05,0.95)))

  (PropTrueCountIn090CI <- mean(transDf$Pred_q_0.05 <= transDf$TrueCount & transDf$Pred_q_0.95 >= transDf$TrueCount))
  set.seed(123)
  (PropTrueCountIn090CI.transBoot <- boot(transDf$Pred_q_0.05 <= transDf$TrueCount & transDf$Pred_q_0.95 >= transDf$TrueCount,s.mean,R=10^4,stype="i"))
  (PropTrueCountIn090CI.uncer <- quantile(PropTrueCountIn090CI.transBoot$t,c(0.05,0.95)))

  transDf$x <- 1:nrow(transDf)

  transDfList[[m]] <- transDf

  covarageList50[[m]] <- c(PropTrueCountIn050CI,PropTrueCountIn050CI.uncer)
  covarageList90[[m]] <- c(PropTrueCountIn090CI,PropTrueCountIn090CI.uncer)

  print(m)
}

### Added

MAE_vec <- rep(NA,4)
for (m in 1:4){
  MAE_vec[m] <- mean(abs(transDfList[[m]]$TrueCount-transDfList[[m]]$Pred_q_0.5))
}

# Using the mean as point estimate
# > MAE_vec
# /home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/REGSpatialLeaveOutAsCV/
#   16.24294
# /home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/GAMSpatialLeaveOutAsCV/
#   13.27737
# /home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/REGSpatialLeaveOutAsCV/
#   695.78428
# /home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/GAMSpatialLeaveOutAsCV/
#   142.04719

# Using the median as point estimate
# /home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/REGSpatialLeaveOutAsCV/
#   13.81481
# /home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHooded/GAMSpatialLeaveOutAsCV/
#   12.44444
# /home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/REGSpatialLeaveOutAsCV/
#   179.18519
# /home/jullum/nr/project_stat/PointProcess/Martin/Prosjekter/FRINATEK-Point/PointProcessesInINLA/Results/finalHarps/GAMSpatialLeaveOutAsCV/
#   130.22222



names(MAE_vec) <- resultsBaseFolderVec
names(transDfList) <- resultsBaseFolderVec
names(covarageList50) <- resultsBaseFolderVec
names(covarageList90) <- resultsBaseFolderVec

(covarageList50)
(covarageList90)

coverage.hooded <- data.frame("50"=c(covarageList50[[1]][1],covarageList50[[2]][1]),
                              "90"=c(covarageList90[[1]][1],covarageList90[[2]][1]))
rownames(coverage.hooded) = c("LGCP","GAM NB")
colnames(coverage.hooded) = c("50%","90%")

coverage.harps <- data.frame("50"=c(covarageList50[[3]][1],covarageList50[[4]][1]),
                              "90"=c(covarageList90[[3]][1],covarageList90[[4]][1]))
rownames(coverage.harps) = c("LGCP","GAM NB")
colnames(coverage.harps) = c("50%","90%")

(coverage.hooded)
(coverage.harps)

transDfHooded <- as.data.frame(rbind(cbind(transDfList[[1]],ID="REG"),cbind(transDfList[[2]],ID="GAM NB")))
transDfHarps <- as.data.frame(rbind(cbind(transDfList[[3]],ID="REG"),cbind(transDfList[[4]],ID="GAM NB")))

# transDfHooded.Lower <- c(transDfHooded$Pred_q_0.05,transDfHooded$Pred_q_0.25)
# transDfHooded.Upper <- c(transDfHooded$Pred_q_0.95,transDfHooded$Pred_q_0.75)
# transDfHooded.Interval <- rep(c("90%","50%"),each=nrow(transDfHooded))
#
# transDfHooded <- rbind(transDfHooded,transDfHooded)
# transDfHooded$Lower <- transDfHooded.Lower
# transDfHooded$Upper <- transDfHooded.Upper
# transDfHooded$Interval <- transDfHooded.Interval

add=2

logbreaks <-c(0,3,10,30,100,400)+add
logbreaks.label <- logbreaks-add


transDfHooded2 = transDfHooded
transDfHooded2$TrueCount=transDfHooded2$TrueCount+add
transDfHooded2$Pred_q_0.05=transDfHooded2$Pred_q_0.05+add
transDfHooded2$Pred_q_0.25=transDfHooded2$Pred_q_0.25+add
transDfHooded2$Pred_q_0.75=transDfHooded2$Pred_q_0.75+add
transDfHooded2$Pred_q_0.95=transDfHooded2$Pred_q_0.95+add
transDfHooded2$Pred_q_0.5=transDfHooded2$Pred_q_0.5+add


gg <- ggplot(transDfHooded2,aes(x=x,y = Pred_q_0.5,color=ID))
gg <- gg + geom_crossbar(aes(ymin=Pred_q_0.05,ymax = Pred_q_0.95),position=position_dodge(width = 0.7),width=0.6,alpha=0.4)
gg <- gg + geom_crossbar(aes(ymin=Pred_q_0.25,ymax = Pred_q_0.75,fill=ID),position=position_dodge(width = 0.7),width=0.6,alpha=0.4)
dfdf <- data.frame(xmin=transDfHooded2$x-0.8/2, xmax=transDfHooded2$x+0.8/2,y= transDfHooded2$TrueCount,ID="Truth")
gg <- gg + geom_segment(aes(x=xmin,xend=xmax,y=y,yend=y,linetype=ID),data=dfdf,size=1,color=gg_color_hue(3)[2])
gg <- gg + theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16)) +
  xlab("Transect number") +
  ylab("Prediction")
#gg <- gg + scale_color_discrete(name="Model")
gg <- gg + scale_color_manual(values=c(gg_color_hue(3)[c(1,3,2)]),guide=FALSE)
#gg <- gg + scale_color_manual(name="Type",labels = c("90% CI","50% CI","Median"),values=c(gg_color_hue(2),"black"),
#                              guide = guide_legend(keyheight = 1, override.aes = list(color="black",linetype=c(1,1,1),size=c(0.4,0.4,0.4),fill=c(1,1,1),alpha=c(0.01,0.3,1)),order = 3))
gg <- gg + scale_fill_manual(name="Model",values=c(gg_color_hue(3)[c(1,3)]),guide = guide_legend(override.aes = list(color=gg_color_hue(3)[c(1,3)],fill=gg_color_hue(3)[c(1,3)],alpha=c(1,1)),order = 1),
                             labels=c("LGCP","GAM NB"))
gg <- gg + scale_linetype_manual(name="",values=1)#,guide_legend(order=2))#values=rev(c(gg_color_hue(2))),guide = guide_legend(override.aes = list(color=rev(gg_color_hue(2)),fill=rev(gg_color_hue(2)),alpha=1),order = 3))

gg


gglog <- gg + scale_y_log10(breaks = logbreaks,
                            labels = logbreaks.label)

gglog

## Use onyl with add=0
##pdf(file=file.path(paperFigureFolder,"transectRes_hooded.pdf"),width=12,height=6)
##gg
##dev.off()

pdf(file=file.path(paperFigureFolder,"transectRes_hooded_log3.pdf"),width=12,height=6)
gglog
dev.off()



table(transDfHooded2$ID[as.logical((transDfHooded2$TrueCount >= transDfHooded2$Pred_q_0.05)*(transDfHooded2$TrueCount <= transDfHooded2$Pred_q_0.95))])
table(transDfHooded2$ID[as.logical((transDfHooded2$TrueCount >= transDfHooded2$Pred_q_0.25)*(transDfHooded2$TrueCount <= transDfHooded2$Pred_q_0.75))])

#> table(transDfHooded2$ID[as.logical((transDfHooded2$TrueCount >= transDfHooded2$Pred_q_0.05)*(transDfHooded2$TrueCount <= transDfHooded2$Pred_q_0.95))])
##REG GAM
#26  18
#> table(transDfHooded2$ID[as.logical((transDfHooded2$TrueCount >= transDfHooded2$Pred_q_0.25)*(transDfHooded2$TrueCount <= transDfHooded2$Pred_q_0.75))])
#
#REG GAM
#16  11

#########

add=2

logbreaks <-c(0,10,100,1000,10000)+add
logbreaks.label <- logbreaks-add


transDfHarps2 = transDfHarps
transDfHarps2$TrueCount=transDfHarps2$TrueCount+add
transDfHarps2$Pred_q_0.05=transDfHarps2$Pred_q_0.05+add
transDfHarps2$Pred_q_0.25=transDfHarps2$Pred_q_0.25+add
transDfHarps2$Pred_q_0.75=transDfHarps2$Pred_q_0.75+add
transDfHarps2$Pred_q_0.95=transDfHarps2$Pred_q_0.95+add
transDfHarps2$Pred_q_0.5=transDfHarps2$Pred_q_0.5+add


gg <- ggplot(transDfHarps2,aes(x=x,y = Pred_q_0.5,color=ID))
gg <- gg + geom_crossbar(aes(ymin=Pred_q_0.05,ymax = Pred_q_0.95),position=position_dodge(width = 0.7),width=0.6,alpha=0.4)
gg <- gg + geom_crossbar(aes(ymin=Pred_q_0.25,ymax = Pred_q_0.75,fill=ID),position=position_dodge(width = 0.7),width=0.6,alpha=0.4)
dfdf <- data.frame(xmin=transDfHarps2$x-0.8/2, xmax=transDfHarps2$x+0.8/2,y= transDfHarps2$TrueCount,ID="Truth")
gg <- gg + geom_segment(aes(x=xmin,xend=xmax,y=y,yend=y,linetype=ID),data=dfdf,size=1,color=gg_color_hue(3)[2])
gg <- gg + theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16)) +
  xlab("Transect number") +
  ylab("Prediction")
#gg <- gg + scale_color_discrete(name="Model")
gg <- gg + scale_color_manual(values=c(gg_color_hue(3)[c(1,3,2)]),guide=FALSE)
#gg <- gg + scale_color_manual(name="Type",labels = c("90% CI","50% CI","Median"),values=c(gg_color_hue(2),"black"),
#                              guide = guide_legend(keyheight = 1, override.aes = list(color="black",linetype=c(1,1,1),size=c(0.4,0.4,0.4),fill=c(1,1,1),alpha=c(0.01,0.3,1)),order = 3))
gg <- gg + scale_fill_manual(name="Model",values=c(gg_color_hue(3)[c(1,3)]),guide = guide_legend(override.aes = list(color=gg_color_hue(3)[c(1,3)],fill=gg_color_hue(3)[c(1,3)],alpha=c(1,1)),order = 1),
                             labels=c("LGCP","GAM NB"))
gg <- gg + scale_linetype_manual(name="",values=1)#,guide_legend(order=2))#values=rev(c(gg_color_hue(2))),guide = guide_legend(override.aes = list(color=rev(gg_color_hue(2)),fill=rev(gg_color_hue(2)),alpha=1),order = 3))

gg


gglog <- gg + scale_y_log10(breaks = logbreaks,
                            labels = logbreaks.label)

gglog

## Use only with add=0
##pdf(file=file.path(paperFigureFolder,"transectRes_Harps.pdf"),width=12,height=6)
##gg
##dev.off()

pdf(file=file.path(paperFigureFolder,"transectRes_harps_log.pdf"),width=12,height=6)
gglog
dev.off()

table(transDfHarps2$ID[as.logical((transDfHarps2$TrueCount >= transDfHarps2$Pred_q_0.05)*(transDfHarps2$TrueCount <= transDfHarps2$Pred_q_0.95))])
table(transDfHarps2$ID[as.logical((transDfHarps2$TrueCount >= transDfHarps2$Pred_q_0.25)*(transDfHarps2$TrueCount <= transDfHarps2$Pred_q_0.75))])
table(transDfHarps2$ID[as.logical((transDfHarps2$TrueCount >= transDfHarps2$Pred_q_0.05)*(transDfHarps2$TrueCount <= transDfHarps2$Pred_q_0.95))])

#REG GAM
#24  18
#> table(transDfHarps2$ID[as.logical((transDfHarps2$TrueCount >= transDfHarps2$Pred_q_0.25)*(transDfHarps2$TrueCount <= transDfHarps2$Pred_q_0.75))])#
#
#REG GAM
#14  11



#####################################################################################################################################################

