source('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/scripts/read.eo.R')
source('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/scripts/find.transects.R')
source('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/scripts/image.dim.R')
source('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/scripts/sum.image.R')
source('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/scripts/read.counts.R')
source('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/scripts/kingsley.R')
source('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/scripts/stageCorrection.R')
source('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/scripts/degree.R')

source('R/input_functions.R')
source('R/processing_functions.R')


require(shapefiles)
require(rgdal)
library(XLConnect)


#############################################
# This could go into a loadData-function....
#############################################

#load('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/data/eo.35002.RData')
#load('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/data/eo.35003.RData')
#load('H:/Work/HI/MarineMammals/PupProduction/files4TA/files4TA/data/eo.35004.RData')

# READ RAW DATA FILES FROM TERRATEC

eo.35002 = readeo("data/FlightPhotos/EO_35002C_C098_20180327A.txt")
eo.35003 = readeo("data/FlightPhotos/EO_35003C_C098_20180327B.txt")
eo.35004 = readeo("data/FlightPhotos/EO_35004_C098_20170328.txt")

areas <- read.csv('data/FlightPhotos/photo.areas.csv', stringsAsFactors=F)
areas$IMAGE_ID <- gsub('.tif', '', areas$IMAGE_ID, fixed=T)

eo.35002$data <- merge(eo.35002$data, areas, by.x='Event.ID', by.y='IMAGE_ID', all.x=T, all.y=F)

eo.35003$data <- merge(eo.35003$data, areas, by.x='Event.ID', by.y='IMAGE_ID', all.x=T, all.y=F)

eo.35004$data <- merge(eo.35004$data, areas, by.x='Event.ID', by.y='IMAGE_ID', all.x=T, all.y=F)

eo.all <- rbind(eo.35002$data, eo.35003$data, eo.35004$data)

batches <- unique(unlist(lapply(eo.all$transect, function(x) unlist(strsplit(x, '_'))[1])))
batches <- batches[which(!is.na(batches))]
t.num <- c(1:100)
tlevels <- expand.grid(t.num, batches)[,c(2,1)]
tlevels <- apply(tlevels, 1, function(x) paste(x[1], x[2], sep='_'))
tlevels <- gsub(' ', '', tlevels)
eo.all$transect <- factor(eo.all$transect, levels=tlevels)
eo.all$transect <- droplevels(eo.all$transect)

############################
# End of load data function
############################

# READ PUP COUNTS FROM FILE
mp <- readCounts("data/CountFiles/MP_MP_count.shp")
ll <- readCounts("data/CountFiles/LottaL_LL_count.shp")
lld <- readCounts("data/CountFiles/LottaL_start_LL_count.shp")
mpd <- readCounts("data/CountFiles/MP_MP_double.shp")

names(mp)[length(names(mp))] <- 'NOTES'
names(mpd)[length(names(mpd))] <- 'NOTES'


###

all.ind <- rbind(ll, mp)
all.ind <- merge(all.ind, eo.all[,c(1,22,24)], by.x='PHOTO_ID', by.y='Event.ID', all.x=T, all.y=F)
all.ind$Date <- as.POSIXct(strptime(format(all.ind$UTC, '%Y-%m-%d'), '%Y-%m-%d'), tz='UTC')

mp.counts <- sumImageList(ims=eo.all$Event.ID)
ll.counts <- sumImageList(ll, ims=eo.all$Event.ID)

pup.counts <- data.frame(PHOTO_ID=mp.counts$PHOTO_ID,
                         HarpPup=mp.counts$HarpPup+ll.counts$HarpPup,
                         HoodPup=mp.counts$HoodPup+ll.counts$HoodPup)

pup.counts$PHOTO_ID <- as.character(pup.counts$PHOTO_ID)

## Merge with image data:
pup.counts <- merge(pup.counts, eo.all[,which(names(eo.all) %in% c("Event.ID", "Grid.Easting", "Grid.Northing", "Latitude", "Longitude", "UTC", "lt", "H", "W", "transect", "area"))], by.x='PHOTO_ID', by.y='Event.ID', all.x=T, all.y=F)


####


double <- XLConnect::readWorksheetFromFile('Q:/ishavssel/PhotoReading/DoubleCounting.xlsx', 1)

mp.sum <- sum.image.list()
mpd.sum <- sum.image.list(mpd)
ll.sum <- sum.image.list(ll)
lld.sum <- sum.image.list(lld)

harp.double <- data.frame(PHOTO_ID=double$Image_ID,
                          LL=ll.sum$HarpPup+lld.sum$HarpPup,
                          MP=mp.sum$HarpPup+mpd.sum$HarpPup)

hood.double <- data.frame(PHOTO_ID=double$Image_ID,
                          LL=ll.sum$HoodPup+lld.sum$HoodPup,
                          MP=mp.sum$HoodPup+mpd.sum$HoodPup)
lm.harp <- lm(LL~MP, data=harp.double)
pred.harp <- predict(lm.harp, newdata=data.frame(MP=c(0:60)),
                     se.fit=T, interval='prediction')

lm.hood <- lm(LL~MP, data=hood.double)
pred.hood <- predict(lm.hood, newdata=data.frame(MP=c(0:60)),
                     se.fit=T, interval='prediction')


######

pup.counts$survey <- rep(1, nrow(pup.counts))
pup.counts$survey[grep('35004', pup.counts$PHOTO_ID)] <- 2
which.batch <- grep('35004', pup.counts$PHOTO_ID)
which.transects <- which(as.numeric(unlist(lapply(as.character(pup.counts$transect), function(x) unlist(strsplit(x, '_'))[2])))>23)
pup.counts$survey[intersect(which.batch, which.transects)] <- 3

survey <- pup.counts[which(!is.na(pup.counts$transect)),]


#########

#Staging data
harpCorrect <- stageCorrection(hoodStage=NA)
harpTab <- harpCorrect$harpStage[,c(1:9)]
names(harpTab)[c(2, 4,5,6,7,8)] <- c('Newborn', 'Thin white', 'Fat white', 'Grey coat', 'Ragged jacket', 'Beater')
harpTab$Date <- format(harpTab$Date, '%B %d')

hoodCorrect <- stageCorrection(harpStage=NA, repfile='correctHoodStageProps2012.txt')
hoodTab <- hoodCorrect$hoodStage[,c(1:7)]
names(hoodTab)[c(2:6)] <- c('Parturient female', 'Newborn', 'Thin blueback', 'Fat blueback', 'Solitary blueback')
hoodTab$Date <- format(hoodTab$Date, '%B %d')

flightStages <- table(all.ind$Date, all.ind$STATUS, all.ind$SPECIES_ID)

hoodFlight <- as.data.frame.matrix(flightStages[,,2])
hoodFlight <- rbind(hoodFlight, apply(hoodFlight, 2, mean))
hoodFlight$FmMpNb <-hoodFlight$Fm+hoodFlight$Mp+hoodFlight$Nb

hoodFl <- data.frame(Date=paste0('March ', paste(substr(dimnames(flightStages)[[1]], 9,10), collapse='-'), '*'),
                     Parturient=NA, Newborn=hoodFlight$FmMpNb[nrow(hoodFlight)],
                     Thin=0, Fat=0, Solitary=hoodFlight$Bb[nrow(hoodFlight)])
hoodFl[,-1] <- as.integer(hoodFl[,-1])
hoodFl$Total <- apply(hoodFl[,-1], 1, sum, na.rm=T)


##############

tmp <- stageCorrection(hoodStage=NA, plotting=T)


###########

########

tmp <- stageCorrection(harpStage=NA, repfile='correctHoodStageProps2012.txt', plotting=T)

######

transSum <- aggregate(survey$UTC, list(survey$transect), min)
names(transSum) <- c('Transect', 'StartTime')
transSum$EndTime <- aggregate(survey$UTC, list(survey$transect), max)$x
transSum$meanLat <- aggregate(survey$Latitude, list(survey$transect), mean, na.rm=T)$x
transSum$minLon <- aggregate(survey$Longitude, list(survey$transect), min, na.rm=T)$x
transSum$maxLon <- aggregate(survey$Longitude, list(survey$transect), max, na.rm=T)$x
transSum$nHarp <- aggregate(survey$HarpPup, list(survey$transect), sum)$x
transSum$nHood <- aggregate(survey$HoodPup, list(survey$transect), sum)$x
transSum$nPhotos <- aggregate(survey$HarpPup, list(survey$transect), length)$x
transSum$survey <- aggregate(survey$survey, list(survey$transect), mean, na.rm=T)$x
##transSum$stratum <- aggregate(survey$stratum, list(survey$transect), mean, na.rm=T)$x
##transSum$stratum4 <- aggregate(survey$stratum4, list(survey$transect), function(x) x[1])$x

## Mean transect spacing for each survey (in NM):
meanDists <- aggregate(transSum$meanLat, list(transSum$survey), function(x) mean(abs(diff(x)))*60)
names(meanDists) <- c('survey', 'Spacing')

transTab <- transSum
transTab$Date <- format(transTab$StartTime, '%B %d')

latrange <- aggregate(transTab$meanLat, list(transTab$Date), min)
names(latrange) <- c('Date', 'min')
latrange$max <- aggregate(transTab$meanLat, list(transTab$Date), max)$x


meanLat <- degree(transTab$meanLat, transTab$minLon, digits=0)$lat
minLon <-  degree(transTab$meanLat, transTab$minLon, digits=0)$lon
maxLon <-  degree(transTab$meanLat, transTab$maxLon, digits=0)$lon
##meanLat <- gsub('°', '\\degree', meanLat, fixed=T)
transTab$StartTime <- format(transTab$StartTime, '%H:%M')
transTab$EndTime <- format(transTab$EndTime, '%H:%M')
transTab$meanLat <- meanLat
transTab$minLon <-  minLon
transTab$maxLon <-  maxLon
##transTab <- transTab[,-c(10:12)]
names(transTab)[c(2:9)] <- c('Start', 'End', 'Latitude', 'West', 'East', 'Harps', 'Hoods', 'nphotos')
transTab <- transTab[,c(1,11,c(2:10))]


#############

#########

reader2 <- match(survey$PHOTO_ID, mp$PHOTO_ID)
reader2[which(!is.na(reader2))] <- '2'
reader2[which(is.na(reader2))] <- ''
reader1 <- match(survey$PHOTO_ID, ll$PHOTO_ID)
reader1[which(!is.na(reader1))] <- '1'
reader1[which(is.na(reader1))] <- ''

reader <- paste0(reader1, reader2)
reader[which(reader=='')] <- '1'
survey$reader <- reader

survey$HarpPup.rc <- survey$HarpPup
survey$HarpPup.rc[which(survey$reader=='2')] <- predict(lm.harp, newdata=data.frame(MP=survey$HarpPup[which(survey$reader=='2')]))

survey$HoodPup.rc <- survey$HoodPup
survey$HoodPup.rc[which(survey$reader=='2')] <- predict(lm.hood, newdata=data.frame(MP=survey$HoodPup[which(survey$reader=='2')]))

###########


#######

meanDists <- aggregate(transSum$meanLat, list(transSum$survey), function(x) mean(abs(diff(x)))*60)
names(meanDists) <- c('survey', 'Spacing')


## Survey 1 (March 27, 2018):
survey1 <- survey[which(survey$survey==1),]
EstHarp.1 <- kingsley((survey1$Grid.Easting-mean(survey1$Grid.Easting))/1852,
                      survey1$HarpPup,
                      survey1$area,
                      survey1$transect,
                      meanDists$Spacing[1], 1)

EstHood.1 <- kingsley((survey1$Grid.Easting-mean(survey1$Grid.Easting))/1852,
                      survey1$HoodPup,
                      survey1$area,
                      survey1$transect,
                      meanDists$Spacing[1], 1)


## Survey 2 (March 28, 2018, northward leg):
survey2 <- survey[which(pup.counts$survey==2),]
EstHarp.2 <- kingsley((survey2$Grid.Easting-mean(survey2$Grid.Easting))/1852,
                      survey2$HarpPup,
                      survey2$area,
                      survey2$transect,
                      meanDists$Spacing[2], 1)

EstHood.2 <- kingsley((survey2$Grid.Easting-mean(survey2$Grid.Easting))/1852,
                      survey2$HoodPup,
                      survey2$area,
                      survey2$transect,
                      meanDists$Spacing[2], 1)


## Survey 3 (March 28, 2018, southward leg):
survey3 <- survey[which(survey$survey==3),]
EstHarp.3 <- kingsley((survey3$Grid.Easting-mean(survey3$Grid.Easting))/1852,
                      survey3$HarpPup,
                      survey3$area,
                      survey3$transect,
                      meanDists$Spacing[3], 1)

EstHood.3 <- kingsley((survey3$Grid.Easting-mean(survey3$Grid.Easting))/1852,
                      survey3$HoodPup,
                      survey3$area,
                      survey3$transect,
                      meanDists$Spacing[3], 1)

allEst <- rbind(EstHarp.1, EstHarp.2, EstHarp.3,
                EstHood.1, EstHood.2, EstHood.3)
allEst$Species <- c(rep('Harp', 3), rep('Hood', 3))
allEst$Survey <- rep(c(c(1:3)), 2)
allEst$lowerCI <- allEst$N-allEst$SE
allEst$upperCI <- allEst$N+allEst$SE
allEst <- allEst[,c(4,5,1,2,6,7,3)]
allEst$N <- as.integer(round(allEst$N))
allEst$SE <- round(allEst$SE, 1)
allEst$CV <- round(allEst$CV, 1)
allEst$lowerCI <- as.integer(round(allEst$lowerCI))
allEst$upperCI <- as.integer(round(allEst$upperCI))


#######

survey$merged.transect <- survey$transect
survey$stratum <- rep(NA, nrow(survey))
survey$stratum[which(survey$survey==1 & survey$Latitude<71.84)] <- 1
survey$stratum[which(survey$survey==2 & survey$Latitude>71.84)] <- 2
survey$stratum[which(survey$transect=='35004_22')] <- NA
survey$stratum[which(survey$transect=='35004_26')] <- 2
survey$stratum[which(survey$transect=='35004_30')] <- 2
survey$stratum[which(survey$transect=='35004_32')] <- 2

survey$merged.transect[which(survey$transect=='35004_26')] <- '35004_19'
survey$merged.transect[which(survey$transect=='35004_30')] <- '35004_16'
survey$merged.transect[which(survey$transect=='35004_32')] <- '35004_15'

# Omit overlapping images
survey$stratum[which(survey$transect=='35004_19' & survey$Longitude>(-18.886))] <- NA

survey$stratum[which(survey$transect=='35004_16' & survey$Longitude>(-19.02))] <- NA

survey$stratum[which(survey$transect=='35004_15' & survey$Longitude>(-19.4))] <- NA

survey$stratum[which(survey$survey==3 &
                       (survey$transect!='35004_26' &
                          survey$transect!='35004_30' &
                          survey$transect!='35004_32'))] <- 3

survey$stratum[which(survey$transect=='35004_22')] <- 3
survey$stratum4 <- survey$stratum>=2

transSum$stratum <- aggregate(survey$stratum, list(survey$transect), mean, na.rm=T)$x

transSum$stratum[which(is.nan(transSum$stratum))] <- NA


strataSum <- aggregate(survey$Latitude, list(survey$stratum), min)
names(strataSum) <- c('Stratum', 'minLat')
strataSum$maxLat <- aggregate(survey$Latitude, list(survey$stratum), max)$x

## Mean transect spacing for each stratum (in NM):
strataSum$spacing <- aggregate(transSum$meanLat, list(transSum$stratum), function(x) mean(abs(diff(x)))*60)$x


######

par(mfrow=c(2,2))
plot(Latitude~Longitude, main='Stratum 1 (March 27)', data=survey, pch=15, col='lightgrey')
points(Latitude~Longitude, survey[which(survey$stratum==1),], pch=15, cex=0.5, col=rgb(0,0,1,0.3))

plot(Latitude~Longitude, main='Stratum 2 (March 28, northward (eastward extensions))', data=survey, pch=15, col='lightgrey')
points(Latitude~Longitude, survey[which(survey$stratum==2),], pch=15, cex=0.5, col=rgb(0,0,1,0.3))

plot(Latitude~Longitude, main='Stratum 3 (March 28, southward)', data=survey, pch=15, col='lightgrey')
points(Latitude~Longitude, survey[which(survey$stratum==3),], pch=15, cex=0.5, col=rgb(0,0,1,0.3))

plot(Latitude~Longitude, main='Stratum 4 (March 28, 1 nm spacing)', data=survey, pch=15, col='lightgrey')
points(Latitude~Longitude, survey[which(survey$stratum4),], pch=15, cex=0.5, col=rgb(0,0,1,0.3))


#######

