#' Load data
#'
#' Loads all data used for pup production
#'
#' @param filename Shape file with pup counts
#' @param lt ?
#' @return
#' @keywords
#' @export
#' @examples
#' loadData()

loadData <- function() {

  #Flight data
  #find names in flight data folder
  #eo = data.frame()
  #fnames = list.files("data/FlightPhotos")
  #for(i in 1:length(fnames)){
  #  if(grepl(".txt",fnames[i])){
  #    eotemp = readeo(paste0("data/FlightPhotos/",fnames[i]))
  #    eo = rbind(eo,eotemp$data)
  #  }
  #names(eo) = names(eotemp$data)

#  }
  #eo.35002 = readeo("data/FlightPhotos/EO_35002C_C098_20180327A.txt")
  #eo.35003 = readeo("data/FlightPhotos/EO_35003C_C098_20180327B.txt")
  #eo.35004 = readeo("data/FlightPhotos/EO_35004_C098_20170328.txt")

  #load("H:/Work/HI/MarineMammals/rPupps/data/FlightPhotos/eo.35002.RData")
  #load("H:/Work/HI/MarineMammals/rPupps/data/FlightPhotos/eo.35003.RData")
  #load("H:/Work/HI/MarineMammals/rPupps/data/FlightPhotos/eo.35003.RData")

  load('data/FlightPhotos/eo.35002.RData')
  load('data/FlightPhotos/eo.35003.RData')
  load('data/FlightPhotos/eo.35004.RData')

  areas <- read.csv('data/FlightPhotos/photo.areas.csv', stringsAsFactors=F)
  areas$IMAGE_ID <- gsub('.tif', '', areas$IMAGE_ID, fixed=T)

  #eo = merge(eo, areas, by.x='Event.ID', by.y='IMAGE_ID', all.x=T, all.y=F)


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

  # READ PUP COUNTS FROM FILE
  mp <- readCounts("data/CountFiles/MP_MP_count.shp")
  ll <- readCounts("data/CountFiles/LottaL_LL_count.shp")
  lld <- readCounts("data/CountFiles/LottaL_start_LL_count.shp")
  mpd <- readCounts("data/CountFiles/MP_MP_double.shp")

  names(mp)[length(names(mp))] <- 'NOTES'
  names(mpd)[length(names(mpd))] <- 'NOTES'

  all.ind <- rbind(ll, mp)
  all.ind <- merge(all.ind, eo.all[,c(1,22,24)], by.x='PHOTO_ID', by.y='Event.ID', all.x=T, all.y=F)
  all.ind$Date <- as.POSIXct(strptime(format(all.ind$UTC, '%Y-%m-%d'), '%Y-%m-%d'), tz='UTC')

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

  double <- XLConnect::readWorksheetFromFile('data/countfiles/DoubleCounting.xlsx', 1)

  mp.sum <- sumImageList()
  mpd.sum <- sumImageList(mpd)
  ll.sum <- sumImageList(ll)
  lld.sum <- sumImageList(lld)

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


}


#' Read the pup count shape files
#'
#' Performs an inverse logit transforms the input number
#' @param countfile Shape file with pup counts
#' @param to.spatial ?
#' @param to.ll ?
#' @return
#' @keywords
#' @export
#' @examples
#' readCounts()

readCounts <- function(countfile=file.choose(),
                        to.spatial=T,
                        to.ll=F) {
  if(to.spatial) {
    shp <- rgdal::readOGR(countfile, verbose=F, stringsAsFactors=F)
  } else {
    shp <- shapefiles::read.shapefile(gsub('.shp', '', countfile))
    if(to.ll) {
      shp <- rgdal::spTransform(shp, CRS('+proj=longlat +ellps=WGS84'))
    }
  }
  shp
}

#' Read the pup count shape files
#'
#' Performs an inverse logit transforms the input number
#' @param counts Shape file with pup counts
#' @param eo ?
#' @param long.form ?
#' @param outform ?
#' @return
#' @keywords
#' @export
#' @examples
#' sumCounts()

sumCounts <- function(counts=cnt,
                       eo=eo.dat,
                       long.form=F,
                       outfile=NA) {
  counts$SPECIES_ID <- factor(counts$SPECIES_ID, levels=c('Harp', 'Hood'))
  counts$STATUS <- factor(counts$STATUS, levels=c('Nb', 'TF', 'U', 'Bb', 'Mp', 'Fm', 'Ad'))
  cnt.eo <- eo[unique(match(counts$PHOTO_ID, eo$Event.ID)),]
  by.spp.status <- aggregate(counts$PHOTO_ID, list(counts$PHOTO_ID, counts$SPECIES_ID, counts$STATUS), length, drop=F)
  names(by.spp.status) <- c('Event.ID', 'SPECIES', 'STATUS', 'COUNT')
  by.spp.status$Event.ID <- as.character(by.spp.status$Event.ID)
  by.spp.status$SPECIES <- as.character(by.spp.status$SPECIES)
  by.spp.status$STATUS <- as.character(by.spp.status$STATUS)
  cnt.eo$READER <- rep(counts$READER[1], nrow(cnt.eo))
  cnt.eo$HoodAd <- cnt.eo$HoodUnk <- cnt.eo$HoodFam <- cnt.eo$HoodMP <- cnt.eo$HoodSol <-
    cnt.eo$HarpAd <- cnt.eo$HarpUnk <- cnt.eo$HarpTF <- cnt.eo$HarpNb <- rep(0, nrow=cnt.eo)
  HarpNb <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Harp' & by.spp.status$STATUS=='Nb')]
  HarpTF <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Harp' & by.spp.status$STATUS=='TF')]
  HarpUnk <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Harp' & by.spp.status$STATUS=='U')]
  HarpAd <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Harp' & by.spp.status$STATUS=='Ad')]
  HoodSol <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Hood' & by.spp.status$STATUS=='Bb')]
  HoodMP <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Hood' & by.spp.status$STATUS=='Mp')]
  HoodFam <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Hood' & by.spp.status$STATUS=='Fm')]
  HoodUnk <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Hood' & by.spp.status$STATUS=='U')]
  HoodAd <- by.spp.status$COUNT[which(by.spp.status$SPECIES=='Hood' & by.spp.status$STATUS=='Ad')]
  if(length(HarpNb)!=0) cnt.eo$HarpNb <- HarpNb
  if(length(HarpTF)!=0) cnt.eo$HarpTF <- HarpTF
  if(length(HarpUnk)!=0) cnt.eo$HarpUnk <- HarpUnk
  if(length(HarpAd)!=0) cnt.eo$HarpAd <- HarpAd
  if(length(HoodSol)!=0) cnt.eo$HoodSol <- HoodSol
  if(length(HoodMP)!=0) cnt.eo$HoodMP <- HoodMP
  if(length(HoodFam)!=0) cnt.eo$HoodFam <- HoodFam
  if(length(HoodUnk)!=0) cnt.eo$HoodUnk <- HoodUnk
  if(length(HoodAd)!=0) cnt.eo$HoodAd <- HoodAd
  cnt.eo$HarpPup <- cnt.eo$HarpNb+cnt.eo$HarpTF+cnt.eo$HarpUnk
  cnt.eo$HoodPup <- cnt.eo$HoodSol+cnt.eo$HoodMP+cnt.eo$HoodFam+cnt.eo$HoodUnk
  if(long.form) {
    if(is.na(outfile)) {
      cnt.eo
    } else {
      write.table(cnt.eo, file=outfile, sep=',', dec='.', col.names=T, row.names=F)
    }
  } else {
    if(is.na(outfile)) {
      cnt.eo[,c(c(1:3), c(19:20), c(22:36))]
    } else {
      write.table(cnt.eo[,c(c(1:3), c(19:20), c(22:36))], file=outfile, sep=',', dec='.', col.names=T, row.names=F)
    }
  }
}


#' TITLE?
#'
#' WHAT DOES IT DO?
#'
#' @param filename Shape file with pup counts
#' @param lt ?
#' @return
#' @keywords
#' @export
#' @examples
#' readeo()

readeo <- function(filename=file.choose(), lt=+1) {
  dat <- readLines(filename)
  header <- dat[c(1:72)]
  dat <- dat[-c(1:72)]
  dat.ll <- lapply(dat, function(x) {
    dd <- unlist(strsplit(x, ' '))
    dd[which(nchar(dd)>0)]
  })
  dat.df <- as.data.frame(do.call('rbind', dat.ll), stringsAsFactors=F)
  names(dat.df) <- header[grep('Output:', header)+c(1:length(dat.df))]
  names(dat.df) <- lapply(names(dat.df), function(x) unlist(strsplit(x, '[(]'))[1])
  names(dat.df) <- gsub('# ', '', names(dat.df))
  names(dat.df) <- trimws(names(dat.df))
  names(dat.df) <- gsub(' ', '.', names(dat.df), fixed=T)
  names(dat.df) <- gsub('-', '.', names(dat.df), fixed=T)
  for(i in 2:length(dat.df)) dat.df[,i] <- as.numeric(dat.df[,i])
  dat.df$UTC <- as.POSIXct((dat.df$Week*(7*86400))+dat.df$ToW, origin='1980-01-06 00:00:00', tz='UTC')
  dat.df$lt <- dat.df$UTC+(lt*3600)
  list(meta=header, data=dat.df)
}

