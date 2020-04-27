#' WHAT DOES THIS FUNCTION DO?
#'
#' WHAT DOES THIS FUNCTION DO?
#' @param countfile ?
#' @param im ?
#' @return
#' @keywords
#' @export
#' @examples
#' sumImages()

sumImage <- function(countfile='Q:/ishavssel/PhotoReading/CountFiles/MP_MP_count.shp', im='35004_55549') {

  if(class(countfile)=='character') {
    cnt <- rgdal::readOGR(countfile)
  } else {
    cnt <- countfile
  }

  cnt$PHOTO_ID <- as.character(cnt$PHOTO_ID)
  cnt <- cnt[grep(im, cnt$PHOTO_ID),]

  if(nrow(cnt)>0) {
    harp <- length(which(cnt$SPECIES_ID=='Harp'))
    hood <- length(which(cnt$SPECIES_ID=='Hood'))
    Ad <- length(which(cnt$STATUS=='Ad'))
    Bb <- length(which(cnt$STATUS=='Bb'))
    Fm <- length(which(cnt$STATUS=='Fm'))
    Mp <- length(which(cnt$STATUS=='Mp'))
    Nb <- length(which(cnt$STATUS=='Nb'))
    TF <- length(which(cnt$STATUS=='TF'))
    U <- length(which(cnt$STATUS=='U'))
    HarpPup <- length(which(cnt$SPECIES_ID=='Harp' &
                              (cnt$STATUS=='Mp' |
                                 cnt$STATUS=='Nb' |
                                 cnt$STATUS=='TF')))
    HoodPup <- length(which(cnt$SPECIES_ID=='Hood' &
                              (cnt$STATUS=='Mp' |
                                 cnt$STATUS=='Nb' |
                                 cnt$STATUS=='TF' |
                                 cnt$STATUS=='Bb' |
                                 cnt$STATUS=='Fm')))
    data.frame(PHOTO_ID=cnt$PHOTO_ID[1], Harp=harp,Hood=hood,Ad=Ad, Bb=Bb,
               Fm=Fm,Mp=Mp,Nb=Nb,TF=TF,U=U,
               HarpPup=HarpPup,HoodPup=HoodPup)
  } else {
    data.frame(PHOTO_ID=im, Harp=0,Hood=0,Ad=0, Bb=0,
               Fm=0,Mp=0,Nb=0,TF=0,U=0,
               HarpPup=0,HoodPup=0)
  }
}


#' WHAT DOES THIS FUNCTION DO?
#'
#' WHAT DOES THIS FUNCTION DO?
#' @param countfile ?
#' @param im ?
#' @return
#' @keywords
#' @export
#' @examples
#' sumImageList()

sumImageList <- function(countfile=mp, ims=double$Image_ID) {
  require(pbapply)
  do.call('rbind', pblapply(ims, function(x) sumImage(countfile, x)))
}


stageCorrection <- function(harpStage='data/StagingData/HarpStaging.txt',
                            hoodStage='data/StagingData/HoodStaging.txt',
                            repfolder='data/StagingData/',
                            repfile='correctStageProps2012.txt',
                            photoStage=all.ind,
                            DatePhoto = c(27,28),
                            plotting=F, return.phst=F) {

  REP <- read.csv(paste(repfolder, repfile, sep='/'))

  if(!is.na(harpStage)) {
    harpStage <- read.csv(harpStage)

    harpStage$Total <- apply(harpStage[,-1], 1, sum)
    harpStage$nn1 <- harpStage$New+harpStage$Yellow
    harpStage$nn2 <- harpStage$Thin
    harpStage$nn3 <- harpStage$Fat+harpStage$Grey
    harpStage$nn1.2 <- harpStage$nn1+harpStage$nn2

    harpStage$p1 <- harpStage$nn1/harpStage$Total
    harpStage$p2 <- harpStage$nn2/harpStage$Total
    harpStage$p3 <- harpStage$nn3/harpStage$Total
    harpStage$p1.2 <- harpStage$nn1.2/harpStage$Total
    harpStage$Date <- as.POSIXct(strptime(as.character(harpStage$Date), '%Y-%m-%d'), tz='UTC')
    harpStage$numDate <- as.numeric(format(harpStage$Date, '%d'))

    minDiff <- which.min(abs(harpStage$p1-REP$p1)+abs(harpStage$p2-REP$p2)+abs(harpStage$p3-REP$p3))

    dayAdj <- REP$t_tot[minDiff]-harpStage$numDate
    flightAdj <- DatePhoto+dayAdj
  }

  if(!is.na(hoodStage)) {
    hoodStage <- read.csv(hoodStage)

    hoodStage$Total <- apply(hoodStage[,-1], 1, sum)
    hoodStage$nn1 <- hoodStage$New+hoodStage$Thin
    hoodStage$nn2 <- hoodStage$Fat
    hoodStage$nn3 <- hoodStage$Solitary
    hoodStage$nn1.2 <- hoodStage$nn1+hoodStage$nn2

    hoodStage$p1 <- hoodStage$nn1/hoodStage$Total
    hoodStage$p2 <- hoodStage$nn2/hoodStage$Total
    hoodStage$p3 <- hoodStage$nn3/hoodStage$Total
    hoodStage$p1.2 <- hoodStage$nn1.2/hoodStage$Total
    hoodStage$Date <- as.POSIXct(strptime(as.character(hoodStage$Date), '%Y-%m-%d'), tz='UTC')
    hoodStage$numDate <- as.numeric(format(hoodStage$Date, '%d'))

    minDiffSt <- which.min(abs(hoodStage$p1-REP$p1)+abs(hoodStage$p2-REP$p2)+abs(hoodStage$p3-REP$p3))

    dayAdjSt <- REP$t_tot[minDiffSt]-hoodStage$numDate

    photoStage <- table(photoStage$Date, photoStage$STATUS, photoStage$SPECIES_ID)
    hood.phst <- photoStage[,,2]

    hood.phst <- as.data.frame.matrix(hood.phst)
    hood.phst <- cbind(Date=sort(unique(all.ind$Date)), hood.phst)
    hood.phst$Date <- as.POSIXct(strptime(row.names(hood.phst), '%Y-%m-%d'), tz='UTC')
    hood.phst <- hood.phst[,match(c('Date', 'Nb', 'Mp', 'Fm', 'Bb'), names(hood.phst))]
    hood.phst$NbMpFm <- hood.phst$Nb+hood.phst$Mp+hood.phst$Fm

    hood.phst <- rbind(hood.phst, data.frame(Date=mean(hood.phst$Date),
                                             Nb=mean(hood.phst$Nb),
                                             Mp=mean(hood.phst$Mp),
                                             Fm=mean(hood.phst$Fm),
                                             Bb=mean(hood.phst$Bb),
                                             NbMpFm=mean(hood.phst$NbMpFm)))
    hood.phst$nn1 <- hood.phst$Nb
    hood.phst$nn1.2 <- hood.phst$NbMpFm
    hood.phst$nn3 <- hood.phst$Bb
    hood.phst$nn <- hood.phst$nn1.2+hood.phst$nn3
    hood.phst$p1.2 <- hood.phst$nn1.2/hood.phst$nn
    hood.phst$p3 <- hood.phst$nn3/hood.phst$nn

    minDiff <- which.min(abs(hood.phst$p1.2[3]-(REP$p1+REP$p2))+abs(hood.phst$p3[3]-REP$p3))

    flightAdj <- REP$t_tot[minDiff]-DatePhoto
    flightAdj <- DatePhoto+mean(flightAdj)

    flightAdjSt <- REP$t_tot[minDiffSt]-DatePhoto
    flightAdjSt <- DatePhoto+mean(flightAdjSt)

  }


  if(plotting) {
    par(mar=c(5, 4, 12, 2) + 0.1, mfrow=c(1,2), cex=0.8)
    matplot(REP[,1], REP[,c(2:4)], ylim=c(0,1),
            type='l', lty=1, lwd=2,
            col=c('red', 'blue', 'green'),
            xlab='Days since March 1st (2012)',
            ylab='Proportion', axes=F)
    axis(1)
    axis(2)
    abline(v=par('usr')[1])
    abline(h=par('usr')[3])

    if(class(harpStage)!='logical') {
      abline(v=c(harpStage$numDate, DatePhoto), lty=2, col='slategrey')
      abline(v=REP$t_tot[minDiff])
      points(rep(REP$t_tot[minDiff], 3),
             harpStage[,match(c('p1', 'p2', 'p3'), names(harpStage))],
             pch=21, bg=c(2,4,3), cex=1.5)
    }

    if(class(hoodStage)!='logical') {
      abline(v=c(hoodStage$numDate, DatePhoto), lty=2, col='slategrey')
      abline(v=REP$t_tot[minDiff], col='grey')
      lines(REP[,1], REP[,2]+REP[,3], lwd=1, col='purple')
      points(rep(REP$t_tot[minDiff], 2),
             hood.phst[3,match(c('p1.2', 'p3'), names(hood.phst))],
             pch=21, bg=c('purple','green'), cex=1.5)
    }

    if(length(flightAdj)>1) {
      rect(flightAdj[-length(flightAdj)], rep(-1, length(flightAdj)-1),
           flightAdj[-1], rep(2, length(flightAdj)-1), col=rgb(0.5,0.5,0,0.3),
           border=NA)
    } else {
      abline(v=flightAdj, col=rgb(0.5,0.5,0, 1))
    }

    if(class(harpStage)!='logical') {
      text(harpStage$numDate, par('usr')[4], adj=c(0, -0.5), srt=45, col='slategrey', 'Staging 2018', cex=0.7, xpd=NA)
      text(REP$t_tot[minDiff], par('usr')[4], adj=c(0, -0.5), srt=45, 'Best fitting 2018 proportions', cex=0.7, xpd=NA)
      text(mean(DatePhoto), par('usr')[4], adj=c(0, -0.5), srt=45, col='slategrey', 'Photo surveys 2018', cex=0.7, xpd=NA)
      text(mean(flightAdj), par('usr')[4], adj=c(0, -0.5), srt=45, 'Corrected photo survey dates', cex=0.7, xpd=NA)
    }

    if(class(hoodStage)!='logical') {
      text(hoodStage$numDate, par('usr')[4], adj=c(0,-0.5), srt=45, col='slategrey', 'Staging 2018', cex=0.7, xpd=NA)
      text(REP$t_tot[minDiff]+3, par('usr')[4], adj=c(0,0), srt=45, 'Best fitting 2018 proportions & corrected photo survey dates
', cex=0.7, xpd=NA)
      text(mean(DatePhoto), par('usr')[4], adj=c(0,-0.5), srt=45, col='slategrey', 'Photo surveys 2018', cex=0.7, xpd=NA)
    }
    text(par('usr')[1], par('usr')[4], pos=3, '(A)', xpd=NA)
    if(class(harpStage)!='logical') {
      legend('right', lwd=2, col=c(2,4,3), cex=0.5, c('Newborn/Yellow', 'Thin', 'Fat/Grey'), bty='n')
    }

    if(class(hoodStage)!='logical') {
      legend('left', lwd=c(rep(2, 3), 1),
             col=c('red','blue','green', 'purple'),
             cex=0.5,
             c('Newborn/Thin', 'Fat', 'Solitary', 'Newborn/Thin + Fat'),
             bty='n')
    }

    matplot(REP[,1], REP[,c(5:7)],
            type='l', lty=1, lwd=2, ylim=c(0,1),
            col=c('red', 'blue', 'green'),
            xlab='Days since March 1st (2012)',
            ylab='Proportion seals on ice', axes=F)
    axis(1)
    axis(2)
    abline(v=par('usr')[1])
    abline(h=par('usr')[3])

    if(length(flightAdj)>1) {
      rect(flightAdj[-length(flightAdj)], rep(-1, length(flightAdj)-1),
           flightAdj[-1], rep(2, length(flightAdj)-1), col=rgb(0.5,0.5,0,0.3),
           border=NA)
    } else {
      abline(v=flightAdj, col=rgb(0.5,0.5,0, 1))
    }

    text(mean(flightAdj),  par('usr')[4], adj=c(0, -0.5), srt=45, 'Corrected photo survey dates', cex=0.7, xpd=NA)
    text(par('usr')[1], par('usr')[4], pos=3, '(B)', xpd=NA)
    if(class(harpStage)!='logical') {
      legend('topright', lwd=2, col=c(2,4,3), cex=0.5, c('Newborn/Yellow', 'Thin', 'Fat/Grey'), bty='n')
    }
    if(class(hoodStage)!='logical') {
      legend('topleft', lwd=2,
             col=c('red','blue','green'),
             cex=0.5,
             c('Newborn/Thin', 'Fat', 'Solitary'),
             bty='n')
    }
  }

  which.prop <- unlist(lapply(flightAdj, function(x) which.min(abs(x-REP$t_tot))))
  which.prop <- c(which.prop[1]:which.prop[2])

  if(class(harpStage)!='logical') {
    list(harpStage=harpStage, dayAdj=dayAdj, flightAdj=flightAdj, prop=REP$pr1.2.3[which.prop])
  } else {

    if(class(hoodStage)!='logical') {
      which.prop.st <- unlist(lapply(flightAdjSt, function(x) which.min(abs(x-REP$t_tot))))
      which.prop.st <- c(which.prop.st[1]:which.prop.st[2])
      if(return.phst) {
        list(hoodStage=hoodStage, phst=hood.phst, dayAdjSt=dayAdjSt, flightAdjSt=flightAdjSt, dayAdj=unique(flightAdj-DatePhoto), flightAdj=flightAdj,
             prop=REP$pr1.2.3[which.prop], propSt=REP$pr1.2.3[which.prop.st])
      } else {
        list(hoodStage=hoodStage, dayAdjSt=dayAdjSt, flightAdjSt=flightAdjSt, dayAdj=unique(flightAdj-DatePhoto), flightAdj=flightAdj,
             prop=REP$pr1.2.3[which.prop], propSt=REP$pr1.2.3[which.prop.st])
      }
    }
  }
}

