stageCorrection <- function(path='C:/Users/a5406/Documents/Populasjonsmodellering/WestIce2012 Martin/BirthDist',
                            datafile='HarpStages2012_1.txt',
                            repfile='birthdist3.rep',
                            stdfile='birthdist3.std',
                            harpStage='C:/Users/a5406/Documents/Telletokt2018/HarpStaging.txt',
                            hoodStage='C:/Users/a5406/Documents/Telletokt2018/HoodStaging.txt',
                            photoStage=all.ind,
                            DatePhoto = 28) {
  
  ## 2012 data:
  data <- read.table(paste(path, datafile, sep='/'),sep = "",header = TRUE)
  days = data$Date; ndays = length(days)
  
  length_stage1 = 2.4
  length_stage2 = 4.42
  length_stage3 = 11.39
  
  kappa = 12.4

  
  rho1 = length_stage1 / kappa;
  rho2 = length_stage2 / kappa;
  rho3 = length_stage3 / kappa;
  
  
  staging =  as.matrix(data[1:ndays,2:8])
  
  #Combine newborn/yellow and fat/grey stages
  staging[,2] = staging[,1] + staging[,2]
  staging[,4] = staging[,4] + staging[,5]
  
  staging = staging[,c(-1,-5,-6,-7)] 
  
  xmax = 45
  xmin = 15 # West ice surveys

  DimStaging = dim(staging); nstages = DimStaging[2]
  stageprop = staging/rowSums(staging)
  
  spacing = 1/24; spacing = 0.05
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
  
  #Not used any more
  priors = matrix(NA,2,2)
  priors[1,] = c(21,5); priors[2,] = c(1,1)
  
  datoind = round(mean(which((t_tot>DatePhoto) & (t_tot<(DatePhoto+1)))))
  
  
  REP <- read.table(paste(path, repfile, sep='/'),header = F, sep = "")
  TMP <- read.table(paste(path, stdfile, sep='/'),col.names = c("Indx","Parameter","Estimate","SE"),skip = 1)
  names(REP) <- c('nn1', 'nn2', 'nn3', 'nn')
  REP$nn1.2 <- REP$nn1+REP$nn2
  REP$nn2.3 <- REP$nn2+REP$nn3
  REP$p1 <- REP$nn1/REP$nn
  REP$p2 <- REP$nn2/REP$nn
  REP$p3 <- REP$nn3/REP$nn
  REP$p1.2 <- REP$nn1.2/REP$nn
  REP$p2.3 <- REP$nn2.3/REP$nn
  
  REP$p1[which(is.nan(REP$p1))] <- NA
  REP$p2[which(is.nan(REP$p2))] <- NA
  REP$p3[which(is.nan(REP$p3))] <- NA
  REP$p1.2[which(is.nan(REP$p1.2))] <- NA
  REP$p2.3[which(is.nan(REP$p2.3))] <- NA
  
  harpStage <- read.csv(harpStage, dec='.')
  harpStage$Total <- apply(harpStage[,-1], 1, sum)
  harpStage$nn1 <- harpStage$New+harpStage$Yellow
  harpStage$nn2 <- harpStage$Thin
  harpStage$nn3 <- harpStage$Fat+harpStage$Grey
  harpStage$nn1.2 <- harpStage$nn1+harpStage$nn2
  
  harpStage$p1 <- harpStage$nn1/harpStage$Total
  harpStage$p2 <- harpStage$nn2/harpStage$Total
  harpStage$p3 <- harpStage$nn3/harpStage$Total
  harpStage$p1.2 <- harpStage$nn1.2/harpStage$Total
  
  ## minDiff <- which.min(abs(harpStage$p1-REP$p1)+abs(harpStage$p2-REP$p2)+abs(harpStage$p3-REP$p3))
  minDiff <- which.min(abs(harpStage$p1.2-REP$p1.2)+abs(harpStage$p3-REP$p3))
  
  ## Add staging from photo surveys in 2018:
  
  photoStage <- table(photoStage$Date, photoStage$STATUS, photoStage$SPECIES_ID)
  harp.phst <- photoStage[,,1]
  hood.phst <- photoStage[,,2]
  
  harp.phst <- as.data.frame.matrix(harp.phst)
  harp.phst <- cbind(Date=sort(unique(all.ind$Date)), harp.phst)
  harp.phst$Date <- as.POSIXct(strptime(row.names(harp.phst), '%Y-%m-%d'), tz='UTC')
  harp.phst <- harp.phst[,match(c('Date', 'Nb', 'TF'), names(harp.phst))]
  
  harp.phst <- rbind(harp.phst, data.frame(Date=mean(harp.phst$Date), 
                                           Nb=mean(harp.phst$Nb),
                                           TF=mean(harp.phst$TF)))
  harp.phst$nn1 <- harp.phst$Nb
  harp.phst$nn2.3 <- harp.phst$TF
  harp.phst$nn <- harp.phst$nn1+harp.phst$nn2.3
  harp.phst$p1 <- harp.phst$nn1/harp.phst$nn
  harp.phst$p2.3 <- harp.phst$nn2.3/harp.phst$nn
  
  minDiffPh <- which.min(abs(harp.phst$p1[3]-REP$p1)+abs(harp.phst$p2.3[3]-REP$p2.3))
  
  palette(RColorBrewer::brewer.pal(8, 'Set1'))
  ##hood.phst <- as.data.frame.matrix(hood.phst)
  ##hood.phst <- cbind(Date=sort(unique(all.ind$Date)), hood.phst)
  ##hood.phst$Date <- as.POSIXct(strptime(row.names(hood.phst), '%Y-%m-%d'), tz='UTC')
  ##hood.phst <- hood.phst[,match(c('Date', 'Nb', 'Mp', 'Fm', 'Bb'), names(hood.phst))]
  ##hood.phst$MpFm <- hood.phst$Mp+hood.phst$Fm
  
  
  
  
  plot(t_tot, REP$nn1, type='l', lwd=2, col=2, 
       xlim = c(xmin,xmax),ylim = c(0,1),
       xlab = "Days since 1. March 2012",
       ylab = "Proportion",cex.lab = 1.5,cex.main = 1.5)
  lines(t_tot, REP$nn1+REP$nn2, lwd=2, col=3)
  lines(t_tot, REP$nn1+REP$nn2+REP$nn3, lwd=2, col=4)
  legend(33,0.6,legend = c("Newborn/Yellow","Thin","Fat/Grey"),lwd = rep(2, 3),col = c(2:4),bty = "n",cex = 1.2)
  
}

## This is with correct curves (as 'included in'stolen' from 2016 report):

stageCorrection <- function(harpStage='C:/Users/a5406/Documents/Telletokt2018/HarpStaging.txt',
                            hoodStage='C:/Users/a5406/Documents/Telletokt2018/HoodStaging.txt',
                            repfolder='C:/Users/a5406/Documents/Populasjonsmodellering/WestIce2012 Martin/BirthDist',
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