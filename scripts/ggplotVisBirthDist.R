library(ggplot2)

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
