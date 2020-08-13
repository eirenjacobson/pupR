
staging = matrix(NA,4,3)
staging[1,1:3]=c(72,2,2)
staging[2,1:3]=c(314,25,83)
staging[3,1:3]=c(107,57,160)
staging[4,1:3]=c(32,26,150)

days = c(22, 24, 27, 29)



phat <- read.table("HoodedBirthDistOrigParam.csv",header = FALSE, sep = ",")
padj <- read.table("HoodedBirthDistAdjParam.csv",header = FALSE, sep = ",")



X11("",9,7)
par(mar=c(6,5,4,5),bg = "white")
plot(phat[,1],phat[,2]/phat[,5],type = "l",col = "darkred",lwd = 4,xlim = c(10,40),ylim = c(0,1),xlab = "Days in March 2012",ylab = "Proportion",cex.lab = 1.5,cex.main = 1.5,bty = "l")

lines(phat[,1],phat[,3]/phat[,5],col = "darkblue",lwd = 4)
lines(phat[,1],phat[,4]/phat[,5],col = "darkgreen",lwd = 4)

lines(days,staging[,1]/rowSums(staging),type = "p",col = "darkred",pch = 15,lwd = 5)
#lines(days,staging[,2]/rowSums(staging),type = "p",col = "darkgreen",pch = 2,lwd = 5)
lines(days,staging[,2]/rowSums(staging),type = "p",col = "darkblue",pch = 19,lwd = 5)
lines(days,staging[,3]/rowSums(staging),type = "p",col = "darkgreen",pch = 17,lwd = 5)

#To versjoner en med og uten disse.
lines(padj[,1],padj[,2]/padj[,5],col = "darkred",lwd = 2,lty = 2)
lines(padj[,1],padj[,3]/padj[,5],col = "darkblue",lwd = 2,lty = 2)
lines(padj[,1],padj[,4]/padj[,5],col = "darkgreen",lwd = 2,lty = 2)

legend(33,0.6,legend = c("Newborn/Thin","Fat","Solitary"),lwd = c(4,4,4),col = c("darkred","darkblue","darkgreen"),bty = "n")

X11("",9,7)
par(mar=c(6,5,4,5),bg = "white")
ind = which((phat[,1]>28) & (phat[,1]<29))

plot(phat[,1],phat[,2],type = "n",col = "darkred",lwd = 4,xlim = c(10,40),ylim = c(0,1),xlab = "Days in March 2012",ylab = "Proportion on ice",cex.lab = 1.5,cex.main = 1.5,bty = "l")
polygon(c(phat[min(ind),1],phat[max(ind),1],phat[max(ind),1],phat[min(ind),1]),c(0,0,1,1),col = "wheat",border = NA)
lines(phat[,1],phat[,2],col = "darkred",lwd = 4)
lines(phat[,1],phat[,2]+phat[,3],col = "darkblue",lwd = 4)
lines(phat[,1],phat[,2]+phat[,3]+phat[,4],col = "darkgreen",lwd = 4)
lines(phat[mean(ind),1]*array(1,10),seq(0,1,length.out = 10),col = "black",lwd = 1,lty = 3)

#To versjoner en med og uten disse.
lines(padj[,1],padj[,2],col = "darkred",lwd = 2,lty = 2)
lines(padj[,1],padj[,2]+padj[,3],col = "darkblue",lwd = 2,lty = 2)
lines(padj[,1],padj[,2]+padj[,3]+padj[,4],col = "darkgreen",lwd = 2,lty = 2)

legend("topleft",legend = c("Newborn/Thin","Fat","Solitary"),lwd = c(4,4,4),col = c("darkred","darkblue","darkgreen"),bty = "n")



#figure(2)
#hold on
#plot(t_tot,nn1,'r',t_tot,nn1+nn2,'g',t_tot,nn1+nn2+nn3,'b')
#axis([1 40 0 1.5])
#xlabel('date (March 2007)');
#ylabel('Probability');
#legend('Stage 1','Stage 1 + 2','Stage 1 + 2 + 3')
#grid on
