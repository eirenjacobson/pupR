library(TMB)

compile("birthDist.cpp","-O1 -g",DLLFLAGS="")
dyn.load(dynlib("birthDist"))

load("birthidstData.RDat")

obj <- MakeADFun(data,parameters,DLL="birthDist",checkParameterOrder = FALSE)

obj$fn()
obj$gr()
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr,control = list(eval.max = 1e6,maxit = 1e6)))

rep<-sdreport(obj, getJointPrecision=TRUE)
rep.matrix <- summary(rep)
rep.rnames = rownames(rep.matrix)

#Extract estimates
indQ = which(rep.rnames == "PropEst")
indmub = which(rep.rnames == "mub")
indsigmab = which(rep.rnames == "sigmab")

Q = rep.matrix[indQ,1]
Qsd = rep.matrix[indQ,2]
mub = rep.matrix[indmub,1]
mubsd = rep.matrix[indmub,2]
sigmab = rep.matrix[indsigmab,1]
sigmabsd = rep.matrix[indsigmab,2]

xax = seq(mub-3*sigmab,mub+3*sigmab,by = 0.1)
bdist = dnorm(xax,mub,sigmab)
#bdistmin = dnorm(xax,mybmin,sigmab)
#bdistmax = dnorm(xax,mybmax,sigmab)

############################
# Brukes ikke
############################

windows(width = 9,height = 6)
plot(xax,bdist,type = "l",
     lwd = 4,
     col = "royalblue",
     xlab = "Dates in March",
     ylab = "Density",
     bty = "l")

mubmin = mub - 1.96*mubsd
mubmax = mub + 1.96*mubsd

sigmabmin = sigmab - 1.96*sigmabsd
sigmabmax = sigmab + 1.96*sigmabsd

bdistmin = dnorm(xax,mubmin,sigmab)
bdistmax = dnorm(xax,mubmax,sigmab)

bdistminsig = dnorm(xax,mub,sigmabmin)
bdistmaxsig = dnorm(xax,mub,sigmabmax)

bdistminnew = dnorm(xax,mubmin,sigmabmin)
bdistmaxnew = dnorm(xax,mubmax,sigmabmax)

windows(width = 9,height = 6)
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(xax,bdist,type = "n",
     ylim = c(0,0.5),
     lwd = 4,
     col = "royalblue",
     xlab = "Dates in March",
     ylab = "Density",
     bty = "l",
     main = "Estimated birth distribution",
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(x = c(xax,rev(xax)),c(bdistminnew,rev(bdistmaxnew)),border = NA,
        col = "lightblue")
polygon(x = c(xax,rev(xax)),c(bdistmin,rev(bdistmax)),
        border = NA,
        col = "lightblue")
polygon(x = c(xax,rev(xax)),c(bdistminsig,rev(bdistmaxsig)),
        border = NA,
        col = "lightblue")

lines(xax,bdist,col = "royalblue",lwd = 4)
lines(mub*rep(1,10),seq(0,0.5,length.out = 10),
      lwd = 4,lty = 2,col = "black")

##################
# Not used
lines(xax,bdistmin,lwd = 4,lty = 2,col = "royalblue")
lines(xax,bdistmax,lwd = 4,lty = 2,col = "royalblue")

# Not used
lines(xax,bdistminsig,lwd = 4,lty = 2,col = "red")
lines(xax,bdistmaxsig,lwd = 4,lty = 2,col = "red")

lines(xax,bdistminnew,lwd = 4,lty = 2,col = "green")
lines(xax,bdistmaxnew,lwd = 4,lty = 2,col = "green")
###################


# Monte Carlo simulations

Nsim = 10000
muv = rnorm(Nsim,mub,mubsd)
sigmav = rnorm(Nsim,sigmab,sigmabsd)

BirthDistCurves = matrix(0,nrow = Nsim,ncol = length(xax))

windows(width = 9,height = 6)
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(xax,bdist,type = "n",
     ylim = c(0,0.5),
     lwd = 4,
     col = "royalblue",
     xlab = "Dates in March",
     ylab = "Density",
     bty = "l",
     main = "Estimated birth distribution",
     cex.axis = 1.5,
     cex.lab = 1.5)

for(i in 1:Nsim){
        BirthDistCurves[i,] = dnorm(xax,muv[i],sigmav[i])
        #points(xax,BirthDistCurves[i,],lwd = 1,col = rgb(red = 0, green = 1, blue = 0, alpha = 0.1))

}

Bdistmin = rep(0,length(xax))
Bdistmax = Bdistmin
for(i in 1:length(xax)){
        Bdistmin[i] = quantile(BirthDistCurves[,i],0.025)
        Bdistmax[i] = quantile(BirthDistCurves[,i],0.975)
}


windows(width = 9,height = 6)
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(xax,bdist,type = "n",
     ylim = c(0,0.5),
     lwd = 4,
     col = "royalblue",
     xlab = "Dates in March",
     ylab = "Density",
     bty = "l",
     main = "Estimated birth distribution",
     cex.axis = 1.5,
     cex.lab = 1.5)
polygon(x = c(xax,rev(xax)),c(Bdistmin,rev(Bdistmax)),border = NA,
        col = "lightblue")

lines(xax,bdist,col = "royalblue",lwd = 4)
lines(mub*rep(1,10),seq(0,0.5,length.out = 10),
      lwd = 4,lty = 2,col = "black")


Report = obj$report()
nn1 = Report$Nout[,1]
nn2 = Report$Nout[,2]
nn3 = Report$Nout[,3]
nn = Report$Nout[,4]

BirthDistCurves = matrix(0,nrow = Nsim,ncol = length(xax))
nn1min = rep(0,length(nn1))
nn1max = nnmin
nn2min = nnmin
nn2max = nnmin
nn3min = nnmin
nn3max = nnmin
nnmin = nnmin
nnmax = nnmin

for(i in 1:length(xax)){
  Bdistmin[i] = quantile(BirthDistCurves[,i],0.025)
  Bdistmax[i] = quantile(BirthDistCurves[,i],0.975)
}



t_tot = data$ttot
days = data$days
staging = data$staging

windows(height = 7,width = 9)
par(mar=c(6,5,4,5),bg = "white")
plot(t_tot,nn1/nn,type = "l",col = "red",lwd = 4,xlim = c(min(days)-5,max(days)+5),ylim = c(0,1),xlab = "Days since 1. March 2012",ylab = "Proportion",cex.lab = 1.5,cex.main = 1.5,bty = "l")
lines(t_tot,nn2/nn,col = "blue",lwd = 4)
lines(t_tot,nn3/nn,col= "green",lwd = 4)
points(days,staging[,1]/rowSums(staging),bg = "red",pch = 21,cex = 1.5)
points(days,staging[,2]/rowSums(staging),bg = "blue",pch = 22,cex = 1.5)
points(days,staging[,3]/rowSums(staging),bg = "green",pch = 24,cex = 1.5)

legend('right', lwd=2, col=c("red","blue","green"), cex=1.0, c('Newborn/Yellow', 'Thin', 'Fat/Grey'), bty='n')

