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

Report = obj$report()
nn1 = Report$Nout[,1]
nn2 = Report$Nout[,2]
nn3 = Report$Nout[,3]
nn = Report$Nout[,4]


colvcol = c("darkred","darkblue","darkgreen")
#colvgr = c("grey60","grey40","grey80","black")
colvgr = c("grey80","grey60","grey40","black")
#colvgr = colvcol
palette(RColorBrewer::brewer.pal(8, 'Set1'))

t_tot = data$ttot
days = data$days
staging = data$staging

windows(height = 7,width = 9)
par(mar=c(6,5,4,5),bg = "white")
plot(t_tot,nn1/nn,type = "l",col = colvcol[1],lwd = 4,xlim = c(xmin,xmax),ylim = c(0,1),xlab = "Days since 1. March 2012",ylab = "Proportion",cex.lab = 1.5,cex.main = 1.5,bty = "l")
lines(t_tot,nn2/nn,col = colvcol[2],lwd = 4)
lines(t_tot,nn3/nn,col= colvcol[3],lwd = 4)
lines(days,staging[,1]/rowSums(staging),type = "p",col = colvcol[1],pch = 15,lwd = 5)
lines(days,staging[,2]/rowSums(staging),type = "p",col = colvcol[2],pch = 19,lwd = 5)
lines(days,staging[,3]/rowSums(staging),type = "p",col = colvcol[3],pch = 17,lwd = 5)
