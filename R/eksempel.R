library(TMB)

compile("eksempel.CPP")
dyn.load(dynlib("eksempel"))

data = list()
data$x = c(1,2,3,4)
data$h = c(4,5,6,7)

