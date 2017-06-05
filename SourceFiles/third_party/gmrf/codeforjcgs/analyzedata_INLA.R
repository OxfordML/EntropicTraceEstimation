setwd("S:/Documents/research/carlik/code/codeforjcgs/")

#install.packages("sp")
#install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/stable")
library(sp)
library(INLA)


# red sea data
data = read.table('readseadatavec.txt',header=T)
coords = as.matrix(data[,1:2])
redsea.dom = cbind(c(20,100,310,310,-10,-10), c(-20,-20,300,480,220,0))

t1 = proc.time()[3]
#mesh2 = inla.mesh.2d(coords, max.edge = 1)
mesh2 = inla.mesh.2d(, redsea.dom, max.edge=6)
#plot(mesh2,asp=1,xlim=c(100,120),ylim=c(100,120))
#plot(mesh2,asp=1)
spde2 = inla.spde2.matern(mesh2,alpha=2)
A2 = inla.spde.make.A(mesh2, loc=coords)
stk2 = inla.stack(data=list(resp=data$aot), A=list(A2,1),
                  effects=list(i=1:spde2$n.spde,
                               m=rep(1,nrow(data))), tag='est')
dim(inla.stack.A(stk2))
res2 <- inla(resp ~ 0 + m + f(i, model=spde2),
             data=inla.stack.data(stk2),
             control.predictor=list(A=inla.stack.A(stk2)),verbose=TRUE)
t2 = proc.time()[3]
(t2-t1)/60

summary(res2)
res2.field <- inla.spde2.result(res2, 'i', spde2, do.transf=TRUE)
muhat = res2$summary.fixed[[6]] # muhat
tausqhat = inla.mmarginal(res2.field$marginals.tau[[1]] )^2 # tau
kappasqhat = inla.mmarginal(res2.field$marginals.kappa[[1]] )^2 # kappa
nughat = 1/res2$summary.hy[1,6]*tausqhat # 1/\sigma_\varepsilon^2
muhat
tausqhat
kappasqhat
nughat









