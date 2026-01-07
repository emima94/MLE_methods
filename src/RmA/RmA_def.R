
## RmA/Bazykin-model (2D model) in natural coordinates
fN <- function(n,p,par=p0) with(par,r*n*(1-n/K) - beta*n*p/(1+beta*n/Cmax))
fP <- function(n,p,par=p0) with(par,epsilon*beta*n*p/(1+beta*n/Cmax) - mu*p)


gN <- function(n,par=p0) with(par,sN*n) 
gP <- function(p,par=p0) with(par,sP*p) 

## In transformed coordinates X=h(N), Y=h(P)
h <- log
hi <- exp
hd <- function(n) 1/n
hdd <- function(n) -1/n^2

fX <- function(x,y,par=p0)
    fN(hi(x),hi(y),par) * hd(hi(x)) + 0.5*gN(hi(x),par)^2*hdd(hi(x))

fY <- function(x,y,par=p0)
    fP(hi(x),hi(y),par) * hd(hi(y)) + 0.5*gP(hi(y),par)^2*hdd(hi(y))

gX <- function(x,par=p0) with(par,gN(hi(x),par)*hd(hi(x)))
gY <- function(y,par=p0) with(par,gP(hi(y),par)*hd(hi(y)))

## The model for simulation: dX = f(X)*dt + g(X)*dB, X=(log N,logP)
fsim <- function(x,par=p0) c(fX(x[1],x[2],par),fY(x[1],x[2],par))
gsim <- function(x,par=p0) diag(c(gX(x[1],par),gY(x[2],par)))
