## Fit model with Laplace approximation and "tiny" approach


## Core function: The joint negative log density of all random variables
nloglik_laplace_tiny <- function(p, YN_gauss, iobs, dt, obs_sd)
{
    tiny <- 1e-3
    
    dX <- diff(p$X)
    dY <- diff(p$Y)

    Xi <- head(p$X,-1)
    Yi <- head(p$Y,-1)

    dXpred <- fX(Xi,Yi,p)*dt + gX(Xi,p)*p$dBN
    dYpred <- fY(Xi,Yi,p)*dt + gY(Yi,p)*p$dBP

    ## Use the Euler approximation for the SDE and the "tiny" approach
    nll <- - sum(dnorm(dX,
                       dXpred,
                       tiny,
                       log=TRUE))

    nll <- nll - sum(dnorm(dY,
                       dYpred,
                       tiny,
                       log=TRUE))
    #print(nll)

    ## The contribution from the Brownian motion 
    nll <- nll - sum(dnorm(p$dBN,0,sqrt(dt),log=TRUE))
    #print(nll)
    nll <- nll - sum(dnorm(p$dBP,0,sqrt(dt),log=TRUE))
    #print(nll)
    
    ## Contribution from observations
    #nll <- nll - sum(dpois(YN,hi(p$X[iobs])*vN,log=TRUE))
    obs_sd <- 0.1
    nll <- nll - sum(dnorm(YN_gauss,
                        mean = hi(p$X[iobs]),
                        sd = obs_sd,
                        log=TRUE))
    #print(nll)
    #nll <- nll - sum(dpois(YP,hi(p$Y[iobs])*vP,log=TRUE))
    #nll <- nll - sum(dnorm(YP_gauss,
    #                    mean = hi(p$Y[iobs]),
    #                    sd = 0.02,
    #                    log=TRUE))
    #print(nll)

    return(nll)
}

fit_laplace_tiny_RmA <- function(p0, Xsim, Ysim, BN, BP, YN_gauss, iobs, dt, obs_sd) {

    ## Initial guess on all parameters
    pinit <- c(p0,
        list(X=0*Xsim,
            dBN=0*diff(BN),
            Y=0*Ysim,
            dBP=0*diff(BP)))

    map <- list(epsilon=factor(NA),Cmax=factor(NA))

    nll <- function(p) {
        nloglik_laplace_tiny(p,
                             YN_gauss=YN_gauss,
                             iobs=iobs,
                             dt=dt,
                             obs_sd=obs_sd)
    }

    obj <- MakeADFun(nll,pinit,random=c("X","Y","dBN","dBP"),
                 map=map)        

    ## Fit parameters
    comp.time <- system.time(fit <- nlminb(obj$par, obj$fn, obj$gr))

    ## Obtain random effects, variances, etc.
    rep <- sdreport(obj)

    ## Extract states: Estimates and marginal posterior variances, confidence limits
    match.name <- function(x,name) x[names(x)==name]
    find.name <- function(x,name) names(x)==name

    Xhat <- match.name(rep$par.random,"X")
    Yhat <- match.name(rep$par.random,"Y")

    pred <- list(Xhat=Xhat,
                 Yhat=Yhat)

    return(list(obj=obj, rep=rep, comp.time=comp.time, pred=pred))
    
}