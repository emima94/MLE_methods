# Generate simulated data for the RmA 2D model

RmA_sim <- function(p0, dt, T, obs_sd) {
    
   
    t <- seq(0,T,dt)

    x0 <- log(c(0.1,0.1))

    ## Simulate using Euler-Maruyama
    BN <- rBM(t)
    BP <- rBM(t)
    sim <- euler(function(x)fsim(x,p0),
                function(x)gsim(x,p0),
                t,x0,B=cbind(BN,BP))

    Nsim <- hi(sim$X[,1])
    Psim <- hi(sim$X[,2])

    ## Generate Poisson observations at these time indeces
    n <- length(t)
    tsample <- 1 # 5 
    iobs <- round(seq(1,n,tsample/dt))
    vP <- 1e-4
    vN <- 10 # 8 # 10 
    YN <- rpois(length(iobs),vN*Nsim[iobs])
    YP <- rpois(length(iobs),vP*Psim[iobs])

    ## Generate Gaussian observations
    YN_gauss <- Nsim[iobs] + rnorm(length(iobs),0,obs_sd)
    YP_gauss <- Psim[iobs] + rnorm(length(iobs),0,obs_sd)

    return(list(t=t, Xsim = sim$X[,1], Ysim = sim$X[,2], 
                Nsim=Nsim, Psim=Psim,
                BN=BN, BP=BP,
                iobs=iobs, YN=YN, YP=YP,
                YN_gauss=YN_gauss, YP_gauss=YP_gauss))

}