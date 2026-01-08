## Generate data

generate_data <- function(fsim, gsim, t, x0, p0, dt, N) {

    # Gaussian observations
    n <- length(t)
    tsample <- 1 # 5 
    iobs <- round(seq(1,n,tsample/dt))

    ## Simulate using Euler-Maruyama
    Xsim <- matrix(0, nrow=length(t), ncol=N)
    Ysim <- matrix(0, nrow=length(iobs), ncol=N)

    for (i in 1:N) {
        B <- rBM(t)

        sim <- euler(function(x)fsim(x,p0),
                    function(x)gsim(x,p0),
                    t,x0[i],B=B)
        Xsim[,i] = sim$X

        Ysim[,i] = Xsim[iobs,i] + rnorm(length(iobs),0,p0$obs_sd) 
    }

    return(list(Xsim=Xsim, Ysim=Ysim, iobs=iobs))
}
