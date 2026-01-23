## Generate data

generate_data_RmA <- function(fsim, gsim, hsim, t, x0, p0, dt, N, tsample=1) {

    # Gaussian observations
    n <- length(t)
    iobs <- round(seq(1,n,tsample/dt))

    # State and observation dimension length
    nx <- ncol(x0)
    ny <- length(hsim(x0[1, , drop = FALSE], p0))
    ## Simulate using Euler-Maruyama
    #Xsim <- matrix(0, nrow=length(t), ncol=N)
    #Ysim <- matrix(0, nrow=length(iobs), ncol=N)

    Xsim <- vector("list", N)
    Ysim <- vector("list", N)

    for (i in 1:N) {
        B <- rvBM(t, nx)

        sim <- euler(function(x)fsim(x,p0),
                    function(x)gsim(x,p0),
                    t,x0[i,],B=B)

        Xsim[[i]] = sim$X
        Ysim[[i]] = matrix(pmax(0,
                       hsim(Xsim[[i]][iobs,], p0) 
                       + rnorm(length(iobs),0,p0$obs_sd)), 
                       nrow=length(iobs), ncol=ny)
        # Poisson observations
        # Nsim[[i]] = exp(Xsim[[i]][,1])
        # vN <- 10
        # Ysim[[i]] = rpois(length(iobs), vN * Nsim[[i]][iobs])
  
                    
    
    }

    # Are any observations negative?
    is_Y_neg <- any(sapply(Ysim, function(y) any(y < 0)))
    if (is_Y_neg) {
        warning("Some generated observations are negative!")
    }


    return(list(Xsim=Xsim, Ysim=Ysim, iobs=iobs))
}
