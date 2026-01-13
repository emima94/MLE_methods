# Heat equation model def for simulation


# Model definitions for simulation
fsim <- function(x, par=p0) {
    K <- exp(par$log_K)
    dx <- par$dx
    dX <- numeric(length(x))
    # First FV
    dX[1] <- K * (x[2] - 2*x[1] + par$T0)/(dx^2)
    # Middle FVs
    for (i in 2:(length(x)-1)) {
        dX[i] <- K * (x[i+1] - 2*x[i] + x[i-1])/(dx^2)
    }
    # Last FV
    dX[length(x)] <- K * (par$TL - 2*x[length(x)] + x[length(x)-1])/(dx^2)
    return(dX)
}

gsim <- function(x, par=p0) {
    sT <- exp(par$log_sT)
    return(diag(rep(sT, length(x))))
}

hsim <- function(x, par=p0) {
    return(x[,par$sensor_pos])
}
