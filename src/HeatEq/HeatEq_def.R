# Heat equation model def for simulation


# Model definitions for simulation
fsim <- function(x, par=p0, nonlinear = FALSE) {
    if (!nonlinear) {
        K <- exp(par$log_K)
    } else {
        K0 <- exp(par$log_K0)
        beta <- exp(par$log_beta)
        K <- K0 * exp(-beta * x)
    }
    L <- par$L
    Nx <- length(x)
    dx <- L / Nx
    dX <- numeric(length(x))
    # First FV
    dX[1] <- K[1] * (x[2] - 2*x[1] + par$T0)/(dx^2)
    # Middle FVs
    for (i in 2:(length(x)-1)) {
        dX[i] <- K[i] * (x[i+1] - 2*x[i] + x[i-1])/(dx^2)
    }
    # Last FV
    dX[length(x)] <- K[length(x)] * (par$TL - 2*x[length(x)] + x[length(x)-1])/(dx^2)
    return(dX)
}

gsim <- function(x, par=p0, state_dependent = FALSE) {
    sT <- exp(par$log_sT)

    if (!state_dependent) {
        g <- sT * diag(length(x))
    } else {
        #g <- sT * diag(x * (par$TL - x))
        g <- sT * diag(sqrt(x))
    }

    return(g)
}

hsim <- function(x, par=p0) {

    Nx <- ncol(x)
    sensor_pos_idx <- ceiling(par$sensor_pos/par$L * Nx)

    return(x[,sensor_pos_idx])
}
