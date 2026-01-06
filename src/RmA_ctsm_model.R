## Source functions for fitting the RmA model with EKF

create_RmA_model_ctsmTMB <- function() {
    model <- ctsmTMB$new()

    # Add system dynamics
    model$addSystem(
        dX ~ (fN * hd_N + 0.5 * gN^2 * hdd_N) * dt + gN * hd_N * dw1,
        dY ~ (fP * hd_P + 0.5 * gP^2 * hdd_P) * dt + gP * hd_P * dw2
    )

    # Add observations
    model$addObs(
        YN_gauss ~ N
    )    

    model$setVariance(
        YN_gauss ~ obs_sd^2
    )

    model$setAlgebraics(
        N ~ exp(X),
        P ~ exp(Y),
        fN ~ r * N * (1 - N / K) - beta * N * P / (1 + beta * N / Cmax),
        fP ~ epsilon * beta * N * P / (1 + beta * N / Cmax) - mu * P,
        gN ~ sN * N,
        gP ~ sP * P,
        hd_N ~ 1 / N,
        hd_P ~ 1 / P,
        hdd_N ~ -1 / N^2,
        hdd_P ~ -1 / P^2
    )

    model$setParameter(
        r = c(init=1, lower=0, upper = 10),
        K = c(init=1, lower=0, upper = 10),
        #epsilon = c(init=3, lower=0, upper = 10),
        epsilon = 3.0,
        beta = c(init=3, lower=0, upper = 10),
        #Cmax = c(init=1, lower=0, upper = 10),
        Cmax = 1.0,
        mu = c(init=1, lower=0, upper = 10),
        sN = c(init=0.2, lower=0, upper = 10),
        sP = c(init=0.1, lower=0, upper = 10),
        obs_sd = 0.1
    )

    model$setInitialState(
        list(x0 = log(c(0.1, 0.1)), p0 = diag(c(0.1, 0.1)))
    )

    return(model)

}
