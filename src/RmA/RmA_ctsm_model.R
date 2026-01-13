## Source functions for fitting the RmA model with EKF

create_RmA_model <- function(p, x0, P0) {
    model <- ctsmTMB$new()

    # Add system dynamics
    model$addSystem(
        dX1 ~ (fN * hd_N + 0.5 * gN^2 * hdd_N) * dt + gN * hd_N * dw1,
        dX2 ~ (fP * hd_P + 0.5 * gP^2 * hdd_P) * dt + gP * hd_P * dw2
    )

    # Add observations
    model$addObs(
        Y ~ N
    )    

    model$setVariance(
        Y ~ obs_sd^2
    )

    model$setAlgebraics(
        N ~ exp(X1),
        P ~ exp(X2),
        r ~ exp(log_r),
        K ~ exp(log_K),
        beta ~ exp(log_beta),
        mu ~ exp(log_mu),
        sN ~ exp(log_sN),
        sP ~ exp(log_sP),
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
        log_r = log(c(init=1.0, lower=0, upper = 10)),
        log_K = log(c(init=1, lower=0, upper = 10)),
        #epsilon = c(init=3, lower=0, upper = 10),
        epsilon = 3.0,
        log_beta = log(c(init=3.0, lower=0, upper = 10)),
        #Cmax = c(init=1, lower=0, upper = 10),
        Cmax = 1.0,
        log_mu = log(c(init=1.0, lower=0, upper = 10)),
        log_sN = log(c(init=0.2, lower=0, upper = 10)),
        #log_sP = log(c(init=0.1, lower=0, upper = 10)),
        log_sP = p$log_sP,
        obs_sd = p$obs_sd
    )

    model$setInitialState(
        list(x0 = x0, p0 = P0)
    )

    return(model)

}
