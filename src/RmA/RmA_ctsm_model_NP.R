## Source functions for fitting the RmA model with EKF

create_RmA_NP_model <- function(p, x0, P0) {
    model <- ctsmTMB$new()

    # Add system dynamics
    model$addSystem(
        dN ~ fN * dt + gN * dw1,
        dP ~ fP * dt + gP * dw2
    )

    # Add observations
    model$addObs(
        Y ~ N
    )    

    model$setVariance(
        Y ~ obs_sd^2
    )

    model$setAlgebraics(
        fN ~ r * N * (1 - N / K) - beta * N * P / (1 + beta * N / Cmax),
        fP ~ epsilon * beta * N * P / (1 + beta * N / Cmax) - mu * P,
        gN ~ sN * N,
        gP ~ sP * P
    )

    model$setParameter(
        r = c(init = 0.5, lower=0, upper = 10),
        K = c(init = 0.5, lower=0, upper = 10),
        #epsilon = c(init=3, lower=0, upper = 10),
        epsilon = p0$epsilon,
        beta = c(init=1.0, lower=0, upper = 10),
        Cmax = c(init=1, lower=0, upper = 10),
        mu = c(init=0.5, lower=0, upper = 10),
        sN = c(init=0.4, lower=0, upper = 10),
        #sP = c(init=0.1, lower=0, upper = 10),
        sP = p$log_sP,
        obs_sd = p$obs_sd
    )

    model$setInitialState(
        list(x0 = x0, p0 = P0)
    )

    return(model)

}
