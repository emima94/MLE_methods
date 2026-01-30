## Source functions for fitting the RmA model with EKF

create_RmA_model_general <- function(p, x0, P0, log_transform_par = FALSE, log_transform_state = FALSE, obs_partial = TRUE) {
    model <- ctsmTMB$new()

    # Add system dynamics
    if (log_transform_state) {
        model$addSystem(
            dX1 ~ (fN * hd_N + 0.5 * gN^2 * hdd_N) * dt + gN * hd_N * dw1,
            dX2 ~ (fP * hd_P + 0.5 * gP^2 * hdd_P) * dt + gP * hd_P * dw2
        )
        model$setAlgebraics(
            N ~ exp(X1),
            P ~ exp(X2),
            hd_N ~ 1 / N,
            hd_P ~ 1 / P,
            hdd_N ~ -1 / N^2,
            hdd_P ~ -1 / P^2
        )
    } else {
        model$addSystem(
            dN ~ fN * dt + gN * dw1,
            dP ~ fP * dt + gP * dw2
        )
    }
    
    # Add observations
    if (obs_partial) {
        model$addObs(
            Y ~ N
        )    
        model$setVariance(
            Y ~ obs_sd^2
        )
    } else {
        model$addObs(
            Y1 ~ N,
            Y2 ~ P
        )   
        model$setVariance(
            Y1 ~ obs_sd^2,
            Y2 ~ obs_sd^2
        ) 
    }
    
    model$setAlgebraics(
        fN ~ r * N * (1 - N / K) - beta * N * P / (1 + beta * N / Cmax),
        fP ~ epsilon * beta * N * P / (1 + beta * N / Cmax) - mu * P,
        gN ~ sN * N,
        gP ~ sP * P
    )

    if (log_transform_par) {
        model$setAlgebraics(
            r ~ exp(log_r),
            K ~ exp(log_K),
            beta ~ exp(log_beta),
            mu ~ exp(log_mu),
            sN ~ exp(log_sN)
        )
        model$setParameter(
            log_r = c(init = log(0.5), lower=log(0.01), upper = log(10)),
            log_K = c(init = log(0.5), lower=log(0.01), upper = log(10)),
            #epsilon = c(init=3, lower=0, upper = 10),
            epsilon = p$epsilon,
            log_beta = c(init=log(1.0), lower=log(0.01), upper = log(10)),
            Cmax = p$Cmax,
            log_mu = c(init=log(0.5), lower=log(0.01), upper = log(10)),
            log_sN = c(init=log(0.4), lower=log(0.01), upper = log(10)),
            #log_sP = c(init=log(0.1), lower=log(0.01), upper = log(10)),
            sP = p$sP,
            obs_sd = p$obs_sd
        )
    } else {

    model$setParameter(
        r = c(init = 0.5, lower=0.01, upper = 10),
        K = c(init = 0.5, lower=0.01, upper = 10),
        #epsilon = c(init=3, lower=0, upper = 10),
        epsilon = p$epsilon,
        beta = c(init=1.0, lower=0.01, upper = 10),
        Cmax = p$Cmax,
        mu = c(init=0.5, lower=0.01, upper = 10),
        sN = c(init=0.4, lower=0.01, upper = 10),
        #sP = c(init=0.1, lower=0, upper = 10),
        sP = p$sP,
        obs_sd = p$obs_sd
    )
    }

    model$setInitialState(
        list(x0 = x0, p0 = P0)
    )

    return(model)

}
