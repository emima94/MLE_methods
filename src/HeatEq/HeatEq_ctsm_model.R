## Source functions for fitting the RmA model with EKF

create_HeatEq_model <- function(obs_sd, x0, P0, Nx, dx, T0, TL, sensor_pos) {
    model <- ctsmTMB$new()

    # Add system dynamics
    # First FV:
    model$addSystem(
        dT1 ~ K * (T2 - 2*T1 + T0)/(dx^2) * dt + sT * dw1
    )

    # Middle FVs:
    for (i in 2:(Nx-1)) {
        model$addSystem(
            as.formula(paste0("dT",i," ~ K * (T",i+1," - 2*T",i," + T",i-1,")/(dx^2) * dt + sT * dw",i))
        )
    }
    # Last FV:
    model$addSystem(
        as.formula(paste0("dT",Nx," ~ K * (TL - 2*T",Nx," + T",Nx-1,")/(dx^2) * dt + sT * dw",Nx))
    )

    # Add observations as sensor positions
    for (j in 1:length(sensor_pos)) {
       model$addObs(
           as.formula(paste0("Y",j," ~ T",sensor_pos[j]))
       )
    } 

    # Add observation variance
    for (j in 1:length(sensor_pos)) {
       model$setVariance(
           as.formula(paste0("Y",j," ~ obs_sd^2"))
       )
    }

    model$setAlgebraics(
        K ~ exp(log_K),
        sT ~ exp(log_sT)
    )

    model$setParameter(
        log_K = log(c(init=1.0, lower=0, upper = 10)),
        log_sT = log(c(init=1.0, lower=0, upper = 10)),
        dx = dx,       # spatial step size
        T0 = T0,        # boundary condition at left end
        TL = TL,       # boundary condition at right end
        obs_sd = obs_sd
    )

    model$setInitialState(
        list(x0 = x0, p0 = P0)
    )

    return(model)

}
