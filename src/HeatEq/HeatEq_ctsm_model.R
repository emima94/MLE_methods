## Source functions for fitting the RmA model with EKF

create_HeatEq_model <- function(obs_sd, x0, P0, Nx, L, T0, TL, 
    sensor_pos_idx, 
    nonlinear = FALSE,
    state_dep_diffusion = FALSE
    ) {

    model <- ctsmTMB$new()

    # Add system dynamics
    # First FV:
    model$addSystem(
        dT1 ~ K1 * (T2 - 2*T1 + T0)/(dx^2) * dt + g1 * dw1
    )

    # Middle FVs:
    for (i in 2:(Nx-1)) {
        model$addSystem(
            as.formula(paste0("dT",i," ~ K",i," * (T",i+1," - 2*T",i," + T",i-1,")/(dx^2) * dt + g",i," * dw",i))
        )
    }
    # Last FV:
    model$addSystem(
        as.formula(paste0("dT",Nx," ~ K",Nx," * (TL - 2*T",Nx," + T",Nx-1,")/(dx^2) * dt + g",Nx," * dw",Nx))
    )

    # Add observations as sensor positions
    for (j in 1:length(sensor_pos_idx)) {
       model$addObs(
           as.formula(paste0("Y",j," ~ T",sensor_pos_idx[j]))
       )
    } 

    # Add observation variance
    for (j in 1:length(sensor_pos_idx)) {
       model$setVariance(
           as.formula(paste0("Y",j," ~ obs_sd^2"))
       )
    }

    # Set algebraic equations
    model$setAlgebraics(
        K0 ~ exp(log_K0),
        beta ~ exp(log_beta),
        sT ~ exp(log_sT),
        dx ~ L / Nx
    )

    # Set algebraic equations for K depending on nonlinear flag
    if (!nonlinear) {
        # Set algebraic equations for constant K
        for (i in 1:Nx) {
            model$setAlgebraics(
                as.formula(paste0("K",i," ~ exp(log_K)"))
            )
        }
    } else {
        # Set algebraic equations for temperature-dependent K
        for (i in 1:Nx) {
            model$setAlgebraics(
                as.formula(paste0("K",i," ~ K0 * exp(-beta * T",i,")"))
            )
        }
    }

    # Set algebraic equations for diffusion term g depending on state_dep_diffusion flag
    if (!state_dep_diffusion) {
        for (i in 1:Nx) {
            model$setAlgebraics(
                as.formula(paste0("g",i," ~ sT"))
            )
        }
    } else {
        for (i in 1:Nx) {
            model$setAlgebraics(
                as.formula(paste0("g",i," ~ sT * sqrt(T",i,")"))
            )
        }
    }

    # Set parameters
    if (!nonlinear) {
        model$setParameter(
            log_K = log(c(init=1.0, lower=0, upper = 10))
        )
    } else {
        model$setParameter(
            log_K0 = log(c(init=1.0, lower=0, upper = 10)),
            log_beta = log(c(init=0.1, lower=0, upper = 10))
        )
    }

 
    model$setParameter(
        log_sT = log(c(init=1.0, lower=0, upper = 10)),
        L = L,         # Length of the rod
        Nx = Nx,       # Number of spatial finite volumes
        T0 = T0,        # boundary condition at left end
        TL = TL,       # boundary condition at right end
        obs_sd = obs_sd
    )

    model$setInitialState(
        list(x0 = x0, p0 = P0)
    )

    return(model)

}
