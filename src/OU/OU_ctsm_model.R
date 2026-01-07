## O-U 1D linear process model definition for ctsmTMB

create_OU_model <- function(obs_sd, x0, p0) {
    model <- ctsmTMB$new()

    # Add system dynamics
    model$addSystem(
        dX ~ (theta * (mu - X)) * dt + sigma * dw1
    )

    # Add observations
    model$addObs(
        Y ~ X
    )
    model$setVariance(
        Y ~ obs_sd^2
    )

    model$setParameter(
        theta = c(init=1.0, lower=0, upper = 10),
        mu = c(init=1.0, lower=-10, upper = 10),
        sigma = c(init=1.0, lower=0, upper = 10),
        obs_sd = obs_sd
    )

    model$setInitialState(
        list(x0 = x0, p0 = p0)
    )

    return(model)

}