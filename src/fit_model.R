## Fit model to data

fit_model <- function(methods, model, t, Ysim, iobs, dt, N) {
    fit <- vector("list", length(methods))
    names(fit) <- methods

    for (m in methods) {
        message("Fitting method: ", m, sep="")
        
        #fit_sub <- future_lapply(1:N, function(i) {
        fit_sub <- lapply(1:N, function(i) {
            
            message(" Fitting dataset ", i, " of ", N, " with ", m)

            df <- data.frame(
                t = t[iobs],
                Y = Ysim[,i]
            )

            model$estimate(
                    data = df,
                    method = m,
                    ode.solver = "rk4",
                    ode.timestep = dt,
                    silent = TRUE,
                    control = list(trace = 0)
            )
            

        })

        fit[[m]] <- fit_sub
    }

    return(fit)
}