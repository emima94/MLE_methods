## Fit model to data

fit_model <- function(methods, model, t, Ysim, iobs, dt, N, df_fun = NULL) {
  
  fit <- vector("list", length(methods))
  names(fit) <- methods
  
  for (m in methods) {
    message("Fitting method: ", m)
    
    fit_sub <- lapply(1:N, function(i) {
      
      message(" Fitting dataset ", i, " of ", N, " with ", m)
      
      if (is.null(df_fun)) {
        # Default data frame
        df <- data.frame(
          t = t[iobs],
          Y = Ysim[[i]]
        )
      } else {
        # Custom data frame
        df <- df_fun(t, Ysim[[i]], iobs)
      }
      
      timing <- system.time({
        res <- model$estimate(
          data = df,
          method = m,
          ode.solver = "rk4",
          ode.timestep = dt,
          silent = TRUE,
          control = list(trace = 0)
        )
      })
      
      list(
        fit = res,
        time = timing["elapsed"]
      )
    })
    
    fit[[m]] <- fit_sub
  }
  
  return(fit)
}
