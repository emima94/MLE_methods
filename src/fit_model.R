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

# New Fit model in parallel
fit_model_par <- function(n_workers, out_dir, methods, model, t, Ysim, iobs, dt, N, df_fun = NULL) {

  for (m in methods) {
      message("Fitting method: ", m)
  
      cl <- makeCluster(n_workers)
      
      clusterExport(
        cl,
        varlist = c("model", "m", "t", "Ysim", "iobs", "dt", "N", "df_fun", "out_dir"),
        envir = environment()
      )

      parLapply(cl, seq_len(N), function(i) {
        
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

        #res$nll_funcs <- model$getLikelihood()
        
        out <- list(
          fit = res,
          model = model,
          time = timing["elapsed"]
        )
        saveRDS(out, file.path(out_dir, sprintf("%s_%04d.rds", m, i)))
        return(NULL)
      })
      stopCluster(cl)
  }
}

