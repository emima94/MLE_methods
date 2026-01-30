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

      # Create method subfolder for results
      out_dir_method_light <- file.path(out_dir$light, m)
      if (!dir.exists(out_dir_method_light)) {
          dir.create(out_dir_method_light, recursive = TRUE)
      }
      out_dir_method_heavy <- file.path(out_dir$heavy, m)
      if (!dir.exists(out_dir_method_heavy)) {
          dir.create(out_dir_method_heavy, recursive = TRUE)
      }
      out_dir_method <- list(
        light = out_dir_method_light,
        heavy = out_dir_method_heavy
      )
  
      cl <- makeCluster(n_workers)
      
      clusterExport(
        cl,
        varlist = c("model", "m", "t", "Ysim", "iobs", "dt", "N", "df_fun", "out_dir_method"),
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
        out_light <- list(
          par.fixed = res$par.fixed,
          nll = res$nll,
          sd.fixed = res$sd.fixed,
          converged = res$converged,
          nll.gradient = res$nll.gradient,
          opt = res$private$opt,
          time.elapsed = res$private$timer_estimation,
          time.compile = res$private$timer_construct_adfun,
          time = timing["elapsed"]
        )

        saveRDS(out, file.path(out_dir_method$heavy, sprintf("%04d.rds", i)))
        saveRDS(out_light, file.path(out_dir_method$light, sprintf("%04d.rds", i)))
        return(NULL)
      })
      stopCluster(cl)
  }
}

