
likelihood_profile_2D_RmA <- function(p0, T, pg_name, n_grid, delta_list) {
    delta_plus <- delta_list$delta_plus_list
    delta_minus <- delta_list$delta_minus_list
    N <- 1
    x0_bar = c(0.5, 0.5)
    P0 <- diag(c(0.01,0.01))
    x0 <- log(matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar)))
    #T <- 200
    dt <- 0.1
    t <- seq(0,T,dt)
    tsample <- 0.2
    sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0[1,,drop=FALSE], p0, dt, N, tsample)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    df <- data.frame(t = t[iobs], Y = Ysim[[1]])

    model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)   
    pg <- c(paste0("log_", pg_name[1]), paste0("log_", pg_name[2]))

    var1 <- seq(exp(p0[[pg[1]]]) - delta_minus[[pg[1]]], exp(p0[[pg[1]]]) + delta_plus[[pg[1]]], length.out=n_grid)
    var2 <- seq(exp(p0[[pg[2]]]) - delta_minus[[pg[2]]], exp(p0[[pg[2]]]) + delta_plus[[pg[2]]], length.out=n_grid)

    # Simulate edge cases
    #pdf(paste0("figures/RmA_likelihood_profile_2D_", pg_name[1], "_", pg_name[2], "_edge_cases.pdf"), width = 8, height = 6)
    par(mfrow = c(2,2))
    for (edge_1 in c(min(var1), max(var1)) ) {
        for (edge_2 in c(min(var2), max(var2))) {
            p0_edge <- p0
            p0_edge[[pg[1]]] <- log(edge_1)
            p0_edge[[pg[2]]] <- log(edge_2)
            sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0[1,,drop=FALSE], p0_edge, dt, 1, tsample)
            Xsim <- sim_data$Xsim
            Ysim <- sim_data$Ysim
            iobs <- sim_data$iobs
            plot(t, exp(Xsim[[1]][,1]), type="l", main=sprintf("%s=%.2f, %s=%.2f", pg_name[1], edge_1, pg_name[2], edge_2))
            points(t[iobs], Ysim[[1]], col="red")
        }
    }
    plt_edge <- recordPlot()
    #dev.off()

    model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)

    if (!("log_r" %in% pg)) {
        model$setParameter(
            log_r = p0$log_r
        )
    }
    if (!("log_mu" %in% pg)) {
        model$setParameter(
            log_mu = p0$log_mu
        )
    }
    if (!("log_beta" %in% pg)) {
        model$setParameter(
            log_beta = p0$log_beta
        )
    }
    if (!("log_K" %in% pg)) {
        model$setParameter(
            log_K = p0$log_K
        )
    }
    if (!("log_sN" %in% pg)) {
        model$setParameter(
            log_sN = p0$log_sN
        )
    }
    model$getParameters()
    nll_func <- model$likelihood(df, method = "laplace",
                    ode.timestep = dt,
                    ode.solver = "rk4",
                    silent = FALSE)


    # Compute likelihood surface
    nll_matrix <- matrix(0, nrow=length(var1), ncol=length(var2))
    for (i in 1:length(var1)) {
        for (j in 1:length(var2)) {
            
            params <- c(
                log(var1[i]),
                log(var2[j])
            )
            nll_value <- nll_func$fn(params)[1]
            nll_matrix[i,j] <- nll_value
        
            #message(sprintf("r=%.2f, beta=%.2f, NLL=%.2f", var1[i], var2[j], nll_value))
        }
    }
    # Plot likelihood surface
    i <- which(nll_matrix == min(nll_matrix), arr.ind = TRUE)
    #pdf(paste0("figures/RmA_likelihood_profile_2D_", pg_name[1], "_", pg_name[2], ".pdf"), width = 8, height = 6)
      filled.contour(var1, var2, nll_matrix,
        xlab = pg_name[1],
        ylab = pg_name[2],
        main = "Likelihood surface",
        plot.axes = {
            axis(1); axis(2)
            points(exp(p0[[pg[1]]]), exp(p0[[pg[2]]]),
                pch = 19, cex = 1.5, col = "green")
            points(var1[i[1]], var2[i[2]], pch = 1, col = "black", cex = 2, lwd = 3)       
        }
        )
        plt_contour <- recordPlot()
    #dev.off()
    return(list(var1 = var1, var2 = var2, nll_matrix = nll_matrix, plt_contour = plt_contour, plt_edge = plt_edge))
}

likelihood_profile_estimate_2D_RmA <- function(df, method, pg_name, delta_list, p0, x0_bar, P0, n_grid, dt, save_fig=TRUE, save_vars = FALSE) {
    # Now, try to estimate also log_sN and log_beta while profiling log_r and log_mu


    pg <- c(paste0("log_", pg_name[1]), paste0("log_", pg_name[2]))

    var1_vals <- seq(exp(p0[[pg[1]]]) - delta_list$delta_minus_list[[pg[1]]], exp(p0[[pg[1]]]) + delta_list$delta_plus_list[[pg[1]]], length.out=n_grid)
    var2_vals <- seq(exp(p0[[pg[2]]]) - delta_list$delta_minus_list[[pg[2]]], exp(p0[[pg[2]]]) + delta_list$delta_plus_list[[pg[2]]], length.out=n_grid)

    nll_matrix_full <- matrix(NA, nrow=n_grid, ncol=n_grid)
    nll_max_grad_matrix <- array(NA, dim = c(n_grid, n_grid))
    k <- 0
    par_est <- list()
    nlminb_opt <- list()
    for (i in 1:length(var1_vals)) {
        for (j in 1:length(var2_vals)) {
            k <- k + 1
            var1_val <- var1_vals[i]
            var2_val <- var2_vals[j]
            
            model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)

            params <- list(
                    log_K = p0$log_K)

            params[[pg[1]]] <- log(var1_val)
            params[[pg[2]]] <- log(var2_val)

            do.call(model$setParameter, params)
            
            fit <- model$estimate(df, method = method,
                        ode.timestep = dt,
                        ode.solver = "rk4",
                        silent = TRUE,
                        control = list(trace = 0)
                        )

            if (!(is.null(fit))) {
                nll_value <- fit$nll
                nll__max_grad <- max(fit$nll.gradient)
                par_est[[k]] <- fit$par.fixed


                nll_matrix_full[i,j] <- nll_value
                nll_max_grad_matrix[i,j] <- nll__max_grad
                par_est[[k]] <- fit$par.fixed
                nlminb_opt[[k]] <- fit$private$opt
            } else {
                nll_value <- NA
                nll_matrix_full[i,j] <- NA
                nll_max_grad_matrix[i,j] <- NA
                par_est[[k]] <- NA
                nlminb_opt[[k]] <- NA
            }

            message(sprintf("Method %s: %s=%.2f, %s=%.2f, NLL=%.2f", method, pg_name[1], var1_val, pg_name[2], var2_val, nll_value))
        }
    }   

    # Plot likelihood surface
    i <- which(nll_matrix_full == min(nll_matrix_full), arr.ind = TRUE)
    # clip extreme values for better plotting
    nll_matrix_full_clip <- nll_matrix_full
    threshold <- min(nll_matrix_full) + 400
    nll_matrix_full_clip[nll_matrix_full > threshold] <- threshold

    if (save_fig) {
        pdf(paste0("figures/RmA_likelihood_profile_full_2D_", method, "_", pg_name[1], "_", pg_name[2], ".pdf"), width = 8, height = 6)
        filled.contour(
        var1_vals, var2_vals, nll_matrix_full_clip,
        xlab = pg_name[1],
        ylab = pg_name[2],
        main = paste0("Likelihood surface, in the ", pg_name[2], " vs. ", pg_name[1], " plane"),

        # ---- Color bar label ----
        key.title = title(main = "NLL"),

        plot.axes = {
            axis(1); axis(2)

            # True value
            points(exp(p0[[pg[1]]]), exp(p0[[pg[2]]]),
                pch = 19, cex = 1.5, col = "green")

            # Optimal (minimum NLL)
            points(var1_vals[i[1]], var2_vals[i[2]],
                pch = 1, col = "black", cex = 2, lwd = 3)

            # ---- Legend ----
            legend("topright",
                legend = c("True value", "Optimal (MLE)"),
                pch = c(19, 1),
                col = c("green", "black"),
                pt.cex = c(1.2, 1.5),
                lwd = c(NA, NA),
                bg = "white")
        }
        )
        dev.off()
    }

    if (save_vars) {
        saveRDS(list(
            var1_vals = var1_vals,
            var2_vals = var2_vals,
            nll_matrix_full = nll_matrix_full,
            nll_max_grad_matrix = nll_max_grad_matrix,
            par_est = par_est,
            nlminb_opt = nlminb_opt
        ), file = paste0("results/likelihood_profiles/RmA_likelihood_profile_full_2D_", method, "_", pg_name[1], "_", pg_name[2], ".rds"))
    }

    return(list(var1 = var1_vals, var2 = var2_vals, nll_matrix = nll_matrix_full, nll_max_grad_matrix = nll_max_grad_matrix, par_est = par_est, nlminb_opt = nlminb_opt))
}