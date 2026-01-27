## Eval results (build plots)

#### Computation time vs var
create_computation_time_plots <- function(fit, vals, methods, var_name, file_prefix) {
    " Create computation time plots for estimation, compilation, and total time vs a sweeped variable.

    Args:
        fit: List of fits for different values of the variable.
        vals: Vector of variable values sweeped.
        methods: Vector of method names.
        var_name: Name of the sweeped variable.
        file_prefix: Prefix for the output file names.
    "

    est_time_mat <- matrix(0, nrow=length(vals), ncol=length(methods))
    rownames(est_time_mat) <- paste0(var_name, "_", vals)
    colnames(est_time_mat) <- methods
    compile_time_mat <- matrix(0, nrow=length(vals), ncol=length(methods))
    rownames(compile_time_mat) <- paste0(var_name, "_", vals)
    colnames(compile_time_mat) <- methods

    null_count <- matrix(0, nrow=length(vals), ncol=length(methods))
    colnames(null_count) <- methods

    for (i in 1:length(vals)) {
        for (m in methods) {

            fit_im <- fit[[i]][[m]]
            is_valid <- !sapply(fit_im, function(fit) is.null(fit$fit))
            null_count[i,m] <- sum(!is_valid)

            valid_fits <- fit_im[is_valid]

            if (length(valid_fits) == 0) {
                est_time_mat[i,m] <- NA
                compile_time_mat[i,m] <- NA
                next
            }

            est_time_mat[i,m] <- mean(sapply(valid_fits, function(res) res$fit$private$timer_estimation))
            compile_time_mat[i,m] <- mean(sapply(valid_fits, function(res) res$fit$private$timer_construct_adfun))


            #est_time_mat[i,m] <- mean(sapply(fit[[i]][[m]], function(res) res$fit$private$timer_estimation))
            #compile_time_mat[i,m] <- mean(sapply(fit[[i]][[m]], function(res) res$fit$private$timer_construct_adfun))
        }
    }
    
    ylims <- c(0, max(est_time_mat+compile_time_mat, na.rm=TRUE))
    pdf(paste0("figures/", file_prefix, "_computation_time", "_", var_name, ".pdf"), width=8, height=3)
    par(mfrow=c(1,3))
    matplot(vals, est_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab=paste0("Values of ", var_name),
        ylab="Average estimation time (s)",
        ylim = ylims,
        main=paste0("Estimation time vs ", var_name))
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    matplot(vals, compile_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab=paste0("Values of ", var_name),
        ylab="Average compilation time (s)",
        ylim = ylims,
        main=paste0("Compilation time vs ", var_name))
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    matplot(vals, est_time_mat + compile_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab=paste0("Values of ", var_name),
        ylab="Average total time (s)",
        ylim = ylims,
        main=paste0("Total computation time vs ", var_name))
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    dev.off()

    return(null_count)
}

#### Parameter estimation distributions ####
create_parameter_estimation_plots <- function(fit_summary_var, vals, methods, var_name, file_prefix, par_names, p0, xlims, var_name_title, var_title_func=NULL) {

    " Create parameter estimation distribution plots for each parameter and method and variable value.
    The variable value is the sweeped variable

    Args:
        fit_summary_var: List of fit summaries for different values of the variable.
        vals: Vector of variable values sweeped.
        methods: Vector of method names.
        var_name: Name of the sweeped variable OBS: should match parameter name if it is a parameter.
        file_prefix: Prefix for the output file names.
        par_names: Vector of FREE parameter names.
        p0: Named vector of true parameter values (both free and fixed).
        xlims: Named list of x-axis limits for each parameter.
        var_name_title: Title for the sweeped variable to be used in plots.
    "

    n_par <- length(par_names)
    N <- nrow(fit_summary_var[[1]][[1]]$est)
    for (i in 1:length(vals)) {
        fit_summary <- fit_summary_var[[i]]

        p0_val <- p0
                # If var_name is a parameter, set p0_val[var_name] to vals[i]
                if (var_name %in% par_names) {
                    p0_val[[var_name]] <- vals[i]
                }

        if (is.null(var_title_func)) {
                val_title <- vals[i]
            } else {
                val_title <- var_title_func(vals[i])
            }   
        #print(p0_val[[var_name]])
        #print(vals[i])
        
        # Plot histograms of parameter estimates
        pdf(paste0("figures/", file_prefix, "_parameter_estimates_", var_name, "_", i, ".pdf"), width=12, height=8)
        par(mfrow=c(length(methods), n_par))
        for (j in seq_along(methods)) {
            for (k in seq_along(par_names)) {
                m <- methods[j]
                par <- par_names[k]
                if (all(is.nan(fit_summary[[m]]$est[,par]))) {
                    next
                }

                

                if (j == 1 && k == 1) {
                    main <- paste0("N=", N, ", ", var_name_title, "=", val_title)
                } else {
                    main <- ""
                }
                if (j == length(methods)) {
                    xlab <- par
                } else {
                    xlab <- ""
                }
                if (k == 1) {
                    ylab <- paste0("density (", m, ")")
                } else {
                    ylab <- ""
                }

                hist(fit_summary[[m]]$est[,par],
                     main = main,
                     xlab = par,
                     ylab = ylab,
                     xlim = xlims[[par]]
                     #ylim = NULL
                     )
                abline(v = p0_val[[par]], col="red", lwd=2)
                est_par_mean <- mean(fit_summary[[m]]$est[,par], na.rm=TRUE)
                sd_par_mean <- mean(fit_summary[[m]]$se[,par], na.rm=TRUE)
                # Plot mean estimation and estimation of 95 % CI with mean sd and mean estimate
                abline(v = est_par_mean, col="green", lwd=1, lty=2)
                abline(v = est_par_mean - 1.96 * sd_par_mean, col="blue", lwd=1, lty=2)
                abline(v = est_par_mean + 1.96 * sd_par_mean, col="blue", lwd=1, lty=2)
                legend("topright",
                       legend = c("True value", "Mean estimate", "95% CI (mean sd)"),
                       col = c("red", "green", "blue"),
                       lty = c(1,2,2),
                       lwd = 2)
            }
        }
        dev.off()
    }

    
}

#### Plot RMSE of parameter estimates vs sweeped variable ####
create_perf_metric_plot <- function(fit_summary_var, vals, methods, var_name, file_prefix, par_names, metric_name) {
    " Create performance metric plots for each parameter and method vs a sweeped variable.

    Args:
        fit_summary_var: List of fit summaries for different values of the variable.
        vals: Vector of variable values sweeped.
        methods: Vector of method names.
        var_name: Name of the sweeped variable.
        file_prefix: Prefix for the output file names.
        par_names: Vector of FREE parameter names.
        metric_name: Name of the performance metric to plot available metrics: rmse, bias, coverage
    "

    if (!(metric_name %in% c("rmse", "bias", "coverage"))) {
        stop("Invalid metric_name. Choose from: rmse, bias, coverage")
    }

    if (metric_name == "coverage") {
        ylims <- c(0,1)
    } else {
        ylims <- NULL
    }

    for (i in 1:length(par_names)) {
        par <- par_names[i]
        pdf(paste0("figures/", file_prefix, "_", metric_name, "_vs_", var_name, "_", i, ".pdf"), width=6, height=4)
        par(mfrow=c(1,1))
        metric_par <- sapply(vals, function(val) {
            fit_summary <- fit_summary_var[[which(vals == val)]]
            sapply(methods, function(m) fit_summary[[m]]$perf[par, metric_name])
        })
        matplot(vals, t(metric_par), type="b", pch=1, lty=1, col=1:length(methods),
            xlab=paste("Values of ", var_name, sep=""),
            ylab=paste(toupper(metric_name), " of ", par, " estimates", sep=""),
            ylim = ylims,
            main=paste(toupper(metric_name), " of ", par, " estimates vs ", var_name, sep=""))
        legend("topright", legend=methods, col=1:length(methods), pch=1, lty=1)
        dev.off()
    }
}


# Plot histogram of likelihood for all methods and values in sweeped variable
create_likelihood_histograms <- function(fit_var, vals, methods, var_name, theta_true, file_prefix, var_name_title, var_title_func=NULL, plot_dev=FALSE) {
    " Create likelihood histograms for each method and variable value.

    Args:
        fit_var: List of fits for different values of the variable.
        vals: Vector of variable values sweeped.
        methods: Vector of method names.
        var_name: Name of the sweeped variable.
        file_prefix: Prefix for the output file names.
        var_name_title: Title for the sweeped variable to use in plots.
        var_title_func: Optional function to transform variable values for display.
    "

    n_vals <- length(vals)
    n_methods <- length(methods)

    if (plot_dev) {
        file_prefix <- paste0(file_prefix, "_dev")
    }

    pdf(paste0("figures/", file_prefix, "_likelihood_histograms_", var_name_title, ".pdf"), width=12, height=8)
    par(mfrow=c(n_methods, n_vals))
    for (i in 1:n_methods) {
        m <- methods[i]
        for (j in 1:n_vals) {

            val <- vals[j]
            fit_ij <- fit_var[[j]][[m]]
            loglik_vals <- sapply(fit_ij, function(res) {
                if (is.null(res$fit)) {
                    return(NA)
                } else {
                    return(res$fit$nll)
                }
            })

            # Get all true likelihoods
            if (plot_dev) {
                nll_true_vals <- sapply(fit_ij, function(res) {
                    if (is.null(res$fit)) {
                        return(NA)
                    } else {
                        return(res$fit$nll_funcs$fn(theta_true))
                    }
                })
                y <- loglik_vals - nll_true_vals
            } else {
                y <- loglik_vals
            }

            if (is.null(var_title_func)) {
                val_title <- vals[j]
            } else {
                val_title <- var_title_func(vals[j])
            }   

            if (i == 1) {
                main <- paste0(var_name_title, "=", val_title)
            } else {
                main <- ""
            }
            if (j == 1) {
                ylab <- paste0("Density (", m, ")")
            } else {
                ylab <- ""
            }
            xlab = "nll"

            hist(y,
                 xlab=xlab, ylab = ylab, main = main)
        }
    }
    dev.off()
}
