compute_fit_summary <- function(fit, p0, methods, n_par, par_names, metric_names) {
    fit_summary <- vector("list", length(methods))
    names(fit_summary) <- methods
    for (m in methods) {


        par_est_m <- matrix(0, nrow=N, ncol=n_par)
        colnames(par_est_m) <- par_names

        par_se_m <- matrix(0, nrow=N, ncol=n_par)
        colnames(par_se_m) <- par_names

        for (par in par_names) {
            message(" Processing parameter ", par, " for method ", m, sep="")
            for (i in 1:N) {
                
                message("  Dataset ", i, sep="")
                x <- fit[[m]][[i]]$fit  

                par_est_m[i,par] <- if (is.null(x)) {
                    NaN
                } else {
                    x$par.fixed[par]
                }

                par_se_m[i,par] <- if (is.null(x)) {
                    NaN
                } else {
                    cov <- x$cov.fixed
                    colnames(cov) <- rownames(cov) <- par_names
                    v <- cov[par, par]
                    

                    if (is.na(v) || !is.finite(v) || v < 0) {
                        NaN
                    } else {
                        sqrt(v)
                    }
                }
                
            }
        }

        # Performance metrics for method m: coverage, bias, rmse, store as a named matrix
        perf_metric_m <- matrix(0, nrow=n_par, ncol=3)
        rownames(perf_metric_m) <- par_names
        colnames(perf_metric_m) <- metric_names
        for (par in par_names) {
            ## Compute coverage
            perf_metric_m[par,"coverage"] <- compute_coverage(
                estimates = par_est_m[,par],
                ses = par_se_m[,par],
                true_value = p0[[par]],
                alpha = 0.05
            )
            ## Compute bias and RMSE
            perf_metric_m[par,"bias"] <- mean(par_est_m[,par] - p0[[par]], na.rm=TRUE)
            perf_metric_m[par,"rmse"] <- sqrt(mean((par_est_m[,par] - p0[[par]])^2, na.rm=TRUE))
        }
        
        fit_summary[[m]] <- list(est = par_est_m, se = par_se_m, perf = perf_metric_m)
    }

    return(fit_summary)

}