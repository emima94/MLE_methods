require(ctsmTMB)
require(SDEtools)
require(RTMB)
library(parallel)


# Source all functions in src/
sofun <- function() {
    src_files <- list.files("src/RmA",pattern="*.R",full.names=TRUE)
    for (f in src_files) {
        source(f)
    }
    source("src/plots.R")
    source("src/generate_data.R")
    source("src/generate_data_RmA.R")
    source("src/fit_model.R")
    source("src/compute_coverage.R")
    source("src/compute_fit_summary.R")
    source("src/eval_results.R")
}
sofun()

# Results output directory
#out_dir <- "results/driver_RmA_sweep_T"
#out_dir <- "hpc/results_and_figs/driver_RmA_sweep_T_260129_1017"
#out_dir <- "results/RmA_model_NP_sweep_T"
out_dir <- "results/RmA_model_log_state_trans_TRUE_partial_obs_FALSE_sweep_T_20260130_102557"


# Load settings
set <- readRDS(file.path(out_dir, "T_settings.rds"))
T_vals <- set$T_vals
methods <- set$methods
n_par <- set$n_par
par_names <- set$par_names
N <- set$N
obs_partial <- set$obs_partial
log_transform_state <- set$log_transform_state

# Load key results for each T value for each method for each dataset
# Read parameter function
get_file_path <- function(T_val, method, i, out_dir) {
    #fit_file <- file.path(out_dir, paste0("T_", T_val), sprintf("%s_%04d.rds", method, i))
    fit_file <- file.path(out_dir, "light", paste0("T_", T_val), method, sprintf("%04d.rds", i))
    return(fit_file)
}

# read_res <- function(T_val, method, i, out_dir) {
#     fit <- readRDS(get_file_path(T_val, method, i, out_dir))

#     res <- list(
#         par = fit$fit$par.fixed,
#         nll = fit$fit$nll,
#         nll.gradient = fit$fit$nll.gradient,
#         sd = fit$fit$sd.fixed,
#         opt = fit$fit$private$opt,
#         time = c(estimation = fit$fit$private$timer_estimation, 
#                  compile = fit$fit$private$timer_construct_adfun)
#     )

#     return(res)
# }

N_test <- N

# Load results into list (testing)
res <- lapply(T_vals, function(T_val) {
    res_method <- lapply(methods, function(m) {
        message("Loading results for T = ", T_val, ", method = ", m)
        res_datasets <- lapply(1:N_test, function(i) {
            readRDS(get_file_path(T_val, m, i, out_dir))
        })
        return(res_datasets)
    })
    names(res_method) <- methods
    return(res_method)
})

res

# Parameter histograms for each T value and method
xlims = list(
            log_r = c(-1,1),
            log_K = c(-1,1),
            log_beta = c(-2,2),
            log_mu = c(-2,2),
            log_sN = c(-16,2)
)
xlims = list(
            r = c(0.3,2.0),
            K = c(0.3,2.0),
            beta = c(0.1,6),
            mu = c(0.1,3),
            sN = c(0.01,2)
)


for (i in seq_along(T_vals)) {
    T_val <- T_vals[i]
    par(mfrow = c(length(methods), n_par))
    for (m in methods) {
        par_mat <- matrix(NA, nrow=N_test, ncol=n_par)
        colnames(par_mat) <- par_names
        sd_mat <- matrix(NA, nrow=N_test, ncol=n_par)
        colnames(sd_mat) <- par_names
        for (j in seq_len(N_test)) {
            par_mat[j, ] <- res[[i]][[m]][[j]]$par
            sd_mat[j, ] <- res[[i]][[m]][[j]]$sd
        } 
        for (k in seq_len(n_par)) {
            p_true <- set$p0[[par_names[k]]]
            p_est_mean <- mean(par_mat[, k])
            p_sd_mean <- mean(sd_mat[ ,k])
            hist(par_mat[, k], 
                    main = paste0("T=", T_val, ", method=", m, ", param=", par_names[k], "\n obs_partial=", obs_partial, ", log_transform_state=", log_transform_state),
                    xlab = "Parameter estimate", 
                    xlim = xlims[[par_names[k]]],
                    breaks = 20)
            abline(v = p_true, col="red", lwd=2)
            abline(v = p_est_mean, col="green", lwd=1, lty=2)
            abline(v = p_est_mean - 1.96 * p_sd_mean, col="blue", lwd=1, lty=2)
            abline(v = p_est_mean + 1.96 * p_sd_mean, col="blue", lwd=1, lty=2)
        }
    }
}

# Plot nll histograms for each T value and method
{
    par(mfrow = c(length(methods), length(T_vals)))
    for (m in methods) {
        for (i in seq_along(T_vals)) {
            T_val <- T_vals[i]
            nll_vals <- sapply(res[[i]][[m]], function(r) r$nll)
            hist(nll_vals, 
                main = paste0("T=", T_val, ", method=", m, ", NLL", "\n obs_partial=", obs_partial, ", log_transform_state=", log_transform_state),
                xlab = "Negative log-likelihood", 
                breaks = 20)
        }
    }
}

res[[1]][[1]][[1]]$nll.gradient

# Plot computation time vs T_vals
{
    est_time_mat <- matrix(NA, nrow=length(T_vals), ncol=length(methods))
    colnames(est_time_mat) <- methods
    rownames(est_time_mat) <- paste0("T_", T_vals)
    compile_time_mat <- matrix(NA, nrow=length(T_vals), ncol=length(methods))
    colnames(compile_time_mat) <- methods
    rownames(compile_time_mat) <- paste0("T_", T_vals)
    for (i in seq_along(T_vals)) {
        T_val <- T_vals[i]
        for (m in methods) {
            #est_times <- sapply(res[[i]][[m]], function(r) r$time["estimation.elapsed"])
            #compile_times <- sapply(res[[i]][[m]], function(r) r$time["compile.elapsed"])
            est_times <- sapply(res[[i]][[m]], function(r) r$time.elapsed)
            compile_times <- sapply(res[[i]][[m]], function(r) r$time.compile)
            est_time_mat[i, m] <- mean(est_times)
            compile_time_mat[i, m] <- mean(compile_times)
        }
    }
    var_name <- "T"
    ylims <- c(0, max(est_time_mat+compile_time_mat, na.rm=TRUE))
    par(mfrow=c(1,3))
    matplot(T_vals, est_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab=paste0("Values of ", var_name),
        ylab="Average estimation time (s)",
        ylim = ylims,
        main=paste0("Estimation time vs ", var_name))
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    matplot(T_vals, compile_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab=paste0("Values of ", var_name),
        ylab="Average compilation time (s)",
        ylim = ylims,
        main=paste0("Compilation time vs ", var_name))
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    matplot(T_vals, est_time_mat + compile_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab=paste0("Values of ", var_name),
        ylab="Average total time (s)",
        ylim = ylims,
        main=paste0("Total computation time vs ", var_name))
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
}

# Plot RMSE and coverage vs. T_vals for each parameter for each method
{
    RMSE <- list()
    Coverage <- list()
    for (par_name in par_names) {
        RMSE[[par_name]] <- matrix(NA, nrow=length(T_vals), ncol=length(methods))
        colnames(RMSE[[par_name]]) <- methods
        rownames(RMSE[[par_name]]) <- paste0("T_", T_vals)
        Coverage[[par_name]] <- matrix(NA, nrow=length(T_vals), ncol=length(methods))
        colnames(Coverage[[par_name]]) <- methods
        rownames(Coverage[[par_name]]) <- paste0("T_", T_vals)

        for (i in seq_along(T_vals)) {
            T_val <- T_vals[i]
            for (m in methods) {
                par_estimates <- sapply(res[[i]][[m]], function(r) r$par[par_name])
                par_sds <- sapply(res[[i]][[m]], function(r) r$sd[par_name])
                p0_value <- set$p0[[par_name]]
                # RMSE
                RMSE[[par_name]][i, m] <- sqrt(mean((par_estimates - p0_value)^2, na.rm=TRUE))
                # Coverage
                lower_bounds <- par_estimates - 1.96 * par_sds
                upper_bounds <- par_estimates + 1.96 * par_sds
                coverage_vals <- (p0_value >= lower_bounds) & (p0_value <= upper_bounds)
                Coverage[[par_name]][i, m] <- mean(coverage_vals, na.rm=TRUE)
            }
        }
    }
}
RMSE
Coverage

# Plot RMSE vs T_vals
{
    par(mfrow=c(n_par,1))
    for (par_name in par_names) {
        matplot(T_vals, RMSE[[par_name]], type="b", pch=1, lty=1, col=1:length(methods),
            xlab=paste0("Values of ", var_name),
            ylab=paste0("RMSE of ", par_name),
            main=paste0("RMSE of ", par_name, " vs ", var_name))
        legend("topright", legend=methods, col=1:length(methods), pch=1, lty=1)
    }
}
# Plot Coverage vs T_vals
{
    par(mfrow=c(n_par,1))
    for (par_name in par_names) {
        matplot(T_vals, Coverage[[par_name]], type="b", pch=1, lty=1, col=1:length(methods),
            xlab=paste0("Values of ", var_name),
            ylab=paste0("Coverage of ", par_name),
            main=paste0("Coverage of ", par_name, " vs ", var_name),
            ylim=c(0,1))
        abline(h=0.95, col="red", lty=2)
        legend("topright", legend=methods, col=1:length(methods), pch=1, lty=1)
    }
}

# Did the estimations converge (in percentage) for each T value and method?
{
    conv_counts <- matrix(0, nrow=length(T_vals), ncol=length(methods))
    colnames(conv_counts) <- methods
    rownames(conv_counts) <- paste0("T_", T_vals)
    for (i in seq_along(T_vals)) {
        T_val <- T_vals[i]
        for (m in methods) {
            conv_vals <- sapply(res[[i]][[m]], function(r) r$opt$convergence == 0)
            conv_counts[i, m] <- mean(conv_vals, na.rm=TRUE) * 100
        }
    }
    conv_counts
}
# What is the max gradient nll for each T value and method?
{
    max_grad_mat <- matrix(NA, nrow=length(T_vals), ncol=length(methods))
    colnames(max_grad_mat) <- methods
    rownames(max_grad_mat) <- paste0("T_", T_vals)
    for (i in seq_along(T_vals)) {
        T_val <- T_vals[i]
        for (m in methods) {
            max_grads <- sapply(res[[i]][[m]], function(r) max(r$nll.gradient))
            max_grad_mat[i, m] <- max(max_grads, na.rm=TRUE)
        }
    }
    max_grad_mat
}
