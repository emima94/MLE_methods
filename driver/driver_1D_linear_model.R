## Driver 2D model

require(SDEtools)
require(ctsmTMB)
require(RTMB)
library(future)
library(future.apply)

set.seed(123456)

#plan(multicore) # use multicore processing

# Source all functions in src/
sofun <- function() {
    src_files <- list.files("src/OU",pattern="*.R",full.names=TRUE)
    for (f in src_files) {
        source(f)
    }
}
sofun()

source("src/plots.R")
source("src/generate_data.R")
source("src/fit_model.R")
source("src/compute_coverage.R")

#### Generate data ####

# Model defintion of the Ornstein-Uhlenbeck (O-U) process
fsim <- function(x,par=p0) with(par, c(theta*(mu - x)))
gsim <- function(x,par=p0) with(par, sigma)

T <- 100
dt <- 0.1
t <- seq(0,T,dt)

# Number of datasets to simulate
N <- 40

x0_bar = 0.3
P0 <- 0.1
x0 <- rnorm(N, x0_bar, P0)

# True parameters
p0 <- list(
    theta = 0.2,
    mu = 0.5,
    sigma = 0.2,
    obs_sd = 0.02
)

sim_data <- generate_data(fsim, gsim, t, x0, p0, dt, N)
Xsim <- sim_data$Xsim
Ysim <- sim_data$Ysim
iobs <- sim_data$iobs

matplot(t, Xsim,
        type = "l",
        lty = 1,
        col = rgb(0, 0, 0, 0.2),
        ylim = c(0, 1),
        xlab = "Time",
        ylab = "X")


#### Fit models ####
# Get ctsmTMB model
sofun()
model <- create_OU_model(obs_sd = p0$obs_sd, x0 = x0_bar, p0 = P0)

par_names <- rownames(model$getParameters(type = "free"))
n_par <- length(par_names)

# fit to each dataset with each method
methods <- c("ekf", "laplace", "laplace2")

# Test over a range of noise lvls:
sigma_vals <- c(0.05, 0.1, 0.2)
fit_sigma <- vector("list", length(sigma_vals))
names(fit_sigma) <- paste0("sigma_", sigma_vals)

for (sigma_val in sigma_vals) {

    message("Fitting for sigma = ", sigma_val, sep="")
    p0$sigma <- sigma_val

    sim_data <- generate_data(fsim, gsim, t, x0, p0, dt, N)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_sigma[[paste0("sigma_", sigma_val)]] <- fit_model(methods, model, t, Ysim, iobs, dt, N)
}

# Summarize fits
metric_names <- c("coverage", "bias", "rmse")

fit_summary_sigma <- vector("list", length(sigma_vals))
names(fit_summary_sigma) <- names(fit_sigma)
for (sigma_val in sigma_vals) {
    fit <- fit_sigma[[paste0("sigma_", sigma_val)]]
    fit_summary <- vector("list", length(methods))
    names(fit_summary) <- methods
    for (m in methods) {

        par_est_m <- matrix(0, nrow=N, ncol=n_par)
        colnames(par_est_m) <- par_names

        par_se_m <- matrix(0, nrow=N, ncol=n_par)
        colnames(par_se_m) <- par_names

        for (par in par_names) {
            par_est_m[,par] <- sapply(fit[[m]], function(x) x$par.fixed[par])
            par_se_m[,par] <- sapply(fit[[m]], function(x) sqrt(x$cov.fixed[par,par]))
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
            perf_metric_m[par,"bias"] <- mean(par_est_m[,par] - p0[[par]])
            perf_metric_m[par,"rmse"] <- sqrt(mean((par_est_m[,par] - p0[[par]])^2))
        }
        
        fit_summary[[m]] <- list(est = par_est_m, se = par_se_m, perf = perf_metric_m)
    }
    fit_summary_sigma[[paste0("sigma_", sigma_val)]] <- fit_summary
}
# Plot metric vs noise level for each method for all parameters
ylims <- list(c(0,1), c(-0.2,0.2), c(0,0.3))
names(ylims) <- metric_names
plot_metric <- function(metric_name) {
    par(mfrow = c(n_par, 1))

    cols <- seq_along(methods)  # consistent colors

    for (par in par_names) {

        plot(sigma_vals,
             sapply(sigma_vals, function(sigma_val)
                 fit_summary_sigma[[paste0("sigma_", sigma_val)]][[methods[1]]]$perf[par, metric_name]),
             type = "b",
             col  = cols[1],
             ylim = ylims[[metric_name]],
             xlab = "Sigma",
             ylab = metric_name,
             main = paste("Parameter:", par))

        for (m in methods[-1]) {
            points(sigma_vals,
                   sapply(sigma_vals, function(sigma_val)
                       fit_summary_sigma[[paste0("sigma_", sigma_val)]][[m]]$perf[par, metric_name]),
                   type = "b",
                   col  = cols[methods == m])
        }

        legend("topright",
               legend = methods,
               col    = cols,
               lty    = 1,
               pch    = 1,
               bty    = "n")
    }
}
plot_metric("coverage")
plot_metric("bias")
plot_metric("rmse")
# For a given noise level, show parameter estimate distributions
sigma_val <- 0.05
fit_summary <- fit_summary_sigma[[paste0("sigma_", sigma_val)]]
# Plot distribution of parameter estimates for every method for every parameter (n_par x n_methods plots)
par(mfrow=c(n_par, length(methods)))
for (par in par_names) {
    for (m in methods) {
        hist(fit_summary[[m]]$est[,par],
             main = paste("Estimates of ", par, " (", m, ")", sep=""),
             xlab = par)
        abline(v = p0[[par]], col="red", lwd=2)
    }
}
# Plot distribution of standard errors for every method for every parameter (n_par x n_methods plots)
par(mfrow=c(n_par, length(methods)))
for (par in par_names) {
    for (m in methods) {
        hist(fit_summary[[m]]$se[,par],
             main = paste("SE of ", par, " (", m, ")", sep=""),
             xlab = par)
    }
}

fit_summary$ekf$perf
fit_summary$laplace$perf
fit_summary$laplace2$perf



