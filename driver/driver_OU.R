## Driver 1D model: The Ornstein-Uhlenbeck process

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
source("src/compute_fit_summary.R")

#### Generate data ####

# Model defintion of the Ornstein-Uhlenbeck (O-U) process
fsim <- function(x,par=p0) with(par, c(theta*(mu - x)))
gsim <- function(x,par=p0) with(par, sigma)
hsim <- function(x,par=p0) x # identity observation function


T <- 100
dt <- 0.1
t <- seq(0,T,dt)

# Number of datasets to simulate
N <- 3

x0_bar = 0.3
P0 <- 0.1
x0 <- matrix(rnorm(N, x0_bar, P0), ncol=1)

# True parameters
p0 <- list(
    theta = 0.2,
    mu = 0.5,
    sigma = 0.1,
    obs_sd = 0.02
)

sim_data <- generate_data(fsim, gsim, hsim, t, x0, p0, dt, N)
Xsim <- sim_data$Xsim
Ysim <- sim_data$Ysim
iobs <- sim_data$iobs

matplot(t, do.call(cbind, Xsim),
        type = "l",
        lty = 1,
        col = rgb(0, 0, 0, 0.2),
        ylim = c(0, 1),
        xlab = "Time",
        ylab = "X")

matpoints(t[iobs], do.call(cbind, Ysim),
           pch = 1,
           col = rgb(1, 0, 0, 0.2),
           cex = 0.5)


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
    p0_k <- p0
    p0_k$sigma <- sigma_val

    sim_data <- generate_data(fsim, gsim, hsim, t, x0, p0_k, dt, N)
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
    # True parameters
    p0_k <- p0
    p0_k$sigma <- sigma_val
    fit_summary <- compute_fit_summary(fit, p0_k, methods, n_par, par_names, metric_names)
    
    fit_summary_sigma[[paste0("sigma_", sigma_val)]] <- fit_summary
}


# Plot metric vs noise level for each method for all parameters
ylims <- list(c(0,1), c(-0.2,0.2), c(0,0.3))
names(ylims) <- metric_names

# Plot metrics vs noise level
plot_metric("coverage", fit_summary_sigma, sigma_vals, "sigma", ylims, methods, par_names)
plot_metric("bias", fit_summary_sigma, sigma_vals, "sigma", ylims, methods, par_names)
plot_metric("rmse", fit_summary_sigma, sigma_vals, "sigma", ylims, methods, par_names)

# For a given noise level, show parameter estimate distributions
sigma_val <- sigma_vals[2]
p0_k <- p0
p0_k$sigma <- sigma_val
fit_summary <- fit_summary_sigma[[paste0("sigma_", sigma_val)]]
# Plot distribution of parameter estimates for every method for every parameter (n_par x n_methods plots)
par(mfrow=c(n_par, length(methods)))
for (par in par_names) {
    for (m in methods) {
        hist(fit_summary[[m]]$est[,par],
             main = paste("Estimates of ", par, " (", m, ")", sep=""),
             xlab = par)
        abline(v = p0_k[[par]], col="red", lwd=2)
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

# Test different lengths of time series
T_vals <- c(20, 50, 100, 200)
fit_T <- vector("list", length(T_vals))
names(fit_T) <- paste0("T_", T_vals)

for (T_val in T_vals) {

    message("Fitting for T = ", T_val, sep="")
    t <- seq(0,T_val,dt)

    sim_data <- generate_data(fsim, gsim, hsim, t, x0, p0, dt, N)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_T[[paste0("T_", T_val)]] <- fit_model(methods, model, t, Ysim, iobs, dt, N)
}

# Summarize fits
fit_summary_T <- vector("list", length(T_vals))
names(fit_summary_T) <- names(fit_T)

for (T_val in T_vals) {
    fit <- fit_T[[paste0("T_", T_val)]]

    fit_summary <- compute_fit_summary(fit, p0, methods, n_par, par_names, metric_names)
    
    fit_summary_T[[paste0("T_", T_val)]] <- fit_summary
}
fit_summary_T

 
plot_metric("coverage", fit_summary_T, T_vals, "T", ylims, methods, par_names)
plot_metric("bias", fit_summary_T, T_vals, "T", ylims, methods, par_names)
plot_metric("rmse", fit_summary_T, T_vals, "T", ylims, methods, par_names)
