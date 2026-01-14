## Driver Heat equtaion model (1D diffusion equation)

require(SDEtools)
require(ctsmTMB)
require(RTMB)
library(future)
library(future.apply)

set.seed(123456)

#plan(multicore) # use multicore processing

# Source all functions in src/
sofun <- function() {
    src_files <- list.files("src/HeatEq",pattern="*.R",full.names=TRUE)
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
source("src/util.R")

#### Generate data ####

## Heat equation model in 1D
T <- 100
dt <- 0.1
t <- seq(0,T,dt)

# Number of datasets to simulate
N <- 20

# Number of spatial finite volumes
Nx <- 10

x0_bar = runif(Nx, min=2, max=8)
#x0_bar = seq(8,2,length.out=Nx)
#x0_bar <- rep(5.0, Nx)
#P0 <- diag(rep(1, Nx))
P0 <- diag(rep(0.001, Nx))
x0 <- matrix(rnorm(N * Nx, x0_bar, diag(sqrt(P0))), nrow=N, ncol=Nx, byrow=TRUE)

# True parameters
p0 <- list(
    log_K = log(0.1),
    log_sT = log(0.05), #0.1
    T0 = 0.0,
    TL = 10.0,
    dx = 1.0,
    obs_sd = 0.05,
    sensor_pos = c(1,3,7,9)
)

source("src/generate_data.R")
sofun()
sim_data <- generate_data(fsim, gsim, hsim, t, x0, p0, dt, N)
Xsim <- sim_data$Xsim
Ysim <- sim_data$Ysim
iobs <- sim_data$iobs

# Plot simulated data
#k <- 2
#matplot(t, Xsim[[k]], type="l")
#matpoints(t[iobs], Ysim[[k]], col="red", pch=1)

# Fit model to data
model <- create_HeatEq_model(
    obs_sd = p0$obs_sd,
    x0 = x0_bar,
    P0 = P0,
    Nx = Nx,
    dx = p0$dx,
    T0 = p0$T0,
    TL = p0$TL,
    sensor_pos = p0$sensor_pos
)

par_names <- rownames(model$getParameters(type = "free"))
n_par <- length(par_names)

methods = c("ekf", "laplace", "laplace.thygesen")

metric_names <- c("coverage", "bias", "rmse")
source("src/compute_fit_summary.R")
source("src/util.R")

# Test parallel par estimation
library(parallel)
# Use all available cores
n_clusters <- detectCores() - 2
n_clusters
  source("src/fit_model.R")
time_par <- system.time({

    fit_par <- fit_model_par(n_clusters, methods, model, t, Ysim, iobs, dt, N, df_fun_HeatEq)
})

print(time_par)

print_computation_time(fit_par, methods)

fit_summary <- compute_fit_summary(fit_par, p0, methods, n_par, par_names, metric_names)

n_par_plot <- c(1,2)
# Define xlims for parameters in histograms
xlims <- list(
    log_K = c(-2.8,-2),
    log_sT = c(-3.5,-2.5)
)
# Compute mean of sds for each parameter
sd <- mean(fit_summary[[1]]$se[,"log_K"], na.rm=TRUE)
# Plot histograms of parameter estimates
par(mfrow=c(length(n_par_plot), length(methods)))
for (par in par_names[n_par_plot]) {
    for (m in methods) {
        hist(fit_summary[[m]]$est[,par],
             main = paste("Estimates of ", par, " (", m, ")", sep=""),
             xlab = par,
             xlim = xlims[[par]])
        abline(v = p0[[par]], col="red", lwd=2)
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

# Compute coverage
coverage_ekf <- fit_summary$ekf$perf[,"coverage"]
coverage_laplace <- fit_summary$laplace$perf[,"coverage"]
coverage_laplace_thygesen <- fit_summary$laplace.thygesen$perf[,"coverage"]

coverage_matrix <- cbind(
    coverage_ekf,
    coverage_laplace,
    coverage_laplace_thygesen
)

colnames(coverage_matrix) <- methods
coverage_matrix


# Plot histograms of computation times
par(mfrow=c(1,length(methods)))
for (m in methods) {
    hist(sapply(fit_par[[m]], function(res) res$time),
         main = paste("Computation times (", m, ")", sep=""),
         xlab = "Time (s)")
}

# Plot histograms of likelihoods
par(mfrow=c(1,length(methods)))
for (m in methods) {
    hist(sapply(fit_par[[m]], function(res) res$fit$nll),
         main = paste("Negative log-likelihoods (", m, ")", sep=""),
         xlab = "Negative log-likelihood")
}

