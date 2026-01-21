## Driver 2D model: The RmA model

require(SDEtools)
require(ctsmTMB)
require(RTMB)
library(parallel)
set.seed(123456)

#plan(multicore) # use multicore processing

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
}
sofun()

# Set font in plots to Times New Roman 

#### Generate data ####

## RmA/Bazykin-model in natural coordinates
T <- 100
dt <- 0.1
t <- seq(0,T,dt)

# Number of datasets to simulate
N <- 5

x0_bar = c(0.5, 0.5)
P0 <- diag(c(0.001,0.001))
x0 <- log(matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar)))
N0 <- exp(x0)
N0

# True parameters
# p0 <- list(
#     r = 2.0,
#     K = 100,
#     epsilon = 0.3,
#     beta = 0.06,
#     Cmax = 3.0,
#     mu = 0.2,
#     sN = 0.2, #0.1
#     sP = 0.1,  #0.1
#     obs_sd = 0.1
# )
p0 <- list(
    log_r = log(1.0),
    log_K = log(1.0),
    epsilon = 3.0,
    log_beta = log(3.0),
    Cmax = 1.0,
    log_mu = log(1.0),
    log_sN = log(0.2), #0.1
    log_sP = log(0.1),  #0.1
    obs_sd = 0.05
)
source("src/generate_data_RmA.R")
sofun()
sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0, p0, dt, N)
Xsim <- sim_data$Xsim
Ysim <- sim_data$Ysim
iobs <- sim_data$iobs
Nsim <- sapply(Xsim, function(x) exp(x[,1]), simplify=FALSE)


# Plot example data
scale <- 0.7
ratio <- 0.3
k <- c(1,2,3)
for (k in 1:3){
    pdf(sprintf("figures/RmA-sim-obs_prey_%d.pdf", k), width = 10*scale, height = 10*ratio*scale)
    par(mar = c(4.0, 4.0, 1.0, 1.0))
    plot(t,exp(Xsim[[k]][,1]), type="l", col="black", xlab = "Time", ylab = "N",
         family = "serif") 
    points(t[iobs], Ysim[[k]], col="red")
    dev.off()
    pdf(sprintf("figures/RmA-sim-obs_predator_%d.pdf", k), width = 10*scale, height = 10*ratio*scale)
    par(mar = c(4.0, 4.0, 1.0, 1.0))   
    plot(t,exp(Xsim[[k]][,2]), type="l", col="black", xlab = "Time", ylab = "P", family = "serif")
    dev.off()
}

matplot(t, do.call(cbind, Nsim),
        type = "l",
        lty = 1,
        col = rgb(0, 0, 0, 0.2),
        xlab = "Time",
        ylab = "N"
        )

matpoints(t[iobs], do.call(cbind, Ysim),
          pch = 16, col = rgb(1, 0, 0, 0.5))

# Are any observations negative?
is_Y_neg <- any(sapply(Ysim, function(y) any(y < 0)))
if (is_Y_neg) {
    message("Some generated observations are negative!")
    } else {
       message("No negative observations generated.")
    }



#### Fit models ####
# Get ctsmTMB model
sofun()
model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)

# # Test the first 10 datasets
# for (k in 1:3) {
#     message("Fitting dataset ", k, sep="")
# df <- data.frame(
#     t = t[iobs],
#     Y = Ysim[[k]]

# )
# message("Fitting with EKF...")
# fit_ekf <- model$estimate(
#     data = df,
#     method = "ekf",
#     ode.solver = "rk4",
#     ode.timestep = dt
#     )
# message("Fitting with Laplace...")
# fit_laplace <- model$estimate(
#     data = df,
#     method = "laplace",
#     ode.solver = "rk4",
#     ode.timestep = dt
#     )
# }   
par_names <- rownames(model$getParameters(type = "free"))
n_par <- length(par_names)

# fit to each dataset with each method
methods = c("ekf", "laplace", "laplace.thygesen")

# Test over a range of noise lvls:
#sigma_vals <- c(0.01, 0.05, 0.1, 0.2, 0.4)
sigma_vals <- c(0.05, 0.1, 0.2)
fit_sigma <- vector("list", length(sigma_vals))
names(fit_sigma) <- paste0("sigma_", sigma_vals)
comp_time <- matrix(NA, nrow=length(sigma_vals), ncol=length(methods))
rownames(comp_time) <- paste0("sigma_", sigma_vals)
colnames(comp_time) <- methods

# Number of clusters for parallel processing
n_clusters <- 4
sofun()

for (sigma_val in sigma_vals) {

    message("Fitting for sigma = ", sigma_val, sep="")
    p0_k <- p0
    p0_k$log_sN <- log(sigma_val)

    sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0, p0_k, dt, N)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_sigma[[paste0("sigma_", sigma_val)]] <- fit_model_par(n_clusters, methods, model, 
                                            t, Ysim, iobs, dt, N, df_fun = NULL)
}

# Save fits with date and time in file name
#filename <- paste0("results/RmA_fits_sigma_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
#saveRDS(fit_sigma, file=filename)

#fit_sigma <- readRDS("results/RmA_fits_sigma_20260113_092925.rds")

fit_sigma[[1]]$ekf[[1]]



# Summarize fits
source("src/compute_coverage.R")
source("src/compute_fit_summary.R")
metric_names <- c("coverage", "bias", "rmse")

fit_summary_sigma <- vector("list", length(sigma_vals))
names(fit_summary_sigma) <- names(fit_sigma)

for (sigma_val in sigma_vals) {
    fit <- fit_sigma[[paste0("sigma_", sigma_val)]]
    message("Computing summary for sigma = ", sigma_val, sep="")
    # True parameters
    p0_k <- p0
    p0_k$log_sN <- log(sigma_val)
    fit_summary <- compute_fit_summary(fit, p0_k, methods, n_par, par_names, metric_names)
    
    fit_summary_sigma[[paste0("sigma_", sigma_val)]] <- fit_summary
}
fit_summary_sigma

# Convergence summary
compute_convergence <- function(fit_sigma, sigma_vals, methods, N) {
    convergence <- matrix(NA, nrow = length(sigma_vals), ncol = length(methods))
    rownames(convergence) <- paste0("sigma_", sigma_vals)
    colnames(convergence) <- methods

    for (i in 1:length(sigma_vals)) {
        sigma_val <- sigma_vals[i]
        fit_i <- fit_sigma[[i]]
        for (m in methods) {
            convergence[i,m] <- sum(sapply(fit_i[[m]], function(fit) {
                !is.null(fit$fit)
            }))/N
        }
        
    }
    return(convergence)
}

convergence <- compute_convergence(fit_sigma, sigma_vals, methods, N)
convergence

# Plot metric vs noise level for each method for all parameters
ylims <- list(c(0,1), c(-0.2,0.2), c(0,0.5))
names(ylims) <- metric_names

# Plot metrics vs noise level
plot_metric("coverage", fit_summary_sigma, sigma_vals, "sigma", ylims, methods, par_names)
plot_metric("bias", fit_summary_sigma, sigma_vals, "sigma", ylims, methods, par_names)
plot_metric("rmse", fit_summary_sigma, sigma_vals, "sigma", ylims, methods, par_names)

fit_sigma[[1]]$ekf[[3]]$fit

# Computation time summary
for (i in 1:length(sigma_vals)) {
    sigma_val <- sigma_vals[i]
    fit_i <- fit_sigma[[paste0("sigma_", sigma_val)]]
    for (m in methods) {
        times <- sapply(fit_i[[m]], function(fit) {
            fit$time
        })
        comp_time[i,m] <- mean(times, na.rm=TRUE)
    }
}
comp_time

# For a given noise level, show parameter estimate distributions
sigma_val <- sigma_vals[4]
p0_k <- p0
p0_k$log_sN <- log(sigma_val)
fit_summary <- fit_summary_sigma[[paste0("sigma_", sigma_val)]]

# Plot distribution of parameter estimates for every method for every parameter (n_par x n_methods plots)
n_par_plot <- c(1,2,3,4,5)
# Define xlims for parameters in histograms
xlims <- list(
    log_r = c(-2,2),
    log_K = c(-2,2),
    log_beta = c(-2,6),
    log_mu = c(-2,3),
    log_sN = c(-3,0)
)
par(mfrow=c(length(n_par_plot), length(methods)))
for (par in par_names[n_par_plot]) {
    for (m in methods) {
        hist(fit_summary[[m]]$est[,par],
             main = paste("Estimates of ", par, " (", m, ")", sep=""),
             xlab = par,
             xlim = xlims[[par]])
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

#### Test different lengths of time series ####
N <- 100
dt <- 0.1
x0_bar = c(0.5, 0.5)
P0 <- diag(c(0.01,0.01))
x0 <- log(matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar)))

methods <- c("ekf", "laplace", "laplace.thygesen")

# Test different lengths of time series
T_vals <- c(20, 50, 100, 200)
fit_T <- vector("list", length(T_vals))
names(fit_T) <- paste0("T_", T_vals)
for (T_val in T_vals) {

    message("Fitting for T = ", T_val, sep="")
    t <- seq(0,T_val,dt)

    sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0, p0, dt, N)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_T[[paste0("T_", T_val)]] <- fit_model(methods, model, t, Ysim, iobs, dt, N)
}

# Save fits with date and time in file name
filename_T <- paste0("results/RmA_fits_T_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
saveRDS(fit_T, file=filename_T)

# Summarize fits
fit_summary_T <- vector("list", length(T_vals))
names(fit_summary_T) <- names(fit_T)

for (T_val in T_vals) {
    fit <- fit_T[[paste0("T_", T_val)]]

    fit_summary <- compute_fit_summary(fit, p0, methods, n_par, par_names, metric_names)
    
    fit_summary_T[[paste0("T_", T_val)]] <- fit_summary
}
fit_summary_T

fit_T
 
ylims <- list(c(0,1), c(-0.2,0.2), c(0,0.5))
names(ylims) <- metric_names
plot_metric("coverage", fit_summary_T, T_vals, "T", ylims, methods, par_names)
plot_metric("bias", fit_summary_T, T_vals, "T", ylims, methods, par_names)
plot_metric("rmse", fit_summary_T, T_vals, "T", ylims, methods, par_names)


# For a given time length, show parameter estimate distributions
T_val <- T_vals[1]
fit_summary <- fit_summary_T[[paste0("T_", T_val)]]

# Plot distribution of parameter estimates for every method for every parameter (n_par x n_methods plots)
n_par_plot <- c(1,2,3,4,5)
# Define xlims for parameters in histograms
xlims <- list(
    log_r = c(-2,2),
    log_K = c(-2,2),
    log_beta = c(-2,6),
    log_mu = c(-2,3),
    log_sN = c(-3,0)
)
par(mfrow=c(length(n_par_plot), length(methods)))
for (par in par_names[n_par_plot]) {
    for (m in methods) {
        hist(fit_summary[[m]]$est[,par],
             main = paste("Estimates of ", par, " (", m, ")", sep=""),
             xlab = par,
             xlim = xlims[[par]])
        abline(v = p0[[par]], col="red", lwd=2)
    }
}
