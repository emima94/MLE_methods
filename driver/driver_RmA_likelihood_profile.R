# Likelihood profile for RmA model - 2D contour plot
# Author:
# Emil Skov Martinsen
# Date: January 2026


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
    source("src/eval_results.R")
}
sofun()

# Model parameters (true)
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

# Define the ranges for the likelihood profiling
delta_plus_list <- list(
        log_r = 0.6,
        log_mu = 0.4,
        log_beta = 1.0
    )

delta_minus_list <- list(
    log_r = 0.6,
    log_mu = 0.7,
    log_beta = 2.0
)
delta_list <- list(delta_plus_list = delta_plus_list, delta_minus_list = delta_minus_list)

# Construct the dataset
N <- 1
x0_bar = log(c(0.5, 0.5))
#P0 <- diag(c(0.01,0.01))
P0 <- diag(c(0.001,0.001))
x0 <- matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar))
T <- 200
dt <- 0.1
t <- seq(0,T,dt)
tsample <- 0.2
sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0[1,,drop=FALSE], p0, dt, N, tsample)
Ysim <- sim_data$Ysim
iobs <- sim_data$iobs
df <- data.frame(t = t[iobs], Y = Ysim[[1]])

# Variables to profile pair-wise
pg_name_list <- list(
    c("r", "mu"),
    c("r", "beta"),
    c("mu", "beta")
)

# Methods to profile
method_list <- c("ekf", "laplace")

# # Test ekf only (to make it work..)
# pg_name_list <- list(
#     c("r", "beta"),
#     c("mu", "beta")
# )
# method_list <- c("ekf")


for (method in method_list) {
    for (pg_name in pg_name_list) {
        cat(sprintf("Profiling %s using %s method\n", paste(pg_name, collapse=", "), method))
        profile_res <- likelihood_profile_estimate_2D_RmA(
            df = df,
            method = method,
            pg_name = pg_name,
            delta_list = delta_list,
            p0 = p0,
            x0_bar = x0_bar,
            P0 = P0,
            n_grid = 10,
            dt = dt,
            save_fig = TRUE,
            save_vars = TRUE
        )
    }
}

# files <- c("results/likelihood_profiles/RmA_likelihood_profile_full_2D_ekf_r_mu.rds",
#            "results/likelihood_profiles/RmA_likelihood_profile_full_2D_ekf_r_beta.rds",
#            "results/likelihood_profiles/RmA_likelihood_profile_full_2D_ekf_mu_beta.rds"
# )

# pg_names_list <- list(
#     c("r", "mu"),
#     c("r", "beta"),
#     c("mu", "beta")
# )
# for (i in 1:length(files)) {
#     out <- readRDS(files[i])
#     pg_name <- pg_names_list[[i]]
#     plot_likelihood_profile(
#         pg_name = pg_name,
#         var1_vals = out$var1_vals,
#         var2_vals = out$var2_vals,
#         nll_matrix_full = out$nll_matrix_full,
#         p0 = p0,
#         df = df,
#         dt = dt,
#         x0_bar = x0_bar,
#         P0 = P0,
#         method = "ekf",
#         save_fig=TRUE
#     )
# }