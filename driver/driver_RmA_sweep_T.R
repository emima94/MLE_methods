
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

sofun()
# Number of datasets to simulate
N <- 100

x0_bar = log(c(0.5, 0.5))
P0 <- diag(c(0.001,0.001))
x0 <- matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar))

model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)
#model$setParameter(
#    log_r = p0$log_r#,
#    log_mu = p0$log_mu,
    #log_K = p0$log_K
#)

par_names <- rownames(model$getParameters(type = "free"))
n_par <- length(par_names)
par_names

# fit to each dataset with each method
methods = c("ekf", "laplace", "laplace.thygesen")
#methods = c("ekf", "laplace")


# Observation frequency
tsample = 0.2

# Time step size
dt <- 0.1


var_name <- "T"
T_vals <- c(100, 200, 400, 800, 1600, 2400)
#T_vals <- c(100, 200)
fit_T <- vector("list", length(T_vals))
names(fit_T) <- paste0(var_name, "_", 1:length(T_vals))
n_clusters <- 8
for (i in 1:length(T_vals)) {
    T <- T_vals[i]

    message("Fitting for T = ", T, sep="")
    t <- seq(0,T,dt)

    sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0, p0, dt, N, tsample)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_T[[i]] <- fit_model(methods, model, 
                                            t, Ysim, iobs, dt, N, df_fun = NULL)
}

saveRDS(fit_T, file = "results/driver_RmA_sweep_T_fit_T.rds")

# # Sweep with model in N and P space
# sofun()

# p0_NP <- list(
#     r = 1.0,
#     K = 1.0,
#     epsilon = 3.0,
#     beta = 3.0,
#     Cmax = 1.0,
#     mu = 1.0,
#     sN = 0.2, #0.1
#     sP = 0.1,  #0.1
#     obs_sd = 0.05
# )
# fNPsim <- function(x,par=p0_NP) c(fN(x[1],x[2],par),fP(x[1],x[2],par))
# gNPsim <- function(x,par=p0_NP) diag(c(gN(x[1],par),gP(x[2],par)))
# hNPsim <- function(x,par=p0_NP) x[,1]  # observe N only


# N0_bar = c(0.5, 0.5)
# P0_N <- diag(c(0.001,0.001))
# N0 <- matrix(rnorm(N * length(N0_bar), N0_bar, diag(sqrt(P0_N))), nrow=N, ncol=length(N0_bar)))

# model <- create_RmA_NP_model(p = p0_NP, x0 = N0_bar, P0 = P0_N)
# fit_T_NP <- vector("list", length(T_vals))
# names(fit_T_NP) <- paste0(var_name, "_", 1:length(T_vals))
# for (i in 1:length(T_vals)) {
#     T <- T_vals[i]
#     message("Fitting for T = ", T, sep="")
#     t <- seq(0,T,dt)

#     sim_data <- generate_data_RmA(fNPsim, gNPsim, hNPsim, t, N0_bar, p0_NP, dt, N, tsample)
#     Ysim <- sim_data$Ysim
#     iobs <- sim_data$iobs
#     fit_T_NP[[i]] <- fit_model_par(n_clusters, methods, model, 
#                                             t, Ysim, iobs, dt, N, df_fun = NULL)
# }


# Summarize fits
metric_names <- c("coverage", "bias", "rmse")
fit_summary_T <- vector("list", length(T_vals))
names(fit_summary_T) <- names(fit_T)
sofun()
for (i in 1:length(T_vals)) {
    T <- T_vals[i]
    fit <- fit_T[[i]]
    message("Computing summary for T = ", T, sep="")
    fit_summary <- compute_fit_summary(fit, p0, methods, n_par, par_names, metric_names)
    
    fit_summary_T[[i]] <- fit_summary
}
fit_summary_T
sofun()
# Computation time vs T_vals #
null_count <- create_computation_time_plots(fit_T, T_vals, methods, "T", "RmA")
null_count  

# Parameter estimation distributions #
sofun()
xlims = list(
            log_r = c(-1,1),
            log_K = c(-1,1),
            log_beta = c(-2,2),
            log_mu = c(-2,2),
            log_sN = c(-16,2)
)
create_parameter_estimation_plots(fit_summary_T, T_vals, methods, 
                        "T", "RmA", par_names, p0, xlims, var_name_title =  "T", var_title_func = function(x) x
                                 )


# Coverage and RMSE plots # 
create_perf_metric_plot(fit_summary_T, T_vals, methods, 
                        "T", "RmA", par_names, "coverage"
                                 )

create_perf_metric_plot(fit_summary_T, T_vals, methods, 
                        "T", "RmA", par_names, "rmse"
                                 )

# Likelihood distribution plots
sofun()
create_likelihood_histograms(fit_T, T_vals, methods, 
                        "T", NA, "RmA", "T", var_title_func = function(x) x
                                 )
# True parameters
theta_true <- sapply(par_names, function(pn) p0[[pn]])
create_likelihood_histograms(fit_T, T_vals, methods, 
                        "T", theta_true, "RmA", "T", var_title_func = function(x) x,
                                 plot_dev = TRUE)
