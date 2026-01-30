
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

log_transform_state_vector <- c(TRUE, TRUE, FALSE, FALSE)
obs_partial_vector <- c(TRUE, FALSE, TRUE, FALSE)

# Load parameter index from HPC arguments
i <- as.integer(commandArgs(trailingOnly=TRUE)[1])

### Inputs ###
log_transform_state = log_transform_state_vector[i]
obs_partial = obs_partial_vector[i]
###



# Default settings
p0 <- list(
    r = 1.0,
    K = 1.0,
    epsilon = 3.0,
    beta = 3.0,
    Cmax = 1.0,
    mu = 1.0,
    sN = 0.2, #0.1
    sP = 0.1,  #0.1
    obs_sd = 0.005
)

# RmA def, no par log transforms, no state log transforms
fN_1 <- function(N,P,par) {with(par, {
    r * N * (1 - N / K) - beta * N * P / (1 + beta * N / Cmax)
})}
fP_1 <- function(N,P,par) {with(par, {
    epsilon * beta * N * P / (1 + beta * N / Cmax) - mu * P
})}
gN_1 <- function(N,par) {with(par, {
    sN * N
})}
gP_1 <- function(P,par) {with(par, {
    sP * P
})}

# Simulation functions
fNPsim <- function(x,par=p0) c(fN_1(x[1],x[2],par),fP_1(x[1],x[2],par))
gNPsim <- function(x,par=p0) diag(c(gN_1(x[1],par),gP_1(x[2],par)))
if (obs_partial) {
    hNPsim <- function(x,par=p0) x[,1]  # observe N only
} else {
    hNPsim <- function(x,par=p0) x      # observe both N and P
}

# Dataframe function for fitting
if (obs_partial) {
    df_fun <- NULL
} else {
    df_fun <- function(t, Ysim_i, iobs) {
        df <- data.frame(t = t[iobs], 
                        Y1 = Ysim_i[,1], 
                        Y2 = Ysim_i[,2])
        return(df)
    }
}

# Number of datasets to simulate
N <- 2

# Initial state for simulation and model
x0_model <- c(0.5, 0.5)
P0_model <- diag(c(1,1)*0.1)
x0_sim <- matrix(rep(c(0.5, 0.5), each=N), nrow=N, ncol=2)

log_transform_par = FALSE

model <- create_RmA_model_general(p = p0, x0 = x0_model, P0 = P0_model, log_transform_par = log_transform_par, log_transform_state = log_transform_state, obs_partial = obs_partial)

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

# Directory to save results
date_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir_base <- file.path("results", paste0("RmA_model_log_state_trans_", log_transform_state, "_partial_obs_", obs_partial, "_sweep_T_", date_time))
# Create output directory
if (!dir.exists(out_dir_base)) {
    dir.create(out_dir_base, recursive = TRUE)
}
# Create directory for heavy and light results
if (!dir.exists(file.path(out_dir_base, "light"))) {
    dir.create(file.path(out_dir_base, "light"), recursive = TRUE)
}
if (!dir.exists(file.path(out_dir_base, "heavy"))) {
    dir.create(file.path(out_dir_base, "heavy"), recursive = TRUE)
}

var_name <- "T"
#T_vals <- c(100, 200, 400, 800, 1600, 2400)
T_vals <- c(50, 400)

# Save parameters and settings
settings <- list(
    p0 = p0,
    n_par = n_par,
    par_names = par_names,
    methods = methods,
    model = model,
    log_transform_par = log_transform_par,
    log_transform_state = log_transform_state,
    obs_partial = obs_partial,
    N = N,
    x0_model = x0_model,
    x0_sim = x0_sim,
    P0 = P0_model,
    dt = dt,
    tsample = tsample,
    T_vals = T_vals,
    var_name = var_name
)

saveRDS(settings, file = file.path(out_dir_base, paste0(var_name, "_settings.rds")))

# t_test <- seq(0,100,dt)
# sim_data_test <- generate_data_RmA(fNPsim, gNPsim, hNPsim, t_test, x0_sim, p0, dt, N, tsample)
# Ysim_test <- sim_data_test$Ysim
# iobs_test <- sim_data_test$iobs
# par(mfrow=c(2,1))
# plot(t_test[iobs_test], Ysim_test[[1]][,1], type='o', col='blue')
# plot(t_test[iobs_test], Ysim_test[[1]][,2], type='o', col='red')

# df <- df_fun(t_test, Ysim_test, iobs_test)
# df
# fit <- model$estimate(
#     data = df,
#     method = "ekf",
#     ode.solver = "rk4",
#     ode.timestep = dt,
#     silent = FALSE,
#     control = list(trace = 1)
# )


#T_vals <- c(100, 200)
n_clusters <- 4
for (i in 1:length(T_vals)) {

    # Create output directory
    out_dir <- list(
        light = file.path(out_dir_base, "light", paste0(var_name, "_", T_vals[i])),
        heavy = file.path(out_dir_base, "heavy", paste0(var_name, "_", T_vals[i]))
    )
    if (!dir.exists(out_dir$light)) {
        dir.create(out_dir$light, recursive = TRUE)
    }
    if (!dir.exists(out_dir$heavy)) {
        dir.create(out_dir$heavy, recursive = TRUE)
    }

    T <- T_vals[i]

    message("Fitting for T = ", T, sep="")
    t <- seq(0,T,dt)

    sim_data <- generate_data_RmA(fNPsim, gNPsim, hNPsim, t, x0_sim, p0, dt, N, tsample)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_model_par(n_clusters, out_dir, methods, model, 
                                            t, Ysim, iobs, dt, N, df_fun = df_fun)
}

fit <- readRDS("results/RmA_model_log_state_trans_TRUE_partial_obs_FALSE_sweep_T_20260130_102557/light/T_400/laplace/0002.rds")
fit
#saveRDS(fit_T, file = "results/driver_RmA_sweep_T_fit_T.rds")

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

