
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
N <- 10

x0_bar = log(c(0.5, 0.5))
P0 <- diag(c(0.001,0.001))
x0 <- matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar))

model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)
model$setParameter(
    log_r = p0$log_r#,
#    log_mu = p0$log_mu,
    #log_K = p0$log_K
)

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
out_dir_base <- "results/RmA_sweep_T"
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
T_vals <- c(100, 200, 400, 800, 1600)

# Save parameters and settings
settings <- list(
    p0 = p0,
    n_par = n_par,
    par_names = par_names,
    methods = methods,
    N = N,
    x0_bar = x0_bar,
    x0 = x0,
    P0 = P0,
    dt = dt,
    tsample = tsample,
    T_vals = T_vals,
    var_name = var_name
)

saveRDS(settings, file = file.path(out_dir_base, paste0(var_name, "_settings.rds")))

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

    sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0, p0, dt, N, tsample)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_model_par(n_clusters, out_dir, methods, model, 
                                            t, Ysim, iobs, dt, N, df_fun = NULL)
}


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

