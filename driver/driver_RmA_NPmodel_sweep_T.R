
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

# RmA def, no par log transforms
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

fNPsim <- function(x,par=p0) c(fN_1(x[1],x[2],par),fP_1(x[1],x[2],par))
gNPsim <- function(x,par=p0) diag(c(gN_1(x[1],par),gP_1(x[2],par)))
hNPsim <- function(x,par=p0) x[,1]  # observe N only

sofun()
# Number of datasets to simulate
N <- 2

{
    set.seed(123456)
    #x0_model <- c(0.001, 0.001)
    x0_model <- c(0.5, 0.5)
    P0_model <- diag(c(1,1)*0.1)

    #x0_model <- c(0.5,0.5)
    #P0_model <- diag(c(0.01,0.01) * 0.01)

    #P0_sim <- diag(c(0.01,0.01) * 100)
    #x0_sim <- matrix(rep(x0_model, each=N), nrow=N, ncol=2)
    #x0_sim <- matrix(pmax(0,rnorm(N * length(x0_model), x0_model, diag(sqrt(P0_sim)))), nrow=N, ncol=length(x0_model))

    # TRUE initial states for simulation
    #x0_sim <- matrix(rep(c(0.0001, 0.0001), each=N), nrow=N, ncol=2)
    #x0_sim <- matrix(rep(c(0.5, 0.5), each=N), nrow=N, ncol=2)
    #x0_sim <- matrix(rep(x0_model, each=N), nrow=N, ncol=2)
    x0_sim <- matrix(rep(c(0.5, 0.5), each=N), nrow=N, ncol=2)
    log_transform_par = TRUE
    model <- create_RmA_NP_model(p = p0, x0 = x0_model, P0 = P0_model, log_transform_par = log_transform_par)
    if (log_transform_par) {
        #model$setParameter(log_r = log(p0$r))
        #model$setParameter(log_beta = log(p0$beta))
    } else {
        #model$setParameter(r = p0$r)
        #model$setParameter(beta = p0$beta)
    }
    par_names <- rownames(model$getParameters(type = "free"))
    n_par <- length(par_names)
    par_names
    #model$setParameter(
    #    r = p0$r#,
    #    log_mu = p0$log_mu,
        #log_K = p0$log_K
    #

    #### TESTING #### 
    sofun()
    # Simulate data
    T <- 1600
    dt <- 0.1
    t <- seq(0,T,dt)
    tsample <- 0.2

    sim_data <- generate_data_RmA(fNPsim, gNPsim, hNPsim, t, x0_sim, p0, dt, N, tsample)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs

    i <- 1

    plot(t, sim_data$Xsim[[i]][,1], type="l", main=sprintf("Dataset %d", i))
    points(t[iobs], Ysim[[i]], col="red")

    # Remove zeros
    #Y_df <- pmax(0, Ysim[[i]])
    Y_df <- Ysim[[i]]

    df <- data.frame(t = t[iobs], Y = Y_df)
    fit <- model$estimate(
        data = df,
        ode.solver = "rk4",
        ode.timestep = dt)
    fit_laplace <- model$estimate(
        data = df,
        method = "laplace",
        ode.solver = "rk4",
        ode.timestep = dt)
    fit_laplace.thygesen <- model$estimate(
        data = df,
        method = "laplace.thygesen",
        ode.solver = "rk4",
        ode.timestep = dt)

    fit
    fit_laplace
    fit_laplace.thygesen

    par_mat <- matrix(NA, nrow=4, ncol=length(rownames(model$getParameters(type = "free"))))
    colnames(par_mat) <- rownames(model$getParameters(type = "free"))
    rownames(par_mat) <- c("ekf", "laplace", "laplace.thygesen", "true")
    if (log_transform_par) {
        par_mat[1, ] <- exp(fit$par.fixed)
        par_mat[2, ] <- if (is.null(fit_laplace)) NA else exp(fit_laplace$par.fixed)
        par_mat[3, ] <- exp(fit_laplace.thygesen$par.fixed)
    } else {
        par_mat[1, ] <- fit$par.fixed
        par_mat[2, ] <- if (is.null(fit_laplace)) NA else fit_laplace$par.fixed
        par_mat[3, ] <- fit_laplace.thygesen$par.fixed
    }
    par_mat[4, ] <- unlist(p0[c("r","K", "beta","mu","sN")])


    # Make predictions with fitted model
    pred_ekf <- model$predict(
        data = df,
        par = fit$par.fixed,
        ode.solver = "rk4",
        ode.timestep = dt)
    pred_laplace <- model$predict(
        data = df,
        par = fit_laplace$par.fixed,
        ode.solver = "rk4",
        ode.timestep = dt)
    pred_laplace.thygesen <- model$predict(
        data = df,
        par = fit_laplace.thygesen$par.fixed,
        ode.solver = "rk4",
        ode.timestep = dt)

    plot(pred_ekf$observations$t.j, pred_ekf$observations$Y, type="l", col="blue", pch=16, ylim = c(0,2))
    lines(pred_laplace$observations$t.j, pred_laplace$observations$Y, col="green", pch=16)
    lines(pred_laplace.thygesen$observations$t.j, pred_laplace.thygesen$observations$Y, col="purple", pch=16)
    points(t[iobs], Ysim[[i]], col="red", pch=1)
    legend("topright", legend=c("EKF", "Laplace", "Laplace Thygesen", "Data"),
        col=c("blue", "green", "purple", "red"), pch=c(16,16,16,1))

    par_mat
}

#-----

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
out_dir_base <- file.path("results", paste0("RmA_model_NP_sweep_T_", date_time))
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

    sim_data <- generate_data_RmA(fNPsim, gNPsim, hNPsim, t, x0_sim, p0, dt, N, tsample)
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

