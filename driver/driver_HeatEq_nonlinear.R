## Driver Heat equtaion model (1D diffusion equation)

require(SDEtools)
require(ctsmTMB)
require(RTMB)
library(parallel)

set.seed(123456)

#plan(multicore) # use multicore processing

# Source all functions in src/
sofun <- function() {
    src_files <- list.files("src/HeatEq",pattern="*.R",full.names=TRUE)
    for (f in src_files) {
        source(f)
    }
    source("src/plots.R")
    source("src/generate_data.R")
    source("src/fit_model.R")
    source("src/compute_coverage.R")
    source("src/compute_fit_summary.R")
    source("src/util.R")
}
sofun()

#### Generate data ####

## Heat equation model in 1D
T <- 100
dt <- 0.1
t <- seq(0,T,dt)

# Number of datasets to simulate
N <- 100

# True parameters
p0 <- list(
    log_sT = log(0.2),
    log_K0 = log(1),
    log_beta = log(0.6),
    T0 = 0.0,
    TL = 10.0,
    L = 10.0,
    obs_sd = 0.05,
    sensor_pos = c(1,3,7,9)
)

fsim_bar <- function(x, par=p0) {
    fsim(x, par, nonlinear=TRUE)
}
gsim_bar <- function(x, par=p0) {
    gsim(x, par, state_dependent=FALSE)
}

methods = c("ekf", "laplace", "laplace.thygesen")

metric_names <- c("coverage", "bias", "rmse")


# Test parallel par estimation

# Number of clusters for parallel processing
n_clusters <- 2


## Sweep over Nx
Nx_vec <- c(5, 10, 15, 20)
fit_Nx <- vector("list", length(Nx_vec))

for (Nx in Nx_vec) {
    message("Fitting for Nx = ", Nx)
    # Set initial conditions
    x0_bar = rep(5.0, Nx)
    P0 <- diag(rep(0.5, Nx))
    x0 <- matrix(runif(N * Nx, min=2, max=8), nrow=N, ncol=Nx, byrow=TRUE)
    # Simulate data
    p0_Nx <- p0
    sim_data <- generate_data(fsim_bar, gsim_bar, hsim, t, x0, p0, dt, N)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs

    # Create model
    model <- create_HeatEq_model(
        obs_sd = p0$obs_sd,
        x0 = x0_bar,
        P0 = P0,
        Nx = Nx,
        L = p0$L,
        T0 = p0$T0,
        TL = p0$TL,
        sensor_pos_idx = ceiling(p0$sensor_pos/p0$L * Nx),
        nonlinear = TRUE,
        state_dep_diffusion = FALSE
    )

    par_names <- rownames(model$getParameters(type = "free"))
    n_par <- length(par_names)

    fit_par <- fit_model_par(n_clusters, methods, model, 
        t, Ysim, iobs, dt, N, df_fun_HeatEq)   

    fit_Nx[[paste0("Nx_", Nx)]] <- fit_par   

}

# Save results
time_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
saveRDS(fit_Nx, file=paste0("results/HeatEq_nonlinear_fit_Nx_", time_str, ".rds"))



fit_Nx_summary <- list()
for (Nx in Nx_vec) {
    fit_par <- fit_Nx[[paste0("Nx_", Nx)]]
    fit_summary <- compute_fit_summary(fit_par, p0, methods, n_par, par_names, metric_names)
    fit_Nx_summary[[paste0("Nx_", Nx)]] <- fit_summary
}   

#### Computation time vs Nx ####
{
    est_time_mat <- matrix(0, nrow=length(Nx_vec), ncol=length(methods))
    rownames(est_time_mat) <- paste0("Nx_", Nx_vec)
    colnames(est_time_mat) <- methods
    compile_time_mat <- matrix(0, nrow=length(Nx_vec), ncol=length(methods))
    rownames(compile_time_mat) <- paste0("Nx_", Nx_vec)
    colnames(compile_time_mat) <- methods

    for (i in 1:length(Nx_vec)) {
        Nx <- Nx_vec[i]
        fit_par <- fit_Nx[[paste0("Nx_", Nx)]]
        for (m in methods) {
            est_time_mat[i,m] <- mean(sapply(fit_par[[m]], function(res) res$fit$private$timer_estimation))
            compile_time_mat[i,m] <- mean(sapply(fit_par[[m]], function(res) res$fit$private$timer_construct_adfun))
        }
    }

    pdf("figures/HeatEq_Nx_computation_time.pdf", width=8, height=3)
    par(mfrow=c(1,3))
    matplot(Nx_vec, est_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab="Number of spatial finite volumes (Nx)",
        ylab="Average estimation time (s)",
        main="Estimation time vs Nx")
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    matplot(Nx_vec, compile_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab="Number of spatial finite volumes (Nx)",
        ylab="Average compilation time (s)",
        main="Compilation time vs Nx")
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    matplot(Nx_vec, est_time_mat + compile_time_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab="Number of spatial finite volumes (Nx)",
        ylab="Average total time (s)",
        main="Total computation time vs Nx")
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    dev.off()
}

#### Iterations vs Nx ####
{
    iter_mat <- matrix(0, nrow=length(Nx_vec), ncol=length(methods))
    rownames(iter_mat) <- paste0("Nx_", Nx_vec)
    colnames(iter_mat) <- methods  
    for (i in 1:length(Nx_vec)) {
        Nx <- Nx_vec[i]
        fit_par <- fit_Nx[[paste0("Nx_", Nx)]]
        for (m in methods) {
            iter_mat[i,m] <- mean(sapply(fit_par[[m]], function(res) res$fit$private$opt$iterations))
        }
    }
    pdf("figures/HeatEq_Nx_iterations.pdf", width=6, height=4)
    par(mfrow=c(1,1))
    matplot(Nx_vec, iter_mat, type="b", pch=1, lty=1, col=1:length(methods),
        xlab="Number of spatial finite volumes (Nx)",
        ylab="Average number of outer iterations",
        main="Number of outer iterations vs Nx")
    legend("topleft", legend=methods, col=1:length(methods), pch=1, lty=1)
    dev.off()
}

#### Parameter estimation distributions ####
# Define xlims for parameters in histograms
xlims <- list(
    log_K0 = c(-2,3.5),
    log_beta = c(-2,1),
    log_sT = c(-2.5,-1)
)
for (Nx in Nx_vec) {
    fit_summary <- fit_Nx_summary[[paste0("Nx_", Nx)]]
    
    # Plot histograms of parameter estimates
    pdf(paste0("figures/HeatEq_parameter_estimates_Nx_", Nx, ".pdf"), width=12, height=8)
    par(mfrow=c(n_par, length(methods)))
    for (par in par_names) {
        for (m in methods) {
            if (all(is.nan(fit_summary[[m]]$est[,par]))) {
                next
            }
    
            hist(fit_summary[[m]]$est[,par],
                 main = paste("Estimates of ", par, " (", m, ", Nx=", Nx, ")", sep=""),
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
    dev.off()
}

#### Plot RMSE of parameter estimates vs Nx ####
for (par in par_names) {
    pdf(paste0("figures/HeatEq_RMSE_vs_Nx_", par, ".pdf"), width=6, height=4)
    par(mfrow=c(1,1))
    rmse_par <- sapply(Nx_vec, function(Nx) {
        fit_summary <- fit_Nx_summary[[paste0("Nx_", Nx)]]
        sapply(methods, function(m) fit_summary[[m]]$perf[par, "rmse"])
    })
    matplot(Nx_vec, t(rmse_par), type="b", pch=1, lty=1, col=1:length(methods),
        xlab="Number of spatial finite volumes (Nx)",
        ylab=paste("RMSE of ", par, " estimates", sep=""),
        main=paste("RMSE of ", par, " estimates vs Nx", sep=""))
    legend("topright", legend=methods, col=1:length(methods), pch=1, lty=1)
    dev.off()
}


#### Plot example data ####
{
    N <- 3
    Nx <- 10
    # Set initial conditions
    x0_bar = rep(5.0, Nx)
    P0 <- diag(rep(0.5, Nx))
    x0 <- matrix(runif(N * Nx, min=2, max=8), nrow=N, ncol=Nx, byrow=TRUE)
    # Simulate data
    sim_data <- generate_data(fsim_nonlinear, gsim_nonlinear, hsim, t, x0, p0, dt, N)
    Xsim <- sim_data$Xsim
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs

    par(mfrow = c(N,1))
    for (i in 1:N) {
        pdf(paste0("figures/HeatEq_example_data_dataset_", i, ".pdf"), width=12, height=8)
        matplot(t, Xsim[[i]] , type="l",
            xlab="Time", ylab="Temperature",
            main=paste("True states, dataset ", i, " Nx =", Nx, sep=""))
        matpoints(t[iobs], Ysim[[i]], pch=1, col="red")
        legend("topright", legend=c("True states", "Observations"), col=c("black", "red"), pch=c(NA,1), lty=c(1,NA))
        dev.off()
    }
}



# Plot diffusion coefficient K(T)
T_seq <- seq(0, 10, length.out=100)
beta <- exp(p0$log_beta)
K0 <- exp(p0$log_K0)
K_T <- K0 * exp(-beta * T_seq)
plot(T_seq, K_T, type="l", main="Diffusion Coefficient K(T)", xlab="Temperature T", ylab="K(T)")
