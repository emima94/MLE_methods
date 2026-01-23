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
    source("src/eval_results.R")
}
sofun()

# Set font in plots to Times New Roman 

#### Generate data ####

## RmA/Bazykin-model in natural coordinates
T <- 100
dt <- 0.1
t <- seq(0,T,dt)

# Number of datasets to simulate
N <- 100

x0_bar = c(0.5, 0.5)
P0 <- diag(c(0.01,0.01))
x0 <- log(matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar)))
N0 <- exp(x0)
N0

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
p0$log_mu <- log(4)
sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0, p0, dt, N)
Xsim <- sim_data$Xsim
Ysim <- sim_data$Ysim
iobs <- sim_data$iobs
Nsim <- sapply(Xsim, function(x) exp(x[,1]), simplify=FALSE)
plot(t, exp(Xsim[[1]][,1]), type="l")


matplot(t, do.call(cbind, Nsim),
        type = "l",
        lty = 1,
        col = rgb(0, 0, 0, 0.2),
        xlab = "Time",
        ylab = "N"
        )
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

par_names <- rownames(model$getParameters(type = "free"))
n_par <- length(par_names)

# fit to each dataset with each method
methods = c("ekf", "laplace", "laplace.thygesen")

# Test over a range of noise lvls:
#sigma_vals <- c(0.01, 0.05, 0.1, 0.2, 0.4)
var_name <- par_names[5]
log_sN_vals <- log(c(0.01, 0.05, 0.1, 0.2, 0.4))
fit_sigma <- vector("list", length(log_sN_vals))
names(fit_sigma) <- paste0(var_name, "_", 1:length(log_sN_vals))

# Number of clusters for parallel processing
n_clusters <- 4

sofun()

for (i in 1:length(log_sN_vals)) {
    log_sN_val <- log_sN_vals[i]

    message("Fitting for sigma = ", exp(log_sN_val), sep="")
    p0_k <- p0
    p0_k$log_sN <- log_sN_val

    sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0, p0_k, dt, N)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_sigma[[i]] <- fit_model_par(n_clusters, methods, model, 
                                            t, Ysim, iobs, dt, N, df_fun = NULL)
}

# Save fits with date and time in file name
#filename <- paste0("results/RmA_fits_sigma_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
#saveRDS(fit_sigma, file=filename)

fit_sigma <- readRDS("results/RmA_fits_sigma_20260121_100920.rds")

#fit_sigma





# Summarize fits
metric_names <- c("coverage", "bias", "rmse")

fit_summary_sigma <- vector("list", length(log_sN_vals))
names(fit_summary_sigma) <- names(fit_sigma)

for (i in 1:length(log_sN_vals)) {
    log_sN_val <- log_sN_vals[i]
    fit <- fit_sigma[[i]]
    message("Computing summary for sigma = ", exp(log_sN_val), sep="")
    # True parameters
    p0_k <- p0
    p0_k$log_sN <- log_sN_val
    fit_summary <- compute_fit_summary(fit, p0_k, methods, n_par, par_names, metric_names)
    
    fit_summary_sigma[[i]] <- fit_summary
}
fit_summary_sigma

#### Evaluation plots ####
sofun()
# Computation time vs sigma_vals #
null_count <- create_computation_time_plots(fit_sigma, log_sN_vals, methods, "log_sN", "RmA")
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
create_parameter_estimation_plots(fit_summary_sigma, log_sN_vals, methods, 
                        "log_sN", "RmA", par_names, p0, xlims, var_name_title =  "sigma_N", var_title_func = function(x) exp(x)
                                 )


# Coverage and RMSE plots # 
create_perf_metric_plot(fit_summary_sigma, log_sN_vals, methods, 
                        "log_sN", "RmA", par_names, "coverage"
                                 )

create_perf_metric_plot(fit_summary_sigma, log_sN_vals, methods, 
                        "log_sN", "RmA", par_names, "rmse"
                                 )

# Likelihood distribution plots
create_likelihood_histograms(fit_sigma, log_sN_vals, methods, 
                        "log_sN", "RmA", "sigma_N", var_title_func = function(x) exp(x)
                                 )

#### ---------------------------------- ####
#### Sweep lengths of time series, T    ####
#### ---------------------------------- ####


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
x0_bar = c(0.5, 0.5)
P0 <- diag(c(0.01,0.01))
x0 <- log(matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar)))

model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)
model$setParameter(
    log_sP = log(c(init = 0.05, lower = 0, upper = 10))
)

par_names <- rownames(model$getParameters(type = "free"))
n_par <- length(par_names)

# fit to each dataset with each method
methods = c("ekf", "laplace", "laplace.thygesen")

# Number of datasets to simulate
N <- 20

# Observation frequency
tsample = 0.1





var_name <- "T"
T_vals <- c(400, 800)
#T_vals <- c(100, 200)
fit_T <- vector("list", length(T_vals))
names(fit_T) <- paste0(var_name, "_", 1:length(T_vals))
n_clusters <- 1
for (i in 1:length(T_vals)) {
    T <- T_vals[i]

    message("Fitting for T = ", T, sep="")
    t <- seq(0,T,dt)

    sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0, p0, dt, N, tsample)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    fit_T[[i]] <- fit_model_par(n_clusters, methods, model, 
                                            t, Ysim, iobs, dt, N, df_fun = NULL)
}

# Summarize fits
metric_names <- c("coverage", "bias", "rmse")
fit_summary_T <- vector("list", length(T_vals))
names(fit_summary_T) <- names(fit_T)

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


# Make a likelihood profile in 2D for a single dataset
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
   
pg_name <- c("r", "beta")
pg <- c(paste0("log_", pg_name[1]), paste0("log_", pg_name[2]))

n_grid = 10
sofun()
delta_list <- list(
    delta_plus_list = delta_plus_list,
    delta_minus_list = delta_minus_list
)
out_2D_nll <- likelihood_profile_2D_RmA(p0, T, pg_name, n_grid, delta_list)
out_2D_nll$nll_matrix
replayPlot(out_2D_nll$plt_edge)
replayPlot(out_2D_nll$plt_contour)

# Now, try to estimate also log_sN and log_beta while profiling log_r and log_mu
n_grid = 10
dt <- 0.1
sofun()
pg_name <- c("mu", "beta")
pg <- c(paste0("log_", pg_name[1]), paste0("log_", pg_name[2]))

var1_vals <- seq(exp(p0[[pg[1]]]) - delta_list$delta_minus_list[[pg[1]]], exp(p0[[pg[1]]]) + delta_list$delta_plus_list[[pg[1]]], length.out=n_grid)
var2_vals <- seq(exp(p0[[pg[2]]]) - delta_list$delta_minus_list[[pg[2]]], exp(p0[[pg[2]]]) + delta_list$delta_plus_list[[pg[2]]], length.out=n_grid)
sofun()
nll_matrix_full <- matrix(NA, nrow=n_grid, ncol=n_grid)
nll_max_grad_matrix <- array(NA, dim = c(n_grid, n_grid))
k <- 0
par_est <- list()
for (i in 1:length(var1_vals)) {
    for (j in 1:length(var2_vals)) {
        k <- k + 1
        var1_val <- var1_vals[i]
        var2_val <- var2_vals[j]
        
        model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)

        params <- list(
                log_K = p0$log_K)

        params[[pg[1]]] <- log(var1_val)
        params[[pg[2]]] <- log(var2_val)

        do.call(model$setParameter, params)

        fit <- model$estimate(df, method = "laplace",
                    ode.timestep = dt,
                    ode.solver = "rk4",
                    silent = TRUE,
                    control = list(trace = 0))

        nll_value <- fit$nll
        nll__max_grad <- max(fit$nll.gradient)
        par_est[[k]] <- fit$par.fixed


        

        nll_matrix_full[i,j] <- nll_value
        nll_max_grad_matrix[i,j] <- nll__max_grad
        par_est[[k]] <- fit$par.fixed
        message(sprintf("%s=%.2f, %s=%.2f, NLL=%.2f", pg_name[1], var1_val, pg_name[2], var2_val, nll_value))
    }
}   

# Plot likelihood surface
i <- which(nll_matrix_full == min(nll_matrix_full), arr.ind = TRUE)
# clip extreme values for better plotting
nll_matrix_full_clip <- nll_matrix_full
threshold <- min(nll_matrix_full) + 400
nll_matrix_full_clip[nll_matrix_full > threshold] <- threshold

pdf(paste0("figures/RmA_likelihood_profile_full_2D_", pg_name[1], "_", pg_name[2], ".pdf"), width = 8, height = 6)
filled.contour(
  var1_vals, var2_vals, nll_matrix_full_clip,
  xlab = pg_name[1],
  ylab = pg_name[2],
  main = paste0("Likelihood surface, in the ", pg_name[2], " vs. ", pg_name[1], " plane"),

  # ---- Color bar label ----
  key.title = title(main = "NLL"),

  plot.axes = {
    axis(1); axis(2)

    # True value
    points(exp(p0[[pg[1]]]), exp(p0[[pg[2]]]),
           pch = 19, cex = 1.5, col = "green")

    # Optimal (minimum NLL)
    points(var1_vals[i[1]], var2_vals[i[2]],
           pch = 1, col = "black", cex = 2, lwd = 3)

    # ---- Legend ----
    legend("topright",
           legend = c("True value", "Optimal (MLE)"),
           pch = c(19, 1),
           col = c("green", "black"),
           pt.cex = c(1.2, 1.5),
           lwd = c(NA, NA),
           bg = "white")
  }
)
dev.off()





filled.contour(var1_vals, var2_vals, log10(abs(nll_max_grad_matrix)),
xlab = pg_name[1],
ylab = pg_name[2],
main = "Max gradient surface (log10)",
plot.axes = {
    axis(1); axis(2)
    points(exp(p0[[pg[1]]]), exp(p0[[pg[2]]]),
        pch = 19, cex = 1.5, col = "green")
    points(var1_vals[i[1]], var2_vals[i[2]], pch = 1, col = "black", cex = 2, lwd = 3)       
}
)
log_mu_est <- sapply(par_est, function(p) p["log_beta"])

filled.contour(var1_vals, var2_vals, matrix(exp(log_mu_est), nrow=n_grid, ncol=n_grid, byrow = TRUE),
xlab = pg_name[1],
ylab = pg_name[2],
main = "Estimated log_mu surface",
plot.axes = {
    axis(1); axis(2)
    points(exp(p0[[pg[1]]]), exp(p0[[pg[2]]]),
        pch = 19, cex = 1.5, col = "green")
    points(var1_vals[i[1]], var2_vals[i[2]], pch = 1, col = "black", cex = 2, lwd = 3)       
}
)

log_sN_est <- sapply(par_est, function(p) p["log_sN"])
filled.contour(var1_vals, var2_vals, matrix(exp(log_sN_est), nrow=n_grid, ncol=n_grid, byrow = TRUE),
xlab = pg_name[1],
ylab = pg_name[2],
main = "Estimated log_sN surface",
plot.axes = {
    axis(1); axis(2)
    points(exp(p0[[pg[1]]]), exp(p0[[pg[2]]]),
        pch = 19, cex = 1.5, col = "green")
    points(var1_vals[i[1]], var2_vals[i[2]], pch = 1, col = "black", cex = 2, lwd = 3)       
}
)
