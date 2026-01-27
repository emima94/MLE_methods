# Driver CIR model

require(SDEtools)
require(ctsmTMB)


# Source all functions in src/
sofun <- function() {
    src_files <- list.files("src/CIR",pattern="*.R",full.names=TRUE)
    for (f in src_files) {
        source(f)
    }
    source("src/plots.R")
    source("src/generate_data.R")
    source("src/fit_model.R")
    source("src/compute_coverage.R")
    source("src/compute_fit_summary.R")
    source("src/eval_results.R")
}
sofun()

# Define CIR process parameters (true)
par_true <- list(
    kappa = 0.5,
    theta = 2.0,
    sigma = 0.3,
    obs_sd = 1e-6
)

# Initial state
x0_bar = 2.0
P0 <- matrix(0.01)

# Number of datasets to simulate
N <- 50
x0 <- matrix(rnorm(N, x0_bar, sqrt(P0)), nrow=N, ncol=1)

fsim <- function(x, p) p$kappa * (p$theta - x)
gsim <- function(x, p) p$sigma * sqrt(x)
hsim <- function(x, p) x

T <- 200
dt <- 0.1
t <- seq(0,T,dt)

sim <- generate_data(fsim, gsim, hsim, t, x0, par_true, dt, N) 
sim

# Plot simulations
matplot(t, sapply(sim$Xsim, function(x) x[,1]), type="l", lty=1, col=rainbow(N), ylab="X", xlab="Time", main="CIR process simulations")
matpoints(t[sim$iobs], sapply(sim$Ysim, function(y) y[,1]), pch=16, col=rainbow(N), cex=0.7)

# Estimate parameters using true likelihood function
par0 <- list(
    kappa = 0.1, 
    theta = 1.0,
    sigma = 0.1
)

# Initial parameter vector
par_init <- unlist(par0[c("kappa", "theta", "sigma")])
par_init

sofun()

fit_true <- list()
for (i in 1:N) {
    message(sprintf("Fitting dataset %d/%d", i, N))
    fit_true[[i]] <- nlminb(
        start = par_init,
        objective = function(par) -loglik_cir(par, sim$Xsim[[i]], dt),
        lower = c(1e-6, 1e-6, 1e-6),
        upper = c(10, 10, 10)   
    )
}
# Parameter estimation distribution
par(mfrow=c(1,3))
for (param_name in names(par0)) {
    par_vals <- sapply(fit_true, function(fit) fit$par[param_name])
    hist(par_vals, main=paste("Estimated", param_name), xlab=param_name)
    abline(v = par_true[[param_name]], col="red", lwd=2)
}

# Likelihood distribution
nll_vals <- sapply(fit_true, function(fit) fit$objective)
hist(nll_vals, main="Negative log-likelihood values", xlab="NLL")
# True negative log-likelihoods
nll_true <- sapply(sim$Xsim, function(x) -loglik_cir(unlist(par_true[c("kappa", "theta", "sigma")]), x, dt))
nll_true

nll_residuals <- nll_vals - nll_true
hist(nll_residuals, main="NLL residuals (est - true)", xlab="NLL residual")

# Define ctsmTMB model for CIR (Lamperti transformed)
model_ctsmTMB <- ctsmTMB$new()

model_ctsmTMB$addSystem(
    dY1 ~ (-kappa / 2 * Y1 + (nu - 1) / (2 * Y1)) * dt + dw
)

# Inverse Lamperti transform
model_ctsmTMB$addObs(
    X ~ (sigma * Y1 / 2)^2
)

model_ctsmTMB$setVariance(
    X ~ obs_sd^2
)

model_ctsmTMB$setAlgebraics(
    nu ~ 4 * kappa * theta / (sigma^2)
)

model_ctsmTMB$setParameter(
    kappa = c(init = 0.1, lower = 1e-6, upper = 10),
    theta = c(init = 1.0, lower = 1e-6, upper = 10),
    sigma = c(init = 0.1, lower = 1e-6, upper = 10),
    obs_sd = p0$obs_sd
)

# Estimate parameters using ctsmTMB model with EKF and laplace
fit_ctsmTMB_laplace <- list()
for (i in 1:N) {
    message(sprintf("Fitting dataset %d/%d with ctsmTMB", i, N))
    df <- data.frame(t = t[sim$iobs], X = sim$Ysim[[i]])
    Y0 = 2/p0$sigma * sqrt(x0[i])
    model_ctsmTMB$setInitialState(list(Y0, P0))
    fit_ctsmTMB_laplace[[i]] <- model_ctsmTMB$estimate(
        df,
        ode.solver = "rk4",
        ode.timestep = dt,
        method = "laplace",
        silent = TRUE
    )
}

nll_laplace <- sapply(fit_ctsmTMB_laplace, function(fit) fit$nll)
hist(nll_laplace, main="ctsmTMB Laplace NLL", xlab="NLL")
# Parameter estimation distribution
par(mfrow=c(1,3))
for (param_name in names(par0)) {
    par_vals <- sapply(fit_ctsmTMB_laplace, function(fit) fit$par.fixed[param_name])
    hist(par_vals, main=paste("Estimated", param_name), xlab=param_name)
    abline(v = par_true[[param_name]], col="red", lwd=2)
}


nll_residuals_laplace <- nll_laplace - nll_true
hist(nll_residuals_laplace, main="ctsmTMB Laplace NLL residuals (est - true)", xlab="NLL residual")



# Estimate parameters using ctsmTMB model with EKF and laplace
fit_ctsmTMB_ekf <- list()
for (i in 1:N) {
    message(sprintf("Fitting dataset %d/%d with ctsmTMB", i, N))
    df <- data.frame(t = t[sim$iobs], X = sim$Ysim[[i]])
    Y0 = 2/p0$sigma * sqrt(x0[i])
    model_ctsmTMB$setInitialState(list(Y0, P0))
    fit_ctsmTMB_ekf[[i]] <- model_ctsmTMB$estimate(
        df,
        ode.solver = "rk4",
        ode.timestep = dt,
        method = "ekf",
        silent = TRUE
    )
}

# Plot chi-squared distribution with 3 degrees of freedom

# # Lamperti transformed system
# fbarsim <- function(z, p) with(p, {

#     nu <- 4 * kappa * theta / (sigma^2)
#     dz <- -kappa / 2 * z + (nu - 1) / (2 * z)
#     return(dz)
# })
# gbarsim <- function(z, p) 1

# z0_bar <- 2 * sqrt(x0_bar) / p0$sigma
# sim_bars <- euler(function(z)fbarsim(z,p0),
#                 function(z)gbarsim(z,p0),
#                 t,z0_bar, B=B)

# x_bars <- (p0$sigma * sim_bars$X / 2)^2

# plot(t, sim$X, type="l")
# points(t, x_bars, col="red")

# sofun()

# # Compute CIR likelihood
# p1 <- list(
#     kappa = 0.5,
#     theta = 1.0,
#     sigma = 0.3
# )
# nll0 <- -loglik_cir(p0, sim$X, dt)
# nll1 <- -loglik_cir(p1, sim$X, dt)
# nll0
# nll1
