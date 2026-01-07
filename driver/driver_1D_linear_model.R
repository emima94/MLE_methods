## Driver 2D model

require(SDEtools)
require(ctsmTMB)
require(RTMB)

set.seed(123456)

# Source all functions in src/
sofun <- function() {
    src_files <- list.files("src/OU",pattern="*.R",full.names=TRUE)
    for (f in src_files) {
        source(f)
    }
}
sofun()

source("src/plots.R")

#### Generate data ####

# Model defintion of the Ornstein-Uhlenbeck (O-U) process
fsim <- function(x,par=p0) with(par, c(theta*(mu - x)))
gsim <- function(x,par=p0) with(par, sigma)

T <- 100
dt <- 0.1
t <- seq(0,T,dt)

# Number of datasets to simulate
N <- 20

P0 <- 0.1
x0 <- rnorm(N, 0.3, P0)

# True parameters
p0 <- list(
    theta = 0.2,
    mu = 0.5,
    sigma = 0.1,
    obs_sd = 0.02
)

# Gaussian observations
n <- length(t)
tsample <- 1 # 5 
iobs <- round(seq(1,n,tsample/dt))


## Simulate using Euler-Maruyama
Xsim <- matrix(0, nrow=length(t), ncol=N)
Xsim
Ysim <- matrix(0, nrow=length(iobs), ncol=N)
Ysim

for (i in 1:N) {
    B <- rBM(t)

    sim <- euler(function(x)fsim(x,p0),
                function(x)gsim(x,p0),
                t,x0[i],B=B)
    Xsim[,i] = sim$X

    Ysim[,i] = Xsim[iobs,i] + rnorm(length(iobs),0,p0$obs_sd) 
}
matplot(t, Xsim,
        type = "l",
        lty = 1,
        col = rgb(0, 0, 0, 0.2),
        ylim = c(0, 1),
        xlab = "Time",
        ylab = "X")


#### Fit models ####
# Get ctsmTMB model
sofun()
model <- create_OU_model(obs_sd = p0$obs_sd, x0 = x0, p0 = P0)

# fit to each dataset

fit_EKF <- vector("list", N)

for (i in 1:N) {
    print(paste("Fitting dataset ", i, " of ", N, sep=""))

    df <- data.frame(
        t = t[iobs],
        Y = Ysim[,i]
    )

    fit_EKF[[i]] <- model$estimate(
        data = df,
        method = "ekf",
        ode.solver = "rk4",
        ode.timestep = dt,
        silent = TRUE
    )

}

# Plot distribution of parameter estimates
theta_est <- sapply(fit_EKF, function(fit) fit$par.fixed["theta"])
mu_est <- sapply(fit_EKF, function(fit) fit$par.fixed["mu"])
sigma_est <- sapply(fit_EKF, function(fit) fit$par.fixed["sigma"])

par(mfrow=c(3,1))
hist(theta_est, main=expression("Estimates of " ~ theta), xlab=expression(theta))
abline(v=p0$theta, col="red", lwd=2)
hist(mu_est, main=expression("Estimates of " ~ mu), xlab=expression(mu))
abline(v=p0$mu, col="red", lwd=2)
hist(sigma_est, main=expression("Estimates of " ~ sigma), xlab=expression(sigma))
abline(v=p0$sigma, col="red", lwd=2)


fit_EKF$nll
fit_EKF

fit_laplace_X <- model$estimate(
    data = df,
    method = "laplace",
    ode.timestep = dt
)

fit_laplace_X$nll
fit_laplace_X
fit_laplace_XdB <- model$estimate(
    data = df,
    method = "laplace2",
    ode.timestep = dt
)
fit_laplace_XdB$nll
fit_laplace_XdB



