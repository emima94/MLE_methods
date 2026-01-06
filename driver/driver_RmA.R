## Driver 2D model

require(SDEtools)
require(ctsmTMB)
require(RTMB)

set.seed(123456)

# Source all functions in src/
sofun <- function() {
    src_files <- list.files("src",pattern="*.R",full.names=TRUE)
    for (f in src_files) {
        source(f)
    }
}
sofun()
#### Generate data ####
## Simulation control
dt <- 0.1
T <- 200

## System parameters
p0 <- list(
    r = 1,
    K = 1,
    epsilon = 3,
    beta = 3,
    Cmax = 1,
    mu = 1,
    sN = 0.2, #0.1
    sP = 0.1  #0.1
)

obs_sd <- 0.1
sim <- RmA_sim(p0, dt, T, obs_sd)

# Plot sim data
plot_RmA_sim_obs(sim$t, sim$Nsim, sim$Psim,
                 sim$iobs, sim$YN_gauss)

#### Fit models ####
# Fit using Laplace approximation with tiny approach
#fit_laplace_tiny <- fit_laplace_tiny_RmA(p0, sim$Xsim, sim$Ysim,
#                                         sim$BN, sim$BP,
#                                         sim$YN_gauss, sim$iobs, dt, obs_sd)

#fit_laplace_tiny$rep

# Fit using EKF approach
# Get ctsmTMB model
sofun()
ctsm_model <- create_RmA_model_ctsmTMB()

df <- data.frame(
    t = sim$t[sim$iobs],
    YN_gauss = sim$YN_gauss
)

fit_EKF <- ctsm_model$estimate(
    data = df,
    method = "ekf",
    ode.solver = "rk4",
    ode.timestep = dt
)

fit_EKF

fit_laplace_X <- ctsm_model$estimate(
    data = df,
    method = "laplace",
    ode.timestep = dt
)

fit_laplace_XdB <- ctsm_model$estimate(
    data = df,
    method = "laplace2",
    ode.timestep = dt
)

fit_EKF
fit_laplace_X    
fit_laplace_XdB
