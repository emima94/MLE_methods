## CIR helper functions ##

cir_params <- function(x_prev, kappa, theta, sigma, dt) {
  c  <- sigma^2 * (1 - exp(-kappa * dt)) / (4 * kappa)
  nu <- 4 * kappa * theta / sigma^2
  lambda <- 4 * kappa * exp(-kappa * dt) * x_prev / (sigma^2 * (1 - exp(-kappa * dt)))
  list(c = c, nu = nu, lambda = lambda)
}

logdens_cir <- function(x, x_prev, kappa, theta, sigma, dt) {
  p <- cir_params(x_prev, kappa, theta, sigma, dt)
  c  <- p$c
  nu <- p$nu
  lam <- p$lambda

  z <- x / c
  logpdf <- dchisq(z, df = nu, ncp = lam, log = TRUE) - log(c)
  return(logpdf)
}

loglik_cir <- function(par, x, dt) {
  kappa <- par[1]
  theta <- par[2]
  sigma <- par[3]

  # enforce positivity
  if (kappa <= 0 || theta <= 0 || sigma <= 0) return(-Inf)

  ll <- 0
  for (i in 2:length(x)) {
    ll <- ll + logdens_cir(x[i], x[i-1], kappa, theta, sigma, dt)
  }
  return(ll)
}
