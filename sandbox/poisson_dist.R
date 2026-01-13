# Sandbox to  understand the poisson distribution

my_dpois <- function(k, lambda) {
    return(lambda^k * exp(-lambda) / factorial(k))
}

# Test the function
k_values <- 0:20
lambda <- c(0.5, 1, 2.5, 4, 10)

# plot the PMF for different lambda values, lines and circles
plot(k_values, my_dpois(k_values, lambda[1]), type = "b", lwd=2, col="blue",
     xlab="k", ylab="P(X=k)", main="Poisson PMF for different lambda values",
     ylim=c(0, max(my_dpois(k_values, lambda))))
for (i in 2:length(lambda)) {
    points(k_values, my_dpois(k_values, lambda[i]), type="b", lwd=2, col=i+1)
}
legend("topright", legend=paste("lambda =", lambda), col=2:(length(lambda)+1), lwd=2)

# Compare with Gaussian distribution for lambda = 4 and lambda = 10
lines(k_values, dnorm(k_values, mean=4, sd=sqrt(4)), type="l", lty = 2, lwd=2, col="red")
lines(k_values, dnorm(k_values, mean=10, sd=sqrt(10)), type="l", lty = 2, lwd=2, col="green")
legend("topright", legend=c("Poisson PMF", "Gaussian Approximation"), col=c("blue", "red"), lwd=2)

lambda_large <- 100
k_large <- 80:120
poisson_pmf <- my_dpois(k_large, lambda_large)
gaussian_approx <- dnorm(k_large, mean=lambda_large, sd=sqrt(lambda_large))

plot(k_large, poisson_pmf, type="b", lwd=2, col="blue",
     xlab="k", ylab="Probability", main="Poisson vs Gaussian Approximation (lambda=100)",
     ylim=c(0, max(poisson_pmf, gaussian_approx)))
points(k_large, gaussian_approx, type="b", lwd=2, col="red")
legend("topright", legend=c("Poisson PMF", "Gaussian Approximation"), col=c("blue", "red"), lwd=2)
