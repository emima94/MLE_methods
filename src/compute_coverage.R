# Compute coverage for a parameter:
compute_coverage <- function(estimates, ses, true_value, alpha = 0.05) {
   
   z <- qnorm(1 - alpha / 2)

    coverage <- mean(abs(estimates - true_value) <= z * ses, na.rm=TRUE)
    return(coverage)
}