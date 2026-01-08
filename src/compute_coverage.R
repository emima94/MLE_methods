# Compute coverage for a parameter:
compute_coverage <- function(estimates, ses, true_value, alpha=0.05) {
    z <- qnorm(1 - alpha/2)
    
    coverage <- sum(abs(estimates - true_value) <= z * ses) / length(estimates)
    
    return(coverage)
}
