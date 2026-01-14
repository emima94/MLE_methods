# Util

print_computation_time <- function(fit, methods) {
    for (m in methods) {
        N <- length(fit[[m]])
        times_m <- sapply(1:N, function(i) fit[[m]][[i]]$time)
        cat("Method:", m, "\n")
        cat(" Average time (s):", mean(times_m), "\n")
        cat(" SD time (s):", sd(times_m), "\n\n")
    }
}