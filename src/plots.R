# Plot RmA simulated data and observations
plot_RmA_sim_obs <- function(t, Nsim, Psim, iobs, YN_gauss) {

    # Plot simulated states and observations
    scale = 1.0
    width = 8 * scale
    height = 2 * scale
    # Save as PDF
    pdf(file="figures/RmA-sim-obs_prey.pdf",width=width,height=height)
    #par(mar=c(3,3,2,2))
    par(mar=c(4,4,1,1))
    plot(
    t, Nsim,
    xlab = expression("Time," ~ italic(t)),
    ylab = expression("Prey, " ~ italic(N)),
    type = "l",
    #main = "Prey population with Gaussian observations"
    )
    points(t[iobs],YN_gauss,col="red",pch=1)
    dev.off()
    pdf(file="figures/RmA-sim-obs_predator.pdf",width=width,height=height)
    par(mar=c(4,4,1,1))
    plot(
        t, Psim, 
        xlab = expression("Time," ~ italic(t)),
        ylab = expression("Predator, " ~ italic(P)), 
        type="l"
        )
    dev.off()

}


plot_metric <- function(metric_name, fit_summary, vals, xlab, ylims, methods, par_names) {
    n_par <- length(par_names)
    par(mfrow = c(n_par, 1))

    cols <- seq_along(methods)  # consistent colors

    for (par in par_names) {

        plot(vals,
             sapply(vals, function(val)
                 fit_summary[[paste0(xlab, "_", val)]][[methods[1]]]$perf[par, metric_name]),
             type = "b",
             col  = cols[1],
             ylim = ylims[[metric_name]],
             xlab = xlab,
             ylab = metric_name,
             main = paste("Parameter:", par))

        for (m in methods[-1]) {
            points(vals,
                   sapply(vals, function(val)
                       fit_summary[[paste0(xlab, "_", val)]][[m]]$perf[par, metric_name]),
                   type = "b",
                   col  = cols[methods == m])
        }

        legend("topright",
               legend = methods,
               col    = cols,
               lty    = 1,
               pch    = 1,
               bty    = "n")
    }
}