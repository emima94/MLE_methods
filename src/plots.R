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