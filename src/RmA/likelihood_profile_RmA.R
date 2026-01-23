
likelihood_profile_2D_RmA <- function(p0, T, pg_name, n_grid, delta_list) {
    delta_plus <- delta_list$delta_plus_list
    delta_minus <- delta_list$delta_minus_list
    N <- 1
    x0_bar = c(0.5, 0.5)
    P0 <- diag(c(0.01,0.01))
    x0 <- log(matrix(rnorm(N * length(x0_bar), x0_bar, diag(sqrt(P0))), nrow=N, ncol=length(x0_bar)))
    #T <- 200
    dt <- 0.1
    t <- seq(0,T,dt)
    tsample <- 0.2
    sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0[1,,drop=FALSE], p0, dt, N, tsample)
    Ysim <- sim_data$Ysim
    iobs <- sim_data$iobs
    df <- data.frame(t = t[iobs], Y = Ysim[[1]])

    model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)   
    pg <- c(paste0("log_", pg_name[1]), paste0("log_", pg_name[2]))

    var1 <- seq(exp(p0[[pg[1]]]) - delta_minus[[pg[1]]], exp(p0[[pg[1]]]) + delta_plus[[pg[1]]], length.out=n_grid)
    var2 <- seq(exp(p0[[pg[2]]]) - delta_minus[[pg[2]]], exp(p0[[pg[2]]]) + delta_plus[[pg[2]]], length.out=n_grid)

    # Simulate edge cases
    #pdf(paste0("figures/RmA_likelihood_profile_2D_", pg_name[1], "_", pg_name[2], "_edge_cases.pdf"), width = 8, height = 6)
    par(mfrow = c(2,2))
    for (edge_1 in c(min(var1), max(var1)) ) {
        for (edge_2 in c(min(var2), max(var2))) {
            p0_edge <- p0
            p0_edge[[pg[1]]] <- log(edge_1)
            p0_edge[[pg[2]]] <- log(edge_2)
            sim_data <- generate_data_RmA(fsim, gsim, hsim, t, x0[1,,drop=FALSE], p0_edge, dt, 1, tsample)
            Xsim <- sim_data$Xsim
            Ysim <- sim_data$Ysim
            iobs <- sim_data$iobs
            plot(t, exp(Xsim[[1]][,1]), type="l", main=sprintf("%s=%.2f, %s=%.2f", pg_name[1], edge_1, pg_name[2], edge_2))
            points(t[iobs], Ysim[[1]], col="red")
        }
    }
    plt_edge <- recordPlot()
    #dev.off()

    model <- create_RmA_model(p = p0, x0 = x0_bar, P0 = P0)

    if (!("log_r" %in% pg)) {
        model$setParameter(
            log_r = p0$log_r
        )
    }
    if (!("log_mu" %in% pg)) {
        model$setParameter(
            log_mu = p0$log_mu
        )
    }
    if (!("log_beta" %in% pg)) {
        model$setParameter(
            log_beta = p0$log_beta
        )
    }
    if (!("log_K" %in% pg)) {
        model$setParameter(
            log_K = p0$log_K
        )
    }
    if (!("log_sN" %in% pg)) {
        model$setParameter(
            log_sN = p0$log_sN
        )
    }
    model$getParameters()
    nll_func <- model$likelihood(df, method = "laplace",
                    ode.timestep = dt,
                    ode.solver = "rk4",
                    silent = FALSE)


    # Compute likelihood surface
    nll_matrix <- matrix(0, nrow=length(var1), ncol=length(var2))
    for (i in 1:length(var1)) {
        for (j in 1:length(var2)) {
            
            params <- c(
                log(var1[i]),
                log(var2[j])
            )
            nll_value <- nll_func$fn(params)[1]
            nll_matrix[i,j] <- nll_value
        
            #message(sprintf("r=%.2f, beta=%.2f, NLL=%.2f", var1[i], var2[j], nll_value))
        }
    }
    # Plot likelihood surface
    i <- which(nll_matrix == min(nll_matrix), arr.ind = TRUE)
    #pdf(paste0("figures/RmA_likelihood_profile_2D_", pg_name[1], "_", pg_name[2], ".pdf"), width = 8, height = 6)
      filled.contour(var1, var2, nll_matrix,
        xlab = pg_name[1],
        ylab = pg_name[2],
        main = "Likelihood surface",
        plot.axes = {
            axis(1); axis(2)
            points(exp(p0[[pg[1]]]), exp(p0[[pg[2]]]),
                pch = 19, cex = 1.5, col = "green")
            points(var1[i[1]], var2[i[2]], pch = 1, col = "black", cex = 2, lwd = 3)       
        }
        )
        plt_contour <- recordPlot()
    #dev.off()
    return(list(var1 = var1, var2 = var2, nll_matrix = nll_matrix, plt_contour = plt_contour, plt_edge = plt_edge))
}

