plot_likelihood_profile <- function(pg_name, var1_vals, var2_vals, nll_matrix_full, p0, df, dt, x0_bar, P0, method, save_fig=TRUE) {
    # Plot likelihood surface
    i <- which(nll_matrix_full == min(nll_matrix_full, na.rm = TRUE), arr.ind = TRUE)
    # clip extreme values for better plotting
    nll_matrix_full_clip <- nll_matrix_full
    threshold <- min(nll_matrix_full, na.rm = TRUE) + 400
    nll_matrix_full_clip[nll_matrix_full > threshold] <- threshold
    pg <- c(paste0("log_", pg_name[1]), paste0("log_", pg_name[2]))

    if (save_fig) {
        pdf(paste0("figures/RmA_likelihood_profile_full_2D_", method, "_", pg_name[1], "_", pg_name[2], ".pdf"), width = 8, height = 6)
        filled.contour(
        var1_vals, var2_vals, nll_matrix_full_clip,
        xlab = pg_name[1],
        ylab = pg_name[2],
        main = paste0("Likelihood surface, in the ", pg_name[2], " vs. ", pg_name[1], " plane"),

        # ---- Color bar label ----
        key.title = title(main = "NLL"),

        plot.axes = {
            axis(1); axis(2)

            # True value
            points(exp(p0[[pg[1]]]), exp(p0[[pg[2]]]),
                pch = 19, cex = 1.5, col = "green")

            # Optimal (minimum NLL)
            points(var1_vals[i[1]], var2_vals[i[2]],
                pch = 1, col = "black", cex = 2, lwd = 3)

            # ---- Legend ----
            legend("topright",
                legend = c("True value", "Optimal (MLE)"),
                pch = c(19, 1),
                col = c("green", "black"),
                pt.cex = c(1.2, 1.5),
                lwd = c(NA, NA),
                bg = "white")
        }
        )
        dev.off()
    }

}