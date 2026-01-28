
# Summarize fits
metric_names <- c("coverage", "bias", "rmse")
fit_summary_T <- vector("list", length(T_vals))
names(fit_summary_T) <- names(fit_T)
sofun()
for (i in 1:length(T_vals)) {
    T <- T_vals[i]
    fit <- fit_T[[i]]
    message("Computing summary for T = ", T, sep="")
    fit_summary <- compute_fit_summary(fit, p0, methods, n_par, par_names, metric_names)
    
    fit_summary_T[[i]] <- fit_summary
}
fit_summary_T
sofun()
# Computation time vs T_vals #
null_count <- create_computation_time_plots(fit_T, T_vals, methods, "T", "RmA")
null_count  

# Parameter estimation distributions #
sofun()
xlims = list(
            log_r = c(-1,1),
            log_K = c(-1,1),
            log_beta = c(-2,2),
            log_mu = c(-2,2),
            log_sN = c(-16,2)
)
create_parameter_estimation_plots(fit_summary_T, T_vals, methods, 
                        "T", "RmA", par_names, p0, xlims, var_name_title =  "T", var_title_func = function(x) x
                                 )


# Coverage and RMSE plots # 
create_perf_metric_plot(fit_summary_T, T_vals, methods, 
                        "T", "RmA", par_names, "coverage"
                                 )

create_perf_metric_plot(fit_summary_T, T_vals, methods, 
                        "T", "RmA", par_names, "rmse"
                                 )

# Likelihood distribution plots
sofun()
create_likelihood_histograms(fit_T, T_vals, methods, 
                        "T", NA, "RmA", "T", var_title_func = function(x) x
                                 )
# True parameters
theta_true <- sapply(par_names, function(pn) p0[[pn]])
create_likelihood_histograms(fit_T, T_vals, methods, 
                        "T", theta_true, "RmA", "T", var_title_func = function(x) x,
                                 plot_dev = TRUE)
