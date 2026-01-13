df_fun_HeatEq <- function(t, Ysim, iobs) {
    

    df <- data.frame(
        t = t[iobs]
    )

    for (j in 1:ncol(Ysim)) {
        df[[paste0("Y", j)]] <- Ysim[ , j]
    }

    return(df)

}
