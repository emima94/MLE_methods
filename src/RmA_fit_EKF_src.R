## Source functions for fitting the RmA model with EKF

create_RmA_model_ctsmTMB <- function() {
model <- ctsmTMB$new()

# Add system dynamics
model$addSystem(
    dX ~ fN * hd_N + 0.5 * gN^2 * hdd_N + gN * hd_N * dW1,
    dY ~ fP * hd_P + 0.5 * gP^2 * hdd_P + gP * hd_P * dW2
)

# Add observations
model$addObs(
    

model$setAlgebraics(
    N = exp(X),
    P = exp(Y),
    fN = r * N * (1 - N / K) - beta * N * P / (1 + beta * N / Cmax),
    fP = epsilon * beta * N * P / (1 + beta * N / Cmax) - mu * P,
    gN = sN * N,
    gP = sP * P,
    hd_N = 1 / N,
    hd_P = 1 / P,
    hdd_N = -1 / N^2,
    hdd_P = -1 / P^2
)



return (model)

}