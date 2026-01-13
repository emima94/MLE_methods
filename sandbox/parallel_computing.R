library(parallel)

# Create a cluster with the number of available cores
num_cores <- detectCores() - 1  # Leave one core free


# Define a slow function to simulate a time-consuming task
square <- function(x) {

  t0 <- Sys.time()  

  Sys.sleep(0.1)  # Simulate a time-consuming computation
  
  y <- x * x

    t1 <- Sys.time()

    out <- list(res = y, time = t1 - t0)
  
  return(out)
}

input <- 1:10

time_seq <- system.time({
  result_seq <- lapply(input, square)
})

result_seq
time_seq

cl <- makeCluster(6)
time_par <- system.time({
  result_par <- parLapply(cl, input, square)
})
stopCluster(cl)

result_par
time_par


# Make an overview of computation times
time_summary <- data.frame(
  Method = c("Sequential", "Parallel"),
  UserTime = c(time_seq["user.self"], time_par["user.self"]),
  SysTime = c(time_seq["sys.self"], time_par["sys.self"]),
  ElapsedTime = c(time_seq["elapsed"], time_par["elapsed"])
)
print(time_summary)
