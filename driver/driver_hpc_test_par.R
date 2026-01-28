#!/usr/bin/env Rscript

## ================================
## HPC-safe parallel benchmark
## ================================

library(parallel)
library(future.apply)

## ---- Configuration ----

n_workers <- 4
n_tasks   <- 100
n_obs     <- 5e4        # per task
out_dir   <- "bench_results"

dir.create(out_dir, showWarnings = FALSE)

Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1
)

## ---- Shared data (large, master-only) ----

set.seed(1)
Ysim <- replicate(n_tasks, rnorm(n_obs), simplify = FALSE)

alpha <- 0.3
beta  <- 1.2

task_fun <- function(Yi, alpha, beta) {
  sum(alpha * Yi^2 + beta * Yi)
}

## ================================
## 1. Sequential
## ================================

cat("Running sequential...\n")
t_seq <- system.time({
  for (i in seq_len(n_tasks)) {
    res <- task_fun(Ysim[[i]], alpha, beta)
    saveRDS(res, file.path(out_dir, sprintf("seq_%03d.rds", i)))
  }
})

## ================================
## 2. Base parallel (PSOCK)
## ================================

cat("Running parallel::parLapply...\n")

cl <- makeCluster(n_workers)

clusterExport(
  cl,
  varlist = c("alpha", "beta", "task_fun"),
  envir = environment()
)

t_par <- system.time({
  parLapply(
    cl,
    X = seq_len(n_tasks),
    fun = function(i) {
      Yi <- Ysim[[i]]
      res <- task_fun(Yi, alpha, beta)
      saveRDS(res, file.path(out_dir, sprintf("par_%03d.rds", i)))
      NULL
    }
  )
})

stopCluster(cl)

## ================================
## 3. future (multisession)
## ================================

cat("Running future.apply...\n")

options(future.globals.maxSize = 2 * 1024^3)
plan(multisession, workers = n_workers)

t_fut <- system.time({
  future_lapply(
    seq_len(n_tasks),
    function(i) {
      Yi <- Ysim[[i]]
      res <- task_fun(Yi, alpha, beta)
      saveRDS(res, file.path(out_dir, sprintf("fut_%03d.rds", i)))
      NULL
    },
    future.seed = TRUE
  )
})

plan(sequential)

## ================================
## Summary
## ================================

summary <- rbind(
  sequential = t_seq,
  parallel   = t_par,
  future     = t_fut
)

print(summary)

saveRDS(summary, file.path(out_dir, "timing_summary.rds"))
write.csv(summary, file.path(out_dir, "timing_summary.csv"))

cat("\nBenchmark complete.\n")
