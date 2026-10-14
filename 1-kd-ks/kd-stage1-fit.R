source("00_global.R")

sim_repN  <- 200
expN_list <- c(1000, 10000, 100000)
ncores    <- 32L

# Task table --------------------------------------------------------------
task_table <- rbind(
    data.frame(expN = 1000L,   piece = 1L,  n_piece = 1L),
    data.frame(expN = 10000L,  piece = 1L,  n_piece = 1L),
    data.frame(expN = 100000L, piece = 1:4, n_piece = 4L))

task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
stopifnot(task >= 1, task <= nrow(task_table))

expN    <- task_table$expN[task]
n_piece <- task_table$n_piece[task]
piece   <- task_table$piece[task]
reps    <- split(seq_len(sim_repN),
                 rep(seq_len(n_piece), each = ceiling(sim_repN / n_piece),
                     length.out = sim_repN))[[piece]]

cat(sprintf("task %d / %d: expN = %d, reps %d-%d\n",
            task, nrow(task_table), expN, min(reps), max(reps)))

# Fitting function (Equation 1 only) --------------------------------------
fitfunc1 <- function(plNoisy, paraTrue, paraRatio, repN, l_conc) {
    pldata    <- data.frame(pl = as.vector(plNoisy), tl = rep(l_conc, each = repN))
    paraGuess <- list(p0 = paraTrue$p0 * paraRatio, kd = paraTrue$kd * paraRatio)
    ## A single error must not abort the whole batch: one failure out of 2e7
    ## fits would otherwise throw away hours of work.
    tryCatch({
        cf <- coef(nlsLM(pl ~ formula1(tl, p0, kd), data = pldata, start = paraGuess))
        c(p0 = cf[["p0"]], kd = cf[["kd"]])
    }, error = function(e) c(p0 = NA_real_, kd = NA_real_))
}

# RNG: one non-overlapping L'Ecuyer stream per task -----------------------
RNGkind("L'Ecuyer-CMRG")
set.seed(59403)
.s <- .Random.seed
if (task > 1) for (i in 2:task) .s <- nextRNGStream(.s)
assign(".Random.seed", .s, envir = .GlobalEnv)

# Run ---------------------------------------------------------------------
clust <- makeCluster(ncores)
clusterEvalQ(clust, library(minpack.lm))
clusterExport(clust, "formula1")

res_kd <- matrix(NA_real_, nrow = length(reps), ncol = expN)
res_p0 <- matrix(NA_real_, nrow = length(reps), ncol = expN)
nfail  <- integer(length(reps))

for (k in seq_along(reps)) {
    pl <- gen_plNoisy_array(plTheoretical_num_large, paraTrue = para,
                            repN = repN5, expN = expN, sd_level = "medium")

    fit <- parApply(clust, pl, 3, fitfunc1, paraTrue = para,
                    paraRatio = ratioList[2], repN = repN5, l_conc = tl_num_large)

    res_kd[k, ] <- fit["kd", ]
    res_p0[k, ] <- fit["p0", ]
    nfail[k]    <- sum(is.na(fit["kd", ]))
    cat(sprintf("  rep %3d  nfail = %d\n", reps[k], nfail[k]))
}

stopCluster(clust)

saveRDS(list(expN = expN, reps = reps, kd = res_kd, p0 = res_p0, nfail = nfail),
        sprintf("res_kd_expN%d_part%d.rds", expN, piece))