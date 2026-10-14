source("00_global.R")

sim_repN <- 200
ncores   <- 32L
sd_level <- "medium"

para <- list(p0 = 5, l0 = 50, kd = 10, ki = 100)
cond <- "ki100"

logTi_num_large <- seq(-10, -5, length = num_large)
ti_num_large    <- 10^logTi_num_large * 1e+9

task_table <- rbind(
    data.frame(expN = 1000L,   piece = 1L,  n_piece = 1L),
    data.frame(expN = 10000L,  piece = 1L,  n_piece = 1L),
    data.frame(expN = 100000L, piece = 1:8, n_piece = 8L))

task    <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
stopifnot(task >= 1, task <= nrow(task_table))
expN    <- task_table$expN[task]
n_piece <- task_table$n_piece[task]
piece   <- task_table$piece[task]
reps    <- split(seq_len(sim_repN),
                 rep(seq_len(n_piece), each = ceiling(sim_repN / n_piece),
                     length.out = sim_repN))[[piece]]

cat(sprintf("task %d / %d: ki = %g nM, expN = %d, reps %d-%d\n",
            task, nrow(task_table), para$ki, expN, min(reps), max(reps)))

plTheo <- with(para, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))

fitfunc_comp1 <- function(plNoisy, paraTrue, paraRatio, repN, i_conc, maxtry = 20L) {
    pldata <- data.frame(pl = as.vector(plNoisy), i0 = rep(i_conc, each = repN))
    l0 <- paraTrue$l0; kd <- paraTrue$kd
    out <- c(ki = NA_real_, p0 = NA_real_, ic50 = NA_real_, b = NA_real_)

    ## 1. direct fit of the competitive binding equation -> ki
    f1 <- tryCatch(
        nlsLM(pl ~ formula_competitive_pl(p0, l0 = l0, i0, kd = kd, ki),
              data = pldata,
              start = list(p0 = paraTrue$p0 * paraRatio,
                           ki = paraTrue$ki * paraRatio)),
        error = function(e) NULL)
    if (!is.null(f1)) {
        cf <- coef(f1)
        out[["ki"]] <- cf[["ki"]]; out[["p0"]] <- cf[["p0"]]
    }

    ## 2. four-parameter logistic -> apparent IC50 (parameter c)
    g <- list(a = 0.01, b = -1, c = 100, d = mean(plNoisy[, 1]))
    f2 <- NULL
    for (i in seq_len(maxtry)) {
        f2 <- tryCatch(nlsLM(pl ~ logistic_func(i0, a, b, c, d),
                             data = pldata, start = g),
                       error = function(e) NULL)
        if (!is.null(f2)) break
        g$a <- g$a / 2
    }
    if (!is.null(f2)) {
        cf <- coef(f2)
        out[["ic50"]] <- cf[["c"]]; out[["b"]] <- cf[["b"]]
    }
    out
}

RNGkind("L'Ecuyer-CMRG")
set.seed(59403)
.s <- .Random.seed
if (task > 1) for (i in 2:task) .s <- nextRNGStream(.s)
assign(".Random.seed", .s, envir = .GlobalEnv)

clust <- makeCluster(ncores)
clusterEvalQ(clust, library(minpack.lm))
clusterExport(clust, c("formula_competitive_pl", "logistic_func"))

res_ki   <- matrix(NA_real_, length(reps), expN)
res_ic50 <- matrix(NA_real_, length(reps), expN)
nfail_ki   <- integer(length(reps))
nfail_ic50 <- integer(length(reps))

for (k in seq_along(reps)) {
    pl <- gen_plNoisy_array(plTheo, paraTrue = para, repN = repN5,
                            expN = expN, sd_level = sd_level)

    fit <- parApply(clust, pl, 3, fitfunc_comp1, paraTrue = para,
                    paraRatio = ratioList[2], repN = repN5,
                    i_conc = ti_num_large)
    rm(pl); gc()

    res_ki[k, ]     <- fit["ki", ]
    res_ic50[k, ]   <- fit["ic50", ]
    nfail_ki[k]     <- sum(is.na(fit["ki", ]))
    nfail_ic50[k]   <- sum(is.na(fit["ic50", ]))
    cat(sprintf("  rep %3d  nfail ki = %d, ic50 = %d\n",
                reps[k], nfail_ki[k], nfail_ic50[k]))

    rm(fit); gc()
}

saveRDS(list(cond = cond, para = para, expN = expN, reps = reps,
             ki = res_ki, ic50 = res_ic50,
             nfail_ki = nfail_ki, nfail_ic50 = nfail_ic50,
             ti_conc = ti_num_large, sd_level = sd_level),
        sprintf("res_comp_%s_expN%d_part%d.rds", cond, expN, piece))

try(stopCluster(clust), silent = TRUE)
