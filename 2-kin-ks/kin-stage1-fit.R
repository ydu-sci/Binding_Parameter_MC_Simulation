source("00_global.R")

sim_repN <- 200
ncores   <- 32L

## kd = koff/kon = 100 nM, matching the equilibrium study
kinetics_para <- list(p0 = 10, koff = 1e-2, kon = 1e-4)
sd_level      <- "low"
time_points   <- 10
repN_kin      <- repN  # 1, as in 2-kin.R

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

cat(sprintf("task %d / %d: expN = %d, reps %d-%d\n",
            task, nrow(task_table), expN, min(reps), max(reps)))

# Ligand concentrations and noise-free curves ------------------------------
l0_list <- with(kinetics_para, l0_calc(p0 = p0, kd = koff / kon))
pl_l0_list <- lapply(l0_list, function(x)
    with(kinetics_para, second_order_asso_points(p0 = p0, l0 = x, koff = koff,
                                                 kon = kon, time_points = time_points)))
stopifnot(length(l0_list) > 1)          # both fits need several l0 values

# Fitting function ---------------------------------------------------------
fitfunc_kin1 <- function(plNoisy, paraTrue, paraRatio, l0_list) {
    l0_num <- length(l0_list)
    p0 <- paraTrue$p0; kon0 <- paraTrue$kon; koff0 <- paraTrue$koff
    out <- c(koff_est = NA_real_, kon_est = NA_real_, resi_sd = NA_real_)

    ## dissociation phase -> koff
    dd <- as.data.frame(plNoisy[, 1:3, drop = FALSE])
    colnames(dd) <- c("pl", "t", "g")
    fit_off <- tryCatch(
        nls(pl ~ formula_binding_kinetics_disso_group(t, pl0, g, koff), data = dd,
            start = list(pl0  = rep(p0 * paraRatio, times = l0_num),
                         koff = koff0 * paraRatio)),
        error = function(e) NULL)
    if (is.null(fit_off)) return(out)
    out[["koff_est"]] <- coef(fit_off)[["koff"]]
    out[["resi_sd"]]  <- sd(residuals(fit_off))

    ## association phase, second-order model, koff fixed at this experiment's
    ## own estimate -> kon
    da <- as.data.frame(plNoisy[, 4:7, drop = FALSE])
    colnames(da) <- c("pl", "t", "l0", "g")
    koff_est <- out[["koff_est"]]
    fit_on <- tryCatch(
        nls(pl ~ formula_second_order_kin(t, l0, p0, koff = koff_est, kon), data = da,
            start = list(p0 = p0 * paraRatio, kon = kon0 * paraRatio)),
        error = function(e) NULL)
    if (!is.null(fit_on)) out[["kon_est"]] <- coef(fit_on)[["kon"]]
    out
}

# RNG: one non-overlapping L'Ecuyer stream per task
RNGkind("L'Ecuyer-CMRG")
set.seed(59403)
.s <- .Random.seed
if (task > 1) for (i in 2:task) .s <- nextRNGStream(.s)
assign(".Random.seed", .s, envir = .GlobalEnv)

# Run ----------------------------------------------------------------------
clust <- makeCluster(ncores)
clusterExport(clust, c("formula_binding_kinetics_disso_group",
                       "formula_second_order_kin"))

res_koff   <- matrix(NA_real_, length(reps), expN)
res_kon    <- matrix(NA_real_, length(reps), expN)
res_kinkd  <- matrix(NA_real_, length(reps), expN)
res_rsd    <- matrix(NA_real_, length(reps), expN)
nfail_koff <- integer(length(reps))
nfail_kon  <- integer(length(reps))
nfail_kinkd<- integer(length(reps))

for (k in seq_along(reps)) {
    sim <- gen_plNoisy_array_kin(pl_l0_list, l0_list, paraTrue = kinetics_para,
                                 repN = repN_kin, expN = expN, sd_level = sd_level)

    ## Put both phases in one array so a single parApply covers both fits.
    nrow_kin <- dim(sim$data_array_asso)[1]
    dat <- array(NA_real_, c(nrow_kin, 7, expN))
    dat[, 1:3, ] <- sim$data_array_disso
    dat[, 4:7, ] <- sim$data_array_asso
    rm(sim); gc()

    fit <- parApply(clust, dat, 3, fitfunc_kin1, paraTrue = kinetics_para,
                    paraRatio = ratioList[2], l0_list = l0_list)
    rm(dat); gc()

    res_koff[k, ]  <- fit["koff_est", ]
    res_kon[k, ]   <- fit["kon_est", ]
    res_kinkd[k, ] <- fit["koff_est", ]/fit["kon_est", ]
    res_rsd[k, ]   <- fit["resi_sd", ]
    nfail_koff[k]  <- sum(is.na(fit["koff_est", ]))
    nfail_kon[k]   <- sum(is.na(fit["kon_est", ]))
    nfail_kinkd[k]   <- sum(is.na(res_kinkd[k, ]))
    cat(sprintf("  rep %3d  nfail koff = %d, kon = %d\n, kinkd = %d",
                reps[k], nfail_koff[k], nfail_kon[k], nfail_kinkd[k]))

    rm(fit); gc()
}

saveRDS(list(expN = expN, reps = reps,
             koff = res_koff, kon = res_kon, 
             kinkd =res_kinkd, resi_sd = res_rsd,
             nfail_koff = nfail_koff, nfail_kon = nfail_kon, 
             nfail_kinkd = nfail_kinkd,
             kinetics_para = kinetics_para, sd_level = sd_level),
        sprintf("res_kin_expN%d_part%d.rds", expN, piece))

try(stopCluster(clust), silent = TRUE)
