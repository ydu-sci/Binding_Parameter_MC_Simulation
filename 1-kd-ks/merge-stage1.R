sim_repN   <- 200
expN_list  <- c(1000, 10000, 100000)
task_table <- rbind(
    data.frame(expN = 1000L,   piece = 1L,  n_piece = 1L),
    data.frame(expN = 10000L,  piece = 1L,  n_piece = 1L),
    data.frame(expN = 100000L, piece = 1:4, n_piece = 4L))

for (i in seq_along(expN_list)) {
    N  <- expN_list[i]
    np <- task_table$n_piece[match(N, task_table$expN)]
    fs <- sprintf("res_kd_expN%d_part%d.rds", N, seq_len(np))
    stopifnot(all(file.exists(fs)))
    ps <- lapply(fs, readRDS)

    m  <- matrix(NA_real_, sim_repN, N)
    nf <- integer(sim_repN)
    filled <- logical(sim_repN)
    ## Each piece carries its own row indices, so read order does not matter.
    for (p in ps) {
        m[p$reps, ]   <- p$kd
        nf[p$reps]    <- p$nfail
        filled[p$reps] <- TRUE
    }
    stopifnot(all(filled)) # no repeat left unwritten

    assign(paste0("res_kd", i), m)
    assign(paste0("nfail",  i), nf)
    cat(sprintf("expN = %6d: %d x %d, total nfail = %d\n",
                N, nrow(m), ncol(m), sum(nf)))
}

save(res_kd1, res_kd2, res_kd3, nfail1, nfail2, nfail3, file = "res_kd.rda")
