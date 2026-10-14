sim_repN  <- 200
expN_list <- c(1000, 10000, 100000)

task_table <- rbind(
    data.frame(expN = 1000L,   piece = 1L,  n_piece = 1L),
    data.frame(expN = 10000L,  piece = 1L,  n_piece = 1L),
    data.frame(expN = 100000L, piece = 1:8, n_piece = 8L))

for (par in c("koff", "kon", "kinkd")) {
    for (i in seq_along(expN_list)) {
        N  <- expN_list[i]
        np <- task_table$n_piece[match(N, task_table$expN)]
        fs <- sprintf("res_kin_expN%d_part%d.rds", N, seq_len(np))
        stopifnot(all(file.exists(fs)))
        ps <- lapply(fs, readRDS)

        m <- matrix(NA_real_, sim_repN, N); nf <- integer(sim_repN)
        filled <- logical(sim_repN)
        for (p in ps) {
            m[p$reps, ]    <- p[[par]]
            nf[p$reps]     <- p[[paste0("nfail_", par)]]
            filled[p$reps] <- TRUE
        }
        stopifnot(all(filled))

        assign(paste0("res_kd", i), m); assign(paste0("nfail", i), nf)
        cat(sprintf("%-4s expN = %6d: %d x %d, total nfail = %d\n",
                    par, N, nrow(m), ncol(m), sum(nf)))
    }
    save(res_kd1, res_kd2, res_kd3, nfail1, nfail2, nfail3,
         file = sprintf("res_%s.rda", par))
}
