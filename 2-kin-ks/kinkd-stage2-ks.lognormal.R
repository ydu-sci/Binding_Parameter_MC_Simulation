library(parallel)

load("res_kinkd.rda") # res_kd1/2/3, nfail1/2/3

B      <- 2000
ncores <- 32L
SEED   <- 59403
alpha  <- 0.05

skew <- function(x) { m <- mean(x); mean((x - m)^3) / mean((x - m)^2)^1.5 }
kurt <- function(x) { m <- mean(x); mean((x - m)^4) / mean((x - m)^2)^2 - 3 }

do1test <- function(kd, B) {
    kd <- kd[is.finite(kd) & kd > 0] # drop fits that failed (NA)
    y  <- log10(kd)
    n  <- length(y)

    param <- MASS::fitdistr(y, "normal")$estimate
    stat  <- unname(ks.test(y, function(z) pnorm(z, param[1], param[2]))$statistic)

    stat.b <- double(B)
    for (i in 1:B) {
        y.b       <- rnorm(n, param[1], param[2])
        fitted.b  <- MASS::fitdistr(y.b, "normal")$estimate
        stat.b[i] <- ks.test(y.b, function(z) pnorm(z, fitted.b[1], fitted.b[2]))$statistic
    }

    c(n        = n,
      D        = stat,
      # Methods report p = k/B. The +0.5/(B+1) form avoids p = 0 and gives
      # identical rejection decisions at B = 2000
      p        = (sum(stat.b >= stat) + 0.5) / (B + 1),
      skew     = skew(y),
      kurt     = kurt(y),
      skew.raw = skew(kd))
}

# Flatten into 600 independent tasks 
kdlist <- c(split(res_kd1, row(res_kd1)),
            split(res_kd2, row(res_kd2)),
            split(res_kd3, row(res_kd3)))
expN   <- rep(c(ncol(res_kd1), ncol(res_kd2), ncol(res_kd3)),
              c(nrow(res_kd1), nrow(res_kd2), nrow(res_kd3)))

## Largest tasks first, so no worker is left holding a long job at the end.
ord <- order(expN, decreasing = TRUE)

cl <- makeCluster(ncores)
clusterExport(cl, c("skew", "kurt"))
clusterSetRNGStream(cl, SEED) 

res <- parLapplyLB(cl, kdlist[ord], do1test, B = B)

stopCluster(cl)

out <- data.frame(expN = expN[ord], do.call(rbind, res))
out <- out[order(out$expN), ]
save(out, B, file = "kinkd-ks.lognormal.rda")


# Summary table
tab <- do.call(rbind, lapply(sort(unique(out$expN)), function(N) {
    z <- out[out$expN == N, ]
    data.frame(expN     = N,
               R        = nrow(z),
               D.median = median(z$D), 
               D.q25    = quantile(z$D, 0.25),
               D.q75    = quantile(z$D, 0.75),
               reject   = mean(z$p < alpha),
               skew     = median(z$skew),
               kurt     = median(z$kurt),
               n.used   = median(z$n)) # compare against expN
}))
print(tab, row.names = FALSE)
