suppressPackageStartupMessages(library(ggplot2))
library(svglite)

N    <- 100000L
runs <- 1:2 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
    stop("Usage: Rscript qqplot_n100000.R <file.rda> [object_name]")

infile <- args[1]
if (!file.exists(infile)) stop("File not found: ", infile)

stem <- tools::file_path_sans_ext(basename(infile))
outdir <- dirname(infile) 

e  <- new.env()
nm <- load(infile, envir = e)

if (length(args) >= 2) {
    obj <- args[2]
    if (!obj %in% nm)
        stop("Object '", obj, "' not found in ", infile,
             ". Available: ", paste(nm, collapse = ", "))
    stem <- paste0(stem, "_", obj)   # keep outputs from one file apart
} else {
    cand <- nm[vapply(nm, function(v) is.matrix(e[[v]]) && ncol(e[[v]]) == N,
                      logical(1))]
    if (length(cand) == 0)
        stop("No matrix with ", N, " columns in ", infile,
             ". Objects: ", paste(nm, collapse = ", "))
    if (length(cand) > 1)
        stop("More than one matrix with ", N, " columns in ", infile, ": ",
             paste(cand, collapse = ", "),
             ". Give the object name as the second argument.")
    obj <- cand
}

m <- e[[obj]]
if (!is.matrix(m) || ncol(m) != N || nrow(m) < max(runs))
    stop("Object '", obj, "' is not a matrix with ", N,
         " columns and at least ", max(runs), " rows")
cat(sprintf("%s: using '%s' (%d x %d)\n", infile, obj, nrow(m), ncol(m)))

zstd <- function(y) (y - mean(y)) / sd(y) 

for (irun in runs) {
    x <- m[irun, ]
    x <- x[is.finite(x) & x > 0] 

    qqdat <- rbind(
        data.frame(scale = "Raw scale",   z = zstd(x)),
        data.frame(scale = "log10 scale", z = zstd(log10(x))))
    qqdat$scale <- factor(qqdat$scale, levels = c("Raw scale", "log10 scale"))

    p <- ggplot(qqdat, aes(sample = z)) +
        stat_qq(size = 0.4, shape = 16) +
        geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.5) +
        facet_wrap(~ scale, nrow = 1) +
        labs(title = sprintf("%s, run %d (n = %s)", stem, irun,
                             formatC(length(x), format = "d", big.mark = ",")),
             x = "Theoretical quantiles (standard normal)",
             y = "Standardized sample quantiles") +
        theme_bw()

    outfile <- file.path(outdir, sprintf("%s_n%d_run%d.png", stem, N, irun))
    png(outfile, width = 7, height = 3.6, units = "in", res = 600,
        type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    
    print(p)
    invisible(dev.off())
    cat(sprintf("  wrote %s (n = %d)\n", outfile, length(x)))
}
