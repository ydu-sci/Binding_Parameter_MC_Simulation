# MIT License
#
# Copyright (c) 2024 Yu Du
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# We recommend using RStudio to view the document outline. The chapter 
# numbers in the program correspond to the names of the output folders.

# Since the number of experiments is set to 100,000, the simulation requires 
# substantial computation time. A multi-core computer is recommended. 


# * Libraries and Settings of the Program -----------------------------------
library(minpack.lm)
library(drc)
library(KSgeneral)
library(grid)
library(ggplot2)
library(plyr) # required by GeomSplitViolin
library(RColorBrewer)
library(extrafont) 
library(svglite)
library(scales)
library(parallel)

options(warn=2) # all warnings are turned into errors
#font_import(prompt=FALSE)
loadfonts(device="win")
theme_set(theme_bw(base_size=10, base_family="Arial",))

set.seed(59403)

# the number of experiments
expN <- 100000

# the number of replicates
repN   <- 1
repN5  <- 5
repN10 <- 10

# parallel config
clust <- makeCluster(round(detectCores()*0.8))

# * Colors and Plot Functions -----------------------------------------------

brewer_pal_paired <- brewer.pal(n=12, "Paired")
c1_s1 <- brewer_pal_paired[1];  c1_s2 <- brewer_pal_paired[2]
c2_s1 <- brewer_pal_paired[3];  c2_s2 <- brewer_pal_paired[4]
c3_s1 <- brewer_pal_paired[5];  c3_s2 <- brewer_pal_paired[6]
c4_s1 <- brewer_pal_paired[7];  c4_s2 <- brewer_pal_paired[8]
c5_s1 <- brewer_pal_paired[9];  c5_s2 <- brewer_pal_paired[10]
c6_s1 <- brewer_pal_paired[11]; c6_s2 <- brewer_pal_paired[12]
f1_color <- c(c1_s2, c2_s2, c3_s2, c4_s2, c5_s2, c6_s2)

# plt_color <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728") below is 50% lighten this color
plt_color <- c("#7fbee9", "#ffbf86", "#87de87", "#eb9293") #50% 

legend_theme <-  theme(legend.title=element_blank(), 
                       legend.background=element_blank(),
                       legend.key=element_blank(),
                       legend.position="inside",
                       legend.position.inside=c(1, 1),
                       legend.justification=c(1, 1),
                       legend.key.size=unit(0.2, 'cm'),
                       )

# Adapted violin plot codes from following two sources:
# 1. https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2 @Trang Q. Nguyen's answer
# 2. https://stackoverflow.com/questions/47651868/split-violin-plot-with-ggplot2-with-quantiles @Axeman's answer
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
   draw_group=function(self, data, ..., nudge=0, draw_quantiles=NULL) {
     # Original function by Jan Gleixner (@jan-glx)
     # Adjustments by Wouter van der Bijl (@Axeman)
     data <- transform(data, xminv=x - violinwidth * (x - xmin), xmaxv=x + violinwidth * (xmax - x))
     grp <- data[1, "group"]
     newdata <- plyr::arrange(transform(data, x=if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
     newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
     newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
     # now nudge them apart
     newdata$x <- ifelse(newdata$group %% 2 == 1, newdata$x - nudge, newdata$x + nudge)
     
     if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
       stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
       quantiles <- create_quantile_segment_frame(data, nudge, draw_quantiles, grp=grp)
       aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop=FALSE]
       aesthetics$alpha <- rep(1, nrow(quantiles))
       both <- cbind(quantiles, aesthetics)
       quantile_grob <- GeomPath$draw_panel(both, ...)
       ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
     }
     else {
       ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
     }
   }
)

#override original func in geom-violin.r
create_quantile_segment_frame <- function(data, nudge=0, draw_quantiles, grp=NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y, ties = "ordered")
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x=ggplot2:::interleave(violin.xs + nudge, violin.xmaxvs + nudge),
      y=rep(ys, each=2),
      group=rep(ys, each=2)
    )
  } else {
    data.frame(
      x=ggplot2:::interleave(violin.xminvs - nudge, violin.xs - nudge),
      y=rep(ys, each=2),
      group=rep(ys, each=2)
    )
  }
}

geom_split_violin <- function(mapping=NULL, data=NULL, stat="ydensity", position="identity", ..., 
                              nudge=0, draw_quantiles=NULL, trim=TRUE, scale="area", na.rm=FALSE, 
                              show.legend=NA, inherit.aes=TRUE) {
  layer(data=data, mapping=mapping, stat=stat, geom=GeomSplitViolin, position=position, 
        show.legend=show.legend, inherit.aes=inherit.aes, 
        params=list(trim=trim, scale=scale, nudge=nudge, draw_quantiles=draw_quantiles, na.rm=na.rm, ...))
}


# Part 1. Binding Experiment ----------------------------------------------

# parameters used to simulate, the unit is nM
para <- list(p0=50, kd=100)

## formula I based on p0 and kd
formula1 <- function(x, p0, kd)  ((p0+x+kd)-sqrt((p0+x+kd)^2-4*p0*x))/2

## formula II based on p0 and kd
formula2 <- function(x, p0, kd)  p0*x/(kd+x)

formulaStringList <- c("Equation 1", "Equation 2")

## formula I based on p0 and log10_kd
formula1_log10 <- function(x, p0, log10_kd)  ((p0+x+10^log10_kd)-sqrt((p0+x+10^log10_kd)^2-4*p0*x))/2

## formula II based on p0 and log10_kd
formula2_log10 <- function(x, p0, log10_kd)  p0*x/(10^log10_kd+x)

clusterExport(clust, c("nlsLM", "formula1", "formula2", "formula1_log10", "formula2_log10"))
clusterExport(clust, c("nlsLM", "formula1", "formula2", "formula1_log10", "formula2_log10", "wfct"))

# the number of concentrations of total ligand (free and bound)
num_small   <- 4
num_medium  <- 6
num_large   <- 9
num_Xlarge  <- 15
num_XXlarge <- 30

# Concentrations of total ligand evenly distributed on a logarithmic scale
# NOTE: if maximal end values of the sequence is -4, fitting with formula1 will give:
# Error in nls step factor 0.000488281 reduced below 'minFactor' of 0.000976562, 
logTl_num_small   <- seq(-9, -5, length=num_small)
logTl_num_medium  <- seq(-9, -5, length=num_medium)
logTl_num_large   <- seq(-9, -5, length=num_large)
logTl_num_Xlarge  <- seq(-9, -5, length=num_Xlarge)
logTl_num_XXlarge <- seq(-9, -5, length=num_XXlarge)

# The unit of total ligand is nM
tl_num_small   <- 10^logTl_num_small   * (1e+9)
tl_num_medium  <- 10^logTl_num_medium  * (1e+9)
tl_num_large   <- 10^logTl_num_large   * (1e+9)
tl_num_Xlarge  <- 10^logTl_num_Xlarge  * (1e+9)
tl_num_XXlarge <- 10^logTl_num_XXlarge * (1e+9)

# theoretical concentrations of protein-ligand complex, its unit is nM
plTheoretical_num_small   <- formula1(tl_num_small,   para$p0, para$kd)
plTheoretical_num_medium  <- formula1(tl_num_medium,  para$p0, para$kd)
plTheoretical_num_large   <- formula1(tl_num_large,   para$p0, para$kd)
plTheoretical_num_Xlarge  <- formula1(tl_num_Xlarge,  para$p0, para$kd)
plTheoretical_num_XXlarge <- formula1(tl_num_XXlarge, para$p0, para$kd)


gen_plNoisy_array <- function(plTheoretical, paraTrue, repN, expN, sd_level, plbias=0, distribution="normal") {
  # sd_level: low 1.25%, medium 2.5%, high 5% of p0
  # the sd of normal and uniform distribution is the same, The range 
  # covering 95% of a normal distribution is approximately within [-2*sigma, 2*sigma].
  # both are 2*sigma = para*2.5%, para*5%, or para*10%, so sigma = para*2.5%/2, para*5%/2, or para*10%/2
  # For uniform distribution, sigma^2 = (max-min)^2/12, max=mu+step, min=mu-step, so sigma = step/sqrt(3)
  # 2*sigma = 2*step/sqrt(3) = para*2.5%, para*5%, or para*10%, 
  # step = sqrt(3)*para*2.5%/2, sqrt(3)*para*5%/2, or sqrt(3)*para*10%/2 
  
  plTheoretical_len <- length(plTheoretical)
  p0 <- paraTrue$p0
  error_sd_low      <- p0*0.025/2
  error_sd_medium   <- p0*0.05/2
  error_sd_high     <- p0*0.1/2
  error_sd_increase <- seq(error_sd_low, error_sd_high, length.out=plTheoretical_len)
  error_sd_decrease <- rev(error_sd_increase)
  
  uniform_step_low    <- sqrt(3)*p0*0.025/2
  uniform_step_medium <- sqrt(3)*p0*0.05/2
  uniform_step_high   <- sqrt(3)*p0*0.1/2
  
  # sd_level selection
  if (sd_level=="low"){
    error_sd <- error_sd_low
    uniform_step <- uniform_step_low
  }
  if (sd_level=="medium"){
    error_sd <- error_sd_medium
    uniform_step <- uniform_step_medium
  }
  if (sd_level=="high"){
    error_sd <- error_sd_high
    uniform_step <- uniform_step_high
  }
  if (sd_level=="increase"){
    error_sd <- error_sd_increase
  }
  if (sd_level=="decrease"){
    error_sd <- error_sd_decrease
  }
  
  # distribution selection
  if (distribution=="normal"){
    error_array <- array(rnorm(repN*plTheoretical_len*expN, sd=error_sd), dim=c(repN, plTheoretical_len, expN)) + plbias
  }
  if (distribution=="uniform"){
    error_array <- array(runif(repN*plTheoretical_len*expN, min=-uniform_step, max=uniform_step), dim=c(repN, plTheoretical_len, expN)) + plbias
  }
  rep(rep(plTheoretical, each=repN), times=expN) + error_array
}

fitfunc <- function(plNoisy, paraTrue, paraRatio, repN, l_conc, algorithm="LM", log10=FALSE, wt=FALSE) {
  # algorithm can be Gauss–Newton (GN) or Levenberg-Marquardt (LM)
  # reshape the wide format matrix to long format data frame
  
  pldata <- data.frame(pl=as.vector(plNoisy), tl=rep(l_conc, each=repN))
  
  if (log10) {
    paraGuess <- list(p0=paraTrue$p0*paraRatio, kd=log10(paraTrue$kd*paraRatio))
  } else {
    paraGuess <- list(p0=paraTrue$p0*paraRatio, kd=paraTrue$kd*paraRatio)
  }
  
  if (algorithm=="GN"){
    if (log10) {
      nls.formula1 <- nls(pl ~ formula1_log10(tl, p0, kd), data=pldata, start=paraGuess)
      nls.formula2 <- nls(pl ~ formula2_log10(tl, p0, kd), data=pldata, start=paraGuess)
    } else {
      nls.formula1 <- nls(pl ~ formula1(tl, p0, kd), data=pldata, start=paraGuess)
      nls.formula2 <- nls(pl ~ formula2(tl, p0, kd), data=pldata, start=paraGuess)
    }
  }
  if (algorithm=="LM"){
    if (log10) {
      nls.formula1 <- nlsLM(pl ~ formula1_log10(tl, p0, kd), data=pldata, start=paraGuess)
      nls.formula2 <- nlsLM(pl ~ formula2_log10(tl, p0, kd), data=pldata, start=paraGuess)
    } else {
      if (wt) {
        var.data <- tapply(pldata$pl, pldata$tl, var)
        var.data <- rep(var.data, each=repN)
        nls.formula1 <- nlsLM(pl ~ formula1(tl, p0, kd), data=pldata, start=paraGuess, weights=1/var.data^2)
        nls.formula2 <- nlsLM(pl ~ formula2(tl, p0, kd), data=pldata, start=paraGuess, weights=1/var.data^2)
      } else {
        nls.formula1 <- nlsLM(pl ~ formula1(tl, p0, kd), data=pldata, start=paraGuess)
        nls.formula2 <- nlsLM(pl ~ formula2(tl, p0, kd), data=pldata, start=paraGuess)
      }
    }
  }
  
  f1_coef <- coef(nls.formula1); f2_coef <- coef(nls.formula2)
  
  c(f1_p0=f1_coef[["p0"]], f1_kd=f1_coef[["kd"]], 
    f2_p0=f2_coef[["p0"]], f2_kd=f2_coef[["kd"]])
}

write_table <- function(filename, listname, ...) {
  arg1 <- as.list(environment())
  filename_base <- arg1$filename
  csvfilename <- paste(filename_base, ".csv", sep="")
  stringList <- arg1$listname

  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("[P]0 (SD)", "Kd (SD)", 
                   "KS test P-value of [P]0", "KS test P-value of log10([P]0)", 
                   "KS test P-value of Kd", "KS test P-value of log10(Kd)", 
                   "Pearson's correlation coefficient of [P]0 and Kd", 
                   "Spearman's correlation coefficient of [P]0 and Kd")
  table_df <- data.frame(matrix(ncol=colN*2, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  for (i in 1:colN) {
    f1_p0 <- arg2[[i]]["f1_p0",]; f1_kd <- arg2[[i]]["f1_kd",]
    f2_p0 <- arg2[[i]]["f2_p0",]; f2_kd <- arg2[[i]]["f2_kd",]
    f1_p0_mean <- mean(f1_p0); f1_kd_mean <- mean(f1_kd)
    f2_p0_mean <- mean(f2_p0); f2_kd_mean <- mean(f2_kd)
    f1_p0_sd <- sd(f1_p0); f1_kd_sd <- sd(f1_kd)
    f2_p0_sd <- sd(f2_p0); f2_kd_sd <- sd(f2_kd)
    
    log10_f1_p0 <- log10(f1_p0); log10_f1_kd <- log10(f1_kd)
    log10_f2_p0 <- log10(f2_p0); log10_f2_kd <- log10(f2_kd)
    log10_f1_p0_mean <- mean(log10_f1_p0); log10_f1_kd_mean <- mean(log10_f1_kd)
    log10_f2_p0_mean <- mean(log10_f2_p0); log10_f2_kd_mean <- mean(log10_f2_kd)
    log10_f1_p0_sd <- sd(log10_f1_p0); log10_f1_kd_sd <- sd(log10_f1_kd)
    log10_f2_p0_sd <- sd(log10_f2_p0); log10_f2_kd_sd <- sd(log10_f2_kd)
    
    f1_p0_ks_pvalue <- cont_ks_test(f1_p0, "pnorm", f1_p0_mean, f1_p0_sd)$p.value
    f1_kd_ks_pvalue <- cont_ks_test(f1_kd, "pnorm", f1_kd_mean, f1_kd_sd)$p.value
    f2_p0_ks_pvalue <- cont_ks_test(f2_p0, "pnorm", f2_p0_mean, f2_p0_sd)$p.value
    f2_kd_ks_pvalue <- cont_ks_test(f2_kd, "pnorm", f2_kd_mean, f2_kd_sd)$p.value
    
    log10_f1_p0_ks_pvalue <- cont_ks_test(log10_f1_p0, "pnorm", log10_f1_p0_mean, log10_f1_p0_sd)$p.value
    log10_f1_kd_ks_pvalue <- cont_ks_test(log10_f1_kd, "pnorm", log10_f1_kd_mean, log10_f1_kd_sd)$p.value
    log10_f2_p0_ks_pvalue <- cont_ks_test(log10_f2_p0, "pnorm", log10_f2_p0_mean, log10_f2_p0_sd)$p.value
    log10_f2_kd_ks_pvalue <- cont_ks_test(log10_f2_kd, "pnorm", log10_f2_kd_mean, log10_f2_kd_sd)$p.value
    
    f1_p0_kd_rp <- cor(f1_p0, f1_kd, method="pearson"); f2_p0_kd_rp <- cor(f2_p0, f2_kd, method="pearson")
    f1_p0_kd_rs <- cor(f1_p0, f1_kd, method="spearman"); f2_p0_kd_rs <- cor(f2_p0, f2_kd, method="spearman")
    
    colName <- paste(stringList[i], "f1", sep="-")
    colName2 <- paste(stringList[i], "f2", sep="-")
    idx1 <- i*2 - 1
    idx2 <- i*2
    table_df[,idx1] <- c(sprintf("%.2f (%.2f)", f1_p0_mean, f1_p0_sd), sprintf("%.2f (%.2f)", f1_kd_mean, f1_kd_sd),
                   sprintf("%.3f", f1_p0_ks_pvalue), sprintf("%.3f", log10_f1_p0_ks_pvalue),
                   sprintf("%.3f", f1_kd_ks_pvalue), sprintf("%.3f", log10_f1_kd_ks_pvalue),
                   sprintf("%.3f", f1_p0_kd_rp), sprintf("%.3f", f1_p0_kd_rs)
                   )
    table_df[,idx2] <- c(sprintf("%.2f (%.2f)", f2_p0_mean, f2_p0_sd), sprintf("%.2f (%.2f)", f2_kd_mean, f2_kd_sd),
                   sprintf("%.3f", f2_p0_ks_pvalue), sprintf("%.3f", log10_f2_p0_ks_pvalue),
                   sprintf("%.3f", f2_kd_ks_pvalue), sprintf("%.3f", log10_f2_kd_ks_pvalue),
                   sprintf("%.3f", f2_p0_kd_rp), sprintf("%.3f", f2_p0_kd_rs)
                   )
    colnames(table_df)[idx1] <- colName
    colnames(table_df)[idx2] <- colName2
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

# * 1.1-SI Comparison of Algorithm and Model  ------------------

# 2*(sd of measurements) should be 2.5%, 5%, and 10% of p0 for reasonable experiments
# so sd is 1.25%, 2.5% and 5% of p0. Note: if sd is too high, NAN error can occur when fitting

# Initial values for estimating parameters
# setting ratio to 0.01, 5 or 10 is bad, this will generate following error when fitting 
# error in sqrt((p0 + x + kd)^2 - 4 * p0 * x) :  (converted from warning) NaNs produced
# setting ratio 0.25 and sd_level="high" can also generate the same error
ratioList <- c(0.25, 0.5, 1, 2) 
ratioStringList <- c("Init-0.25", "Init-0.50", "Init-1.00", "Init-2.00")

# simulation of experiments
plNoisy_array <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="medium")

fitResultRatio05  <- parApply(clust, plNoisy_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

#log10
fitResultRatio05_log10  <- parApply(clust, plNoisy_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large, log10=TRUE)

fitResultRatio05_GN  <- parApply(clust, plNoisy_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large, algorithm="GN")

# *- Note-01 ---------------------------------------------------------------
# fitting using different algorithms generates very similar estimates 
sub2 <- fitResultRatio05["f1_p0",] - fitResultRatio05_GN["f1_p0",]
fivenum(sub2)
sub2 <- fitResultRatio05["f1_kd",] - fitResultRatio05_GN["f1_kd",]
fivenum(sub2)
sub2 <- fitResultRatio05["f2_p0",] - fitResultRatio05_GN["f2_p0",]
fivenum(sub2)
sub2 <- fitResultRatio05["f2_kd",] - fitResultRatio05_GN["f2_kd",]
fivenum(sub2)

sub2 <- fitResultRatio05["f1_p0",] - fitResultRatio05_log10["f1_p0",]
fivenum(sub2)
sub2_10 <- fitResultRatio05["f1_kd",] - 10^fitResultRatio05_log10["f1_kd",]
fivenum(sub2_10)
sub2_log10 <- log10(fitResultRatio05["f1_kd",]) - fitResultRatio05_log10["f1_kd",]
fivenum(sub2_log10)
sub2 <- fitResultRatio05["f2_p0",] - fitResultRatio05_log10["f2_p0",]
fivenum(sub2)
sub2_10 <- fitResultRatio05["f2_kd",] - 10^fitResultRatio05_log10["f2_kd",]
fivenum(sub2_10)
sub2_log10 <- log10(fitResultRatio05["f2_kd",]) - fitResultRatio05_log10["f2_kd",]
fivenum(sub2_log10)

# * 1.2-SI Errors and Replication of Experiments ---------------------------

# *- 1.2.1-SI Standard Deviation of Experimental Errors ----------------------

dir.create("Part1.2.1-SI"); setwd("Part1.2.1-SI")

sdStringList <- c("Low SD", "Medium SD", "High SD", "Increasing SD", "Decrease SD")

# simulation of experiments
plNoisy_errorSDLow_array      <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN10, expN=expN, sd_level="low")
plNoisy_errorSDMedium_array   <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN10, expN=expN, sd_level="medium")
plNoisy_errorSDHigh_array     <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN10, expN=expN, sd_level="high")
plNoisy_errorSDIncrease_array <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN10, expN=expN, sd_level="increase")
plNoisy_errorSDDecrease_array <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN10, expN=expN, sd_level="decrease")

# fitting process
fitResultErrorSD_low      <- parApply(clust, plNoisy_errorSDLow_array,      3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN10, l_conc=tl_num_large)
fitResultErrorSD_medium   <- parApply(clust, plNoisy_errorSDMedium_array,   3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN10, l_conc=tl_num_large)
fitResultErrorSD_high     <- parApply(clust, plNoisy_errorSDHigh_array,     3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN10, l_conc=tl_num_large)
fitResultErrorSD_increase <- parApply(clust, plNoisy_errorSDIncrease_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN10, l_conc=tl_num_large, wt=TRUE)
fitResultErrorSD_decrease <- parApply(clust, plNoisy_errorSDDecrease_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN10, l_conc=tl_num_large, wt=TRUE)

# violin plots of different sd comparing formula 1 and 2
fig02_p0_df <- data.frame(
  x=rep(rep(sdStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResultErrorSD_low["f1_p0",], fitResultErrorSD_medium["f1_p0",], fitResultErrorSD_high["f1_p0",], fitResultErrorSD_increase["f1_p0",], fitResultErrorSD_decrease["f1_p0",],
      fitResultErrorSD_low["f2_p0",], fitResultErrorSD_medium["f2_p0",], fitResultErrorSD_high["f2_p0",], fitResultErrorSD_increase["f2_p0",], fitResultErrorSD_decrease["f2_p0",]),
  f=rep(formulaStringList, each=expN*length(sdStringList))
)
fig02_p0_violin <- ggplot(fig02_p0_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated [P]0 (nM)") +
  scale_x_discrete(limits=sdStringList, guide=guide_axis(n.dodge=2)) +
  theme(axis.title.x=element_blank())+
  legend_theme

fig02_kd_df <- data.frame(
  x=rep(rep(sdStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResultErrorSD_low["f1_kd",], fitResultErrorSD_medium["f1_kd",], fitResultErrorSD_high["f1_kd",], fitResultErrorSD_increase["f1_kd",], fitResultErrorSD_decrease["f1_kd",],
      fitResultErrorSD_low["f2_kd",], fitResultErrorSD_medium["f2_kd",], fitResultErrorSD_high["f2_kd",], fitResultErrorSD_increase["f2_kd",], fitResultErrorSD_decrease["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(sdStringList))
)
fig02_kd_violin <- ggplot(fig02_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=sdStringList, guide=guide_axis(n.dodge=2)) +
  theme(axis.title.x=element_blank())+
  legend_theme

ggsave("p0.png", fig02_p0_violin, dpi=800, width=7, height=7, units='cm')
ggsave("p0.svg", fig02_p0_violin, width=7, height=7, units='cm')
ggsave("kd.png", fig02_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig02_kd_violin, width=7, height=7, units='cm')

setwd("..")

# *- 1.2.2-SI Accuracy of Experiments  ---------------------------------------

dir.create("Part1.2.2-SI"); setwd("Part1.2.2-SI")

biasStringList <- c("Positive bias", "No bias", "Negative bias")

# simulation of experiments
plNoisy_posBias_array <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="medium", plbias=para$p0*0.025)
plNoisy_noBias_array  <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="medium")
plNoisy_negBias_array <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="medium", plbias=-para$p0*0.025)

# fitting process
fitResult_posBias <- parApply(clust, plNoisy_posBias_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_noBias  <- parApply(clust, plNoisy_noBias_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_negBias <- parApply(clust, plNoisy_negBias_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

# violin plots of different biases comparing formula 1 and 2
fig03_p0_df <- data.frame(
  x=rep(rep(biasStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_posBias["f1_p0",], fitResult_noBias["f1_p0",], fitResult_negBias["f1_p0",],
      fitResult_posBias["f2_p0",], fitResult_noBias["f2_p0",], fitResult_negBias["f2_p0",]),
  f=rep(formulaStringList, each=expN*length(biasStringList))
)
fig03_p0_violin <- ggplot(fig03_p0_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated [P]0 (nM)") +
  scale_x_discrete(limits=biasStringList) +
  theme(axis.title.x=element_blank())+
  legend_theme

fig03_kd_df <- data.frame(
  x=rep(rep(biasStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_posBias["f1_kd",], fitResult_noBias["f1_kd",], fitResult_negBias["f1_kd",],
      fitResult_posBias["f2_kd",], fitResult_noBias["f2_kd",], fitResult_negBias["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(biasStringList))
)
fig03_kd_violin <- ggplot(fig03_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=biasStringList) +
  theme(axis.title.x=element_blank())+
  legend_theme

ggsave("p0.png", fig03_p0_violin, dpi=800, width=7, height=7, units='cm')
ggsave("p0.svg", fig03_p0_violin, width=7, height=7, units='cm')
ggsave("kd.png", fig03_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig03_kd_violin, width=7, height=7, units='cm')

setwd("..")

# *- 1.2.3-SI Number of Replicates ----------------------------------------------

dir.create("Part1.2.3-SI"); setwd("Part1.2.3-SI")

repNList <- c(1, 3, 6, 9, 20)
repNStringList <- c("1Rep", "3Rep", "6Rep", "9Rep", "20Rep")

# simulation of experiments
plNoisy_repN1_array  <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repNList[1], expN=expN, sd_level="medium")
plNoisy_repN3_array  <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repNList[2], expN=expN, sd_level="medium")
plNoisy_repN6_array  <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repNList[3], expN=expN, sd_level="medium")
plNoisy_repN9_array  <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repNList[4], expN=expN, sd_level="medium")
plNoisy_repN20_array <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repNList[5], expN=expN, sd_level="medium")

# parallel
fitResult_repN1  <- parApply(clust, plNoisy_repN1_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[1], l_conc=tl_num_large)
fitResult_repN3  <- parApply(clust, plNoisy_repN3_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[2], l_conc=tl_num_large)
fitResult_repN6  <- parApply(clust, plNoisy_repN6_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[3], l_conc=tl_num_large)
fitResult_repN9  <- parApply(clust, plNoisy_repN9_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[4], l_conc=tl_num_large)
fitResult_repN20 <- parApply(clust, plNoisy_repN20_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[5], l_conc=tl_num_large)

# violin plots of different numbers of replicates
fig05_p0_df <- data.frame(
  x=rep(rep(repNStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_repN1["f1_p0",], fitResult_repN3["f1_p0",], fitResult_repN6["f1_p0",], fitResult_repN9["f1_p0",], fitResult_repN20["f1_p0",],
      fitResult_repN1["f2_p0",], fitResult_repN3["f2_p0",], fitResult_repN6["f2_p0",], fitResult_repN9["f2_p0",], fitResult_repN20["f2_p0",]),
  f=rep(formulaStringList, each=expN*length(repNStringList))
)
fig05_p0_violin <- ggplot(fig05_p0_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="area") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated [P]0 (nM)") +
  scale_x_discrete(limits=repNStringList) +
  theme(axis.title.x=element_blank())+
  legend_theme

fig05_kd_df <- data.frame(
  x=rep(rep(repNStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_repN1["f1_kd",], fitResult_repN3["f1_kd",], fitResult_repN6["f1_kd",], fitResult_repN9["f1_kd",], fitResult_repN20["f1_kd",],
      fitResult_repN1["f2_kd",], fitResult_repN3["f2_kd",], fitResult_repN6["f2_kd",], fitResult_repN9["f2_kd",], fitResult_repN20["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(repNStringList))
)
fig05_kd_violin <- ggplot(fig05_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="area") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=repNStringList) +
  theme(axis.title.x=element_blank())+
  legend_theme

ggsave("p0.png", fig05_p0_violin, dpi=800, width=7, height=7, units='cm')
ggsave("p0.svg", fig05_p0_violin, width=7, height=7, units='cm')
ggsave("kd.png", fig05_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig05_kd_violin, width=7, height=7, units='cm')

setwd("..")

# 1.3 Experimental Design of Ligand and Protein Conc. --------

# *- 1.3.1 Number of Concentrations of Total Ligand --------------------------------

dir.create("Part1.3.1"); setwd("Part1.3.1")

ligConcStringList <- c("4LigConcs", "6LigConcs", "9LigConcs", "15LigConcs", "30LigConcs")

# simulation of experiments
plNoisy_4LigConcs_array  <- gen_plNoisy_array(plTheoretical_num_small,   paraTrue=para, repN=repN5, expN=expN, sd_level="medium")
plNoisy_6LigConcs_array  <- gen_plNoisy_array(plTheoretical_num_medium,  paraTrue=para, repN=repN5, expN=expN, sd_level="medium")
plNoisy_9LigConcs_array  <- gen_plNoisy_array(plTheoretical_num_large,   paraTrue=para, repN=repN5, expN=expN, sd_level="medium")
plNoisy_15LigConcs_array <- gen_plNoisy_array(plTheoretical_num_Xlarge,  paraTrue=para, repN=repN5, expN=expN, sd_level="medium")
plNoisy_30LigConcs_array <- gen_plNoisy_array(plTheoretical_num_XXlarge, paraTrue=para, repN=repN5, expN=expN, sd_level="medium")

# fitting process
fitResult_4LigConcs  <- parApply(clust, plNoisy_4LigConcs_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_small)
fitResult_6LigConcs  <- parApply(clust, plNoisy_6LigConcs_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_medium)
fitResult_9LigConcs  <- parApply(clust, plNoisy_9LigConcs_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_15LigConcs <- parApply(clust, plNoisy_15LigConcs_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_Xlarge)
fitResult_30LigConcs <- parApply(clust, plNoisy_30LigConcs_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_XXlarge)

# violin plots of different numbers of ligand concentrations comparing formula 1 and 2
fig06_p0_df <- data.frame(
  x=rep(rep(ligConcStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_4LigConcs["f1_p0",], fitResult_6LigConcs["f1_p0",], fitResult_9LigConcs["f1_p0",], fitResult_15LigConcs["f1_p0",], fitResult_30LigConcs["f1_p0",],
      fitResult_4LigConcs["f2_p0",], fitResult_6LigConcs["f2_p0",], fitResult_9LigConcs["f2_p0",], fitResult_15LigConcs["f2_p0",], fitResult_30LigConcs["f2_p0",]),
  f=rep(formulaStringList, each=expN*length(ligConcStringList))
)
fig06_p0_violin <- ggplot(fig06_p0_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated [P]0 (nM)") +
  scale_x_discrete(limits=ligConcStringList, guide=guide_axis(n.dodge=2)) +
  theme(axis.title.x=element_blank()) +
  legend_theme

fig06_kd_df <- data.frame(
  x=rep(rep(ligConcStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_4LigConcs["f1_kd",], fitResult_6LigConcs["f1_kd",], fitResult_9LigConcs["f1_kd",], fitResult_15LigConcs["f1_kd",], fitResult_30LigConcs["f1_kd",],
      fitResult_4LigConcs["f2_kd",], fitResult_6LigConcs["f2_kd",], fitResult_9LigConcs["f2_kd",], fitResult_15LigConcs["f2_kd",], fitResult_30LigConcs["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(ligConcStringList))
)
fig06_kd_violin <- ggplot(fig06_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=ligConcStringList, guide=guide_axis(n.dodge=2)) +
  theme(axis.title.x=element_blank()) +
  legend_theme

ggsave("p0.png", fig06_p0_violin, dpi=800, width=7, height=7, units='cm')
ggsave("p0.svg", fig06_p0_violin, width=7, height=7, units='cm')
ggsave("kd.png", fig06_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig06_kd_violin, width=7, height=7, units='cm')

write_table(filename="Output", listname=ligConcStringList, fitResult_4LigConcs, fitResult_6LigConcs, fitResult_9LigConcs, fitResult_15LigConcs, fitResult_30LigConcs)
setwd("..")

# *- 1.3.2 Concentrations of Total Protein -----------------------------------------------

dir.create("Part1.3.2"); setwd("Part1.3.2")

para1 <- list(p0=10, kd=100)
para2 <- list(p0=20, kd=100)
para3 <- list(p0=100, kd=100)
para4 <- list(p0=200, kd=100) 
# p0=500, will generate error:
# Error in sqrt((p0 + x + kd)^2 - 4 * p0 * x) : (converted from warning) NaNs produced

plTheoretical_para1 <- formula1(tl_num_large, para1$p0, para1$kd)
plTheoretical_para2 <- formula1(tl_num_large, para2$p0, para2$kd)
plTheoretical_para3 <- formula1(tl_num_large, para3$p0, para3$kd)
plTheoretical_para4 <- formula1(tl_num_large, para4$p0, para4$kd)

p0StringList <- c("[P]0 = 0.1Kd", "[P]0 = 0.2Kd", "[P]0 = 1.0Kd", "[P]0 = 2.0Kd")

# simulation of experiments
plNoisy_para1 <- gen_plNoisy_array(plTheoretical_para1, paraTrue=para1, repN=repN5, expN=expN, sd_level="medium")
plNoisy_para2 <- gen_plNoisy_array(plTheoretical_para2, paraTrue=para2, repN=repN5, expN=expN, sd_level="medium")
plNoisy_para3 <- gen_plNoisy_array(plTheoretical_para3, paraTrue=para3, repN=repN5, expN=expN, sd_level="medium")
plNoisy_para4 <- gen_plNoisy_array(plTheoretical_para4, paraTrue=para4, repN=repN5, expN=expN, sd_level="medium")

# fitting process
fitResult_para1 <- parApply(clust, plNoisy_para1, 3, fitfunc, paraTrue=para1, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_para2 <- parApply(clust, plNoisy_para2, 3, fitfunc, paraTrue=para2, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_para3 <- parApply(clust, plNoisy_para3, 3, fitfunc, paraTrue=para3, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_para4 <- parApply(clust, plNoisy_para4, 3, fitfunc, paraTrue=para4, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

# violin plots of different p0 concentrations comparing formula 1 and 2
fig09_p0_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_para1["f1_p0",], fitResult_para2["f1_p0",], fitResult_para3["f1_p0",], fitResult_para4["f1_p0",],
      fitResult_para1["f2_p0",], fitResult_para2["f2_p0",], fitResult_para3["f2_p0",], fitResult_para4["f2_p0",]),
  f=rep(formulaStringList, each=expN*length(p0StringList))
)
fig09_p0_violin <- ggplot(fig09_p0_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated [P]0 (nM)") +
  scale_x_discrete(limits=p0StringList, guide=guide_axis(n.dodge=2)) +
  theme(axis.title.x=element_blank())+
  legend_theme+
  theme(legend.position.inside=c(0, 1), legend.justification=c(0, 1))

fig09_kd_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_para1["f1_kd",], fitResult_para2["f1_kd",], fitResult_para3["f1_kd",], fitResult_para4["f1_kd",],
      fitResult_para1["f2_kd",], fitResult_para2["f2_kd",], fitResult_para3["f2_kd",], fitResult_para4["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(p0StringList))
)
fig09_kd_violin <- ggplot(fig09_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=p0StringList, guide=guide_axis(n.dodge=2)) +
  theme(axis.title.x=element_blank())+
  legend_theme+
  theme(legend.position.inside=c(0, 1), legend.justification=c(0, 1))

ggsave("p0.png", fig09_p0_violin, dpi=800, width=7, height=7, units='cm')
ggsave("p0.svg", fig09_p0_violin, width=7, height=7, units='cm')
ggsave("kd.png", fig09_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig09_kd_violin, width=7, height=7, units='cm')

write_table(filename="Output", listname=p0StringList, fitResult_para1, fitResult_para2, fitResult_para3, fitResult_para4)
setwd("..")

# * 1.4 Sampling Distributions of Estimated [P]0 and Kd -----------------------
# See Part "Random Seeds and expN"

# * 1.5-SI Concentration Errors in [P]0 and [L]0 -----------------------------

# * 1.5.1-SI Estimation of [P]0 and Kd ---------------------------------------

dir.create("Part1.5.1-SI"); setwd("Part1.5.1-SI")

gen_plNoisy_array2 <- function(tl_theo_list, paraTrue, repN, expN, sd_level, tp_cv, tl_cv) {
  
  tl_list_len <- length(tl_theo_list)
  p0 <- paraTrue$p0; kd <- paraTrue$kd
  
  # tp_cv selection
  if (tp_cv=="low"){
    tp_cv_num <- 0.025
  }
  if (tp_cv=="medium"){
    tp_cv_num <- 0.05
  }
  if (tp_cv=="high"){
    tp_cv_num <- 0.1
  }
  # tl_cv selection
  if (tl_cv=="low"){
    tl_cv_num <- 0.025
  }
  if (tl_cv=="medium"){
    tl_cv_num <- 0.05
  }
  if (tl_cv=="high"){
    tl_cv_num <- 0.1
  }
  
  error_sd_low      <- p0*0.025/2
  error_sd_medium   <- p0*0.05/2
  error_sd_high     <- p0*0.1/2
  
  # sd_level selection
  if (sd_level=="low"){
    error_sd <- error_sd_low
  }
  if (sd_level=="medium"){
    error_sd <- error_sd_medium
  }
  if (sd_level=="high"){
    error_sd <- error_sd_high
  }
  
  
  tl_real_list <- as.vector(replicate(expN, sapply(tl_theo_list, function(x) rnorm(repN, x, x*tl_cv_num))))
  tp_real_list <- rnorm(repN*tl_list_len*expN, p0, p0*tp_cv_num)
  pl_tl_real_list <- formula1(tl_real_list, tp_real_list, kd)
  pl_tl_real_array <- array(pl_tl_real_list, dim=c(repN, tl_list_len, expN))
  
  error_array <- array(rnorm(repN*tl_list_len*expN, sd=error_sd), dim=c(repN, tl_list_len, expN))
  pl_tl_real_array + error_array
}

l0NoisyStringList <- c("Noisy [P]0 and [L]0", "Noise-free [P]0 and [L]0")

plNoisy_tlNoisy_array <- gen_plNoisy_array2(tl_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="medium", tp_cv="medium", tl_cv="medium")

fitResult_tlNoisy <- parApply(clust, plNoisy_tlNoisy_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

# violin plots
res2_p0_df <- data.frame(
  x=rep(rep(l0NoisyStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_tlNoisy["f1_p0",], fitResultRatio05["f1_p0",],
      fitResult_tlNoisy["f2_p0",], fitResultRatio05["f2_p0",]),
  f=rep(formulaStringList, each=expN*length(l0NoisyStringList))
)
res2_p0_violin <- ggplot(res2_p0_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated [P]0 (nM)") +
  theme(axis.title.x=element_blank())+
  legend_theme

res2_kd_df <- data.frame(
  x=rep(rep(l0NoisyStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_tlNoisy["f1_kd",], fitResultRatio05["f1_kd",], 
      fitResult_tlNoisy["f2_kd",], fitResultRatio05["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(l0NoisyStringList))
)
res2_kd_violin <- ggplot(res2_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  theme(axis.title.x=element_blank())+
  legend_theme

ggsave("p0.png", res2_p0_violin, dpi=800, width=7, height=7, units='cm')
ggsave("p0.svg", res2_p0_violin, width=7, height=7, units='cm')
ggsave("kd.png", res2_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", res2_kd_violin, width=7, height=7, units='cm')

write_table(filename="Output", listname=l0NoisyStringList, fitResult_tlNoisy, fitResultRatio05)

setwd("..")

# *- 1.5.2-SI Sampling Distributions of Estimated [P]0 and Kd --------------------
# See Part "Random Seeds and expN"

# Part 2. Binding Kinetics Experiment ------------------------------------

# * 2.1 Estimation of koff and kon ---------

dir.create("Part2.1.1"); setwd("Part2.1.1")

time_points <- 10

formula_binding_kinetics <- function(t, l0, pl_eq, koff, kon) pl_eq * (1 - exp(0 - (kon * l0 + koff) * t)) 
formula_binding_kinetics_group <- function(t, l0, pl_eq, g, koff, kon)  pl_eq[g] * (1 - exp(0 - (kon * l0 + koff) * t)) 

formula_binding_kinetics_obs <- function(t, pl_eq, kobs) pl_eq * (1 - exp(0 - kobs * t)) 
formula_binding_kinetics_obs_group <- function(t, pl_eq, g, kobs) pl_eq[g] * (1 - exp(0 - kobs[g] * t)) 

formula_binding_kinetics_disso <- function(t, pl0, koff) pl0 * exp(0 - koff * t) 
formula_binding_kinetics_disso_group <- function(t, pl0, g, koff) pl0[g] * exp(0 - koff * t) 

clusterExport(clust, c("formula_binding_kinetics", "formula_binding_kinetics_group", 
                       "formula_binding_kinetics_obs", "formula_binding_kinetics_obs_group", 
                       "formula_binding_kinetics_disso", "formula_binding_kinetics_disso_group"))

l0_calc <- function(p0, kd) {
  points <- seq(0.2, 0.9, 0.1)
  pl <- points*p0
  l0 <- (p0+kd-pl)*pl/(p0-pl)
  l0
}

# the unit of p0, koff, and kon is nM, s-1, and nM-1 s-1
kinetics_para1 <- list(p0=100, koff=1e-2, kon=1e-4)
kinetics_para2 <- list(p0=50, koff=1e-2, kon=1e-4)
kinetics_para3 <- list(p0=10, koff=1e-2, kon=1e-4)

l0_list1 <- with(kinetics_para1, l0_calc(p0=p0, kd=koff/kon))
l0_list2 <- with(kinetics_para2, l0_calc(p0=p0, kd=koff/kon))
l0_list3 <- with(kinetics_para3, l0_calc(p0=p0, kd=koff/kon))

second_order_asso_points <- function(p0, l0, koff, kon, time_points=10) {
  a <- 1
  b <- 0 - (p0 + l0 + koff / kon)
  c <- p0 * l0
  pl_1 <- (0 - b - sqrt(b * b - 4 * a * c)) / (2 * a)  # pl_1 is [PL]eq, pl_1 < pl_2
  pl_2 <- (0 - b + sqrt(b * b - 4 * a * c)) / (2 * a)
  A <- 1 / (pl_1 - pl_2)
  
  num_array <- seq(0, 0.999, length.out=time_points)
  second_order_pl <- pl_1 * num_array
  second_order_time <- A / kon * log((pl_2 * (pl_1 - second_order_pl)) / (pl_1 * (pl_2 - second_order_pl)))
  
  num_array_disso <- rev(seq(0.001, 1, length.out=time_points))
  disso_pl <- pl_1 * num_array_disso
  disso_t <- (-1/koff)*log(disso_pl/pl_1)
  res <- list(second_order_pl=second_order_pl, second_order_time=second_order_time, 
              disso_pl=disso_pl, disso_t=disso_t)
  res
}

gen_plNoisy_array_kin <- function(pl_l0_list, l0_list, paraTrue, repN, expN, sd_level="low") {
  
  pl_l0_list_len <- length(pl_l0_list)
  t_len <- length(pl_l0_list[[1]]$second_order_time)
  
  p0 <- paraTrue$p0
  error_sd_low      <- p0*0.025/2
  error_sd_medium   <- p0*0.05/2
  error_sd_high     <- p0*0.1/2
  
  # sd_level selection
  if (sd_level=="low"){
    error_sd <- error_sd_low
  }
  if (sd_level=="medium"){
    error_sd <- error_sd_medium
  }
  if (sd_level=="high"){
    error_sd <- error_sd_high
  }
  
  pl1 <- as.vector(sapply(pl_l0_list, function(x) rep(x$second_order_pl, each=repN)))
  t1 <- as.vector(sapply(pl_l0_list, function(x) rep(x$second_order_time, each=repN)))
  l01 <- rep(l0_list, each=t_len*repN)
  group <- rep(1:pl_l0_list_len, each=t_len*repN)
  a1 <- array(c(pl1, t1, l01, group), dim=c(t_len*repN*pl_l0_list_len, 4))
  
  pl_disso <- as.vector(sapply(pl_l0_list, function(x) rep(x$disso_pl, each=repN)))
  t_disso <- as.vector(sapply(pl_l0_list, function(x) rep(x$disso_t, each=repN)))
  a_disso <- array(c(pl_disso, t_disso, group), dim=c(t_len*repN*pl_l0_list_len, 3))
  
  data_array_asso <- replicate(expN, a1)
  data_array_disso <- replicate(expN, a_disso)
  errors <- rnorm(t_len*repN*pl_l0_list_len*expN, sd=error_sd)
  error_disso <- rnorm(t_len*repN*pl_l0_list_len*expN, sd=error_sd)
  
  data_array_asso[,1,] <- data_array_asso[,1,] + errors
  data_array_disso[,1,] <- data_array_disso[,1,] + error_disso
  
  res <- list(data_array_asso=data_array_asso, data_array_disso=data_array_disso)
  res
}

fitfunc_kin <- function(plNoisy, paraTrue, paraRatio, type, l0_list, offres) {
  l0_num <- length(l0_list)
  
  p0 <- paraTrue$p0; kon <- paraTrue$kon; koff <- paraTrue$koff
  if (type=="asso") {
    colnames(plNoisy) <- c("pl", "t", "l0", "g")
    pldata <- as.data.frame(plNoisy)
    koff_est <- offres[["koff_est"]]
    if (l0_num==1) {
      warning("The prediction of koff and kon needs more than one l0")
    } else {
      # group guess for pl_eq is important setting.
      paraGuess <- list(pl_eq=rep(p0*paraRatio, times=l0_num), kon=kon*paraRatio)
      nls.asso <- nls(pl ~ formula_binding_kinetics_group(t, l0, pl_eq, g, koff=koff_est, kon), data=pldata, start=paraGuess)
    }
    coef_res <- coef(nls.asso)
    res <- coef_res[["kon"]]
  }
  if (type=="asso_obs") {
    colnames(plNoisy) <- c("pl", "t", "l0", "g")
    pldata <- as.data.frame(plNoisy)
    if (l0_num==1) {
      warning("The prediction of koff and kon needs more than one l0")
    } else {
      paraGuess <- list(pl_eq=rep(p0*paraRatio, times=l0_num), kobs=sapply(l0_list, function(x) (kon*x+koff)*paraRatio))
      nls.asso_obs <- nls(pl ~ formula_binding_kinetics_obs_group(t, pl_eq, g, kobs), data=pldata, start=paraGuess)
    }
    coef_res <- coef(nls.asso_obs)
    paranames <- sapply(1:l0_num, function(i) gsub(" ", "", paste("kobs", as.character(i))))
    res <- sapply(paranames, function(x) coef_res[[x]])
  }
  if (type=="disso") {
    colnames(plNoisy) <- c("pl", "t", "g")
    pldata <- as.data.frame(plNoisy)
    if (l0_num==1) {
      paraGuess <- list(pl0=p0*paraRatio, koff=koff*paraRatio)
      nls.disso <- nls(pl ~ formula_binding_kinetics_disso(t, pl0, koff), data=pldata, start=paraGuess)
    } else {
      paraGuess <- list(pl0=rep(p0*paraRatio, times=l0_num), koff=koff*paraRatio)
      nls.disso <- nls(pl ~ formula_binding_kinetics_disso_group(t, pl0, g, koff), data=pldata, start=paraGuess)
    }
    coef_res <- coef(nls.disso)
    res <- c(koff_est=coef_res[["koff"]], resi_sd=sd(residuals(nls.disso)) )
  }
  res
}

kon_koff_calc <- function(l0_list, data) {
  df <- data.frame(l0=l0_list, kobs=as.vector(data))
  res <- lm(kobs ~ l0, data=df)
  kon_est <- coef(res)[[2]]
  koff_est <- coef(res)[[1]]
  c(kon_est=kon_est, koff_est=koff_est)
}

kon_calc_koff_corrected <- function(l0_list, obsres, offres) {
  kobs_list <- as.vector(obsres)
  koff_est <- offres[["koff_est"]]
  kon_est <- mean((kobs_list-koff_est)/l0_list)
  kon_est
}

pl_l0_list1 <- lapply(l0_list1, function(x) with(kinetics_para1, second_order_asso_points(p0=p0, l0=x, koff=koff, kon=kon, time_points=time_points)))
pl_l0_list2 <- lapply(l0_list2, function(x) with(kinetics_para2, second_order_asso_points(p0=p0, l0=x, koff=koff, kon=kon, time_points=time_points)))
pl_l0_list3 <- lapply(l0_list3, function(x) with(kinetics_para3, second_order_asso_points(p0=p0, l0=x, koff=koff, kon=kon, time_points=time_points)))

plNoisy_array_kinetics1 <- gen_plNoisy_array_kin(pl_l0_list1, l0_list1, paraTrue=kinetics_para1, repN=repN, expN=expN, sd_level="low")
plNoisy_array_kinetics2 <- gen_plNoisy_array_kin(pl_l0_list2, l0_list2, paraTrue=kinetics_para2, repN=repN, expN=expN, sd_level="low")
plNoisy_array_kinetics3 <- gen_plNoisy_array_kin(pl_l0_list3, l0_list3, paraTrue=kinetics_para3, repN=repN, expN=expN, sd_level="low")

fitResult_kin_asso_obs1  <- parApply(clust, plNoisy_array_kinetics1$data_array_asso, 3, fitfunc_kin, paraTrue=kinetics_para1, paraRatio=ratioList[2], type="asso_obs", l0_list1)
fitResult_kin_disso1  <- parApply(clust, plNoisy_array_kinetics1$data_array_disso, 3, fitfunc_kin, paraTrue=kinetics_para1, paraRatio=ratioList[2], type="disso", l0_list1)

fitResult_kin_asso_obs2  <- parApply(clust, plNoisy_array_kinetics2$data_array_asso, 3, fitfunc_kin, paraTrue=kinetics_para2, paraRatio=ratioList[2], type="asso_obs", l0_list2)
fitResult_kin_disso2  <- parApply(clust, plNoisy_array_kinetics2$data_array_disso, 3, fitfunc_kin, paraTrue=kinetics_para2, paraRatio=ratioList[2], type="disso", l0_list2)

fitResult_kin_asso_obs3  <- parApply(clust, plNoisy_array_kinetics3$data_array_asso, 3, fitfunc_kin, paraTrue=kinetics_para3, paraRatio=ratioList[2], type="asso_obs", l0_list3)
fitResult_kin_disso3  <- parApply(clust, plNoisy_array_kinetics3$data_array_disso, 3, fitfunc_kin, paraTrue=kinetics_para3, paraRatio=ratioList[2], type="disso", l0_list3)

# Because parSapply needs lots of memory, I don't use it below
fitResult_kin_asso_obs_conv1 <- sapply(1:expN, function(x) kon_koff_calc(l0_list1, fitResult_kin_asso_obs1[,x]))
fitResult_kin_asso_obs_conv2 <- sapply(1:expN, function(x) kon_koff_calc(l0_list2, fitResult_kin_asso_obs2[,x]))
fitResult_kin_asso_obs_conv3 <- sapply(1:expN, function(x) kon_koff_calc(l0_list3, fitResult_kin_asso_obs3[,x]))

kon_linear_koff_corrected1 <- sapply(1:expN, function(x) kon_calc_koff_corrected(l0_list1, fitResult_kin_asso_obs1[,x], fitResult_kin_disso1[,x]))
kon_linear_koff_corrected2 <- sapply(1:expN, function(x) kon_calc_koff_corrected(l0_list2, fitResult_kin_asso_obs2[,x], fitResult_kin_disso2[,x]))
kon_linear_koff_corrected3 <- sapply(1:expN, function(x) kon_calc_koff_corrected(l0_list3, fitResult_kin_asso_obs3[,x], fitResult_kin_disso3[,x]))

kon_nonlinear_koff_corrected1 <- sapply(1:expN, function(x) fitfunc_kin(plNoisy_array_kinetics1$data_array_asso[,,x], paraTrue=kinetics_para1, paraRatio=ratioList[2], type="asso", l0_list1, fitResult_kin_disso1[,x]))
kon_nonlinear_koff_corrected2 <- sapply(1:expN, function(x) fitfunc_kin(plNoisy_array_kinetics2$data_array_asso[,,x], paraTrue=kinetics_para2, paraRatio=ratioList[2], type="asso", l0_list2, fitResult_kin_disso2[,x]))
kon_nonlinear_koff_corrected3 <- sapply(1:expN, function(x) fitfunc_kin(plNoisy_array_kinetics3$data_array_asso[,,x], paraTrue=kinetics_para3, paraRatio=ratioList[2], type="asso", l0_list3, fitResult_kin_disso3[,x]))

kd_kin1 <- fitResult_kin_disso1["koff_est",]/kon_nonlinear_koff_corrected1
kd_kin2 <- fitResult_kin_disso2["koff_est",]/kon_nonlinear_koff_corrected2
kd_kin3 <- fitResult_kin_disso3["koff_est",]/kon_nonlinear_koff_corrected3

para_kin_like1 <- list(p0=100, kd=100)
para_kin_like2 <- list(p0=50, kd=100)
para_kin_like3 <- list(p0=10, kd=100)
plTheoretical_num_large_kin_cmp1 <- formula1(tl_num_large, para_kin_like1$p0, para_kin_like1$kd)
plTheoretical_num_large_kin_cmp2 <- formula1(tl_num_large, para_kin_like2$p0, para_kin_like2$kd)
plTheoretical_num_large_kin_cmp3 <- formula1(tl_num_large, para_kin_like3$p0, para_kin_like3$kd)
plNoisy_array_kin_cmp1 <- gen_plNoisy_array(plTheoretical_num_large_kin_cmp1, paraTrue=para_kin_like1, repN=repN5, expN=expN, sd_level="low")
plNoisy_array_kin_cmp2 <- gen_plNoisy_array(plTheoretical_num_large_kin_cmp2, paraTrue=para_kin_like2, repN=repN5, expN=expN, sd_level="low")
plNoisy_array_kin_cmp3 <- gen_plNoisy_array(plTheoretical_num_large_kin_cmp3, paraTrue=para_kin_like3, repN=repN5, expN=expN, sd_level="low")
fitResult_kin_cmp1 <- parApply(clust, plNoisy_array_kin_cmp1, 3, fitfunc, paraTrue=para_kin_like1, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_kin_cmp2 <- parApply(clust, plNoisy_array_kin_cmp2, 3, fitfunc, paraTrue=para_kin_like2, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_kin_cmp3 <- parApply(clust, plNoisy_array_kin_cmp3, 3, fitfunc, paraTrue=para_kin_like3, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

p0StringList <- c("100", "50", "10")
kinStringList <- c("Asso_obs", "Disso")
# violin plot of koff
kin_asso_koff_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(kinStringList)),
  y=c(fitResult_kin_asso_obs_conv1["koff_est",], fitResult_kin_asso_obs_conv2["koff_est",], fitResult_kin_asso_obs_conv3["koff_est",],
      fitResult_kin_disso1["koff_est",], fitResult_kin_disso2["koff_est",], fitResult_kin_disso3["koff_est",]),
  f=rep(kinStringList, each=expN*length(p0StringList))
)
asso_koff_violin <- ggplot(kin_asso_koff_df, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Asso_obs", "Disso"),
                    values=c("Asso_obs"=c1_s1, "Disso"=c2_s1)) +
  labs(x="[P]0 (nM)", y="Estimated koff (s-1)") +
  scale_x_discrete(limits=p0StringList) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  legend_theme

# violin plot of kon
kinStringList <- c("Asso_obs", "Nonlinear (koff corrected)")
kin_asso_kon_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(kinStringList)),
  y=c(fitResult_kin_asso_obs_conv1["kon_est",], fitResult_kin_asso_obs_conv2["kon_est",], fitResult_kin_asso_obs_conv3["kon_est",],
      kon_nonlinear_koff_corrected1, kon_nonlinear_koff_corrected2, kon_nonlinear_koff_corrected3),
  f=rep(kinStringList, each=expN*length(p0StringList))
)
asso_kon_violin <- ggplot(kin_asso_kon_df, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Asso_obs", "Nonlinear (koff corrected)"),
                    values=c("Asso_obs"=c1_s1, "Nonlinear (koff corrected)"=c2_s1)) +
  labs(x="[P]0 (nM)", y="Estimated kon (nM-1 s-1)") +
  scale_x_discrete(limits=p0StringList) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  coord_cartesian(ylim=c(7.8e-5, 1.2e-4)) +
  legend_theme

# violin plot of kon
kinStringList <- c("Linear (koff corrected)", "Nonlinear (koff corrected)")
kin_asso_kon_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(kinStringList)),
  y=c(kon_linear_koff_corrected1, kon_linear_koff_corrected2, kon_linear_koff_corrected3,
      kon_nonlinear_koff_corrected1, kon_nonlinear_koff_corrected2, kon_nonlinear_koff_corrected3),
  f=rep(kinStringList, each=expN*length(p0StringList))
)
asso_kon_violin2 <- ggplot(kin_asso_kon_df, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Linear (koff corrected)", "Nonlinear (koff corrected)"),
                    values=c("Linear (koff corrected)"=c1_s1, "Nonlinear (koff corrected)"=c2_s1)) +
  labs(x="[P]0 (nM)", y="Estimated kon (nM-1 s-1)") +
  scale_x_discrete(limits=p0StringList) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  legend_theme

# violin plot of kd
kinStringList <- c("Equilibrium", "Kinetics")
kin_kd_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(kinStringList)),
  y=c(fitResult_kin_cmp1["f1_kd",], fitResult_kin_cmp2["f1_kd",], fitResult_kin_cmp3["f1_kd",], 
      kd_kin1, kd_kin2, kd_kin3),
  f=rep(kinStringList, each=expN*length(p0StringList))
)
kin_kd_violin <- ggplot(kin_kd_df, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equilibrium", "Kinetics"),
                    values=c("Equilibrium"=c1_s1, "Kinetics"=c2_s1)) +
  labs(x="[P]0 (nM)", y="Estimated Kd (nM)") +
  scale_x_discrete(limits=p0StringList) +
  legend_theme

ggsave("koff.png", asso_koff_violin, dpi=800, width=7, height=7, units='cm')
ggsave("koff.svg", asso_koff_violin, width=7, height=7, units='cm')
ggsave("kon.png", asso_kon_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kon.svg", asso_kon_violin, width=7, height=7, units='cm')
ggsave("kon2.png", asso_kon_violin2, dpi=800, width=7, height=7, units='cm')
ggsave("kon2.svg", asso_kon_violin2, width=7, height=7, units='cm')
ggsave("kd.png", kin_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", kin_kd_violin, width=7, height=7, units='cm')

setwd("..")

dir.create("Part2.1.2"); setwd("Part2.1.2")

write_table_kin <- function(filename_base, stringList, ...) {
#func for fitResult_kin_asso_obs_conv1,2,3
  csvfilename <- paste(filename_base, ".csv", sep="")

  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("koff (SD)", "kon (SD)", 
                   "KS test P-value of koff", "KS test P-value of log10(koff)", 
                   "KS test P-value of kon", "KS test P-value of log10(kon)", 
                   "Pearson's correlation coefficient of koff and kon", "Spearman's correlation coefficient of koff and kon")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    koff_est <- arg2[[i]]["koff_est",]; kon_est <- arg2[[i]]["kon_est",]
    koff_est_mean <- mean(koff_est); kon_est_mean <- mean(kon_est)
    koff_est_sd <- sd(koff_est); kon_est_sd <- sd(kon_est)
    
    log10_koff_est <- log10(koff_est); log10_kon_est <- log10(kon_est)
    log10_koff_est_mean <- mean(log10_koff_est); log10_kon_est_mean <- mean(log10_kon_est)
    log10_koff_est_sd <- sd(log10_koff_est); log10_kon_est_sd <- sd(log10_kon_est)
    
    koff_est_ks_pvalue <- cont_ks_test(koff_est, "pnorm", koff_est_mean, koff_est_sd)$p.value
    kon_est_ks_pvalue <- cont_ks_test(kon_est, "pnorm", kon_est_mean, kon_est_sd)$p.value
    
    log10_koff_est_ks_pvalue <- cont_ks_test(log10_koff_est, "pnorm", log10_koff_est_mean, log10_koff_est_sd)$p.value
    log10_kon_est_ks_pvalue <- cont_ks_test(log10_kon_est, "pnorm", log10_kon_est_mean, log10_kon_est_sd)$p.value
    
    koff_est_kon_rp <- cor(koff_est, kon_est, method="pearson")
    koff_est_kon_rs <- cor(koff_est, kon_est, method="spearman")
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", koff_est_mean, koff_est_sd), sprintf("%.2e (%.2e)", kon_est_mean, kon_est_sd),
                   sprintf("%.3f", koff_est_ks_pvalue), sprintf("%.3f", log10_koff_est_ks_pvalue),
                   sprintf("%.3f", kon_est_ks_pvalue), sprintf("%.3f", log10_kon_est_ks_pvalue),
                   sprintf("%.3f", koff_est_kon_rp), sprintf("%.3f", koff_est_kon_rs)
                   )
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_kin_disso <- function(filename_base, stringList, ...) {
#func for fitResult_kin_disso1,2,3
  csvfilename <- paste(filename_base, ".csv", sep="")

  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("koff (SD)", "KS test P-value of koff", "KS test P-value of log10(koff)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    koff_est <- arg2[[i]]["koff_est",]
    koff_est_mean <- mean(koff_est)
    koff_est_sd <- sd(koff_est); 
    
    log10_koff_est <- log10(koff_est)
    log10_koff_est_mean <- mean(log10_koff_est)
    log10_koff_est_sd <- sd(log10_koff_est)
    
    koff_est_ks_pvalue <- cont_ks_test(koff_est, "pnorm", koff_est_mean, koff_est_sd)$p.value
    
    log10_koff_est_ks_pvalue <- cont_ks_test(log10_koff_est, "pnorm", log10_koff_est_mean, log10_koff_est_sd)$p.value
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", koff_est_mean, koff_est_sd),
                   sprintf("%.3f", koff_est_ks_pvalue), sprintf("%.3f", log10_koff_est_ks_pvalue) )
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_kin_asso <- function(filename_base, stringList, koffList, ...) {
#func for kon_linear_koff_corrected1,2,3 and kon_nonlinear_koff_corrected1,2,3
  csvfilename <- paste(filename_base, ".csv", sep="")
  
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("kon (SD)", "KS test P-value of kon", "KS test P-value of log10(kon)", 
                   "Pearson's correlation coefficient of koff and kon", "Spearman's correlation coefficient of koff and kon")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    koff_est <- koffList[[i]]["koff_est",]; kon_est <- arg2[[i]]
    kon_est_mean <- mean(kon_est); kon_est_sd <- sd(kon_est)
    
    log10_kon_est <- log10(kon_est)
    log10_kon_est_mean <- mean(log10_kon_est); log10_kon_est_sd <- sd(log10_kon_est)
    
    kon_est_ks_pvalue <- cont_ks_test(kon_est, "pnorm", kon_est_mean, kon_est_sd)$p.value
    
    log10_kon_est_ks_pvalue <- cont_ks_test(log10_kon_est, "pnorm", log10_kon_est_mean, log10_kon_est_sd)$p.value
    
    koff_est_kon_rp <- cor(koff_est, kon_est, method="pearson")
    koff_est_kon_rs <- cor(koff_est, kon_est, method="spearman")
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", kon_est_mean, kon_est_sd),
                      sprintf("%.3f", kon_est_ks_pvalue), sprintf("%.3f", log10_kon_est_ks_pvalue),
                      sprintf("%.3f", koff_est_kon_rp), sprintf("%.3f", koff_est_kon_rs)
    )
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_kin_kd <- function(filename_base, stringList, type, ...) {
#func for fitResult_kin_cmp1,2,3 and kd_kin1,2,3
  csvfilename <- paste(filename_base, ".csv", sep="")
  
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("Kd (SD)", "KS test P-value of Kd", "KS test P-value of log10(Kd)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    if (type=="equi") {
      kd_est <- arg2[[i]]["f1_kd",]
    }
    if (type=="kin") {
      kd_est <- arg2[[i]]
    }
    kd_est_mean <- mean(kd_est); kd_est_sd <- sd(kd_est)
    
    log10_kd_est <- log10(kd_est)
    log10_kd_est_mean <- mean(log10_kd_est); log10_kd_est_sd <- sd(log10_kd_est)
    
    kd_est_ks_pvalue <- cont_ks_test(kd_est, "pnorm", kd_est_mean, kd_est_sd)$p.value
    
    log10_kd_est_ks_pvalue <- cont_ks_test(log10_kd_est, "pnorm", log10_kd_est_mean, log10_kd_est_sd)$p.value
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", kd_est_mean, kd_est_sd),
                      sprintf("%.3f", kd_est_ks_pvalue), sprintf("%.3f", log10_kd_est_ks_pvalue))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

kinOutString <- c("[P]0=100nM", "[P]0=50nM", "[P]0=10nM")
kinkoffList <- list(fitResult_kin_disso1, fitResult_kin_disso2, fitResult_kin_disso3)
write_table_kin(filename_base="Output.1.asso_obs", stringList=kinOutString, fitResult_kin_asso_obs_conv1, fitResult_kin_asso_obs_conv2, fitResult_kin_asso_obs_conv3)
write_table_kin_disso(filename_base="Output.2.koff", stringList=kinOutString, fitResult_kin_disso1, fitResult_kin_disso2, fitResult_kin_disso3)
write_table_kin_asso(filename_base="Output.3.kon_linear", stringList=kinOutString, kinkoffList, kon_linear_koff_corrected1, kon_linear_koff_corrected2, kon_linear_koff_corrected3)
write_table_kin_asso(filename_base="Output.4.kon_nonlinear", stringList=kinOutString, kinkoffList, kon_nonlinear_koff_corrected1, kon_nonlinear_koff_corrected2, kon_nonlinear_koff_corrected3)
write_table_kin_kd(filename_base="Output.5.kd-equi", stringList=kinOutString, type="equi", fitResult_kin_cmp1, fitResult_kin_cmp2, fitResult_kin_cmp3)
write_table_kin_kd(filename_base="Output.7.kd-kin", stringList=kinOutString, type="kin", kd_kin1, kd_kin2, kd_kin3)
setwd("..")

# * 2.2 Sampling Distributions of Estimated koff and kon  --------
# See Part "Random Seeds and expN"

# Part 3. Competitive Binding Experiment ----------------------------------

# * 3.1 Estimation of Apparent IC50, Apparent Ki, and Ki ------------------

dir.create("Part3.1-1"); setwd("Part3.1-1")

# logTi_num_large <- seq(-9, -4, length=num_large) 
# if logTi range [-9, -4] is used, ki estimated by direct fitting under the condition of para_comp4
# list(p0=5, l0=c(25,50,75,100),  kd=5, ki=100) will show strong numerical instability
logTi_num_large <- seq(-9, -5, length=num_large) 
ti_num_large <- 10^logTi_num_large * (1e+9)

para_comp1_l1 <- list(p0=50, l0=25,  kd=50, ki=100)
para_comp1_l2 <- list(p0=50, l0=50,  kd=50, ki=100)
para_comp1_l3 <- list(p0=50, l0=75,  kd=50, ki=100)
para_comp1_l4 <- list(p0=50, l0=100, kd=50, ki=100)

para_comp2_l1 <- list(p0=50, l0=25,  kd=5, ki=100)
para_comp2_l2 <- list(p0=50, l0=50,  kd=5, ki=100)
para_comp2_l3 <- list(p0=50, l0=75,  kd=5, ki=100)
para_comp2_l4 <- list(p0=50, l0=100, kd=5, ki=100)

para_comp3_l1 <- list(p0=5, l0=25,  kd=50, ki=100)
para_comp3_l2 <- list(p0=5, l0=50,  kd=50, ki=100)
para_comp3_l3 <- list(p0=5, l0=75,  kd=50, ki=100)
para_comp3_l4 <- list(p0=5, l0=100, kd=50, ki=100)

para_comp4_l1 <- list(p0=5, l0=25,  kd=5, ki=100)
para_comp4_l2 <- list(p0=5, l0=50,  kd=5, ki=100)
para_comp4_l3 <- list(p0=5, l0=75,  kd=5, ki=100)
para_comp4_l4 <- list(p0=5, l0=100, kd=5, ki=100)

inhibition_calc <- function(p0, l0, kd, ki, returnlist=FALSE) {
  # [I] = [PI]*Ki/[P], for the derivation, see paper's SI 
  # Part "Calculation of Theoretical IC50 and Apparent IC50".
  # At 50% inhibition, [PI] and [P] can be determined by Kd, 
  # [PL]eq-50, [P]0, and [L]0. As a result, [PI] and [P] 
  # isn't related to Ki, IC50 is proportional to Ki, 
  # i.e. 0.9Ki => 0.9IC50. But apparent IC50 (i.e., [I]0) is 
  # not proportional to Ki
  res <- list()
  a <- 1
  b <- 0 - (p0 + l0 + kd)
  c <- p0 * l0
  pl_eq_0 <- (0 - b - sqrt(b * b - 4 * a * c)) / (2 * a) 
  # At 50% inhibition
  pl_eq_50 <- pl_eq_0 / 2
  p <- pl_eq_50 * kd / (l0 - pl_eq_50)
  pi_eq <- p0 - p - pl_eq_50
  app_IC50 <- pi_eq * (ki + p) / p # i0
  IC50 <- app_IC50 - pi_eq
  
  if (returnlist) {
    # inhibition list
    inhibition_ratio <- c(seq(0.9, 0.1, -0.1), 0.05)
    pl_eq_list <- pl_eq_0 * inhibition_ratio
    p_list <- pl_eq_list * kd / (l0 - pl_eq_list)
    pi_eq_list <- p0 - p_list - pl_eq_list
    i0_list <- pi_eq_list * (ki + p_list) / p_list
    res$i0_list <- i0_list
  }
  res$app_IC50 <- app_IC50
  res$IC50 <- IC50
  res
}

formula_competitive_pl <- function(p0, l0, i0, kd, ki) {
  # The units of p0, l0, i0, kd and ki are nM
  # Equations from FEBS Lett. 1995, 360 (2), 111-4.
  a <- kd + ki + l0 + i0 - p0
  b <- ki * (l0 - p0) + kd * (i0 - p0) + kd * ki
  c <- -kd * ki * p0
  theta <- acos((-2 * a ** 3 + 9 * a * b - 27 * c) / (2 * sqrt((a ** 2 - 3 * b) ** 3)))
  p_eq <- -a / 3 + (2 / 3) * sqrt(a ** 2 - 3 * b) * cos(theta / 3)
  pl_eq <- p_eq * l0 / (kd + p_eq)
  pl_eq
}

calc_ki_app1 <- function(l0, i0, kd) {
  # i0 approximates the i_eq, l0 approximates the l_eq
  # Cheng-Prusoff equation
  # Biochem. Pharmacol. 1973, 22 (23), 3099-3108.
  ki_app1 <- i0 / (1 + l0 / kd)
  ki_app1
}

calc_ki <- function(p0, l0, i0, kd) {
  # use the solutions in Anal. Biochem. 2004, 332 (2), 261-73.
  # p0 is [P]total, l0 is [L]total, i0 is [I]total
  a <- 1
  b <- 0 - (p0 + l0 + kd)
  c <- p0 * l0
  pl_eq0 <- (0 - b - sqrt(b * b - 4 * a * c)) / (2 * a) 
  
  P0 <- p0 - pl_eq0 # p0 here equals [P]total
  PL50 <- pl_eq0 / 2
  L50 <- l0 - PL50
  I50 <- i0 - p0 + kd * PL50 / L50 + PL50
  ki <- I50 / (L50 / kd + P0 / kd + 1)
  ki
}

logistic_func <- function(i0, a, b, c, d) (a - d) / (1 + (i0 / c) ** b) + d

fitfunc_app_ic50 <- function(plNoisy, paraTrue, repN, i_conc, use_drc=FALSE) {
  # reshape the wide format matrix to long format data frame 
  pldata <- data.frame(pl=as.vector(plNoisy), i0=rep(i_conc, each=repN))
  d_guess <- mean(plNoisy[,1])
  paraGuess <- list(a=0.01, b=-1, c=100, d=d_guess)
  res <- list()
  if (use_drc) {
    nls.logistic_drc <- drm(pl ~ i0, data=pldata, fct=LL.4())
    coef_logistic_drc <- coef(nls.logistic_drc)
    res$a <- coef_logistic_drc[["c:(Intercept)"]]; res$b <- coef_logistic_drc[["b:(Intercept)"]]
    res$c <- coef_logistic_drc[["e:(Intercept)"]]; res$d <- coef_logistic_drc[["d:(Intercept)"]]
  } else {
    nls.logistic<-NULL
    while(TRUE){
      try(nls.logistic <- nlsLM(pl ~ logistic_func(i0, a, b, c, d), data=pldata, start=paraGuess))
      if (!is.null(nls.logistic)) break
      paraGuess$a <- paraGuess$a/2
    }
    coef_logistic <- coef(nls.logistic)
    res$a <- coef_logistic[["a"]]; res$b <- coef_logistic[["b"]]
    res$c <- coef_logistic[["c"]]; res$d <- coef_logistic[["d"]]
  }
  c(a=res$a, b=res$b, c=res$c, d=res$d)
}

fitfunc_ki <- function(plNoisy, paraTrue, paraRatio, repN, i_conc) {
  # reshape the wide format matrix to long format data frame 
  pldata <- data.frame(pl=as.vector(plNoisy), i0=rep(i_conc, each=repN))
  l0 <- paraTrue$l0; kd <- paraTrue$kd
  paraGuess <- list(p0=paraTrue$p0*paraRatio, ki=paraTrue$ki*paraRatio)

  # nls can't use formula_competitive_pl to fit data and produce error:
  # step factor 0.000244141 reduced below 'minFactor' of 0.000488281
  # nls.comp <- nls(pl ~ formula_competitive_pl(p0, l0=l0, i0, kd=kd, ki), data=pldata, start=paraGuess, control=list(minFactor=1/2048))

  nls.comp <- nlsLM(pl ~ formula_competitive_pl(p0, l0=l0, i0, kd=kd, ki), data=pldata, start=paraGuess)

  coef_comp <- coef(nls.comp)
  c(p0=coef_comp[["p0"]], ki=coef_comp[["ki"]])
}

clusterExport(clust, c("formula_competitive_pl", "logistic_func", "drm", "LL.4"))

# theo. data
plTheo_comp1_l1 <- with(para_comp1_l1, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp1_l2 <- with(para_comp1_l2, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp1_l3 <- with(para_comp1_l3, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp1_l4 <- with(para_comp1_l4, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))

plTheo_comp2_l1 <- with(para_comp2_l1, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp2_l2 <- with(para_comp2_l2, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp2_l3 <- with(para_comp2_l3, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp2_l4 <- with(para_comp2_l4, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))

plTheo_comp3_l1 <- with(para_comp3_l1, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp3_l2 <- with(para_comp3_l2, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp3_l3 <- with(para_comp3_l3, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp3_l4 <- with(para_comp3_l4, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))

plTheo_comp4_l1 <- with(para_comp4_l1, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp4_l2 <- with(para_comp4_l2, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp4_l3 <- with(para_comp4_l3, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp4_l4 <- with(para_comp4_l4, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))

# exp. data
plNoisy_comp1_l1 <- gen_plNoisy_array(plTheo_comp1_l1, paraTrue=para_comp1_l1, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp1_l2 <- gen_plNoisy_array(plTheo_comp1_l2, paraTrue=para_comp1_l2, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp1_l3 <- gen_plNoisy_array(plTheo_comp1_l3, paraTrue=para_comp1_l3, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp1_l4 <- gen_plNoisy_array(plTheo_comp1_l4, paraTrue=para_comp1_l4, repN=repN5, expN=expN, sd_level="low")

plNoisy_comp2_l1 <- gen_plNoisy_array(plTheo_comp2_l1, paraTrue=para_comp2_l1, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp2_l2 <- gen_plNoisy_array(plTheo_comp2_l2, paraTrue=para_comp2_l2, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp2_l3 <- gen_plNoisy_array(plTheo_comp2_l3, paraTrue=para_comp2_l3, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp2_l4 <- gen_plNoisy_array(plTheo_comp2_l4, paraTrue=para_comp2_l4, repN=repN5, expN=expN, sd_level="low")

plNoisy_comp3_l1 <- gen_plNoisy_array(plTheo_comp3_l1, paraTrue=para_comp3_l1, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp3_l2 <- gen_plNoisy_array(plTheo_comp3_l2, paraTrue=para_comp3_l2, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp3_l3 <- gen_plNoisy_array(plTheo_comp3_l3, paraTrue=para_comp3_l3, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp3_l4 <- gen_plNoisy_array(plTheo_comp3_l4, paraTrue=para_comp3_l4, repN=repN5, expN=expN, sd_level="low")

plNoisy_comp4_l1 <- gen_plNoisy_array(plTheo_comp4_l1, paraTrue=para_comp4_l1, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp4_l2 <- gen_plNoisy_array(plTheo_comp4_l2, paraTrue=para_comp4_l2, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp4_l3 <- gen_plNoisy_array(plTheo_comp4_l3, paraTrue=para_comp4_l3, repN=repN5, expN=expN, sd_level="low")
plNoisy_comp4_l4 <- gen_plNoisy_array(plTheo_comp4_l4, paraTrue=para_comp4_l4, repN=repN5, expN=expN, sd_level="low")

# fitted ki
fitki_comp1_l1 <- parApply(clust, plNoisy_comp1_l1, 3, fitfunc_ki, paraTrue=para_comp1_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp1_l2 <- parApply(clust, plNoisy_comp1_l2, 3, fitfunc_ki, paraTrue=para_comp1_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp1_l3 <- parApply(clust, plNoisy_comp1_l3, 3, fitfunc_ki, paraTrue=para_comp1_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp1_l4 <- parApply(clust, plNoisy_comp1_l4, 3, fitfunc_ki, paraTrue=para_comp1_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)

fitki_comp2_l1 <- parApply(clust, plNoisy_comp2_l1, 3, fitfunc_ki, paraTrue=para_comp2_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp2_l2 <- parApply(clust, plNoisy_comp2_l2, 3, fitfunc_ki, paraTrue=para_comp2_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp2_l3 <- parApply(clust, plNoisy_comp2_l3, 3, fitfunc_ki, paraTrue=para_comp2_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp2_l4 <- parApply(clust, plNoisy_comp2_l4, 3, fitfunc_ki, paraTrue=para_comp2_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)

fitki_comp3_l1 <- parApply(clust, plNoisy_comp3_l1, 3, fitfunc_ki, paraTrue=para_comp3_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp3_l2 <- parApply(clust, plNoisy_comp3_l2, 3, fitfunc_ki, paraTrue=para_comp3_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp3_l3 <- parApply(clust, plNoisy_comp3_l3, 3, fitfunc_ki, paraTrue=para_comp3_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp3_l4 <- parApply(clust, plNoisy_comp3_l4, 3, fitfunc_ki, paraTrue=para_comp3_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)

fitki_comp4_l1 <- parApply(clust, plNoisy_comp4_l1, 3, fitfunc_ki, paraTrue=para_comp4_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp4_l2 <- parApply(clust, plNoisy_comp4_l2, 3, fitfunc_ki, paraTrue=para_comp4_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp4_l3 <- parApply(clust, plNoisy_comp4_l3, 3, fitfunc_ki, paraTrue=para_comp4_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp4_l4 <- parApply(clust, plNoisy_comp4_l4, 3, fitfunc_ki, paraTrue=para_comp4_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)

# fitted results
fitres_comp1_l1 <- parApply(clust, plNoisy_comp1_l1, 3, fitfunc_app_ic50, paraTrue=para_comp1_l1, repN=repN5, i_conc=ti_num_large)
fitres_comp1_l2 <- parApply(clust, plNoisy_comp1_l2, 3, fitfunc_app_ic50, paraTrue=para_comp1_l2, repN=repN5, i_conc=ti_num_large)
fitres_comp1_l3 <- parApply(clust, plNoisy_comp1_l3, 3, fitfunc_app_ic50, paraTrue=para_comp1_l3, repN=repN5, i_conc=ti_num_large)
fitres_comp1_l4 <- parApply(clust, plNoisy_comp1_l4, 3, fitfunc_app_ic50, paraTrue=para_comp1_l4, repN=repN5, i_conc=ti_num_large)

fitres_comp2_l1 <- parApply(clust, plNoisy_comp2_l1, 3, fitfunc_app_ic50, paraTrue=para_comp2_l1, repN=repN5, i_conc=ti_num_large)
fitres_comp2_l2 <- parApply(clust, plNoisy_comp2_l2, 3, fitfunc_app_ic50, paraTrue=para_comp2_l2, repN=repN5, i_conc=ti_num_large)
fitres_comp2_l3 <- parApply(clust, plNoisy_comp2_l3, 3, fitfunc_app_ic50, paraTrue=para_comp2_l3, repN=repN5, i_conc=ti_num_large)
fitres_comp2_l4 <- parApply(clust, plNoisy_comp2_l4, 3, fitfunc_app_ic50, paraTrue=para_comp2_l4, repN=repN5, i_conc=ti_num_large)

fitres_comp3_l1 <- parApply(clust, plNoisy_comp3_l1, 3, fitfunc_app_ic50, paraTrue=para_comp3_l1, repN=repN5, i_conc=ti_num_large)
fitres_comp3_l2 <- parApply(clust, plNoisy_comp3_l2, 3, fitfunc_app_ic50, paraTrue=para_comp3_l2, repN=repN5, i_conc=ti_num_large)
fitres_comp3_l3 <- parApply(clust, plNoisy_comp3_l3, 3, fitfunc_app_ic50, paraTrue=para_comp3_l3, repN=repN5, i_conc=ti_num_large)
fitres_comp3_l4 <- parApply(clust, plNoisy_comp3_l4, 3, fitfunc_app_ic50, paraTrue=para_comp3_l4, repN=repN5, i_conc=ti_num_large)

fitres_comp4_l1 <- parApply(clust, plNoisy_comp4_l1, 3, fitfunc_app_ic50, paraTrue=para_comp4_l1, repN=repN5, i_conc=ti_num_large)
fitres_comp4_l2 <- parApply(clust, plNoisy_comp4_l2, 3, fitfunc_app_ic50, paraTrue=para_comp4_l2, repN=repN5, i_conc=ti_num_large)
fitres_comp4_l3 <- parApply(clust, plNoisy_comp4_l3, 3, fitfunc_app_ic50, paraTrue=para_comp4_l3, repN=repN5, i_conc=ti_num_large)
fitres_comp4_l4 <- parApply(clust, plNoisy_comp4_l4, 3, fitfunc_app_ic50, paraTrue=para_comp4_l4, repN=repN5, i_conc=ti_num_large)

# fitted apparent IC50
fitappic50_comp1_l1 <- fitres_comp1_l1["c",]
fitappic50_comp1_l2 <- fitres_comp1_l2["c",]
fitappic50_comp1_l3 <- fitres_comp1_l3["c",]
fitappic50_comp1_l4 <- fitres_comp1_l4["c",]

fitappic50_comp2_l1 <- fitres_comp2_l1["c",]
fitappic50_comp2_l2 <- fitres_comp2_l2["c",]
fitappic50_comp2_l3 <- fitres_comp2_l3["c",]
fitappic50_comp2_l4 <- fitres_comp2_l4["c",]

fitappic50_comp3_l1 <- fitres_comp3_l1["c",]
fitappic50_comp3_l2 <- fitres_comp3_l2["c",]
fitappic50_comp3_l3 <- fitres_comp3_l3["c",]
fitappic50_comp3_l4 <- fitres_comp3_l4["c",]

fitappic50_comp4_l1 <- fitres_comp4_l1["c",]
fitappic50_comp4_l2 <- fitres_comp4_l2["c",]
fitappic50_comp4_l3 <- fitres_comp4_l3["c",]
fitappic50_comp4_l4 <- fitres_comp4_l4["c",]

# theoretical apparent IC50
theoappic50_comp1_l1 <- with(para_comp1_l1, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp1_l2 <- with(para_comp1_l2, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp1_l3 <- with(para_comp1_l3, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp1_l4 <- with(para_comp1_l4, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
                        
theoappic50_comp2_l1 <- with(para_comp2_l1, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp2_l2 <- with(para_comp2_l2, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp2_l3 <- with(para_comp2_l3, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp2_l4 <- with(para_comp2_l4, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
                        
theoappic50_comp3_l1 <- with(para_comp3_l1, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp3_l2 <- with(para_comp3_l2, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp3_l3 <- with(para_comp3_l3, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp3_l4 <- with(para_comp3_l4, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
                        
theoappic50_comp4_l1 <- with(para_comp4_l1, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp4_l2 <- with(para_comp4_l2, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp4_l3 <- with(para_comp4_l3, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp4_l4 <- with(para_comp4_l4, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50

# ratio of fitted to theoretical apparent IC50
diffappic50_comp1_l1 <- fitappic50_comp1_l1 / theoappic50_comp1_l1
diffappic50_comp1_l2 <- fitappic50_comp1_l2 / theoappic50_comp1_l2
diffappic50_comp1_l3 <- fitappic50_comp1_l3 / theoappic50_comp1_l3
diffappic50_comp1_l4 <- fitappic50_comp1_l4 / theoappic50_comp1_l4

diffappic50_comp2_l1 <- fitappic50_comp2_l1 / theoappic50_comp2_l1
diffappic50_comp2_l2 <- fitappic50_comp2_l2 / theoappic50_comp2_l2
diffappic50_comp2_l3 <- fitappic50_comp2_l3 / theoappic50_comp2_l3
diffappic50_comp2_l4 <- fitappic50_comp2_l4 / theoappic50_comp2_l4

diffappic50_comp3_l1 <- fitappic50_comp3_l1 / theoappic50_comp3_l1
diffappic50_comp3_l2 <- fitappic50_comp3_l2 / theoappic50_comp3_l2
diffappic50_comp3_l3 <- fitappic50_comp3_l3 / theoappic50_comp3_l3
diffappic50_comp3_l4 <- fitappic50_comp3_l4 / theoappic50_comp3_l4

diffappic50_comp4_l1 <- fitappic50_comp4_l1 / theoappic50_comp4_l1
diffappic50_comp4_l2 <- fitappic50_comp4_l2 / theoappic50_comp4_l2
diffappic50_comp4_l3 <- fitappic50_comp4_l3 / theoappic50_comp4_l3
diffappic50_comp4_l4 <- fitappic50_comp4_l4 / theoappic50_comp4_l4

# ref: https://forum.posit.co/t/how-properly-use-scale-fill-manual-and-its-arguments-labels-breaks-limits-colours-being-messed-up/177618/2
diffappic50_df1 <- data.frame(ki=c(diffappic50_comp1_l1, diffappic50_comp2_l1, diffappic50_comp3_l1, diffappic50_comp4_l1), 
                           type=rep(c("A", "B", "C", "D"), each=expN), 
                           l0=para_comp1_l1$l0)
diffappic50_df2 <- data.frame(ki=c(diffappic50_comp1_l2, diffappic50_comp2_l2, diffappic50_comp3_l2, diffappic50_comp4_l2), 
                           type=rep(c("A", "B", "C", "D"), each=expN), 
                           l0=para_comp1_l2$l0)
diffappic50_df3 <- data.frame(ki=c(diffappic50_comp1_l3, diffappic50_comp2_l3, diffappic50_comp3_l3, diffappic50_comp4_l3), 
                           type=rep(c("A", "B", "C", "D"), each=expN), 
                           l0=para_comp1_l3$l0)
diffappic50_df4 <- data.frame(ki=c(diffappic50_comp1_l4, diffappic50_comp2_l4, diffappic50_comp3_l4, diffappic50_comp4_l4), 
                           type=rep(c("A", "B", "C", "D"), each=expN), 
                           l0=para_comp1_l4$l0)

diffappic50_df <- rbind(diffappic50_df1, diffappic50_df2, diffappic50_df3, diffappic50_df4)
diffappic50_df$l0 <- as.factor(diffappic50_df$l0)
diffappic50_violin <- ggplot(diffappic50_df, aes(x=l0, y=ki, fill=type)) +
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count", position=position_dodge(1)) +
  scale_fill_manual(breaks=c("A", "B", "C", "D"),
                    labels=c("[P]0 = 50 nM, Kd = 50 nM", "[P]0 = 50 nM, Kd = 5 nM", "[P]0 = 5 nM, Kd = 50 nM", "[P]0 = 5 nM, Kd = 5 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], 
                             "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[L]0 (nM)", y="Estimated app. IC50 / Theoretical app. IC50") + 
  legend_theme

ggsave("diffappic50_violin.png", diffappic50_violin, dpi=800, width=7, height=7, units='cm')
ggsave("diffappic50_violin.svg", diffappic50_violin, width=7, height=7, units='cm')

# ki calculated by cheng-prusoff equation
ki_cheng_comp1_l1 <- calc_ki_app1(l0=para_comp1_l1$l0, i0=fitappic50_comp1_l1, kd=para_comp1_l1$kd)
ki_cheng_comp1_l2 <- calc_ki_app1(l0=para_comp1_l2$l0, i0=fitappic50_comp1_l2, kd=para_comp1_l2$kd)
ki_cheng_comp1_l3 <- calc_ki_app1(l0=para_comp1_l3$l0, i0=fitappic50_comp1_l3, kd=para_comp1_l3$kd)
ki_cheng_comp1_l4 <- calc_ki_app1(l0=para_comp1_l4$l0, i0=fitappic50_comp1_l4, kd=para_comp1_l4$kd)

ki_cheng_comp2_l1 <- calc_ki_app1(l0=para_comp2_l1$l0, i0=fitappic50_comp2_l1, kd=para_comp2_l1$kd)
ki_cheng_comp2_l2 <- calc_ki_app1(l0=para_comp2_l2$l0, i0=fitappic50_comp2_l2, kd=para_comp2_l2$kd)
ki_cheng_comp2_l3 <- calc_ki_app1(l0=para_comp2_l3$l0, i0=fitappic50_comp2_l3, kd=para_comp2_l3$kd)
ki_cheng_comp2_l4 <- calc_ki_app1(l0=para_comp2_l4$l0, i0=fitappic50_comp2_l4, kd=para_comp2_l4$kd)

ki_cheng_comp3_l1 <- calc_ki_app1(l0=para_comp3_l1$l0, i0=fitappic50_comp3_l1, kd=para_comp3_l1$kd)
ki_cheng_comp3_l2 <- calc_ki_app1(l0=para_comp3_l2$l0, i0=fitappic50_comp3_l2, kd=para_comp3_l2$kd)
ki_cheng_comp3_l3 <- calc_ki_app1(l0=para_comp3_l3$l0, i0=fitappic50_comp3_l3, kd=para_comp3_l3$kd)
ki_cheng_comp3_l4 <- calc_ki_app1(l0=para_comp3_l4$l0, i0=fitappic50_comp3_l4, kd=para_comp3_l4$kd)

ki_cheng_comp4_l1 <- calc_ki_app1(l0=para_comp4_l1$l0, i0=fitappic50_comp4_l1, kd=para_comp4_l1$kd)
ki_cheng_comp4_l2 <- calc_ki_app1(l0=para_comp4_l2$l0, i0=fitappic50_comp4_l2, kd=para_comp4_l2$kd)
ki_cheng_comp4_l3 <- calc_ki_app1(l0=para_comp4_l3$l0, i0=fitappic50_comp4_l3, kd=para_comp4_l3$kd)
ki_cheng_comp4_l4 <- calc_ki_app1(l0=para_comp4_l4$l0, i0=fitappic50_comp4_l4, kd=para_comp4_l4$kd)

ki_cheng_df1 <- data.frame(ki=c(ki_cheng_comp1_l1, ki_cheng_comp2_l1, ki_cheng_comp3_l1, ki_cheng_comp4_l1), 
                          type=rep(c("A", "B", "C", "D"), each=expN), 
                          l0=para_comp1_l1$l0)
ki_cheng_df2 <- data.frame(ki=c(ki_cheng_comp1_l2, ki_cheng_comp2_l2, ki_cheng_comp3_l2, ki_cheng_comp4_l2), 
                          type=rep(c("A", "B", "C", "D"), each=expN), 
                          l0=para_comp1_l2$l0)
ki_cheng_df3 <- data.frame(ki=c(ki_cheng_comp1_l3, ki_cheng_comp2_l3, ki_cheng_comp3_l3, ki_cheng_comp4_l3), 
                          type=rep(c("A", "B", "C", "D"), each=expN), 
                          l0=para_comp1_l3$l0)
ki_cheng_df4 <- data.frame(ki=c(ki_cheng_comp1_l4, ki_cheng_comp2_l4, ki_cheng_comp3_l4, ki_cheng_comp4_l4), 
                          type=rep(c("A", "B", "C", "D"), each=expN), 
                          l0=para_comp1_l4$l0)

ki_cheng_df <- rbind(ki_cheng_df1, ki_cheng_df2, ki_cheng_df3, ki_cheng_df4)
ki_cheng_df$l0 <- as.factor(ki_cheng_df$l0)
ki_cheng_violin <- ggplot(ki_cheng_df, aes(x=l0, y=ki, fill=type)) +
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count", position=position_dodge(1)) +
  scale_fill_manual(breaks=c("A", "B", "C", "D"),
                    labels=c("[P]0 = 50 nM, Kd = 50 nM", "[P]0 = 50 nM, Kd = 5 nM", "[P]0 = 5 nM, Kd = 50 nM", "[P]0 = 5 nM, Kd = 5 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[L]0 (nM)", y="Estimated app. Ki (nM)") + 
  legend_theme
  
ggsave("ki_cheng_violin.png", ki_cheng_violin, dpi=800, width=7, height=7, units='cm')
ggsave("ki_cheng_violin.svg", ki_cheng_violin, width=7, height=7, units='cm')

# ki calculated by wang's group equation
ki_wanggroup_comp1_l1 <- with(para_comp1_l1, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp1_l1, kd=kd))
ki_wanggroup_comp1_l2 <- with(para_comp1_l2, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp1_l2, kd=kd))
ki_wanggroup_comp1_l3 <- with(para_comp1_l3, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp1_l3, kd=kd))
ki_wanggroup_comp1_l4 <- with(para_comp1_l4, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp1_l4, kd=kd))

ki_wanggroup_comp2_l1 <- with(para_comp2_l1, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp2_l1, kd=kd))
ki_wanggroup_comp2_l2 <- with(para_comp2_l2, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp2_l2, kd=kd))
ki_wanggroup_comp2_l3 <- with(para_comp2_l3, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp2_l3, kd=kd))
ki_wanggroup_comp2_l4 <- with(para_comp2_l4, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp2_l4, kd=kd))

ki_wanggroup_comp3_l1 <- with(para_comp3_l1, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp3_l1, kd=kd))
ki_wanggroup_comp3_l2 <- with(para_comp3_l2, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp3_l2, kd=kd))
ki_wanggroup_comp3_l3 <- with(para_comp3_l3, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp3_l3, kd=kd))
ki_wanggroup_comp3_l4 <- with(para_comp3_l4, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp3_l4, kd=kd))

ki_wanggroup_comp4_l1 <- with(para_comp4_l1, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp4_l1, kd=kd))
ki_wanggroup_comp4_l2 <- with(para_comp4_l2, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp4_l2, kd=kd))
ki_wanggroup_comp4_l3 <- with(para_comp4_l3, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp4_l3, kd=kd))
ki_wanggroup_comp4_l4 <- with(para_comp4_l4, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp4_l4, kd=kd))

ki_wanggroup_df1 <- data.frame(ki=c(ki_wanggroup_comp1_l1, ki_wanggroup_comp2_l1, ki_wanggroup_comp3_l1, ki_wanggroup_comp4_l1), 
                          type=rep(c("A", "B", "C", "D"), each=expN), 
                          l0=para_comp1_l1$l0)
ki_wanggroup_df2 <- data.frame(ki=c(ki_wanggroup_comp1_l2, ki_wanggroup_comp2_l2, ki_wanggroup_comp3_l2, ki_wanggroup_comp4_l2), 
                          type=rep(c("A", "B", "C", "D"), each=expN), 
                          l0=para_comp1_l2$l0)
ki_wanggroup_df3 <- data.frame(ki=c(ki_wanggroup_comp1_l3, ki_wanggroup_comp2_l3, ki_wanggroup_comp3_l3, ki_wanggroup_comp4_l3), 
                          type=rep(c("A", "B", "C", "D"), each=expN), 
                          l0=para_comp1_l3$l0)
ki_wanggroup_df4 <- data.frame(ki=c(ki_wanggroup_comp1_l4, ki_wanggroup_comp2_l4, ki_wanggroup_comp3_l4, ki_wanggroup_comp4_l4), 
                          type=rep(c("A", "B", "C", "D"), each=expN), 
                          l0=para_comp1_l4$l0)

ki_wanggroup_df <- rbind(ki_wanggroup_df1, ki_wanggroup_df2, ki_wanggroup_df3, ki_wanggroup_df4)
ki_wanggroup_df$l0 <- as.factor(ki_wanggroup_df$l0)
ki_wanggroup_violin <- ggplot(ki_wanggroup_df, aes(x=l0, y=ki, fill=type)) +
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count", position=position_dodge(1)) +
  scale_fill_manual(breaks=c("A", "B", "C", "D"),
                    labels=c("[P]0 = 50 nM, Kd = 50 nM", "[P]0 = 50 nM, Kd = 5 nM", "[P]0 = 5 nM, Kd = 50 nM", "[P]0 = 5 nM, Kd = 5 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)") + 
  legend_theme

ggsave("ki_wanggroup_violin.png", ki_wanggroup_violin, dpi=800, width=7, height=7, units='cm')
ggsave("ki_wanggroup_violin.svg", ki_wanggroup_violin, width=7, height=7, units='cm')

# ki directly fitted by wang equation
ki_direct_fit_comp1_l1 <- fitki_comp1_l1["ki",]
ki_direct_fit_comp1_l2 <- fitki_comp1_l2["ki",]
ki_direct_fit_comp1_l3 <- fitki_comp1_l3["ki",]
ki_direct_fit_comp1_l4 <- fitki_comp1_l4["ki",]

ki_direct_fit_comp2_l1 <- fitki_comp2_l1["ki",]
ki_direct_fit_comp2_l2 <- fitki_comp2_l2["ki",]
ki_direct_fit_comp2_l3 <- fitki_comp2_l3["ki",]
ki_direct_fit_comp2_l4 <- fitki_comp2_l4["ki",]

ki_direct_fit_comp3_l1 <- fitki_comp3_l1["ki",]
ki_direct_fit_comp3_l2 <- fitki_comp3_l2["ki",]
ki_direct_fit_comp3_l3 <- fitki_comp3_l3["ki",]
ki_direct_fit_comp3_l4 <- fitki_comp3_l4["ki",]

ki_direct_fit_comp4_l1 <- fitki_comp4_l1["ki",]
ki_direct_fit_comp4_l2 <- fitki_comp4_l2["ki",]
ki_direct_fit_comp4_l3 <- fitki_comp4_l3["ki",]
ki_direct_fit_comp4_l4 <- fitki_comp4_l4["ki",]

ki_direct_fit_df1 <- data.frame(ki=c(ki_direct_fit_comp1_l1, ki_direct_fit_comp2_l1, ki_direct_fit_comp3_l1, ki_direct_fit_comp4_l1), 
                               type=rep(c("A", "B", "C", "D"), each=expN), 
                               l0=para_comp1_l1$l0)
ki_direct_fit_df2 <- data.frame(ki=c(ki_direct_fit_comp1_l2, ki_direct_fit_comp2_l2, ki_direct_fit_comp3_l2, ki_direct_fit_comp4_l2), 
                               type=rep(c("A", "B", "C", "D"), each=expN), 
                               l0=para_comp1_l2$l0)
ki_direct_fit_df3 <- data.frame(ki=c(ki_direct_fit_comp1_l3, ki_direct_fit_comp2_l3, ki_direct_fit_comp3_l3, ki_direct_fit_comp4_l3), 
                               type=rep(c("A", "B", "C", "D"), each=expN), 
                               l0=para_comp1_l3$l0)
ki_direct_fit_df4 <- data.frame(ki=c(ki_direct_fit_comp1_l4, ki_direct_fit_comp2_l4, ki_direct_fit_comp3_l4, ki_direct_fit_comp4_l4), 
                               type=rep(c("A", "B", "C", "D"), each=expN), 
                               l0=para_comp1_l4$l0)

ki_direct_fit_df <- rbind(ki_direct_fit_df1, ki_direct_fit_df2, ki_direct_fit_df3, ki_direct_fit_df4)
ki_direct_fit_df$l0 <- as.factor(ki_direct_fit_df$l0)
ki_direct_fit_violin <- ggplot(ki_direct_fit_df, aes(x=l0, y=ki, fill=type)) +
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count", position=position_dodge(1)) +
  scale_fill_manual(breaks=c("A", "B", "C", "D"),
                    labels=c("[P]0 = 50 nM, Kd = 50 nM", "[P]0 = 50 nM, Kd = 5 nM", "[P]0 = 5 nM, Kd = 50 nM", "[P]0 = 5 nM, Kd = 5 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)") + 
  legend_theme

ggsave("ki_direct_fit_violin.png", ki_direct_fit_violin, dpi=800, width=7, height=7, units='cm')
ggsave("ki_direct_fit_violin.svg", ki_direct_fit_violin, width=7, height=7, units='cm')

dfcomp_comp1_l1 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp1_l1[,,1]), group="A")
dfcomp_comp1_l2 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp1_l2[,,1]), group="B")
dfcomp_comp1_l3 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp1_l3[,,1]), group="C")
dfcomp_comp1_l4 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp1_l4[,,1]), group="D")

dfcomp_comp2_l1 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp2_l1[,,1]), group="A")
dfcomp_comp2_l2 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp2_l2[,,1]), group="B")
dfcomp_comp2_l3 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp2_l3[,,1]), group="C")
dfcomp_comp2_l4 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp2_l4[,,1]), group="D")

dfcomp_comp3_l1 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp3_l1[,,1]), group="A")
dfcomp_comp3_l2 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp3_l2[,,1]), group="B")
dfcomp_comp3_l3 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp3_l3[,,1]), group="C")
dfcomp_comp3_l4 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp3_l4[,,1]), group="D")

dfcomp_comp4_l1 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp4_l1[,,1]), group="A")
dfcomp_comp4_l2 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp4_l2[,,1]), group="B")
dfcomp_comp4_l3 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp4_l3[,,1]), group="C")
dfcomp_comp4_l4 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy_comp4_l4[,,1]), group="D")

dfcomp_comp1 <- rbind(dfcomp_comp1_l1, dfcomp_comp1_l2, dfcomp_comp1_l3, dfcomp_comp1_l4)
dfcomp_comp1$group <- as.factor(dfcomp_comp1$group)
dfcomp_comp1_point <- ggplot(dfcomp_comp1, aes(x=x, y=y, color=group)) + 
  geom_point(size=0.1) +
  scale_color_manual(breaks=c("D", "C", "B", "A"),
                    labels=c("[L]0 = 100 nM", "[L]0 = 75 nM", "[L]0 = 50 nM", "[L]0 = 25 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[I]0 (nM)", y="[PL]obs (nM)", title="[P]0 = 50 nM, Kd = 50 nM") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  legend_theme
ggsave("dfcomp_comp1_point.png", dfcomp_comp1_point, dpi=800, width=7, height=7, units='cm')
ggsave("dfcomp_comp1_point.svg", dfcomp_comp1_point, width=7, height=7, units='cm')

dfcomp_comp2 <- rbind(dfcomp_comp2_l1, dfcomp_comp2_l2, dfcomp_comp2_l3, dfcomp_comp2_l4)
dfcomp_comp2$group <- as.factor(dfcomp_comp2$group)
dfcomp_comp2_point <- ggplot(dfcomp_comp2, aes(x=x, y=y, color=group)) + 
  geom_point(size=0.1) +
  scale_color_manual(breaks=c("D", "C", "B", "A"),
                    labels=c("[L]0 = 100 nM", "[L]0 = 75 nM", "[L]0 = 50 nM", "[L]0 = 25 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[I]0 (nM)", y="[PL]obs (nM)", title="[P]0 = 50 nM, Kd = 5 nM") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  legend_theme
ggsave("dfcomp_comp2_point.png", dfcomp_comp2_point, dpi=800, width=7, height=7, units='cm')
ggsave("dfcomp_comp2_point.svg", dfcomp_comp2_point, width=7, height=7, units='cm')

dfcomp_comp3 <- rbind(dfcomp_comp3_l1, dfcomp_comp3_l2, dfcomp_comp3_l3, dfcomp_comp3_l4)
dfcomp_comp3$group <- as.factor(dfcomp_comp3$group)
dfcomp_comp3_point <- ggplot(dfcomp_comp3, aes(x=x, y=y, color=group)) + 
  geom_point(size=0.1) +
  scale_color_manual(breaks=c("D", "C", "B", "A"),
                    labels=c("[L]0 = 100 nM", "[L]0 = 75 nM", "[L]0 = 50 nM", "[L]0 = 25 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[I]0 (nM)", y="[PL]obs (nM)", title="[P]0 = 5 nM, Kd = 50 nM") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  legend_theme
ggsave("dfcomp_comp3_point.png", dfcomp_comp3_point, dpi=800, width=7, height=7, units='cm')
ggsave("dfcomp_comp3_point.svg", dfcomp_comp3_point, width=7, height=7, units='cm')

dfcomp_comp4 <- rbind(dfcomp_comp4_l1, dfcomp_comp4_l2, dfcomp_comp4_l3, dfcomp_comp4_l4)
dfcomp_comp4$group <- as.factor(dfcomp_comp4$group)
dfcomp_comp4_point <- ggplot(dfcomp_comp4, aes(x=x, y=y, color=group)) + 
  geom_point(size=0.1) +
  scale_color_manual(breaks=c("D", "C", "B", "A"),
                    labels=c("[L]0 = 100 nM", "[L]0 = 75 nM", "[L]0 = 50 nM", "[L]0 = 25 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[I]0 (nM)", y="[PL]obs (nM)", title="[P]0 = 5 nM, Kd = 5 nM") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  legend_theme
ggsave("dfcomp_comp4_point.png", dfcomp_comp4_point, dpi=800, width=7, height=7, units='cm')
ggsave("dfcomp_comp4_point.svg", dfcomp_comp4_point, width=7, height=7, units='cm')

# compare direct fit and wang's group approaches
l0StringList <- c("25", "50", "75", "100")
kiTypeStringList <- c("Ki (direct fit)", "Ki (Wang's group)")

violin_ki_df1 <- data.frame(
  x=rep(rep(l0StringList, each=expN), times=length(kiTypeStringList)),
  y=c(ki_direct_fit_comp1_l1, ki_direct_fit_comp1_l2, ki_direct_fit_comp1_l3, ki_direct_fit_comp1_l4,
      ki_wanggroup_comp1_l1, ki_wanggroup_comp1_l2, ki_wanggroup_comp1_l3, ki_wanggroup_comp1_l4),
  f=rep(kiTypeStringList, each=expN*length(l0StringList))
)
ki_violin1 <- ggplot(violin_ki_df1, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Ki (direct fit)", "Ki (Wang's group)"),
                    values=c("Ki (direct fit)"=c1_s1, "Ki (Wang's group)"=c2_s1)) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)", title="[P]0 = 50 nM, Kd = 50 nM") +
  scale_x_discrete(limits=l0StringList) +
  legend_theme

ggsave("ki_violin1.png", ki_violin1, dpi=800, width=7, height=7, units='cm')
ggsave("ki_violin1.svg", ki_violin1, width=7, height=7, units='cm')

violin_ki_df2 <- data.frame(
  x=rep(rep(l0StringList, each=expN), times=length(kiTypeStringList)),
  y=c(ki_direct_fit_comp2_l1, ki_direct_fit_comp2_l2, ki_direct_fit_comp2_l3, ki_direct_fit_comp2_l4,
      ki_wanggroup_comp2_l1, ki_wanggroup_comp2_l2, ki_wanggroup_comp2_l3, ki_wanggroup_comp2_l4),
  f=rep(kiTypeStringList, each=expN*length(l0StringList))
)
ki_violin2 <- ggplot(violin_ki_df2, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Ki (direct fit)", "Ki (Wang's group)"),
                    values=c("Ki (direct fit)"=c1_s1, "Ki (Wang's group)"=c2_s1)) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)", title="[P]0 = 50 nM, Kd = 5 nM") +
  scale_x_discrete(limits=l0StringList) +
  legend_theme

ggsave("ki_violin2.png", ki_violin2, dpi=800, width=7, height=7, units='cm')
ggsave("ki_violin2.svg", ki_violin2, width=7, height=7, units='cm')

violin_ki_df3 <- data.frame(
  x=rep(rep(l0StringList, each=expN), times=length(kiTypeStringList)),
  y=c(ki_direct_fit_comp3_l1, ki_direct_fit_comp3_l2, ki_direct_fit_comp3_l3, ki_direct_fit_comp3_l4,
      ki_wanggroup_comp3_l1, ki_wanggroup_comp3_l2, ki_wanggroup_comp3_l3, ki_wanggroup_comp3_l4),
  f=rep(kiTypeStringList, each=expN*length(l0StringList))
)
ki_violin3 <- ggplot(violin_ki_df3, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Ki (direct fit)", "Ki (Wang's group)"),
                    values=c("Ki (direct fit)"=c1_s1, "Ki (Wang's group)"=c2_s1)) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)", title="[P]0 = 5 nM, Kd = 50 nM") +
  scale_x_discrete(limits=l0StringList) +
  legend_theme

ggsave("ki_violin3.png", ki_violin3, dpi=800, width=7, height=7, units='cm')
ggsave("ki_violin3.svg", ki_violin3, width=7, height=7, units='cm')

violin_ki_df4 <- data.frame(
  x=rep(rep(l0StringList, each=expN), times=length(kiTypeStringList)),
  y=c(ki_direct_fit_comp4_l1, ki_direct_fit_comp4_l2, ki_direct_fit_comp4_l3, ki_direct_fit_comp4_l4,
      ki_wanggroup_comp4_l1, ki_wanggroup_comp4_l2, ki_wanggroup_comp4_l3, ki_wanggroup_comp4_l4),
  f=rep(kiTypeStringList, each=expN*length(l0StringList))
)
ki_violin4 <- ggplot(violin_ki_df4, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Ki (direct fit)", "Ki (Wang's group)"),
                    values=c("Ki (direct fit)"=c1_s1, "Ki (Wang's group)"=c2_s1)) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)", title="[P]0 = 5 nM, Kd = 5 nM") +
  scale_x_discrete(limits=l0StringList) +
  legend_theme

ggsave("ki_violin4.png", ki_violin4, dpi=800, width=7, height=7, units='cm')
ggsave("ki_violin4.svg", ki_violin4, width=7, height=7, units='cm')
setwd("..")

dir.create("Part3.1-2"); setwd("Part3.1-2")

write_table_comp_appic50 <- function(filename_base, stringList, ...) {
  csvfilename <- paste(filename_base, ".csv", sep="")
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("a (SD)", "KS test P-value of a", 
                   "b (SD)", "KS test P-value of b", "KS test P-value of log10(b)",
                   "c (SD)", "KS test P-value of c", "KS test P-value of log10(c)",
                   "d (SD)", "KS test P-value of d", "KS test P-value of log10(d)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    # a 
    a_est <- arg2[[i]]["a",]
    a_est_median <- median(a_est)
    a_est_mean <- mean(a_est)
    a_est_sd <- sd(a_est); 
    a_est_ks_pvalue <- cont_ks_test(a_est, "pnorm", a_est_mean, a_est_sd)$p.value
    
    # b 
    b_est <- abs(arg2[[i]]["b",]) # abs of b
    b_est_median <- median(b_est)
    b_est_mean <- mean(b_est)
    b_est_sd <- sd(b_est); 
    
    log10_b_est <- log10(b_est)
    log10_b_est_mean <- mean(log10_b_est)
    log10_b_est_sd <- sd(log10_b_est)
    
    b_est_ks_pvalue <- cont_ks_test(b_est, "pnorm", b_est_mean, b_est_sd)$p.value
    log10_b_est_ks_pvalue <- cont_ks_test(log10_b_est, "pnorm", log10_b_est_mean, log10_b_est_sd)$p.value
    
    # c 
    c_est <- arg2[[i]]["c",]
    c_est_median <- median(c_est)
    c_est_mean <- mean(c_est)
    c_est_sd <- sd(c_est); 
    
    log10_c_est <- log10(c_est)
    log10_c_est_mean <- mean(log10_c_est)
    log10_c_est_sd <- sd(log10_c_est)
    
    c_est_ks_pvalue <- cont_ks_test(c_est, "pnorm", c_est_mean, c_est_sd)$p.value
    log10_c_est_ks_pvalue <- cont_ks_test(log10_c_est, "pnorm", log10_c_est_mean, log10_c_est_sd)$p.value
    
    # d
    d_est <- arg2[[i]]["d",]
    d_est_median <- median(d_est)
    d_est_mean <- mean(d_est)
    d_est_sd <- sd(d_est); 
    
    log10_d_est <- log10(d_est)
    log10_d_est_mean <- mean(log10_d_est)
    log10_d_est_sd <- sd(log10_d_est)
    
    d_est_ks_pvalue <- cont_ks_test(d_est, "pnorm", d_est_mean, d_est_sd)$p.value
    log10_d_est_ks_pvalue <- cont_ks_test(log10_d_est, "pnorm", log10_d_est_mean, log10_d_est_sd)$p.value
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", a_est_median, a_est_sd),
                   sprintf("%.3f", a_est_ks_pvalue), 
                   sprintf("%.2e (%.2e)", b_est_median, b_est_sd),
                   sprintf("%.3f", b_est_ks_pvalue), sprintf("%.3f", log10_b_est_ks_pvalue),
                   sprintf("%.2e (%.2e)", c_est_median, c_est_sd),
                   sprintf("%.3f", c_est_ks_pvalue), sprintf("%.3f", log10_c_est_ks_pvalue), 
                   sprintf("%.2e (%.2e)", d_est_median, d_est_sd),
                   sprintf("%.3f", d_est_ks_pvalue), sprintf("%.3f", log10_d_est_ks_pvalue) )
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_comp_cheng <- function(filename_base, stringList, ...) {
  csvfilename <- paste(filename_base, ".csv", sep="")
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("Apparent Ki (SD)", "KS test P-value of apparent Ki", "KS test P-value of log10(apparent Ki)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    c_est <- arg2[[i]]
    c_est_median <- median(c_est)
    c_est_mean <- mean(c_est)
    c_est_sd <- sd(c_est); 
    
    log10_c_est <- log10(c_est)
    log10_c_est_mean <- mean(log10_c_est)
    log10_c_est_sd <- sd(log10_c_est)
    
    c_est_ks_pvalue <- cont_ks_test(c_est, "pnorm", c_est_mean, c_est_sd)$p.value
    log10_c_est_ks_pvalue <- cont_ks_test(log10_c_est, "pnorm", log10_c_est_mean, log10_c_est_sd)$p.value
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", c_est_median, c_est_sd),
                      sprintf("%.3f", c_est_ks_pvalue), sprintf("%.3f", log10_c_est_ks_pvalue))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_comp_ki <- function(filename_base, stringList, ...) {
  csvfilename <- paste(filename_base, ".csv", sep="")
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("Ki (SD)", "KS test P-value of Ki", "KS test P-value of log10(Ki)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    c_est <- arg2[[i]]
    c_est_median <- median(c_est)
    c_est_mean <- mean(c_est)
    c_est_sd <- sd(c_est); 
    
    log10_c_est <- log10(c_est)
    log10_c_est_mean <- mean(log10_c_est)
    log10_c_est_sd <- sd(log10_c_est)
    
    c_est_ks_pvalue <- cont_ks_test(c_est, "pnorm", c_est_mean, c_est_sd)$p.value
    log10_c_est_ks_pvalue <- cont_ks_test(log10_c_est, "pnorm", log10_c_est_mean, log10_c_est_sd)$p.value
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", c_est_median, c_est_sd),
                      sprintf("%.3f", c_est_ks_pvalue), sprintf("%.3f", log10_c_est_ks_pvalue))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}


compOutString <- c("p50kd50l25", "p50kd50l50", "p50kd50l75", "p50kd50l100", 
                   "p50kd5l25", "p50kd5l50", "p50kd5l75", "p50kd5l100", 
                   "p5kd50l25", "p5kd50l50", "p5kd50l75", "p5kd50l100", 
                   "p5kd5l25", "p5kd5l50", "p5kd5l75", "p5kd5l100")
write_table_comp_appic50(filename_base="appic50", stringList=compOutString, 
                              fitres_comp1_l1, fitres_comp1_l2, fitres_comp1_l3, fitres_comp1_l4, 
                              fitres_comp2_l1, fitres_comp2_l2, fitres_comp2_l3, fitres_comp2_l4, 
                              fitres_comp3_l1, fitres_comp3_l2, fitres_comp3_l3, fitres_comp3_l4, 
                              fitres_comp4_l1, fitres_comp4_l2, fitres_comp4_l3, fitres_comp4_l4)

write_table_comp_cheng(filename_base="cheng", stringList=compOutString, 
                       ki_cheng_comp1_l1, ki_cheng_comp1_l2, ki_cheng_comp1_l3, ki_cheng_comp1_l4, 
                       ki_cheng_comp2_l1, ki_cheng_comp2_l2, ki_cheng_comp2_l3, ki_cheng_comp2_l4, 
                       ki_cheng_comp3_l1, ki_cheng_comp3_l2, ki_cheng_comp3_l3, ki_cheng_comp3_l4, 
                       ki_cheng_comp4_l1, ki_cheng_comp4_l2, ki_cheng_comp4_l3, ki_cheng_comp4_l4)

write_table_comp_ki(filename_base="wanggroup", stringList=compOutString, 
                           ki_wanggroup_comp1_l1, ki_wanggroup_comp1_l2, ki_wanggroup_comp1_l3, ki_wanggroup_comp1_l4, 
                           ki_wanggroup_comp2_l1, ki_wanggroup_comp2_l2, ki_wanggroup_comp2_l3, ki_wanggroup_comp2_l4, 
                           ki_wanggroup_comp3_l1, ki_wanggroup_comp3_l2, ki_wanggroup_comp3_l3, ki_wanggroup_comp3_l4, 
                           ki_wanggroup_comp4_l1, ki_wanggroup_comp4_l2, ki_wanggroup_comp4_l3, ki_wanggroup_comp4_l4)

write_table_comp_ki(filename_base="directfit", stringList=compOutString, 
                           ki_direct_fit_comp1_l1, ki_direct_fit_comp1_l2, ki_direct_fit_comp1_l3, ki_direct_fit_comp1_l4, 
                           ki_direct_fit_comp2_l1, ki_direct_fit_comp2_l2, ki_direct_fit_comp2_l3, ki_direct_fit_comp2_l4, 
                           ki_direct_fit_comp3_l1, ki_direct_fit_comp3_l2, ki_direct_fit_comp3_l3, ki_direct_fit_comp3_l4, 
                           ki_direct_fit_comp4_l1, ki_direct_fit_comp4_l2, ki_direct_fit_comp4_l3, ki_direct_fit_comp4_l4)
setwd("..")


# * 3.2 Sampling Distributions of Estimated Apparent IC50 and Ki ------------
# See Part "Random Seeds and expN"

# * 3.3-SI Concentration Errors in [P]0, [L]0, and [I]0 ----------------------

# *- 3.3.1-SI Estimation of Apparent IC50, Apparent Ki, and Ki ---------------

dir.create("Part3.3.1-1-SI"); setwd("Part3.3.1-1-SI")

gen_plNoisy_array3 <- function(ti_theo_list, paraTrue, repN, expN, sd_level, tp_cv, tl_cv, ti_cv) {
  
  ti_list_len <- length(ti_theo_list)
  p0 <- paraTrue$p0; l0 <- paraTrue$l0; kd <- paraTrue$kd; ki <- paraTrue$ki
  
  # tp_cv selection
  if (tp_cv=="low"){
    tp_cv_num <- 0.025
  }
  if (tp_cv=="medium"){
    tp_cv_num <- 0.05
  }
  if (tp_cv=="high"){
    tp_cv_num <- 0.1
  }
  # tl_cv selection
  if (tl_cv=="low"){
    tl_cv_num <- 0.025
  }
  if (tl_cv=="medium"){
    tl_cv_num <- 0.05
  }
  if (tl_cv=="high"){
    tl_cv_num <- 0.1
  }
  # ti_cv selection
  if (ti_cv=="low"){
    ti_cv_num <- 0.025
  }
  if (ti_cv=="medium"){
    ti_cv_num <- 0.05
  }
  if (ti_cv=="high"){
    ti_cv_num <- 0.1
  }
  
  error_sd_low    <- p0*0.025/2
  error_sd_medium <- p0*0.05/2
  error_sd_high   <- p0*0.1/2
  
  # sd_level selection
  if (sd_level=="low"){
    error_sd <- error_sd_low
  }
  if (sd_level=="medium"){
    error_sd <- error_sd_medium
  }
  if (sd_level=="high"){
    error_sd <- error_sd_high
  }
  
  # replicate(2, sapply(c(1,2,3,4), function(x) rnorm(5, x, x*0.1)))
  ti_real_list <- as.vector(replicate(expN, sapply(ti_theo_list, function(x) rnorm(repN, x, x*ti_cv_num))))
  tp_real_list <- rnorm(repN*ti_list_len*expN, p0, p0*tp_cv_num)
  tl_real_list <- rnorm(repN*ti_list_len*expN, l0, l0*tl_cv_num)
  
  pl_real_list <- formula_competitive_pl(tp_real_list, tl_real_list, ti_real_list, kd, ki)
  pl_real_array <- array(pl_real_list, dim=c(repN, ti_list_len, expN))
  
  error_array <- array(rnorm(repN*ti_list_len*expN, sd=error_sd), dim=c(repN, ti_list_len, expN))
  pl_real_array + error_array
}

# exp. data
plNoisy2_comp1_l1 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp1_l1, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp1_l2 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp1_l2, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp1_l3 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp1_l3, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp1_l4 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp1_l4, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")

plNoisy2_comp2_l1 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp2_l1, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp2_l2 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp2_l2, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp2_l3 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp2_l3, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp2_l4 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp2_l4, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")

plNoisy2_comp3_l1 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp3_l1, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp3_l2 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp3_l2, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp3_l3 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp3_l3, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp3_l4 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp3_l4, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")

plNoisy2_comp4_l1 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp4_l1, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp4_l2 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp4_l2, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp4_l3 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp4_l3, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
plNoisy2_comp4_l4 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp4_l4, repN=repN5, expN=expN, sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")

# fitted ki
fitki2_comp1_l1 <- parApply(clust, plNoisy2_comp1_l1, 3, fitfunc_ki, paraTrue=para_comp1_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp1_l2 <- parApply(clust, plNoisy2_comp1_l2, 3, fitfunc_ki, paraTrue=para_comp1_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp1_l3 <- parApply(clust, plNoisy2_comp1_l3, 3, fitfunc_ki, paraTrue=para_comp1_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp1_l4 <- parApply(clust, plNoisy2_comp1_l4, 3, fitfunc_ki, paraTrue=para_comp1_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)

fitki2_comp2_l1 <- parApply(clust, plNoisy2_comp2_l1, 3, fitfunc_ki, paraTrue=para_comp2_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp2_l2 <- parApply(clust, plNoisy2_comp2_l2, 3, fitfunc_ki, paraTrue=para_comp2_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp2_l3 <- parApply(clust, plNoisy2_comp2_l3, 3, fitfunc_ki, paraTrue=para_comp2_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp2_l4 <- parApply(clust, plNoisy2_comp2_l4, 3, fitfunc_ki, paraTrue=para_comp2_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)

fitki2_comp3_l1 <- parApply(clust, plNoisy2_comp3_l1, 3, fitfunc_ki, paraTrue=para_comp3_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp3_l2 <- parApply(clust, plNoisy2_comp3_l2, 3, fitfunc_ki, paraTrue=para_comp3_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp3_l3 <- parApply(clust, plNoisy2_comp3_l3, 3, fitfunc_ki, paraTrue=para_comp3_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp3_l4 <- parApply(clust, plNoisy2_comp3_l4, 3, fitfunc_ki, paraTrue=para_comp3_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)

fitki2_comp4_l1 <- parApply(clust, plNoisy2_comp4_l1, 3, fitfunc_ki, paraTrue=para_comp4_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp4_l2 <- parApply(clust, plNoisy2_comp4_l2, 3, fitfunc_ki, paraTrue=para_comp4_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp4_l3 <- parApply(clust, plNoisy2_comp4_l3, 3, fitfunc_ki, paraTrue=para_comp4_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki2_comp4_l4 <- parApply(clust, plNoisy2_comp4_l4, 3, fitfunc_ki, paraTrue=para_comp4_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)

# fitted results
fitres2_comp1_l1 <- parApply(clust, plNoisy2_comp1_l1, 3, fitfunc_app_ic50, paraTrue=para_comp1_l1, repN=repN5, i_conc=ti_num_large)
fitres2_comp1_l2 <- parApply(clust, plNoisy2_comp1_l2, 3, fitfunc_app_ic50, paraTrue=para_comp1_l2, repN=repN5, i_conc=ti_num_large)
fitres2_comp1_l3 <- parApply(clust, plNoisy2_comp1_l3, 3, fitfunc_app_ic50, paraTrue=para_comp1_l3, repN=repN5, i_conc=ti_num_large)
fitres2_comp1_l4 <- parApply(clust, plNoisy2_comp1_l4, 3, fitfunc_app_ic50, paraTrue=para_comp1_l4, repN=repN5, i_conc=ti_num_large)

fitres2_comp2_l1 <- parApply(clust, plNoisy2_comp2_l1, 3, fitfunc_app_ic50, paraTrue=para_comp2_l1, repN=repN5, i_conc=ti_num_large)
fitres2_comp2_l2 <- parApply(clust, plNoisy2_comp2_l2, 3, fitfunc_app_ic50, paraTrue=para_comp2_l2, repN=repN5, i_conc=ti_num_large)
fitres2_comp2_l3 <- parApply(clust, plNoisy2_comp2_l3, 3, fitfunc_app_ic50, paraTrue=para_comp2_l3, repN=repN5, i_conc=ti_num_large)
fitres2_comp2_l4 <- parApply(clust, plNoisy2_comp2_l4, 3, fitfunc_app_ic50, paraTrue=para_comp2_l4, repN=repN5, i_conc=ti_num_large)

fitres2_comp3_l1 <- parApply(clust, plNoisy2_comp3_l1, 3, fitfunc_app_ic50, paraTrue=para_comp3_l1, repN=repN5, i_conc=ti_num_large)
fitres2_comp3_l2 <- parApply(clust, plNoisy2_comp3_l2, 3, fitfunc_app_ic50, paraTrue=para_comp3_l2, repN=repN5, i_conc=ti_num_large)
fitres2_comp3_l3 <- parApply(clust, plNoisy2_comp3_l3, 3, fitfunc_app_ic50, paraTrue=para_comp3_l3, repN=repN5, i_conc=ti_num_large)
fitres2_comp3_l4 <- parApply(clust, plNoisy2_comp3_l4, 3, fitfunc_app_ic50, paraTrue=para_comp3_l4, repN=repN5, i_conc=ti_num_large)

fitres2_comp4_l1 <- parApply(clust, plNoisy2_comp4_l1, 3, fitfunc_app_ic50, paraTrue=para_comp4_l1, repN=repN5, i_conc=ti_num_large)
fitres2_comp4_l2 <- parApply(clust, plNoisy2_comp4_l2, 3, fitfunc_app_ic50, paraTrue=para_comp4_l2, repN=repN5, i_conc=ti_num_large)
fitres2_comp4_l3 <- parApply(clust, plNoisy2_comp4_l3, 3, fitfunc_app_ic50, paraTrue=para_comp4_l3, repN=repN5, i_conc=ti_num_large)
fitres2_comp4_l4 <- parApply(clust, plNoisy2_comp4_l4, 3, fitfunc_app_ic50, paraTrue=para_comp4_l4, repN=repN5, i_conc=ti_num_large)

# fitted apparent IC50
fitappic502_comp1_l1 <- fitres2_comp1_l1["c",]
fitappic502_comp1_l2 <- fitres2_comp1_l2["c",]
fitappic502_comp1_l3 <- fitres2_comp1_l3["c",]
fitappic502_comp1_l4 <- fitres2_comp1_l4["c",]

fitappic502_comp2_l1 <- fitres2_comp2_l1["c",]
fitappic502_comp2_l2 <- fitres2_comp2_l2["c",]
fitappic502_comp2_l3 <- fitres2_comp2_l3["c",]
fitappic502_comp2_l4 <- fitres2_comp2_l4["c",]

fitappic502_comp3_l1 <- fitres2_comp3_l1["c",]
fitappic502_comp3_l2 <- fitres2_comp3_l2["c",]
fitappic502_comp3_l3 <- fitres2_comp3_l3["c",]
fitappic502_comp3_l4 <- fitres2_comp3_l4["c",]

fitappic502_comp4_l1 <- fitres2_comp4_l1["c",]
fitappic502_comp4_l2 <- fitres2_comp4_l2["c",]
fitappic502_comp4_l3 <- fitres2_comp4_l3["c",]
fitappic502_comp4_l4 <- fitres2_comp4_l4["c",]

# ratio of fitted to theoretical apparent IC50
diffappic502_comp1_l1 <- fitappic502_comp1_l1 / theoappic50_comp1_l1
diffappic502_comp1_l2 <- fitappic502_comp1_l2 / theoappic50_comp1_l2
diffappic502_comp1_l3 <- fitappic502_comp1_l3 / theoappic50_comp1_l3
diffappic502_comp1_l4 <- fitappic502_comp1_l4 / theoappic50_comp1_l4

diffappic502_comp2_l1 <- fitappic502_comp2_l1 / theoappic50_comp2_l1
diffappic502_comp2_l2 <- fitappic502_comp2_l2 / theoappic50_comp2_l2
diffappic502_comp2_l3 <- fitappic502_comp2_l3 / theoappic50_comp2_l3
diffappic502_comp2_l4 <- fitappic502_comp2_l4 / theoappic50_comp2_l4

diffappic502_comp3_l1 <- fitappic502_comp3_l1 / theoappic50_comp3_l1
diffappic502_comp3_l2 <- fitappic502_comp3_l2 / theoappic50_comp3_l2
diffappic502_comp3_l3 <- fitappic502_comp3_l3 / theoappic50_comp3_l3
diffappic502_comp3_l4 <- fitappic502_comp3_l4 / theoappic50_comp3_l4

diffappic502_comp4_l1 <- fitappic502_comp4_l1 / theoappic50_comp4_l1
diffappic502_comp4_l2 <- fitappic502_comp4_l2 / theoappic50_comp4_l2
diffappic502_comp4_l3 <- fitappic502_comp4_l3 / theoappic50_comp4_l3
diffappic502_comp4_l4 <- fitappic502_comp4_l4 / theoappic50_comp4_l4


# ref: https://forum.posit.co/t/how-properly-use-scale-fill-manual-and-its-arguments-labels-breaks-limits-colours-being-messed-up/177618/2
diffappic50_df1 <- data.frame(ki=c(diffappic502_comp1_l1, diffappic502_comp2_l1, diffappic502_comp3_l1, diffappic502_comp4_l1), 
                              type=rep(c("A", "B", "C", "D"), each=expN), 
                              l0=para_comp1_l1$l0)
diffappic50_df2 <- data.frame(ki=c(diffappic502_comp1_l2, diffappic502_comp2_l2, diffappic502_comp3_l2, diffappic502_comp4_l2), 
                              type=rep(c("A", "B", "C", "D"), each=expN), 
                              l0=para_comp1_l2$l0)
diffappic50_df3 <- data.frame(ki=c(diffappic502_comp1_l3, diffappic502_comp2_l3, diffappic502_comp3_l3, diffappic502_comp4_l3), 
                              type=rep(c("A", "B", "C", "D"), each=expN), 
                              l0=para_comp1_l3$l0)
diffappic50_df4 <- data.frame(ki=c(diffappic502_comp1_l4, diffappic502_comp2_l4, diffappic502_comp3_l4, diffappic502_comp4_l4), 
                              type=rep(c("A", "B", "C", "D"), each=expN), 
                              l0=para_comp1_l4$l0)

diffappic50_df <- rbind(diffappic50_df1, diffappic50_df2, diffappic50_df3, diffappic50_df4)
diffappic50_df$l0 <- as.factor(diffappic50_df$l0)
diffappic50_violin <- ggplot(diffappic50_df, aes(x=l0, y=ki, fill=type)) +
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count", position=position_dodge(1)) +
  scale_fill_manual(breaks=c("A", "B", "C", "D"),
                    labels=c("[P]0 = 50 nM, Kd = 50 nM", "[P]0 = 50 nM, Kd = 5 nM", "[P]0 = 5 nM, Kd = 50 nM", "[P]0 = 5 nM, Kd = 5 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], 
                             "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[L]0 (nM)", y="Estimated app. IC50 / Theoretical app. IC50") + 
  legend_theme

ggsave("diffappic50_violin.png", diffappic50_violin, dpi=800, width=7, height=7, units='cm')
ggsave("diffappic50_violin.svg", diffappic50_violin, width=7, height=7, units='cm')

# ki calculated by cheng-prusoff equation
ki_cheng2_comp1_l1 <- calc_ki_app1(l0=para_comp1_l1$l0, i0=fitappic502_comp1_l1, kd=para_comp1_l1$kd)
ki_cheng2_comp1_l2 <- calc_ki_app1(l0=para_comp1_l2$l0, i0=fitappic502_comp1_l2, kd=para_comp1_l2$kd)
ki_cheng2_comp1_l3 <- calc_ki_app1(l0=para_comp1_l3$l0, i0=fitappic502_comp1_l3, kd=para_comp1_l3$kd)
ki_cheng2_comp1_l4 <- calc_ki_app1(l0=para_comp1_l4$l0, i0=fitappic502_comp1_l4, kd=para_comp1_l4$kd)

ki_cheng2_comp2_l1 <- calc_ki_app1(l0=para_comp2_l1$l0, i0=fitappic502_comp2_l1, kd=para_comp2_l1$kd)
ki_cheng2_comp2_l2 <- calc_ki_app1(l0=para_comp2_l2$l0, i0=fitappic502_comp2_l2, kd=para_comp2_l2$kd)
ki_cheng2_comp2_l3 <- calc_ki_app1(l0=para_comp2_l3$l0, i0=fitappic502_comp2_l3, kd=para_comp2_l3$kd)
ki_cheng2_comp2_l4 <- calc_ki_app1(l0=para_comp2_l4$l0, i0=fitappic502_comp2_l4, kd=para_comp2_l4$kd)

ki_cheng2_comp3_l1 <- calc_ki_app1(l0=para_comp3_l1$l0, i0=fitappic502_comp3_l1, kd=para_comp3_l1$kd)
ki_cheng2_comp3_l2 <- calc_ki_app1(l0=para_comp3_l2$l0, i0=fitappic502_comp3_l2, kd=para_comp3_l2$kd)
ki_cheng2_comp3_l3 <- calc_ki_app1(l0=para_comp3_l3$l0, i0=fitappic502_comp3_l3, kd=para_comp3_l3$kd)
ki_cheng2_comp3_l4 <- calc_ki_app1(l0=para_comp3_l4$l0, i0=fitappic502_comp3_l4, kd=para_comp3_l4$kd)

ki_cheng2_comp4_l1 <- calc_ki_app1(l0=para_comp4_l1$l0, i0=fitappic502_comp4_l1, kd=para_comp4_l1$kd)
ki_cheng2_comp4_l2 <- calc_ki_app1(l0=para_comp4_l2$l0, i0=fitappic502_comp4_l2, kd=para_comp4_l2$kd)
ki_cheng2_comp4_l3 <- calc_ki_app1(l0=para_comp4_l3$l0, i0=fitappic502_comp4_l3, kd=para_comp4_l3$kd)
ki_cheng2_comp4_l4 <- calc_ki_app1(l0=para_comp4_l4$l0, i0=fitappic502_comp4_l4, kd=para_comp4_l4$kd)

ki_cheng_df1 <- data.frame(ki=c(ki_cheng2_comp1_l1, ki_cheng2_comp2_l1, ki_cheng2_comp3_l1, ki_cheng2_comp4_l1), 
                           type=rep(c("A", "B", "C", "D"), each=expN), 
                           l0=para_comp1_l1$l0)
ki_cheng_df2 <- data.frame(ki=c(ki_cheng2_comp1_l2, ki_cheng2_comp2_l2, ki_cheng2_comp3_l2, ki_cheng2_comp4_l2), 
                           type=rep(c("A", "B", "C", "D"), each=expN), 
                           l0=para_comp1_l2$l0)
ki_cheng_df3 <- data.frame(ki=c(ki_cheng2_comp1_l3, ki_cheng2_comp2_l3, ki_cheng2_comp3_l3, ki_cheng2_comp4_l3), 
                           type=rep(c("A", "B", "C", "D"), each=expN), 
                           l0=para_comp1_l3$l0)
ki_cheng_df4 <- data.frame(ki=c(ki_cheng2_comp1_l4, ki_cheng2_comp2_l4, ki_cheng2_comp3_l4, ki_cheng2_comp4_l4), 
                           type=rep(c("A", "B", "C", "D"), each=expN), 
                           l0=para_comp1_l4$l0)

ki_cheng_df <- rbind(ki_cheng_df1, ki_cheng_df2, ki_cheng_df3, ki_cheng_df4)
ki_cheng_df$l0 <- as.factor(ki_cheng_df$l0)
ki_cheng_violin <- ggplot(ki_cheng_df, aes(x=l0, y=ki, fill=type)) +
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count", position=position_dodge(1)) +
  scale_fill_manual(breaks=c("A", "B", "C", "D"),
                    labels=c("[P]0 = 50 nM, Kd = 50 nM", "[P]0 = 50 nM, Kd = 5 nM", "[P]0 = 5 nM, Kd = 50 nM", "[P]0 = 5 nM, Kd = 5 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[L]0 (nM)", y="Estimated app. Ki (nM)") + 
  legend_theme

ggsave("ki_cheng_violin.png", ki_cheng_violin, dpi=800, width=7, height=7, units='cm')
ggsave("ki_cheng_violin.svg", ki_cheng_violin, width=7, height=7, units='cm')

# ki calculated by wang's group equation
ki_wanggroup2_comp1_l1 <- with(para_comp1_l1, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp1_l1, kd=kd))
ki_wanggroup2_comp1_l2 <- with(para_comp1_l2, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp1_l2, kd=kd))
ki_wanggroup2_comp1_l3 <- with(para_comp1_l3, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp1_l3, kd=kd))
ki_wanggroup2_comp1_l4 <- with(para_comp1_l4, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp1_l4, kd=kd))

ki_wanggroup2_comp2_l1 <- with(para_comp2_l1, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp2_l1, kd=kd))
ki_wanggroup2_comp2_l2 <- with(para_comp2_l2, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp2_l2, kd=kd))
ki_wanggroup2_comp2_l3 <- with(para_comp2_l3, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp2_l3, kd=kd))
ki_wanggroup2_comp2_l4 <- with(para_comp2_l4, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp2_l4, kd=kd))

ki_wanggroup2_comp3_l1 <- with(para_comp3_l1, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp3_l1, kd=kd))
ki_wanggroup2_comp3_l2 <- with(para_comp3_l2, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp3_l2, kd=kd))
ki_wanggroup2_comp3_l3 <- with(para_comp3_l3, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp3_l3, kd=kd))
ki_wanggroup2_comp3_l4 <- with(para_comp3_l4, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp3_l4, kd=kd))

ki_wanggroup2_comp4_l1 <- with(para_comp4_l1, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp4_l1, kd=kd))
ki_wanggroup2_comp4_l2 <- with(para_comp4_l2, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp4_l2, kd=kd))
ki_wanggroup2_comp4_l3 <- with(para_comp4_l3, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp4_l3, kd=kd))
ki_wanggroup2_comp4_l4 <- with(para_comp4_l4, calc_ki(p0=p0, l0=l0, i0=fitappic502_comp4_l4, kd=kd))

ki_wanggroup_df1 <- data.frame(ki=c(ki_wanggroup2_comp1_l1, ki_wanggroup2_comp2_l1, ki_wanggroup2_comp3_l1, ki_wanggroup2_comp4_l1), 
                               type=rep(c("A", "B", "C", "D"), each=expN), 
                               l0=para_comp1_l1$l0)
ki_wanggroup_df2 <- data.frame(ki=c(ki_wanggroup2_comp1_l2, ki_wanggroup2_comp2_l2, ki_wanggroup2_comp3_l2, ki_wanggroup2_comp4_l2), 
                               type=rep(c("A", "B", "C", "D"), each=expN), 
                               l0=para_comp1_l2$l0)
ki_wanggroup_df3 <- data.frame(ki=c(ki_wanggroup2_comp1_l3, ki_wanggroup2_comp2_l3, ki_wanggroup2_comp3_l3, ki_wanggroup2_comp4_l3), 
                               type=rep(c("A", "B", "C", "D"), each=expN), 
                               l0=para_comp1_l3$l0)
ki_wanggroup_df4 <- data.frame(ki=c(ki_wanggroup2_comp1_l4, ki_wanggroup2_comp2_l4, ki_wanggroup2_comp3_l4, ki_wanggroup2_comp4_l4), 
                               type=rep(c("A", "B", "C", "D"), each=expN), 
                               l0=para_comp1_l4$l0)

ki_wanggroup_df <- rbind(ki_wanggroup_df1, ki_wanggroup_df2, ki_wanggroup_df3, ki_wanggroup_df4)
ki_wanggroup_df$l0 <- as.factor(ki_wanggroup_df$l0)
ki_wanggroup_violin <- ggplot(ki_wanggroup_df, aes(x=l0, y=ki, fill=type)) +
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count", position=position_dodge(1)) +
  scale_fill_manual(breaks=c("A", "B", "C", "D"),
                    labels=c("[P]0 = 50 nM, Kd = 50 nM", "[P]0 = 50 nM, Kd = 5 nM", "[P]0 = 5 nM, Kd = 50 nM", "[P]0 = 5 nM, Kd = 5 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)") + 
  legend_theme

ggsave("ki_wanggroup_violin.png", ki_wanggroup_violin, dpi=800, width=7, height=7, units='cm')
ggsave("ki_wanggroup_violin.svg", ki_wanggroup_violin, width=7, height=7, units='cm')

# ki directly fitted by wang equation
ki_direct_fit2_comp1_l1 <- fitki2_comp1_l1["ki",]
ki_direct_fit2_comp1_l2 <- fitki2_comp1_l2["ki",]
ki_direct_fit2_comp1_l3 <- fitki2_comp1_l3["ki",]
ki_direct_fit2_comp1_l4 <- fitki2_comp1_l4["ki",]

ki_direct_fit2_comp2_l1 <- fitki2_comp2_l1["ki",]
ki_direct_fit2_comp2_l2 <- fitki2_comp2_l2["ki",]
ki_direct_fit2_comp2_l3 <- fitki2_comp2_l3["ki",]
ki_direct_fit2_comp2_l4 <- fitki2_comp2_l4["ki",]

ki_direct_fit2_comp3_l1 <- fitki2_comp3_l1["ki",]
ki_direct_fit2_comp3_l2 <- fitki2_comp3_l2["ki",]
ki_direct_fit2_comp3_l3 <- fitki2_comp3_l3["ki",]
ki_direct_fit2_comp3_l4 <- fitki2_comp3_l4["ki",]

ki_direct_fit2_comp4_l1 <- fitki2_comp4_l1["ki",]
ki_direct_fit2_comp4_l2 <- fitki2_comp4_l2["ki",]
ki_direct_fit2_comp4_l3 <- fitki2_comp4_l3["ki",]
ki_direct_fit2_comp4_l4 <- fitki2_comp4_l4["ki",]

ki_direct_fit_df1 <- data.frame(ki=c(ki_direct_fit2_comp1_l1, ki_direct_fit2_comp2_l1, ki_direct_fit2_comp3_l1, ki_direct_fit2_comp4_l1), 
                                type=rep(c("A", "B", "C", "D"), each=expN), 
                                l0=para_comp1_l1$l0)
ki_direct_fit_df2 <- data.frame(ki=c(ki_direct_fit2_comp1_l2, ki_direct_fit2_comp2_l2, ki_direct_fit2_comp3_l2, ki_direct_fit2_comp4_l2), 
                                type=rep(c("A", "B", "C", "D"), each=expN), 
                                l0=para_comp1_l2$l0)
ki_direct_fit_df3 <- data.frame(ki=c(ki_direct_fit2_comp1_l3, ki_direct_fit2_comp2_l3, ki_direct_fit2_comp3_l3, ki_direct_fit2_comp4_l3), 
                                type=rep(c("A", "B", "C", "D"), each=expN), 
                                l0=para_comp1_l3$l0)
ki_direct_fit_df4 <- data.frame(ki=c(ki_direct_fit2_comp1_l4, ki_direct_fit2_comp2_l4, ki_direct_fit2_comp3_l4, ki_direct_fit2_comp4_l4), 
                                type=rep(c("A", "B", "C", "D"), each=expN), 
                                l0=para_comp1_l4$l0)

ki_direct_fit_df <- rbind(ki_direct_fit_df1, ki_direct_fit_df2, ki_direct_fit_df3, ki_direct_fit_df4)
ki_direct_fit_df$l0 <- as.factor(ki_direct_fit_df$l0)
ki_direct_fit_violin <- ggplot(ki_direct_fit_df, aes(x=l0, y=ki, fill=type)) +
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count", position=position_dodge(1)) +
  scale_fill_manual(breaks=c("A", "B", "C", "D"),
                    labels=c("[P]0 = 50 nM, Kd = 50 nM", "[P]0 = 50 nM, Kd = 5 nM", "[P]0 = 5 nM, Kd = 50 nM", "[P]0 = 5 nM, Kd = 5 nM"),
                    values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)") + 
  legend_theme

ggsave("ki_direct_fit_violin.png", ki_direct_fit_violin, dpi=800, width=7, height=7, units='cm')
ggsave("ki_direct_fit_violin.svg", ki_direct_fit_violin, width=7, height=7, units='cm')

dfcomp_comp1_l1 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp1_l1[,,1]), group="A")
dfcomp_comp1_l2 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp1_l2[,,1]), group="B")
dfcomp_comp1_l3 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp1_l3[,,1]), group="C")
dfcomp_comp1_l4 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp1_l4[,,1]), group="D")

dfcomp_comp2_l1 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp2_l1[,,1]), group="A")
dfcomp_comp2_l2 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp2_l2[,,1]), group="B")
dfcomp_comp2_l3 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp2_l3[,,1]), group="C")
dfcomp_comp2_l4 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp2_l4[,,1]), group="D")

dfcomp_comp3_l1 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp3_l1[,,1]), group="A")
dfcomp_comp3_l2 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp3_l2[,,1]), group="B")
dfcomp_comp3_l3 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp3_l3[,,1]), group="C")
dfcomp_comp3_l4 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp3_l4[,,1]), group="D")

dfcomp_comp4_l1 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp4_l1[,,1]), group="A")
dfcomp_comp4_l2 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp4_l2[,,1]), group="B")
dfcomp_comp4_l3 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp4_l3[,,1]), group="C")
dfcomp_comp4_l4 <- data.frame(x=rep(ti_num_large, each=repN5), y=as.vector(plNoisy2_comp4_l4[,,1]), group="D")

dfcomp_comp1 <- rbind(dfcomp_comp1_l1, dfcomp_comp1_l2, dfcomp_comp1_l3, dfcomp_comp1_l4)
dfcomp_comp1$group <- as.factor(dfcomp_comp1$group)
dfcomp_comp1_point <- ggplot(dfcomp_comp1, aes(x=x, y=y, color=group)) + 
  geom_point(size=0.1) +
  scale_color_manual(breaks=c("D", "C", "B", "A"),
                     labels=c("[L]0 = 100 nM", "[L]0 = 75 nM", "[L]0 = 50 nM", "[L]0 = 25 nM"),
                     values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[I]0 (nM)", y="[PL]obs (nM)", title="[P]0 = 50 nM, Kd = 50 nM") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  legend_theme
ggsave("dfcomp_comp1_point.png", dfcomp_comp1_point, dpi=800, width=7, height=7, units='cm')
ggsave("dfcomp_comp1_point.svg", dfcomp_comp1_point, width=7, height=7, units='cm')

dfcomp_comp2 <- rbind(dfcomp_comp2_l1, dfcomp_comp2_l2, dfcomp_comp2_l3, dfcomp_comp2_l4)
dfcomp_comp2$group <- as.factor(dfcomp_comp2$group)
dfcomp_comp2_point <- ggplot(dfcomp_comp2, aes(x=x, y=y, color=group)) + 
  geom_point(size=0.1) +
  scale_color_manual(breaks=c("D", "C", "B", "A"),
                     labels=c("[L]0 = 100 nM", "[L]0 = 75 nM", "[L]0 = 50 nM", "[L]0 = 25 nM"),
                     values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[I]0 (nM)", y="[PL]obs (nM)", title="[P]0 = 50 nM, Kd = 5 nM") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  legend_theme
ggsave("dfcomp_comp2_point.png", dfcomp_comp2_point, dpi=800, width=7, height=7, units='cm')
ggsave("dfcomp_comp2_point.svg", dfcomp_comp2_point, width=7, height=7, units='cm')

dfcomp_comp3 <- rbind(dfcomp_comp3_l1, dfcomp_comp3_l2, dfcomp_comp3_l3, dfcomp_comp3_l4)
dfcomp_comp3$group <- as.factor(dfcomp_comp3$group)
dfcomp_comp3_point <- ggplot(dfcomp_comp3, aes(x=x, y=y, color=group)) + 
  geom_point(size=0.1) +
  scale_color_manual(breaks=c("D", "C", "B", "A"),
                     labels=c("[L]0 = 100 nM", "[L]0 = 75 nM", "[L]0 = 50 nM", "[L]0 = 25 nM"),
                     values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[I]0 (nM)", y="[PL]obs (nM)", title="[P]0 = 5 nM, Kd = 50 nM") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  legend_theme
ggsave("dfcomp_comp3_point.png", dfcomp_comp3_point, dpi=800, width=7, height=7, units='cm')
ggsave("dfcomp_comp3_point.svg", dfcomp_comp3_point, width=7, height=7, units='cm')

dfcomp_comp4 <- rbind(dfcomp_comp4_l1, dfcomp_comp4_l2, dfcomp_comp4_l3, dfcomp_comp4_l4)
dfcomp_comp4$group <- as.factor(dfcomp_comp4$group)
dfcomp_comp4_point <- ggplot(dfcomp_comp4, aes(x=x, y=y, color=group)) + 
  geom_point(size=0.1) +
  scale_color_manual(breaks=c("D", "C", "B", "A"),
                     labels=c("[L]0 = 100 nM", "[L]0 = 75 nM", "[L]0 = 50 nM", "[L]0 = 25 nM"),
                     values=c("A"=plt_color[1], "B"=plt_color[2], "C"=plt_color[3], "D"=plt_color[4])) +
  labs(x="[I]0 (nM)", y="[PL]obs (nM)", title="[P]0 = 5 nM, Kd = 5 nM") + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  legend_theme
ggsave("dfcomp_comp4_point.png", dfcomp_comp4_point, dpi=800, width=7, height=7, units='cm')
ggsave("dfcomp_comp4_point.svg", dfcomp_comp4_point, width=7, height=7, units='cm')

# compare direct fit and wang's group approaches
l0StringList <- c("25", "50", "75", "100")
kiTypeStringList <- c("Ki (direct fit)", "Ki (Wang's group)")

violin_ki_df1 <- data.frame(
  x=rep(rep(l0StringList, each=expN), times=length(kiTypeStringList)),
  y=c(ki_direct_fit2_comp1_l1, ki_direct_fit2_comp1_l2, ki_direct_fit2_comp1_l3, ki_direct_fit2_comp1_l4,
      ki_wanggroup2_comp1_l1, ki_wanggroup2_comp1_l2, ki_wanggroup2_comp1_l3, ki_wanggroup2_comp1_l4),
  f=rep(kiTypeStringList, each=expN*length(l0StringList))
)
ki_violin1 <- ggplot(violin_ki_df1, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Ki (direct fit)", "Ki (Wang's group)"),
                    values=c("Ki (direct fit)"=c1_s1, "Ki (Wang's group)"=c2_s1)) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)", title="[P]0 = 50 nM, Kd = 50 nM") +
  scale_x_discrete(limits=l0StringList) +
  legend_theme

ggsave("ki_violin1.png", ki_violin1, dpi=800, width=7, height=7, units='cm')
ggsave("ki_violin1.svg", ki_violin1, width=7, height=7, units='cm')

violin_ki_df2 <- data.frame(
  x=rep(rep(l0StringList, each=expN), times=length(kiTypeStringList)),
  y=c(ki_direct_fit2_comp2_l1, ki_direct_fit2_comp2_l2, ki_direct_fit2_comp2_l3, ki_direct_fit2_comp2_l4,
      ki_wanggroup2_comp2_l1, ki_wanggroup2_comp2_l2, ki_wanggroup2_comp2_l3, ki_wanggroup2_comp2_l4),
  f=rep(kiTypeStringList, each=expN*length(l0StringList))
)
ki_violin2 <- ggplot(violin_ki_df2, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Ki (direct fit)", "Ki (Wang's group)"),
                    values=c("Ki (direct fit)"=c1_s1, "Ki (Wang's group)"=c2_s1)) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)", title="[P]0 = 50 nM, Kd = 5 nM") +
  scale_x_discrete(limits=l0StringList) +
  legend_theme

ggsave("ki_violin2.png", ki_violin2, dpi=800, width=7, height=7, units='cm')
ggsave("ki_violin2.svg", ki_violin2, width=7, height=7, units='cm')

violin_ki_df3 <- data.frame(
  x=rep(rep(l0StringList, each=expN), times=length(kiTypeStringList)),
  y=c(ki_direct_fit2_comp3_l1, ki_direct_fit2_comp3_l2, ki_direct_fit2_comp3_l3, ki_direct_fit2_comp3_l4,
      ki_wanggroup2_comp3_l1, ki_wanggroup2_comp3_l2, ki_wanggroup2_comp3_l3, ki_wanggroup2_comp3_l4),
  f=rep(kiTypeStringList, each=expN*length(l0StringList))
)
ki_violin3 <- ggplot(violin_ki_df3, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Ki (direct fit)", "Ki (Wang's group)"),
                    values=c("Ki (direct fit)"=c1_s1, "Ki (Wang's group)"=c2_s1)) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)", title="[P]0 = 5 nM, Kd = 50 nM") +
  scale_x_discrete(limits=l0StringList) +
  legend_theme

ggsave("ki_violin3.png", ki_violin3, dpi=800, width=7, height=7, units='cm')
ggsave("ki_violin3.svg", ki_violin3, width=7, height=7, units='cm')

violin_ki_df4 <- data.frame(
  x=rep(rep(l0StringList, each=expN), times=length(kiTypeStringList)),
  y=c(ki_direct_fit2_comp4_l1, ki_direct_fit2_comp4_l2, ki_direct_fit2_comp4_l3, ki_direct_fit2_comp4_l4,
      ki_wanggroup2_comp4_l1, ki_wanggroup2_comp4_l2, ki_wanggroup2_comp4_l3, ki_wanggroup2_comp4_l4),
  f=rep(kiTypeStringList, each=expN*length(l0StringList))
)
ki_violin4 <- ggplot(violin_ki_df4, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, draw_quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Ki (direct fit)", "Ki (Wang's group)"),
                    values=c("Ki (direct fit)"=c1_s1, "Ki (Wang's group)"=c2_s1)) +
  labs(x="[L]0 (nM)", y="Estimated Ki (nM)", title="[P]0 = 5 nM, Kd = 5 nM") +
  scale_x_discrete(limits=l0StringList) +
  legend_theme

ggsave("ki_violin4.png", ki_violin4, dpi=800, width=7, height=7, units='cm')
ggsave("ki_violin4.svg", ki_violin4, width=7, height=7, units='cm')
setwd("..")

dir.create("Part3.3.1-2-SI"); setwd("Part3.3.1-2-SI")

compOutString <- c("p50kd50l25", "p50kd50l50", "p50kd50l75", "p50kd50l100", 
                   "p50kd5l25", "p50kd5l50", "p50kd5l75", "p50kd5l100", 
                   "p5kd50l25", "p5kd50l50", "p5kd50l75", "p5kd50l100", 
                   "p5kd5l25", "p5kd5l50", "p5kd5l75", "p5kd5l100")
write_table_comp_appic50(filename_base="appic50", stringList=compOutString, 
                              fitres2_comp1_l1, fitres2_comp1_l2, fitres2_comp1_l3, fitres2_comp1_l4, 
                              fitres2_comp2_l1, fitres2_comp2_l2, fitres2_comp2_l3, fitres2_comp2_l4, 
                              fitres2_comp3_l1, fitres2_comp3_l2, fitres2_comp3_l3, fitres2_comp3_l4, 
                              fitres2_comp4_l1, fitres2_comp4_l2, fitres2_comp4_l3, fitres2_comp4_l4)

write_table_comp_cheng(filename_base="cheng", stringList=compOutString, 
                       ki_cheng2_comp1_l1, ki_cheng2_comp1_l2, ki_cheng2_comp1_l3, ki_cheng2_comp1_l4, 
                       ki_cheng2_comp2_l1, ki_cheng2_comp2_l2, ki_cheng2_comp2_l3, ki_cheng2_comp2_l4, 
                       ki_cheng2_comp3_l1, ki_cheng2_comp3_l2, ki_cheng2_comp3_l3, ki_cheng2_comp3_l4, 
                       ki_cheng2_comp4_l1, ki_cheng2_comp4_l2, ki_cheng2_comp4_l3, ki_cheng2_comp4_l4)

write_table_comp_ki(filename_base="wanggroup", stringList=compOutString, 
                           ki_wanggroup2_comp1_l1, ki_wanggroup2_comp1_l2, ki_wanggroup2_comp1_l3, ki_wanggroup2_comp1_l4, 
                           ki_wanggroup2_comp2_l1, ki_wanggroup2_comp2_l2, ki_wanggroup2_comp2_l3, ki_wanggroup2_comp2_l4, 
                           ki_wanggroup2_comp3_l1, ki_wanggroup2_comp3_l2, ki_wanggroup2_comp3_l3, ki_wanggroup2_comp3_l4, 
                           ki_wanggroup2_comp4_l1, ki_wanggroup2_comp4_l2, ki_wanggroup2_comp4_l3, ki_wanggroup2_comp4_l4)

write_table_comp_ki(filename_base="directfit", stringList=compOutString, 
                           ki_direct_fit2_comp1_l1, ki_direct_fit2_comp1_l2, ki_direct_fit2_comp1_l3, ki_direct_fit2_comp1_l4, 
                           ki_direct_fit2_comp2_l1, ki_direct_fit2_comp2_l2, ki_direct_fit2_comp2_l3, ki_direct_fit2_comp2_l4, 
                           ki_direct_fit2_comp3_l1, ki_direct_fit2_comp3_l2, ki_direct_fit2_comp3_l3, ki_direct_fit2_comp3_l4, 
                           ki_direct_fit2_comp4_l1, ki_direct_fit2_comp4_l2, ki_direct_fit2_comp4_l3, ki_direct_fit2_comp4_l4)
setwd("..")

# *- 3.3.2.1-SI Sampling Distributions of Estimated Apparent IC50 and Ki ----------
# See Part "Random Seeds and expN"

# Random Seeds and expN ---------------------------------------------------
# All random seed tests are placed at the end of the program to prevent 
# interference of the main process's random seed setting

seed_list <- 1:100
seed_list_len <- length(seed_list)
expN_kstest_list <- c(1000, 5000, 10000, 100000, 200000)
expN_kstest_list_len <- length(expN_kstest_list)
expN_kstest_stringlist <- c("expN = 1000", "expN = 5000", "expN = 10000", "expN = 100000", "expN = 200000")

ks_parallel <- function(e, paralist, reslist) {
  com_name <- paralist[e]
  dir.create(com_name); setwd(com_name)
  strspl <- strsplit(com_name, "-")
  idx <- as.integer(strspl[[1]][1])
  idx2 <- strspl[[1]][2]
  data <- reslist[[idx]][[idx2]]
  datamean <- mean(data)
  datasd <- sd(data)
  ks_pvalue <- KSgeneral::cont_ks_test(data, "pnorm", datamean, datasd)$p.value
  setwd(".."); unlink(com_name, recursive = TRUE)
  ks_pvalue
}

# * 1.4 Sampling Distributions of Estimated [P]0 and Kd -----------------------

dir.create("Part1.4"); setwd("Part1.4")

kstest_df_f1_kd <- data.frame()
kstest_df_log10_f1_kd <- data.frame()

kstest_df_f1_p0 <- data.frame()
kstest_df_log10_f1_p0 <- data.frame()

control_df <- data.frame()

for (eidx in 1:expN_kstest_list_len) {
  
  reslist <- list()
  for (idx in 1:seed_list_len) {
    set.seed(seed_list[idx])

    plNoisy_array_kstest <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN_kstest_list[eidx], sd_level="medium")
    fitResult_kstest  <- parApply(clust, plNoisy_array_kstest, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

    f1_kd <- fitResult_kstest["f1_kd",]; log10_f1_kd <- log10(f1_kd)
    f1_p0 <- fitResult_kstest["f1_p0",]; log10_f1_p0 <- log10(f1_p0)

    #control
    set.seed(seed_list[idx])
    control_data <- rnorm(expN_kstest_list[eidx])
    
    reslist[[idx]] <- list(f1_kd=f1_kd, log10_f1_kd=log10_f1_kd,
                           f1_p0=f1_p0, log10_f1_p0=log10_f1_p0,
                           control_data=control_data)
    
  }
  colns <- c("f1_kd", "log10_f1_kd", "f1_p0", "log10_f1_p0",
             "control_data")
  tmpdf <- expand.grid(seed_list, colns)
  paralist <- paste(tmpdf$Var1, tmpdf$Var2, sep = "-")
  
  # parSapply preserves order
  parallelres <- parSapply(clust, seq(paralist), ks_parallel, paralist, reslist)
  
  resarray <- array(parallelres, dim=c(seed_list_len, length(colns)))
  colnames(resarray) <- colns
  kspvaluelist_f1_kd <- resarray[,"f1_kd"]
  kspvaluelist_log10_f1_kd <- resarray[,"log10_f1_kd"]
  
  kspvaluelist_f1_p0 <- resarray[,"f1_p0"]
  kspvaluelist_log10_f1_p0 <- resarray[,"log10_f1_p0"]
  
  kspvaluelist_control <- resarray[,"control_data"]
  
  #kd
  new_data_f1_kd <- data.frame(x=seed_list, y=kspvaluelist_f1_kd, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_f1_kd <- data.frame(x=seed_list, y=kspvaluelist_log10_f1_kd, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_f1_kd <- rbind(kstest_df_f1_kd, new_data_f1_kd)
  kstest_df_log10_f1_kd <- rbind(kstest_df_log10_f1_kd, new_data_log10_f1_kd)

  #p0
  new_data_f1_p0 <- data.frame(x=seed_list, y=kspvaluelist_f1_p0, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_f1_p0 <- data.frame(x=seed_list, y=kspvaluelist_log10_f1_p0, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_f1_p0 <- rbind(kstest_df_f1_p0, new_data_f1_p0)
  kstest_df_log10_f1_p0 <- rbind(kstest_df_log10_f1_p0, new_data_log10_f1_p0)

  # control
  new_data_control <- data.frame(x=seed_list, y=kspvaluelist_control, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  control_df <- rbind(control_df, new_data_control)
  
}
#kd
fig_seed_kspvalue_f1_kd <- ggplot(kstest_df_f1_kd, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Kd)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_f1_kd <- ggplot(kstest_df_log10_f1_kd, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Kd)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

#p0
fig_seed_kspvalue_f1_p0 <- ggplot(kstest_df_f1_p0, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value ([P]0)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_f1_p0 <- ggplot(kstest_df_log10_f1_p0, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 [P]0)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_control <- ggplot(control_df, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist,
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (normal control)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme


#kd
ggsave("seed_kstest_pvalue_f1_kd.png", fig_seed_kspvalue_f1_kd, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_f1_kd.svg", fig_seed_kspvalue_f1_kd, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_f1_kd.png", fig_seed_kspvalue_log10_f1_kd, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_f1_kd.svg", fig_seed_kspvalue_log10_f1_kd, width=7, height=7, units='cm')

#p0
ggsave("seed_kstest_pvalue_f1_p0.png", fig_seed_kspvalue_f1_p0, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_f1_p0.svg", fig_seed_kspvalue_f1_p0, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_f1_p0.png", fig_seed_kspvalue_log10_f1_p0, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_f1_p0.svg", fig_seed_kspvalue_log10_f1_p0, width=7, height=7, units='cm')

#control
ggsave("seed_kstest_pvalue_control.png", fig_seed_kspvalue_control, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_control.svg", fig_seed_kspvalue_control, width=7, height=7, units='cm')
setwd("..")

# * 1.5-SI Concentration Errors in [P]0 and [L]0 -----------------------------
# *- 1.5.2-SI Sampling Distributions of Estimated [P]0 and Kd --------------------

dir.create("Part1.5.2-SI"); setwd("Part1.5.2-SI")

kstest_df_f1_kd <- data.frame()
kstest_df_log10_f1_kd <- data.frame()

kstest_df_f1_p0 <- data.frame()
kstest_df_log10_f1_p0 <- data.frame()

for (eidx in 1:expN_kstest_list_len) {
  kspvaluelist_f1_kd <- vector(mode="numeric", length=seed_list_len)
  kspvaluelist_log10_f1_kd <- vector(mode="numeric", length=seed_list_len)
  
  kspvaluelist_f1_p0 <- vector(mode="numeric", length=seed_list_len)
  kspvaluelist_log10_f1_p0 <- vector(mode="numeric", length=seed_list_len)
  
  for (idx in 1:seed_list_len) {
    set.seed(seed_list[idx])

    plNoisy_array_kstest <- gen_plNoisy_array2(tl_num_large, paraTrue=para, repN=repN5, expN=expN_kstest_list[eidx], sd_level="medium", tp_cv="medium", tl_cv="medium")
    fitResult_kstest  <- parApply(clust, plNoisy_array_kstest, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
    
    f1_kd <- fitResult_kstest["f1_kd",]; log10_f1_kd <- log10(f1_kd)
    f1_p0 <- fitResult_kstest["f1_p0",]; log10_f1_p0 <- log10(f1_p0)
    
    reslist[[idx]] <- list(f1_kd=f1_kd, log10_f1_kd=log10_f1_kd,
                           f1_p0=f1_p0, log10_f1_p0=log10_f1_p0)
  }
  colns <- c("f1_kd", "log10_f1_kd", "f1_p0", "log10_f1_p0")
  tmpdf <- expand.grid(seed_list, colns)
  paralist <- paste(tmpdf$Var1, tmpdf$Var2, sep = "-")
  
  # parSapply preserves order
  parallelres <- parSapply(clust, seq(paralist), ks_parallel, paralist, reslist)
  
  resarray <- array(parallelres, dim=c(seed_list_len, length(colns)))
  colnames(resarray) <- colns
  kspvaluelist_f1_kd <- resarray[,"f1_kd"]
  kspvaluelist_log10_f1_kd <- resarray[,"log10_f1_kd"]
  
  kspvaluelist_f1_p0 <- resarray[,"f1_p0"]
  kspvaluelist_log10_f1_p0 <- resarray[,"log10_f1_p0"]
  
  #kd
  new_data_f1_kd <- data.frame(x=seed_list, y=kspvaluelist_f1_kd, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_f1_kd <- data.frame(x=seed_list, y=kspvaluelist_log10_f1_kd, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_f1_kd <- rbind(kstest_df_f1_kd, new_data_f1_kd)
  kstest_df_log10_f1_kd <- rbind(kstest_df_log10_f1_kd, new_data_log10_f1_kd)

  #p0
  new_data_f1_p0 <- data.frame(x=seed_list, y=kspvaluelist_f1_p0, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_f1_p0 <- data.frame(x=seed_list, y=kspvaluelist_log10_f1_p0, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_f1_p0 <- rbind(kstest_df_f1_p0, new_data_f1_p0)
  kstest_df_log10_f1_p0 <- rbind(kstest_df_log10_f1_p0, new_data_log10_f1_p0)
}

#kd
fig_seed_kspvalue_f1_kd <- ggplot(kstest_df_f1_kd, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Kd)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_f1_kd <- ggplot(kstest_df_log10_f1_kd, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Kd)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

#p0
fig_seed_kspvalue_f1_p0 <- ggplot(kstest_df_f1_p0, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value ([P]0)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_f1_p0 <- ggplot(kstest_df_log10_f1_p0, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 [P]0)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

#kd
ggsave("seed_kstest_pvalue_f1_kd.png", fig_seed_kspvalue_f1_kd, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_f1_kd.svg", fig_seed_kspvalue_f1_kd, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_f1_kd.png", fig_seed_kspvalue_log10_f1_kd, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_f1_kd.svg", fig_seed_kspvalue_log10_f1_kd, width=7, height=7, units='cm')

#p0
ggsave("seed_kstest_pvalue_f1_p0.png", fig_seed_kspvalue_f1_p0, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_f1_p0.svg", fig_seed_kspvalue_f1_p0, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_f1_p0.png", fig_seed_kspvalue_log10_f1_p0, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_f1_p0.svg", fig_seed_kspvalue_log10_f1_p0, width=7, height=7, units='cm')

setwd("..")

# * 2.2 Sampling Distributions of Estimated koff and kon --------------------

dir.create("Part2.2"); setwd("Part2.2")

kstest_df_koff_disso <- data.frame()
kstest_df_log10_koff_disso <- data.frame()

kstest_df_kon_asso <- data.frame()
kstest_df_log10_kon_asso <- data.frame()

kstest_df_koff_asso <- data.frame()
kstest_df_log10_koff_asso <- data.frame()

# kon nonlinear koff corrected
kstest_df_kon_nonlinear <- data.frame()
kstest_df_log10_kon_nonlinear <- data.frame()

# kon linear koff corrected
kstest_df_kon_linear <- data.frame()
kstest_df_log10_kon_linear <- data.frame()

kstest_df_kd_kin <- data.frame()
kstest_df_log10_kd_kin <- data.frame()

for (eidx in 1:expN_kstest_list_len) {
  
  reslist <- list()
  for (idx in 1:seed_list_len) {
    set.seed(seed_list[idx])

    plNoisy_array_kinetics2 <- gen_plNoisy_array_kin(pl_l0_list2, l0_list2, paraTrue=kinetics_para2, repN=repN, expN=expN_kstest_list[eidx], sd_level="low")
    fitResult_kin_asso_obs2  <- parApply(clust, plNoisy_array_kinetics2$data_array_asso, 3, fitfunc_kin, paraTrue=kinetics_para2, paraRatio=ratioList[2], type="asso_obs", l0_list2)
    fitResult_kin_disso2  <- parApply(clust, plNoisy_array_kinetics2$data_array_disso, 3, fitfunc_kin, paraTrue=kinetics_para2, paraRatio=ratioList[2], type="disso", l0_list2)
    fitResult_kin_asso_obs_conv2 <- sapply(1:expN_kstest_list[eidx], function(x) kon_koff_calc(l0_list2, fitResult_kin_asso_obs2[,x]))
    kon_linear_koff_corrected2 <- sapply(1:expN_kstest_list[eidx], function(x) kon_calc_koff_corrected(l0_list2, fitResult_kin_asso_obs2[,x], fitResult_kin_disso2[,x]))
    kon_nonlinear_koff_corrected2 <- sapply(1:expN_kstest_list[eidx], function(x) fitfunc_kin(plNoisy_array_kinetics2$data_array_asso[,,x], paraTrue=kinetics_para2, paraRatio=ratioList[2], type="asso", l0_list2, fitResult_kin_disso2[,x]))
    kd_kin2 <- fitResult_kin_disso2["koff_est",]/kon_nonlinear_koff_corrected2
    
    koff_disso <- fitResult_kin_disso2["koff_est",]; log10_koff_disso <- log10(koff_disso)
    kon_asso <- fitResult_kin_asso_obs_conv2["kon_est",]; log10_kon_asso <- log10(kon_asso)
    koff_asso <- fitResult_kin_asso_obs_conv2["koff_est",]; log10_koff_asso <- log10(koff_asso)
    kon_nonlinear <- kon_nonlinear_koff_corrected2; log10_kon_nonlinear <- log10(kon_nonlinear)
    kon_linear <- kon_linear_koff_corrected2; log10_kon_linear <- log10(kon_linear)
    kd_kin <- kd_kin2; log10_kd_kin <- log10(kd_kin)

    reslist[[idx]] <- list(koff_disso=koff_disso, log10_koff_disso=log10_koff_disso, 
                           kon_asso=kon_asso, log10_kon_asso=log10_kon_asso, 
                           koff_asso=koff_asso, log10_koff_asso=log10_koff_asso, 
                           kon_nonlinear=kon_nonlinear, log10_kon_nonlinear=log10_kon_nonlinear,
                           kon_linear=kon_linear, log10_kon_linear=log10_kon_linear,
                           kd_kin=kd_kin, log10_kd_kin=log10_kd_kin)
  }
  colns <- c("koff_disso", "log10_koff_disso", 
             "kon_asso", "log10_kon_asso", 
             "koff_asso", "log10_koff_asso", 
             "kon_nonlinear", "log10_kon_nonlinear", 
             "kon_linear", "log10_kon_linear", 
             "kd_kin", "log10_kd_kin")
  tmpdf <- expand.grid(seed_list, colns)
  paralist <- paste(tmpdf$Var1, tmpdf$Var2, sep = "-")
  
  # parSapply preserves order
  parallelres <- parSapply(clust, seq(paralist), ks_parallel, paralist, reslist)
  
  resarray <- array(parallelres, dim=c(seed_list_len, length(colns)))
  colnames(resarray) <- colns
  
  kspvaluelist_koff_disso <- resarray[,"koff_disso"]
  kspvaluelist_log10_koff_disso <- resarray[,"log10_koff_disso"]
  
  kspvaluelist_kon_asso <- resarray[,"kon_asso"]
  kspvaluelist_log10_kon_asso <- resarray[,"log10_kon_asso"]
  
  kspvaluelist_koff_asso <- resarray[,"koff_asso"]
  kspvaluelist_log10_koff_asso <- resarray[,"log10_koff_asso"]
  
  kspvaluelist_kon_nonlinear <- resarray[,"kon_nonlinear"]
  kspvaluelist_log10_kon_nonlinear <- resarray[,"log10_kon_nonlinear"]
  
  kspvaluelist_kon_linear <- resarray[,"kon_linear"]
  kspvaluelist_log10_kon_linear <- resarray[,"log10_kon_linear"]
  
  kspvaluelist_kd_kin <- resarray[,"kd_kin"]
  kspvaluelist_log10_kd_kin <- resarray[,"log10_kd_kin"]
  

  #koff_disso
  new_data_koff_disso <- data.frame(x=seed_list, y=kspvaluelist_koff_disso, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_koff_disso <- data.frame(x=seed_list, y=kspvaluelist_log10_koff_disso, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_koff_disso <- rbind(kstest_df_koff_disso, new_data_koff_disso)
  kstest_df_log10_koff_disso <- rbind(kstest_df_log10_koff_disso, new_data_log10_koff_disso)

  #kon_asso
  new_data_kon_asso <- data.frame(x=seed_list, y=kspvaluelist_kon_asso, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_kon_asso <- data.frame(x=seed_list, y=kspvaluelist_log10_kon_asso, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_kon_asso <- rbind(kstest_df_kon_asso, new_data_kon_asso)
  kstest_df_log10_kon_asso <- rbind(kstest_df_log10_kon_asso, new_data_log10_kon_asso)

  #koff_asso
  new_data_koff_asso <- data.frame(x=seed_list, y=kspvaluelist_koff_asso, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_koff_asso <- data.frame(x=seed_list, y=kspvaluelist_log10_koff_asso, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_koff_asso <- rbind(kstest_df_koff_asso, new_data_koff_asso)
  kstest_df_log10_koff_asso <- rbind(kstest_df_log10_koff_asso, new_data_log10_koff_asso)

  #kon_nonlinear
  new_data_kon_nonlinear <- data.frame(x=seed_list, y=kspvaluelist_kon_nonlinear, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_kon_nonlinear <- data.frame(x=seed_list, y=kspvaluelist_log10_kon_nonlinear, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_kon_nonlinear <- rbind(kstest_df_kon_nonlinear, new_data_kon_nonlinear)
  kstest_df_log10_kon_nonlinear <- rbind(kstest_df_log10_kon_nonlinear, new_data_log10_kon_nonlinear)

  #kon_linear
  new_data_kon_linear <- data.frame(x=seed_list, y=kspvaluelist_kon_linear, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_kon_linear <- data.frame(x=seed_list, y=kspvaluelist_log10_kon_linear, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_kon_linear <- rbind(kstest_df_kon_linear, new_data_kon_linear)
  kstest_df_log10_kon_linear <- rbind(kstest_df_log10_kon_linear, new_data_log10_kon_linear)

  #kd_kin
  new_data_kd_kin <- data.frame(x=seed_list, y=kspvaluelist_kd_kin, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_kd_kin <- data.frame(x=seed_list, y=kspvaluelist_log10_kd_kin, group=rep(expN_kstest_stringlist[eidx], seed_list_len))

  kstest_df_kd_kin <- rbind(kstest_df_kd_kin, new_data_kd_kin)
  kstest_df_log10_kd_kin <- rbind(kstest_df_log10_kd_kin, new_data_log10_kd_kin)
}

#koff_disso
fig_seed_kspvalue_koff_disso <- ggplot(kstest_df_koff_disso, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (koff)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_koff_disso <- ggplot(kstest_df_log10_koff_disso, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 koff)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

#kon_asso
fig_seed_kspvalue_kon_asso <- ggplot(kstest_df_kon_asso, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (kon)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_kon_asso <- ggplot(kstest_df_log10_kon_asso, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 kon)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

#koff_asso
fig_seed_kspvalue_koff_asso <- ggplot(kstest_df_koff_asso, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (koff)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_koff_asso <- ggplot(kstest_df_log10_koff_asso, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 koff)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

#kon_nonlinear
fig_seed_kspvalue_kon_nonlinear <- ggplot(kstest_df_kon_nonlinear, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (kon)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_kon_nonlinear <- ggplot(kstest_df_log10_kon_nonlinear, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 kon)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

#kon_linear
fig_seed_kspvalue_kon_linear <- ggplot(kstest_df_kon_linear, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (kon)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_kon_linear <- ggplot(kstest_df_log10_kon_linear, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 kon)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

#kd_kin
fig_seed_kspvalue_kd_kin <- ggplot(kstest_df_kd_kin, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Kd)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_kd_kin <- ggplot(kstest_df_log10_kd_kin, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Kd)") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

  
#koff_disso
ggsave("seed_kstest_pvalue_koff_disso.png", fig_seed_kspvalue_koff_disso, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_koff_disso.svg", fig_seed_kspvalue_koff_disso, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_koff_disso.png", fig_seed_kspvalue_log10_koff_disso, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_koff_disso.svg", fig_seed_kspvalue_log10_koff_disso, width=7, height=7, units='cm')

#kon_asso
ggsave("seed_kstest_pvalue_kon_asso.png", fig_seed_kspvalue_kon_asso, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_kon_asso.svg", fig_seed_kspvalue_kon_asso, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_kon_asso.png", fig_seed_kspvalue_log10_kon_asso, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_kon_asso.svg", fig_seed_kspvalue_log10_kon_asso, width=7, height=7, units='cm')

#koff_asso
ggsave("seed_kstest_pvalue_koff_asso.png", fig_seed_kspvalue_koff_asso, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_koff_asso.svg", fig_seed_kspvalue_koff_asso, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_koff_asso.png", fig_seed_kspvalue_log10_koff_asso, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_koff_asso.svg", fig_seed_kspvalue_log10_koff_asso, width=7, height=7, units='cm')

#kon_nonlinear
ggsave("seed_kstest_pvalue_kon_nonlinear.png", fig_seed_kspvalue_kon_nonlinear, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_kon_nonlinear.svg", fig_seed_kspvalue_kon_nonlinear, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_kon_nonlinear.png", fig_seed_kspvalue_log10_kon_nonlinear, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_kon_nonlinear.svg", fig_seed_kspvalue_log10_kon_nonlinear, width=7, height=7, units='cm')

#kon_linear
ggsave("seed_kstest_pvalue_kon_linear.png", fig_seed_kspvalue_kon_linear, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_kon_linear.svg", fig_seed_kspvalue_kon_linear, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_kon_linear.png", fig_seed_kspvalue_log10_kon_linear, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_kon_linear.svg", fig_seed_kspvalue_log10_kon_linear, width=7, height=7, units='cm')

#kd_kin
ggsave("seed_kstest_pvalue_kd_kin.png", fig_seed_kspvalue_kd_kin, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_kd_kin.svg", fig_seed_kspvalue_kd_kin, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_kd_kin.png", fig_seed_kspvalue_log10_kd_kin, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_kd_kin.svg", fig_seed_kspvalue_log10_kd_kin, width=7, height=7, units='cm')
setwd("..")

# * 3.2 Sampling Distributions of Estimated Apparent IC50 and Ki ------------

dir.create("Part3.2"); setwd("Part3.2")

# comp3
kstest_df_appic50_comp3_l1 <- data.frame()
kstest_df_log10_appic50_comp3_l1 <- data.frame()
kstest_df_ki_comp3_l1 <- data.frame()
kstest_df_log10_ki_comp3_l1 <- data.frame()

kstest_df_appic50_comp3_l2 <- data.frame()
kstest_df_log10_appic50_comp3_l2 <- data.frame()
kstest_df_ki_comp3_l2 <- data.frame()
kstest_df_log10_ki_comp3_l2 <- data.frame()

kstest_df_appic50_comp3_l3 <- data.frame()
kstest_df_log10_appic50_comp3_l3 <- data.frame()
kstest_df_ki_comp3_l3 <- data.frame()
kstest_df_log10_ki_comp3_l3 <- data.frame()

kstest_df_appic50_comp3_l4 <- data.frame()
kstest_df_log10_appic50_comp3_l4 <- data.frame()
kstest_df_ki_comp3_l4 <- data.frame()
kstest_df_log10_ki_comp3_l4 <- data.frame()

for (eidx in 1:expN_kstest_list_len) {
  
  reslist <- list()
  for (idx in 1:seed_list_len) {
    set.seed(seed_list[idx])
    
    # comp3
    kstest_plTheo_comp3_l1 <- with(para_comp3_l1, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
    kstest_plNoisy_comp3_l1 <- gen_plNoisy_array(kstest_plTheo_comp3_l1, paraTrue=para_comp3_l1, repN=repN5, expN=expN_kstest_list[eidx], sd_level="low")
    kstest_fitres_comp3_l1 <- parApply(clust, kstest_plNoisy_comp3_l1, 3, fitfunc_app_ic50, paraTrue=para_comp3_l1, repN=repN5, i_conc=ti_num_large)
    kstest_fitki_comp3_l1 <- parApply(clust, kstest_plNoisy_comp3_l1, 3, fitfunc_ki, paraTrue=para_comp3_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
    kstest_fitappic50_comp3_l1 <- kstest_fitres_comp3_l1["c",]
    kstest_ki_direct_fit_comp3_l1 <- kstest_fitki_comp3_l1["ki",]
    log10_appic50_comp3_l1 <- log10(kstest_fitappic50_comp3_l1)
    log10_ki_comp3_l1 <- log10(kstest_ki_direct_fit_comp3_l1)
    
    kstest_plTheo_comp3_l2 <- with(para_comp3_l2, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
    kstest_plNoisy_comp3_l2 <- gen_plNoisy_array(kstest_plTheo_comp3_l2, paraTrue=para_comp3_l2, repN=repN5, expN=expN_kstest_list[eidx], sd_level="low")
    kstest_fitres_comp3_l2 <- parApply(clust, kstest_plNoisy_comp3_l2, 3, fitfunc_app_ic50, paraTrue=para_comp3_l2, repN=repN5, i_conc=ti_num_large)
    kstest_fitki_comp3_l2 <- parApply(clust, kstest_plNoisy_comp3_l2, 3, fitfunc_ki, paraTrue=para_comp3_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
    kstest_fitappic50_comp3_l2 <- kstest_fitres_comp3_l2["c",]
    kstest_ki_direct_fit_comp3_l2 <- kstest_fitki_comp3_l2["ki",]
    log10_appic50_comp3_l2 <- log10(kstest_fitappic50_comp3_l2)
    log10_ki_comp3_l2 <- log10(kstest_ki_direct_fit_comp3_l2)
    
    kstest_plTheo_comp3_l3 <- with(para_comp3_l3, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
    kstest_plNoisy_comp3_l3 <- gen_plNoisy_array(kstest_plTheo_comp3_l3, paraTrue=para_comp3_l3, repN=repN5, expN=expN_kstest_list[eidx], sd_level="low")
    kstest_fitres_comp3_l3 <- parApply(clust, kstest_plNoisy_comp3_l3, 3, fitfunc_app_ic50, paraTrue=para_comp3_l3, repN=repN5, i_conc=ti_num_large)
    kstest_fitki_comp3_l3 <- parApply(clust, kstest_plNoisy_comp3_l3, 3, fitfunc_ki, paraTrue=para_comp3_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
    kstest_fitappic50_comp3_l3 <- kstest_fitres_comp3_l3["c",]
    kstest_ki_direct_fit_comp3_l3 <- kstest_fitki_comp3_l3["ki",]
    log10_appic50_comp3_l3 <- log10(kstest_fitappic50_comp3_l3)
    log10_ki_comp3_l3 <- log10(kstest_ki_direct_fit_comp3_l3)
    
    kstest_plTheo_comp3_l4 <- with(para_comp3_l4, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
    kstest_plNoisy_comp3_l4 <- gen_plNoisy_array(kstest_plTheo_comp3_l4, paraTrue=para_comp3_l4, repN=repN5, expN=expN_kstest_list[eidx], sd_level="low")
    kstest_fitres_comp3_l4 <- parApply(clust, kstest_plNoisy_comp3_l4, 3, fitfunc_app_ic50, paraTrue=para_comp3_l4, repN=repN5, i_conc=ti_num_large)
    kstest_fitki_comp3_l4 <- parApply(clust, kstest_plNoisy_comp3_l4, 3, fitfunc_ki, paraTrue=para_comp3_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
    kstest_fitappic50_comp3_l4 <- kstest_fitres_comp3_l4["c",]
    kstest_ki_direct_fit_comp3_l4 <- kstest_fitki_comp3_l4["ki",]
    log10_appic50_comp3_l4 <- log10(kstest_fitappic50_comp3_l4)
    log10_ki_comp3_l4 <- log10(kstest_ki_direct_fit_comp3_l4)
    
    reslist[[idx]] <- list(# comp3
                           appic50_comp3_l1=kstest_fitappic50_comp3_l1, log10_appic50_comp3_l1=log10_appic50_comp3_l1,
                           appic50_comp3_l2=kstest_fitappic50_comp3_l2, log10_appic50_comp3_l2=log10_appic50_comp3_l2,
                           appic50_comp3_l3=kstest_fitappic50_comp3_l3, log10_appic50_comp3_l3=log10_appic50_comp3_l3,
                           appic50_comp3_l4=kstest_fitappic50_comp3_l4, log10_appic50_comp3_l4=log10_appic50_comp3_l4,
                           ki_comp3_l1=kstest_ki_direct_fit_comp3_l1, log10_ki_comp3_l1=log10_ki_comp3_l1,
                           ki_comp3_l2=kstest_ki_direct_fit_comp3_l2, log10_ki_comp3_l2=log10_ki_comp3_l2,
                           ki_comp3_l3=kstest_ki_direct_fit_comp3_l3, log10_ki_comp3_l3=log10_ki_comp3_l3,
                           ki_comp3_l4=kstest_ki_direct_fit_comp3_l4, log10_ki_comp3_l4=log10_ki_comp3_l4)
  }
  colns <- c("appic50_comp3_l1", "log10_appic50_comp3_l1", "ki_comp3_l1", "log10_ki_comp3_l1",
             "appic50_comp3_l2", "log10_appic50_comp3_l2", "ki_comp3_l2", "log10_ki_comp3_l2",
             "appic50_comp3_l3", "log10_appic50_comp3_l3", "ki_comp3_l3", "log10_ki_comp3_l3",
             "appic50_comp3_l4", "log10_appic50_comp3_l4", "ki_comp3_l4", "log10_ki_comp3_l4")
  tmpdf <- expand.grid(seed_list, colns)
  paralist <- paste(tmpdf$Var1, tmpdf$Var2, sep = "-")
  
  # parSapply preserves order
  parallelres <- parSapply(clust, seq(paralist), ks_parallel, paralist, reslist)
  
  resarray <- array(parallelres, dim=c(seed_list_len, length(colns)))
  colnames(resarray) <- colns
  
  # comp3
  kspvaluelist_appic50_comp3_l1 <- resarray[,"appic50_comp3_l1"]
  kspvaluelist_log10_appic50_comp3_l1 <- resarray[,"log10_appic50_comp3_l1"]
  kspvaluelist_ki_comp3_l1 <- resarray[,"ki_comp3_l1"]
  kspvaluelist_log10_ki_comp3_l1 <- resarray[,"log10_ki_comp3_l1"]
  
  kspvaluelist_appic50_comp3_l2 <- resarray[,"appic50_comp3_l2"]
  kspvaluelist_log10_appic50_comp3_l2 <- resarray[,"log10_appic50_comp3_l2"]
  kspvaluelist_ki_comp3_l2 <- resarray[,"ki_comp3_l2"]
  kspvaluelist_log10_ki_comp3_l2 <- resarray[,"log10_ki_comp3_l2"]
  
  kspvaluelist_appic50_comp3_l3 <- resarray[,"appic50_comp3_l3"]
  kspvaluelist_log10_appic50_comp3_l3 <- resarray[,"log10_appic50_comp3_l3"]
  kspvaluelist_ki_comp3_l3 <- resarray[,"ki_comp3_l3"]
  kspvaluelist_log10_ki_comp3_l3 <- resarray[,"log10_ki_comp3_l3"]
  
  kspvaluelist_appic50_comp3_l4 <- resarray[,"appic50_comp3_l4"]
  kspvaluelist_log10_appic50_comp3_l4 <- resarray[,"log10_appic50_comp3_l4"]
  kspvaluelist_ki_comp3_l4 <- resarray[,"ki_comp3_l4"]
  kspvaluelist_log10_ki_comp3_l4 <- resarray[,"log10_ki_comp3_l4"]
  
  # comp3
  new_data_appic50_comp3_l1 <- data.frame(x=seed_list, y=kspvaluelist_appic50_comp3_l1, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_appic50_comp3_l1 <- data.frame(x=seed_list, y=kspvaluelist_log10_appic50_comp3_l1, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_appic50_comp3_l1 <- rbind(kstest_df_appic50_comp3_l1, new_data_appic50_comp3_l1)
  kstest_df_log10_appic50_comp3_l1 <- rbind(kstest_df_log10_appic50_comp3_l1, new_data_log10_appic50_comp3_l1)
  new_data_ki_comp3_l1 <- data.frame(x=seed_list, y=kspvaluelist_ki_comp3_l1, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_ki_comp3_l1 <- data.frame(x=seed_list, y=kspvaluelist_log10_ki_comp3_l1, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_ki_comp3_l1 <- rbind(kstest_df_ki_comp3_l1, new_data_ki_comp3_l1)
  kstest_df_log10_ki_comp3_l1 <- rbind(kstest_df_log10_ki_comp3_l1, new_data_log10_ki_comp3_l1)
  
  new_data_appic50_comp3_l2 <- data.frame(x=seed_list, y=kspvaluelist_appic50_comp3_l2, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_appic50_comp3_l2 <- data.frame(x=seed_list, y=kspvaluelist_log10_appic50_comp3_l2, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_appic50_comp3_l2 <- rbind(kstest_df_appic50_comp3_l2, new_data_appic50_comp3_l2)
  kstest_df_log10_appic50_comp3_l2 <- rbind(kstest_df_log10_appic50_comp3_l2, new_data_log10_appic50_comp3_l2)
  new_data_ki_comp3_l2 <- data.frame(x=seed_list, y=kspvaluelist_ki_comp3_l2, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_ki_comp3_l2 <- data.frame(x=seed_list, y=kspvaluelist_log10_ki_comp3_l2, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_ki_comp3_l2 <- rbind(kstest_df_ki_comp3_l2, new_data_ki_comp3_l2)
  kstest_df_log10_ki_comp3_l2 <- rbind(kstest_df_log10_ki_comp3_l2, new_data_log10_ki_comp3_l2)
  
  new_data_appic50_comp3_l3 <- data.frame(x=seed_list, y=kspvaluelist_appic50_comp3_l3, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_appic50_comp3_l3 <- data.frame(x=seed_list, y=kspvaluelist_log10_appic50_comp3_l3, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_appic50_comp3_l3 <- rbind(kstest_df_appic50_comp3_l3, new_data_appic50_comp3_l3)
  kstest_df_log10_appic50_comp3_l3 <- rbind(kstest_df_log10_appic50_comp3_l3, new_data_log10_appic50_comp3_l3)
  new_data_ki_comp3_l3 <- data.frame(x=seed_list, y=kspvaluelist_ki_comp3_l3, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_ki_comp3_l3 <- data.frame(x=seed_list, y=kspvaluelist_log10_ki_comp3_l3, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_ki_comp3_l3 <- rbind(kstest_df_ki_comp3_l3, new_data_ki_comp3_l3)
  kstest_df_log10_ki_comp3_l3 <- rbind(kstest_df_log10_ki_comp3_l3, new_data_log10_ki_comp3_l3)
  
  new_data_appic50_comp3_l4 <- data.frame(x=seed_list, y=kspvaluelist_appic50_comp3_l4, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_appic50_comp3_l4 <- data.frame(x=seed_list, y=kspvaluelist_log10_appic50_comp3_l4, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_appic50_comp3_l4 <- rbind(kstest_df_appic50_comp3_l4, new_data_appic50_comp3_l4)
  kstest_df_log10_appic50_comp3_l4 <- rbind(kstest_df_log10_appic50_comp3_l4, new_data_log10_appic50_comp3_l4)
  new_data_ki_comp3_l4 <- data.frame(x=seed_list, y=kspvaluelist_ki_comp3_l4, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_ki_comp3_l4 <- data.frame(x=seed_list, y=kspvaluelist_log10_ki_comp3_l4, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_ki_comp3_l4 <- rbind(kstest_df_ki_comp3_l4, new_data_ki_comp3_l4)
  kstest_df_log10_ki_comp3_l4 <- rbind(kstest_df_log10_ki_comp3_l4, new_data_log10_ki_comp3_l4)
}

# comp3
# l1
fig_seed_kspvalue_appic50_comp3_l1 <- ggplot(kstest_df_appic50_comp3_l1, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (app. IC50)", title="comp3_l1") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_appic50_comp3_l1 <- ggplot(kstest_df_log10_appic50_comp3_l1, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 app. IC50)", title="comp3_l1") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_ki_comp3_l1 <- ggplot(kstest_df_ki_comp3_l1, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Ki)", title="comp3_l1") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_ki_comp3_l1 <- ggplot(kstest_df_log10_ki_comp3_l1, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Ki)", title="comp3_l1") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

# l2
fig_seed_kspvalue_appic50_comp3_l2 <- ggplot(kstest_df_appic50_comp3_l2, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (app. IC50)", title="comp3_l2") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_appic50_comp3_l2 <- ggplot(kstest_df_log10_appic50_comp3_l2, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 app. IC50)", title="comp3_l2") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_ki_comp3_l2 <- ggplot(kstest_df_ki_comp3_l2, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Ki)", title="comp3_l2") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_ki_comp3_l2 <- ggplot(kstest_df_log10_ki_comp3_l2, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Ki)", title="comp3_l2") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

# l3
fig_seed_kspvalue_appic50_comp3_l3 <- ggplot(kstest_df_appic50_comp3_l3, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (app. IC50)", title="comp3_l3") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_appic50_comp3_l3 <- ggplot(kstest_df_log10_appic50_comp3_l3, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 app. IC50)", title="comp3_l3") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_ki_comp3_l3 <- ggplot(kstest_df_ki_comp3_l3, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Ki)", title="comp3_l3") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_ki_comp3_l3 <- ggplot(kstest_df_log10_ki_comp3_l3, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Ki)", title="comp3_l3") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

# l4
fig_seed_kspvalue_appic50_comp3_l4 <- ggplot(kstest_df_appic50_comp3_l4, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (app. IC50)", title="comp3_l4") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_appic50_comp3_l4 <- ggplot(kstest_df_log10_appic50_comp3_l4, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 app. IC50)", title="comp3_l4") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_ki_comp3_l4 <- ggplot(kstest_df_ki_comp3_l4, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Ki)", title="comp3_l4") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_ki_comp3_l4 <- ggplot(kstest_df_log10_ki_comp3_l4, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Ki)", title="comp3_l4") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

# comp3
ggsave("seed_kstest_pvalue_appic50_comp3_l1.png", fig_seed_kspvalue_appic50_comp3_l1, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l1.png", fig_seed_kspvalue_log10_appic50_comp3_l1, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l1.png", fig_seed_kspvalue_ki_comp3_l1, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l1.png", fig_seed_kspvalue_log10_ki_comp3_l1, dpi=800, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l2.png", fig_seed_kspvalue_appic50_comp3_l2, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l2.png", fig_seed_kspvalue_log10_appic50_comp3_l2, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l2.png", fig_seed_kspvalue_ki_comp3_l2, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l2.png", fig_seed_kspvalue_log10_ki_comp3_l2, dpi=800, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l3.png", fig_seed_kspvalue_appic50_comp3_l3, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l3.png", fig_seed_kspvalue_log10_appic50_comp3_l3, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l3.png", fig_seed_kspvalue_ki_comp3_l3, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l3.png", fig_seed_kspvalue_log10_ki_comp3_l3, dpi=800, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l4.png", fig_seed_kspvalue_appic50_comp3_l4, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l4.png", fig_seed_kspvalue_log10_appic50_comp3_l4, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l4.png", fig_seed_kspvalue_ki_comp3_l4, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l4.png", fig_seed_kspvalue_log10_ki_comp3_l4, dpi=800, width=7, height=7, units='cm')

# comp3
ggsave("seed_kstest_pvalue_appic50_comp3_l1.svg", fig_seed_kspvalue_appic50_comp3_l1, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l1.svg", fig_seed_kspvalue_log10_appic50_comp3_l1, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l1.svg", fig_seed_kspvalue_ki_comp3_l1, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l1.svg", fig_seed_kspvalue_log10_ki_comp3_l1, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l2.svg", fig_seed_kspvalue_appic50_comp3_l2, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l2.svg", fig_seed_kspvalue_log10_appic50_comp3_l2, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l2.svg", fig_seed_kspvalue_ki_comp3_l2, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l2.svg", fig_seed_kspvalue_log10_ki_comp3_l2, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l3.svg", fig_seed_kspvalue_appic50_comp3_l3, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l3.svg", fig_seed_kspvalue_log10_appic50_comp3_l3, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l3.svg", fig_seed_kspvalue_ki_comp3_l3, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l3.svg", fig_seed_kspvalue_log10_ki_comp3_l3, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l4.svg", fig_seed_kspvalue_appic50_comp3_l4, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l4.svg", fig_seed_kspvalue_log10_appic50_comp3_l4, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l4.svg", fig_seed_kspvalue_ki_comp3_l4, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l4.svg", fig_seed_kspvalue_log10_ki_comp3_l4, width=7, height=7, units='cm')

setwd("..")

# * 3.3-SI Concentration Errors in [P]0, [L]0, and [I]0 ----------------------
# *- 3.3.2-SI Sampling Distributions of Estimated Apparent IC50 and Ki ----------
# tp_cv="low", tl_cv="low", ti_cv="low"

dir.create("Part3.3.2-SI"); setwd("Part3.3.2-SI")

# comp3
kstest_df_appic50_comp3_l1 <- data.frame()
kstest_df_log10_appic50_comp3_l1 <- data.frame()
kstest_df_ki_comp3_l1 <- data.frame()
kstest_df_log10_ki_comp3_l1 <- data.frame()

kstest_df_appic50_comp3_l2 <- data.frame()
kstest_df_log10_appic50_comp3_l2 <- data.frame()
kstest_df_ki_comp3_l2 <- data.frame()
kstest_df_log10_ki_comp3_l2 <- data.frame()

kstest_df_appic50_comp3_l3 <- data.frame()
kstest_df_log10_appic50_comp3_l3 <- data.frame()
kstest_df_ki_comp3_l3 <- data.frame()
kstest_df_log10_ki_comp3_l3 <- data.frame()

kstest_df_appic50_comp3_l4 <- data.frame()
kstest_df_log10_appic50_comp3_l4 <- data.frame()
kstest_df_ki_comp3_l4 <- data.frame()
kstest_df_log10_ki_comp3_l4 <- data.frame()

for (eidx in 1:expN_kstest_list_len) {
  
  reslist <- list()
  for (idx in 1:seed_list_len) {
    set.seed(seed_list[idx])
    
    # comp3
    kstest_plTheo_comp3_l1 <- with(para_comp3_l1, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
    kstest_plNoisy_comp3_l1 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp3_l1, repN=repN5, expN=expN_kstest_list[eidx], sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
    kstest_fitres_comp3_l1 <- parApply(clust, kstest_plNoisy_comp3_l1, 3, fitfunc_app_ic50, paraTrue=para_comp3_l1, repN=repN5, i_conc=ti_num_large)
    kstest_fitki_comp3_l1 <- parApply(clust, kstest_plNoisy_comp3_l1, 3, fitfunc_ki, paraTrue=para_comp3_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
    kstest_fitappic50_comp3_l1 <- kstest_fitres_comp3_l1["c",]
    kstest_ki_direct_fit_comp3_l1 <- kstest_fitki_comp3_l1["ki",]
    log10_appic50_comp3_l1 <- log10(kstest_fitappic50_comp3_l1)
    log10_ki_comp3_l1 <- log10(kstest_ki_direct_fit_comp3_l1)
    
    kstest_plTheo_comp3_l2 <- with(para_comp3_l2, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
    kstest_plNoisy_comp3_l2 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp3_l2, repN=repN5, expN=expN_kstest_list[eidx], sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
    kstest_fitres_comp3_l2 <- parApply(clust, kstest_plNoisy_comp3_l2, 3, fitfunc_app_ic50, paraTrue=para_comp3_l2, repN=repN5, i_conc=ti_num_large)
    kstest_fitki_comp3_l2 <- parApply(clust, kstest_plNoisy_comp3_l2, 3, fitfunc_ki, paraTrue=para_comp3_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
    kstest_fitappic50_comp3_l2 <- kstest_fitres_comp3_l2["c",]
    kstest_ki_direct_fit_comp3_l2 <- kstest_fitki_comp3_l2["ki",]
    log10_appic50_comp3_l2 <- log10(kstest_fitappic50_comp3_l2)
    log10_ki_comp3_l2 <- log10(kstest_ki_direct_fit_comp3_l2)
    
    kstest_plTheo_comp3_l3 <- with(para_comp3_l3, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
    kstest_plNoisy_comp3_l3 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp3_l3, repN=repN5, expN=expN_kstest_list[eidx], sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
    kstest_fitres_comp3_l3 <- parApply(clust, kstest_plNoisy_comp3_l3, 3, fitfunc_app_ic50, paraTrue=para_comp3_l3, repN=repN5, i_conc=ti_num_large)
    kstest_fitki_comp3_l3 <- parApply(clust, kstest_plNoisy_comp3_l3, 3, fitfunc_ki, paraTrue=para_comp3_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
    kstest_fitappic50_comp3_l3 <- kstest_fitres_comp3_l3["c",]
    kstest_ki_direct_fit_comp3_l3 <- kstest_fitki_comp3_l3["ki",]
    log10_appic50_comp3_l3 <- log10(kstest_fitappic50_comp3_l3)
    log10_ki_comp3_l3 <- log10(kstest_ki_direct_fit_comp3_l3)
    
    kstest_plTheo_comp3_l4 <- with(para_comp3_l4, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
    kstest_plNoisy_comp3_l4 <- gen_plNoisy_array3(ti_num_large, paraTrue=para_comp3_l4, repN=repN5, expN=expN_kstest_list[eidx], sd_level="low", tp_cv="low", tl_cv="low", ti_cv="low")
    kstest_fitres_comp3_l4 <- parApply(clust, kstest_plNoisy_comp3_l4, 3, fitfunc_app_ic50, paraTrue=para_comp3_l4, repN=repN5, i_conc=ti_num_large)
    kstest_fitki_comp3_l4 <- parApply(clust, kstest_plNoisy_comp3_l4, 3, fitfunc_ki, paraTrue=para_comp3_l4, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
    kstest_fitappic50_comp3_l4 <- kstest_fitres_comp3_l4["c",]
    kstest_ki_direct_fit_comp3_l4 <- kstest_fitki_comp3_l4["ki",]
    log10_appic50_comp3_l4 <- log10(kstest_fitappic50_comp3_l4)
    log10_ki_comp3_l4 <- log10(kstest_ki_direct_fit_comp3_l4)
    
    reslist[[idx]] <- list(
      # comp3
      appic50_comp3_l1=kstest_fitappic50_comp3_l1, log10_appic50_comp3_l1=log10_appic50_comp3_l1,
      appic50_comp3_l2=kstest_fitappic50_comp3_l2, log10_appic50_comp3_l2=log10_appic50_comp3_l2,
      appic50_comp3_l3=kstest_fitappic50_comp3_l3, log10_appic50_comp3_l3=log10_appic50_comp3_l3,
      appic50_comp3_l4=kstest_fitappic50_comp3_l4, log10_appic50_comp3_l4=log10_appic50_comp3_l4,
      ki_comp3_l1=kstest_ki_direct_fit_comp3_l1, log10_ki_comp3_l1=log10_ki_comp3_l1,
      ki_comp3_l2=kstest_ki_direct_fit_comp3_l2, log10_ki_comp3_l2=log10_ki_comp3_l2,
      ki_comp3_l3=kstest_ki_direct_fit_comp3_l3, log10_ki_comp3_l3=log10_ki_comp3_l3,
      ki_comp3_l4=kstest_ki_direct_fit_comp3_l4, log10_ki_comp3_l4=log10_ki_comp3_l4)
  }
  colns <- c("appic50_comp3_l1", "log10_appic50_comp3_l1", "ki_comp3_l1", "log10_ki_comp3_l1",
             "appic50_comp3_l2", "log10_appic50_comp3_l2", "ki_comp3_l2", "log10_ki_comp3_l2",
             "appic50_comp3_l3", "log10_appic50_comp3_l3", "ki_comp3_l3", "log10_ki_comp3_l3",
             "appic50_comp3_l4", "log10_appic50_comp3_l4", "ki_comp3_l4", "log10_ki_comp3_l4")
  tmpdf <- expand.grid(seed_list, colns)
  paralist <- paste(tmpdf$Var1, tmpdf$Var2, sep = "-")
  
  # parSapply preserves order
  parallelres <- parSapply(clust, seq(paralist), ks_parallel, paralist, reslist)
  
  resarray <- array(parallelres, dim=c(seed_list_len, length(colns)))
  colnames(resarray) <- colns
  
  # comp3
  kspvaluelist_appic50_comp3_l1 <- resarray[,"appic50_comp3_l1"]
  kspvaluelist_log10_appic50_comp3_l1 <- resarray[,"log10_appic50_comp3_l1"]
  kspvaluelist_ki_comp3_l1 <- resarray[,"ki_comp3_l1"]
  kspvaluelist_log10_ki_comp3_l1 <- resarray[,"log10_ki_comp3_l1"]
  
  kspvaluelist_appic50_comp3_l2 <- resarray[,"appic50_comp3_l2"]
  kspvaluelist_log10_appic50_comp3_l2 <- resarray[,"log10_appic50_comp3_l2"]
  kspvaluelist_ki_comp3_l2 <- resarray[,"ki_comp3_l2"]
  kspvaluelist_log10_ki_comp3_l2 <- resarray[,"log10_ki_comp3_l2"]
  
  kspvaluelist_appic50_comp3_l3 <- resarray[,"appic50_comp3_l3"]
  kspvaluelist_log10_appic50_comp3_l3 <- resarray[,"log10_appic50_comp3_l3"]
  kspvaluelist_ki_comp3_l3 <- resarray[,"ki_comp3_l3"]
  kspvaluelist_log10_ki_comp3_l3 <- resarray[,"log10_ki_comp3_l3"]
  
  kspvaluelist_appic50_comp3_l4 <- resarray[,"appic50_comp3_l4"]
  kspvaluelist_log10_appic50_comp3_l4 <- resarray[,"log10_appic50_comp3_l4"]
  kspvaluelist_ki_comp3_l4 <- resarray[,"ki_comp3_l4"]
  kspvaluelist_log10_ki_comp3_l4 <- resarray[,"log10_ki_comp3_l4"]
  
  # comp3
  new_data_appic50_comp3_l1 <- data.frame(x=seed_list, y=kspvaluelist_appic50_comp3_l1, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_appic50_comp3_l1 <- data.frame(x=seed_list, y=kspvaluelist_log10_appic50_comp3_l1, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_appic50_comp3_l1 <- rbind(kstest_df_appic50_comp3_l1, new_data_appic50_comp3_l1)
  kstest_df_log10_appic50_comp3_l1 <- rbind(kstest_df_log10_appic50_comp3_l1, new_data_log10_appic50_comp3_l1)
  new_data_ki_comp3_l1 <- data.frame(x=seed_list, y=kspvaluelist_ki_comp3_l1, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_ki_comp3_l1 <- data.frame(x=seed_list, y=kspvaluelist_log10_ki_comp3_l1, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_ki_comp3_l1 <- rbind(kstest_df_ki_comp3_l1, new_data_ki_comp3_l1)
  kstest_df_log10_ki_comp3_l1 <- rbind(kstest_df_log10_ki_comp3_l1, new_data_log10_ki_comp3_l1)
  
  new_data_appic50_comp3_l2 <- data.frame(x=seed_list, y=kspvaluelist_appic50_comp3_l2, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_appic50_comp3_l2 <- data.frame(x=seed_list, y=kspvaluelist_log10_appic50_comp3_l2, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_appic50_comp3_l2 <- rbind(kstest_df_appic50_comp3_l2, new_data_appic50_comp3_l2)
  kstest_df_log10_appic50_comp3_l2 <- rbind(kstest_df_log10_appic50_comp3_l2, new_data_log10_appic50_comp3_l2)
  new_data_ki_comp3_l2 <- data.frame(x=seed_list, y=kspvaluelist_ki_comp3_l2, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_ki_comp3_l2 <- data.frame(x=seed_list, y=kspvaluelist_log10_ki_comp3_l2, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_ki_comp3_l2 <- rbind(kstest_df_ki_comp3_l2, new_data_ki_comp3_l2)
  kstest_df_log10_ki_comp3_l2 <- rbind(kstest_df_log10_ki_comp3_l2, new_data_log10_ki_comp3_l2)
  
  new_data_appic50_comp3_l3 <- data.frame(x=seed_list, y=kspvaluelist_appic50_comp3_l3, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_appic50_comp3_l3 <- data.frame(x=seed_list, y=kspvaluelist_log10_appic50_comp3_l3, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_appic50_comp3_l3 <- rbind(kstest_df_appic50_comp3_l3, new_data_appic50_comp3_l3)
  kstest_df_log10_appic50_comp3_l3 <- rbind(kstest_df_log10_appic50_comp3_l3, new_data_log10_appic50_comp3_l3)
  new_data_ki_comp3_l3 <- data.frame(x=seed_list, y=kspvaluelist_ki_comp3_l3, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_ki_comp3_l3 <- data.frame(x=seed_list, y=kspvaluelist_log10_ki_comp3_l3, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_ki_comp3_l3 <- rbind(kstest_df_ki_comp3_l3, new_data_ki_comp3_l3)
  kstest_df_log10_ki_comp3_l3 <- rbind(kstest_df_log10_ki_comp3_l3, new_data_log10_ki_comp3_l3)
  
  new_data_appic50_comp3_l4 <- data.frame(x=seed_list, y=kspvaluelist_appic50_comp3_l4, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_appic50_comp3_l4 <- data.frame(x=seed_list, y=kspvaluelist_log10_appic50_comp3_l4, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_appic50_comp3_l4 <- rbind(kstest_df_appic50_comp3_l4, new_data_appic50_comp3_l4)
  kstest_df_log10_appic50_comp3_l4 <- rbind(kstest_df_log10_appic50_comp3_l4, new_data_log10_appic50_comp3_l4)
  new_data_ki_comp3_l4 <- data.frame(x=seed_list, y=kspvaluelist_ki_comp3_l4, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  new_data_log10_ki_comp3_l4 <- data.frame(x=seed_list, y=kspvaluelist_log10_ki_comp3_l4, group=rep(expN_kstest_stringlist[eidx], seed_list_len))
  kstest_df_ki_comp3_l4 <- rbind(kstest_df_ki_comp3_l4, new_data_ki_comp3_l4)
  kstest_df_log10_ki_comp3_l4 <- rbind(kstest_df_log10_ki_comp3_l4, new_data_log10_ki_comp3_l4)
}

# comp3
# l1
fig_seed_kspvalue_appic50_comp3_l1 <- ggplot(kstest_df_appic50_comp3_l1, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (app. IC50)", title="comp3_l1") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_appic50_comp3_l1 <- ggplot(kstest_df_log10_appic50_comp3_l1, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 app. IC50)", title="comp3_l1") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_ki_comp3_l1 <- ggplot(kstest_df_ki_comp3_l1, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Ki)", title="comp3_l1") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_ki_comp3_l1 <- ggplot(kstest_df_log10_ki_comp3_l1, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Ki)", title="comp3_l1") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

# l2
fig_seed_kspvalue_appic50_comp3_l2 <- ggplot(kstest_df_appic50_comp3_l2, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (app. IC50)", title="comp3_l2") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_appic50_comp3_l2 <- ggplot(kstest_df_log10_appic50_comp3_l2, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 app. IC50)", title="comp3_l2") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_ki_comp3_l2 <- ggplot(kstest_df_ki_comp3_l2, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Ki)", title="comp3_l2") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_ki_comp3_l2 <- ggplot(kstest_df_log10_ki_comp3_l2, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Ki)", title="comp3_l2") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

# l3
fig_seed_kspvalue_appic50_comp3_l3 <- ggplot(kstest_df_appic50_comp3_l3, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (app. IC50)", title="comp3_l3") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_appic50_comp3_l3 <- ggplot(kstest_df_log10_appic50_comp3_l3, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 app. IC50)", title="comp3_l3") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_ki_comp3_l3 <- ggplot(kstest_df_ki_comp3_l3, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Ki)", title="comp3_l3") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_ki_comp3_l3 <- ggplot(kstest_df_log10_ki_comp3_l3, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Ki)", title="comp3_l3") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

# l4
fig_seed_kspvalue_appic50_comp3_l4 <- ggplot(kstest_df_appic50_comp3_l4, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (app. IC50)", title="comp3_l4") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_appic50_comp3_l4 <- ggplot(kstest_df_log10_appic50_comp3_l4, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 app. IC50)", title="comp3_l4") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_ki_comp3_l4 <- ggplot(kstest_df_ki_comp3_l4, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (Ki)", title="comp3_l4") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

fig_seed_kspvalue_log10_ki_comp3_l4 <- ggplot(kstest_df_log10_ki_comp3_l4, aes(x=x, y=y, color=group)) +
  geom_line(aes(group=group), linewidth=0.2) +
  scale_color_manual(breaks=expN_kstest_stringlist, 
                     values=f1_color[1:expN_kstest_list_len]) +
  geom_hline(yintercept=0.05, color="red", linetype="dashed", linewidth=0.2) +
  labs(x="Random seed", y="KS Test P-value (log10 Ki)", title="comp3_l4") +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100),
                     labels = c("1", "25", "50", "75", "100")) +
  coord_cartesian(ylim=c(0, 1.45)) +
  legend_theme

# comp3
ggsave("seed_kstest_pvalue_appic50_comp3_l1.png", fig_seed_kspvalue_appic50_comp3_l1, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l1.png", fig_seed_kspvalue_log10_appic50_comp3_l1, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l1.png", fig_seed_kspvalue_ki_comp3_l1, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l1.png", fig_seed_kspvalue_log10_ki_comp3_l1, dpi=800, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l2.png", fig_seed_kspvalue_appic50_comp3_l2, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l2.png", fig_seed_kspvalue_log10_appic50_comp3_l2, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l2.png", fig_seed_kspvalue_ki_comp3_l2, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l2.png", fig_seed_kspvalue_log10_ki_comp3_l2, dpi=800, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l3.png", fig_seed_kspvalue_appic50_comp3_l3, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l3.png", fig_seed_kspvalue_log10_appic50_comp3_l3, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l3.png", fig_seed_kspvalue_ki_comp3_l3, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l3.png", fig_seed_kspvalue_log10_ki_comp3_l3, dpi=800, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l4.png", fig_seed_kspvalue_appic50_comp3_l4, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l4.png", fig_seed_kspvalue_log10_appic50_comp3_l4, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l4.png", fig_seed_kspvalue_ki_comp3_l4, dpi=800, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l4.png", fig_seed_kspvalue_log10_ki_comp3_l4, dpi=800, width=7, height=7, units='cm')

# comp3
ggsave("seed_kstest_pvalue_appic50_comp3_l1.svg", fig_seed_kspvalue_appic50_comp3_l1, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l1.svg", fig_seed_kspvalue_log10_appic50_comp3_l1, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l1.svg", fig_seed_kspvalue_ki_comp3_l1, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l1.svg", fig_seed_kspvalue_log10_ki_comp3_l1, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l2.svg", fig_seed_kspvalue_appic50_comp3_l2, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l2.svg", fig_seed_kspvalue_log10_appic50_comp3_l2, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l2.svg", fig_seed_kspvalue_ki_comp3_l2, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l2.svg", fig_seed_kspvalue_log10_ki_comp3_l2, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l3.svg", fig_seed_kspvalue_appic50_comp3_l3, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l3.svg", fig_seed_kspvalue_log10_appic50_comp3_l3, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l3.svg", fig_seed_kspvalue_ki_comp3_l3, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l3.svg", fig_seed_kspvalue_log10_ki_comp3_l3, width=7, height=7, units='cm')

ggsave("seed_kstest_pvalue_appic50_comp3_l4.svg", fig_seed_kspvalue_appic50_comp3_l4, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_appic50_comp3_l4.svg", fig_seed_kspvalue_log10_appic50_comp3_l4, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_ki_comp3_l4.svg", fig_seed_kspvalue_ki_comp3_l4, width=7, height=7, units='cm')
ggsave("seed_kstest_pvalue_log10_ki_comp3_l4.svg", fig_seed_kspvalue_log10_ki_comp3_l4, width=7, height=7, units='cm')

setwd("..")

# The Ending Work -----------------------------------------------------------------

stopCluster(clust)
