# * Libraries and Settings of the Program -----------------------------------
library(minpack.lm)
library(drc)
library(grid)
library(ggplot2)
library(plyr) # required by GeomSplitViolin
library(RColorBrewer)
library(extrafont) 
library(svglite)
library(scales)
library(parallel)

options(warn=2) # all warnings are turned into errors
theme_set(theme_bw(base_size=10, base_family="Arial",))

set.seed(59403)

# the number of experiments
expN <- 100000

# the number of replicates
repN   <- 1
repN5  <- 5
repN10 <- 10

# Initial values for estimating parameters
# setting ratio to 0.01, 5 or 10 is bad, this will generate following error when fitting 
# error in sqrt((p0 + x + kd)^2 - 4 * p0 * x) :  (converted from warning) NaNs produced
# setting ratio 0.25 and sd_level="high" can also generate the same error
ratioList <- c(0.25, 0.5, 1, 2)

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
   draw_group=function(self, data, ..., nudge=0, quantiles=NULL) {
     # Original function by Jan Gleixner (@jan-glx)
     # Adjustments by Wouter van der Bijl (@Axeman)
     data <- transform(data, xminv=x - violinwidth * (x - xmin), xmaxv=x + violinwidth * (xmax - x))
     grp <- data[1, "group"]
     newdata <- plyr::arrange(transform(data, x=if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
     newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
     newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
     # now nudge them apart
     newdata$x <- ifelse(newdata$group %% 2 == 1, newdata$x - nudge, newdata$x + nudge)
     
     if (length(quantiles) > 0 & !scales::zero_range(range(data$y))) {
       stopifnot(all(quantiles >= 0), all(quantiles <= 1))
       quantiles <- create_quantile_segment_frame(data, nudge, quantiles, grp=grp)
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

interleave <- function(a, b) {
  as.vector(rbind(a, b))
}
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
      x=interleave(violin.xs + nudge, violin.xmaxvs + nudge),
      y=rep(ys, each=2),
      group=rep(ys, each=2)
    )
  } else {
    data.frame(
      x=interleave(violin.xminvs - nudge, violin.xs - nudge),
      y=rep(ys, each=2),
      group=rep(ys, each=2)
    )
  }
}

geom_split_violin <- function(mapping=NULL, data=NULL, stat="ydensity", position="identity", ..., 
                              nudge=0, quantiles=NULL, trim=TRUE, scale="area", na.rm=FALSE, 
                              show.legend=NA, inherit.aes=TRUE) {
  layer(data=data, mapping=mapping, stat=stat, geom=GeomSplitViolin, position=position, 
        show.legend=show.legend, inherit.aes=inherit.aes, 
        params=list(trim=trim, scale=scale, nudge=nudge, quantiles=quantiles, na.rm=na.rm, ...))
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

# the number of concentrations of total ligand (free and bound)
num_large   <- 10

# Concentrations of total ligand evenly distributed on a logarithmic scale
# NOTE: if maximal end values of the sequence is -4, fitting with formula1 will give:
# Error in nls step factor 0.000488281 reduced below 'minFactor' of 0.000976562, 
logTl_num_large   <- seq(-9, -5, length=num_large)

# The unit of total ligand is nM
tl_num_large   <- 10^logTl_num_large   * (1e+9)

# theoretical concentrations of protein-ligand complex, its unit is nM
plTheoretical_num_large   <- formula1(tl_num_large,   para$p0, para$kd)


gen_plNoisy_array <- function(plTheoretical, paraTrue, repN, expN, sd_level, plbias=0) {
  # 2*(sd of measurements) should be 2.5%, 5%, and 10% of p0 for reasonable experiments
  # so sd is 1.25%, 2.5% and 5% of p0. Note: if sd is too high, NAN error can occur when fitting
  
  # sd_level: low 1.25%, medium 2.5%, high 5% of p0
  # the sd of normal distribution is range covering 95% of a normal 
  # distribution is approximately within [-2*sigma, 2*sigma].
  # 2*sigma = para*2.5%, para*5%, or para*10%, so sigma = para*2.5%/2, para*5%/2, or para*10%/2
  
  plTheoretical_len <- length(plTheoretical)
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
  
  error_array <- array(rnorm(repN*plTheoretical_len*expN, sd=error_sd), dim=c(repN, plTheoretical_len, expN)) + plbias
  rep(rep(plTheoretical, each=repN), times=expN) + error_array
}

fitfunc <- function(plNoisy, paraTrue, paraRatio, repN, l_conc, algorithm="LM", log10=FALSE) {
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
      nls.formula1 <- nlsLM(pl ~ formula1(tl, p0, kd), data=pldata, start=paraGuess)
      nls.formula2 <- nlsLM(pl ~ formula2(tl, p0, kd), data=pldata, start=paraGuess)
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
  rownamelist <- c("Kd (SD)")
  table_df <- data.frame(matrix(ncol=colN*2, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  for (i in 1:colN) {
    f1_kd <- arg2[[i]]["f1_kd",]
    f2_kd <- arg2[[i]]["f2_kd",]
    f1_kd_mean <- mean(f1_kd)
    f2_kd_mean <- mean(f2_kd)
    f1_kd_sd <- sd(f1_kd)
    f2_kd_sd <- sd(f2_kd)
    
    colName <- paste(stringList[i], "f1", sep="-")
    colName2 <- paste(stringList[i], "f2", sep="-")
    idx1 <- i*2 - 1
    idx2 <- i*2
    table_df[,idx1] <- c(sprintf("%.2f (%.2f)", f1_kd_mean, f1_kd_sd))
    table_df[,idx2] <- c(sprintf("%.2f (%.2f)", f2_kd_mean, f2_kd_sd))
    colnames(table_df)[idx1] <- colName
    colnames(table_df)[idx2] <- colName2
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}
write_table2 <- function(filename, listname, ...) {
  arg1 <- as.list(environment())
  filename_base <- arg1$filename
  csvfilename <- paste(filename_base, ".csv", sep="")
  stringList <- arg1$listname
  
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("Kd (SD)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  for (i in 1:colN) {
    f1_kd <- arg2[[i]]["f1_kd",]
    f1_kd_mean <- mean(f1_kd)
    f1_kd_sd <- sd(f1_kd)
    
    colName <- paste(stringList[i], "f1", sep="-")
    table_df[,i] <- c(sprintf("%.2f (%.2f)", f1_kd_mean, f1_kd_sd))
    colnames(table_df)[i] <- colName
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}


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

# Part 2. Binding Kinetics Experiment ------------------------------------

formula_binding_kinetics <- function(t, l0, pl_eq, koff, kon) pl_eq * (1 - exp(0 - (kon * l0 + koff) * t)) 
formula_binding_kinetics_group <- function(t, l0, pl_eq, g, koff, kon)  pl_eq[g] * (1 - exp(0 - (kon * l0 + koff) * t)) 

formula_binding_kinetics_obs <- function(t, pl_eq, kobs) pl_eq * (1 - exp(0 - kobs * t)) 
formula_binding_kinetics_obs_group <- function(t, pl_eq, g, kobs) pl_eq[g] * (1 - exp(0 - kobs[g] * t)) 

formula_binding_kinetics_disso <- function(t, pl0, koff) pl0 * exp(0 - koff * t) 
formula_binding_kinetics_disso_group <- function(t, pl0, g, koff) pl0[g] * exp(0 - koff * t) 

formula_second_order_kin <- function(t, l0, p0, koff, kon) {
  a <- 1
  b <- 0 - (p0 + l0 + koff / kon)
  c <- p0 * l0
  pl_1 <- (0 - b - sqrt(b * b - 4 * a * c)) / (2 * a)  # pl_1 is [PL]eq, pl_1 < pl_2
  pl_2 <- (0 - b + sqrt(b * b - 4 * a * c)) / (2 * a)
  E <- exp(kon * (pl_1 - pl_2) * t)
  pl_1 * pl_2 * (1 - E) / (pl_2 - pl_1 * E)
}

l0_calc <- function(p0, kd) {
  points <- seq(0.5, 0.9, 0.1)
  pl <- points*p0
  l0 <- (p0+kd-pl)*pl/(p0-pl)
  l0
}

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
  list(second_order_pl=second_order_pl, 
       second_order_time=second_order_time,  
       disso_pl=disso_pl, disso_t=disso_t)
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
  
  list(data_array_asso=data_array_asso, 
       data_array_disso=data_array_disso)
}

pl_disso_t <- function(t, p0, l0, koff, kon) {
  b   <- -(p0 + l0 + koff / kon)
  dsc <- sqrt(b * b - 4 * p0 * l0)
  ((-b - dsc) / 2) * exp(-koff * t)
}

fitfunc_kin <- function(plNoisy, paraTrue, paraRatio, type, l0_list, offres) {
# nlsLM cannot use group, so we use nls instead
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
  if (type=="sec_order_asso") {
    colnames(plNoisy) <- c("pl", "t", "l0", "g")
    pldata <- as.data.frame(plNoisy)
    koff_est <- offres[["koff_est"]]
    if (l0_num==1) {
      warning("The prediction of koff and kon needs more than one l0")
    } else {
      paraGuess <- list(p0=p0*paraRatio, kon=kon*paraRatio)
      nls.sec_order_asso <- nls(pl ~ formula_second_order_kin(t, l0, p0, koff=koff_est, kon), data=pldata, start=paraGuess)
    }
    coef_res <- coef(nls.sec_order_asso)
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

write_table_kin <- function(filename_base, stringList, ...) {
#func for fitResult_kin_asso_obs_conv1,2,3
  csvfilename <- paste(filename_base, ".csv", sep="")

  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("koff (SD)", "kon (SD)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    koff_est <- arg2[[i]]["koff_est",]; kon_est <- arg2[[i]]["kon_est",]
    koff_est_mean <- mean(koff_est); kon_est_mean <- mean(kon_est)
    koff_est_sd <- sd(koff_est); kon_est_sd <- sd(kon_est)
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", koff_est_mean, koff_est_sd), sprintf("%.2e (%.2e)", kon_est_mean, kon_est_sd))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_kin_disso <- function(filename_base, stringList, ...) {
#func for fitResult_kin_disso1,2,3
  csvfilename <- paste(filename_base, ".csv", sep="")

  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("koff (SD)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    koff_est <- arg2[[i]]["koff_est",]
    koff_est_mean <- mean(koff_est)
    koff_est_sd <- sd(koff_est); 
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", koff_est_mean, koff_est_sd))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_kin_asso <- function(filename_base, stringList, koffList, ...) {
#func for kon_linear_koff_corrected1,2,3 and kon_nonlinear_koff_corrected1,2,3
  csvfilename <- paste(filename_base, ".csv", sep="")
  
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("kon (SD)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    koff_est <- koffList[[i]]["koff_est",]; kon_est <- arg2[[i]]
    kon_est_mean <- mean(kon_est); kon_est_sd <- sd(kon_est)
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", kon_est_mean, kon_est_sd))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_kin_kd <- function(filename_base, stringList, type, ...) {
#func for fitResult_kin_cmp1,2,3 and kd_kin1,2,3
  csvfilename <- paste(filename_base, ".csv", sep="")
  
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("Kd (SD)")
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
    
    table_df[,i] <- c(sprintf("%.2f (%.2f)", kd_est_mean, kd_est_sd))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

# Part 3. Competitive Binding Experiment ----------------------------------

# logTi_num_large <- seq(-9, -4, length=num_large) 
# if logTi range [-9, -4] is used, ki estimated by direct fitting under the condition of para_comp4
# list(p0=5, l0=c(25,50,75,100),  kd=5, ki=100) will show strong numerical instability
logTi_num_large <- seq(-10, -5, length=num_large) 
ti_num_large <- 10^logTi_num_large * (1e+9)

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
  nls.comp <- nlsLM(pl ~ formula_competitive_pl(p0, l0=l0, i0, kd=kd, ki), data=pldata, start=paraGuess)
  coef_comp <- coef(nls.comp)
  c(p0=coef_comp[["p0"]], ki=coef_comp[["ki"]])
}

write_table_comp_appic50 <- function(filename_base, stringList, ...) {
  csvfilename <- paste(filename_base, ".csv", sep="")
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("a (SD)", "b (SD)", "c (SD)", "d (SD)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    # a 
    a_est <- arg2[[i]]["a",]
    a_est_median <- median(a_est)
    a_est_mean <- mean(a_est)
    a_est_sd <- sd(a_est); 
    
    # b 
    b_est <- abs(arg2[[i]]["b",]) # abs of b
    b_est_median <- median(b_est)
    b_est_mean <- mean(b_est)
    b_est_sd <- sd(b_est); 
    
    # c 
    c_est <- arg2[[i]]["c",]
    c_est_median <- median(c_est)
    c_est_mean <- mean(c_est)
    c_est_sd <- sd(c_est); 
    
    # d
    d_est <- arg2[[i]]["d",]
    d_est_median <- median(d_est)
    d_est_mean <- mean(d_est)
    d_est_sd <- sd(d_est); 
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", a_est_mean, a_est_sd),
                      sprintf("%.2e (%.2e)", b_est_mean, b_est_sd),
                      sprintf("%.2e (%.2e)", c_est_mean, c_est_sd),
                      sprintf("%.2e (%.2e)", d_est_mean, d_est_sd))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_comp_cheng <- function(filename_base, stringList, ...) {
  csvfilename <- paste(filename_base, ".csv", sep="")
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("Apparent Ki (SD)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    c_est <- arg2[[i]]
    c_est_median <- median(c_est)
    c_est_mean <- mean(c_est)
    c_est_sd <- sd(c_est); 
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", c_est_mean, c_est_sd))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

write_table_comp_ki <- function(filename_base, stringList, ...) {
  csvfilename <- paste(filename_base, ".csv", sep="")
  arg2 <- list(...)
  colN <- length(stringList)
  rownamelist <- c("Ki (SD)")
  table_df <- data.frame(matrix(ncol=colN, nrow=length(rownamelist)))
  rownames(table_df) <- rownamelist
  colnames(table_df) <- stringList
  for (i in 1:colN) {
    c_est <- arg2[[i]]
    c_est_median <- median(c_est)
    c_est_mean <- mean(c_est)
    c_est_sd <- sd(c_est); 
    
    table_df[,i] <- c(sprintf("%.2e (%.2e)", c_est_mean, c_est_sd))
  }
  write.csv(table_df, csvfilename, row.names=TRUE)
}

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
