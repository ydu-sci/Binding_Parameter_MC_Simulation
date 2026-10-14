source("00_global.R")

# Part 1. Binding Experiment ----------------------------------------------
clust <- makeCluster(32)
clusterExport(clust, c("nlsLM", "formula1", "formula2", "formula1_log10", "formula2_log10"))

# * 1.1 Comparison of Algorithm and Model  ------------------

# simulation of experiments
plNoisy_array <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="medium")
fitResultRatio05  <- parApply(clust, plNoisy_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

#log10
fitResultRatio05_log10  <- parApply(clust, plNoisy_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large, log10=TRUE)
fitResultRatio05_GN  <- parApply(clust, plNoisy_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large, algorithm="GN")

# *- Note-01 ---------------------------------------------------------------
# fitting using different algorithms generates very similar estimates 
sub2 <- fitResultRatio05["f1_kd",] - fitResultRatio05_GN["f1_kd",]
fivenum(sub2)
sub2 <- fitResultRatio05["f2_kd",] - fitResultRatio05_GN["f2_kd",]
fivenum(sub2)

sub2_10 <- fitResultRatio05["f1_kd",] - 10^fitResultRatio05_log10["f1_kd",]
fivenum(sub2_10)
sub2_log10 <- log10(fitResultRatio05["f1_kd",]) - fitResultRatio05_log10["f1_kd",]
fivenum(sub2_log10)
sub2_10 <- fitResultRatio05["f2_kd",] - 10^fitResultRatio05_log10["f2_kd",]
fivenum(sub2_10)
sub2_log10 <- log10(fitResultRatio05["f2_kd",]) - fitResultRatio05_log10["f2_kd",]
fivenum(sub2_log10)

# * 1.2 Errors and Replicates in Experimental Measurements -------------

# *- 1.2.1 Standard Deviation of Measurement Errors --------------------

dir.create("Part1.2.1"); setwd("Part1.2.1")
sdStringList <- c("Low SD", "Medium SD", "High SD")

# simulation of experiments
plNoisy_errorSDLow_array      <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="low")
plNoisy_errorSDMedium_array   <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="medium")
plNoisy_errorSDHigh_array     <- gen_plNoisy_array(plTheoretical_num_large, paraTrue=para, repN=repN5, expN=expN, sd_level="high")

# fitting process
fitResultErrorSD_low     <- parApply(clust, plNoisy_errorSDLow_array,    3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResultErrorSD_medium  <- parApply(clust, plNoisy_errorSDMedium_array, 3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResultErrorSD_high    <- parApply(clust, plNoisy_errorSDHigh_array,   3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

# violin plots of different sd comparing formula 1 and 2
fig01_kd_df <- data.frame(
  x=rep(rep(sdStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResultErrorSD_low["f1_kd",], fitResultErrorSD_medium["f1_kd",], fitResultErrorSD_high["f1_kd",],
      fitResultErrorSD_low["f2_kd",], fitResultErrorSD_medium["f2_kd",], fitResultErrorSD_high["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(sdStringList))
)
fig01_kd_violin <- ggplot(fig01_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=sdStringList, guide=guide_axis(n.dodge=2)) +
  theme(axis.title.x=element_blank())+
  legend_theme

ggsave("kd.png", fig01_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig01_kd_violin, width=7, height=7, units='cm')

setwd("..")

# *- 1.2.2 Accuracy of Measurements  ---------------------------------------

dir.create("Part1.2.2"); setwd("Part1.2.2")
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
fig02_kd_df <- data.frame(
  x=rep(rep(biasStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_posBias["f1_kd",], fitResult_noBias["f1_kd",], fitResult_negBias["f1_kd",],
      fitResult_posBias["f2_kd",], fitResult_noBias["f2_kd",], fitResult_negBias["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(biasStringList))
)
fig02_kd_violin <- ggplot(fig02_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=biasStringList) +
  theme(axis.title.x=element_blank())+
  legend_theme

ggsave("kd.png", fig02_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig02_kd_violin, width=7, height=7, units='cm')

setwd("..")

# *- 1.2.3 Number of Replicates ----------------------------------------------

dir.create("Part1.2.3"); setwd("Part1.2.3")
repNList <- c(1, 2, 4, 8, 12)
concNlist <- c(48, 24, 12, 6, 4)
repNStringList <- c("1Rep", "2Rep", "4Rep", "8Rep", "12Rep")

logTl_num_1Rep <- seq(-9, -5, length=concNlist[1])
logTl_num_2Rep <- seq(-9, -5, length=concNlist[2])
logTl_num_4Rep <- seq(-9, -5, length=concNlist[3])
logTl_num_8Rep <- seq(-9, -5, length=concNlist[4])
logTl_num_12Rep <- seq(-9, -5, length=concNlist[5])

# The unit of total ligand is nM
tl_num_1Rep <- 10^logTl_num_1Rep * (1e+9)
tl_num_2Rep <- 10^logTl_num_2Rep * (1e+9)
tl_num_4Rep <- 10^logTl_num_4Rep * (1e+9)
tl_num_8Rep <- 10^logTl_num_8Rep * (1e+9)
tl_num_12Rep <- 10^logTl_num_12Rep * (1e+9)

# theoretical concentrations of protein-ligand complex, its unit is nM
plTheoretical_num_1Rep <- formula1(tl_num_1Rep, para$p0, para$kd)
plTheoretical_num_2Rep <- formula1(tl_num_2Rep, para$p0, para$kd)
plTheoretical_num_4Rep <- formula1(tl_num_4Rep, para$p0, para$kd)
plTheoretical_num_8Rep <- formula1(tl_num_8Rep, para$p0, para$kd)
plTheoretical_num_12Rep <- formula1(tl_num_12Rep, para$p0, para$kd)

# simulation of experiments
plNoisy_1Rep_array  <- gen_plNoisy_array(plTheoretical_num_1Rep, paraTrue=para, repN=repNList[1], expN=expN, sd_level="medium")
plNoisy_2Rep_array  <- gen_plNoisy_array(plTheoretical_num_2Rep, paraTrue=para, repN=repNList[2], expN=expN, sd_level="medium")
plNoisy_4Rep_array  <- gen_plNoisy_array(plTheoretical_num_4Rep, paraTrue=para, repN=repNList[3], expN=expN, sd_level="medium")
plNoisy_8Rep_array  <- gen_plNoisy_array(plTheoretical_num_8Rep, paraTrue=para, repN=repNList[4], expN=expN, sd_level="medium")
plNoisy_12Rep_array  <- gen_plNoisy_array(plTheoretical_num_12Rep, paraTrue=para, repN=repNList[5], expN=expN, sd_level="medium")

# parallel
fitResult_1Rep  <- parApply(clust, plNoisy_1Rep_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[1], l_conc=tl_num_1Rep)
fitResult_2Rep  <- parApply(clust, plNoisy_2Rep_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[2], l_conc=tl_num_2Rep)
fitResult_4Rep  <- parApply(clust, plNoisy_4Rep_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[3], l_conc=tl_num_4Rep)
fitResult_8Rep  <- parApply(clust, plNoisy_8Rep_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[4], l_conc=tl_num_8Rep)
fitResult_12Rep  <- parApply(clust, plNoisy_12Rep_array,  3, fitfunc, paraTrue=para, paraRatio=ratioList[2], repN=repNList[5], l_conc=tl_num_12Rep)

# violin plots of different numbers of replicates
fig03_kd_df <- data.frame(
  x=rep(rep(repNStringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_1Rep["f1_kd",], fitResult_2Rep["f1_kd",], fitResult_4Rep["f1_kd",], fitResult_8Rep["f1_kd",], fitResult_12Rep["f1_kd",],
      fitResult_1Rep["f2_kd",], fitResult_2Rep["f2_kd",], fitResult_4Rep["f2_kd",], fitResult_8Rep["f2_kd",], fitResult_12Rep["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(repNStringList))
)
fig03_kd_violin <- ggplot(fig03_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="area") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=repNStringList) +
  theme(axis.title.x=element_blank())+
  legend_theme

ggsave("kd.png", fig03_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig03_kd_violin, width=7, height=7, units='cm')

write_table(filename="Output", listname=repNStringList, fitResult_1Rep, fitResult_2Rep, fitResult_4Rep, fitResult_8Rep, fitResult_12Rep)
setwd("..")

# * 1.3 Concentrations of Total Protein -----------------------------------------------
dir.create("Part1.3"); setwd("Part1.3")

para1 <- list(p0=10, kd=100)
para2 <- list(p0=100, kd=100)
para3 <- list(p0=1000, kd=100)

plTheoretical_para1 <- formula1(tl_num_large, para1$p0, para1$kd)
plTheoretical_para2 <- formula1(tl_num_large, para2$p0, para2$kd)
plTheoretical_para3 <- formula1(tl_num_large, para3$p0, para3$kd)

p0StringList <- c("[P]0 = 0.1Kd", "[P]0 = 1.0Kd", "[P]0 = 10.0Kd")

# simulation of experiments
plNoisy_para1 <- gen_plNoisy_array(plTheoretical_para1, paraTrue=para1, repN=repN5, expN=expN, sd_level="medium")
plNoisy_para2 <- gen_plNoisy_array(plTheoretical_para2, paraTrue=para2, repN=repN5, expN=expN, sd_level="medium")
plNoisy_para3 <- gen_plNoisy_array(plTheoretical_para3, paraTrue=para3, repN=repN5, expN=expN, sd_level="medium")

# fitting process
fitResult_para1 <- parApply(clust, plNoisy_para1, 3, fitfunc, paraTrue=para1, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_para2 <- parApply(clust, plNoisy_para2, 3, fitfunc, paraTrue=para2, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_para3 <- parApply(clust, plNoisy_para3, 3, fitfunc, paraTrue=para3, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

# violin plots of different p0 concentrations comparing formula 1 and 2
fig05_kd_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(formulaStringList)),
  y=c(fitResult_para1["f1_kd",], fitResult_para2["f1_kd",], fitResult_para3["f1_kd",],
      fitResult_para1["f2_kd",], fitResult_para2["f2_kd",], fitResult_para3["f2_kd",]),
  f=rep(formulaStringList, each=expN*length(p0StringList))
)
fig05_kd_violin <- ggplot(fig05_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Equation 1", "Equation 2"),
                    values=c("Equation 1"=c1_s1, "Equation 2"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  scale_x_discrete(limits=p0StringList, guide=guide_axis(n.dodge=2)) +
  theme(axis.title.x=element_blank())+
  legend_theme+
  theme(legend.position.inside=c(0, 1), legend.justification=c(0, 1))

ggsave("kd.png", fig05_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", fig05_kd_violin, width=7, height=7, units='cm')

write_table(filename="Output", listname=p0StringList, fitResult_para1, fitResult_para2, fitResult_para3)
setwd("..")

# * 1.4 Concentration Errors in [P]0 and [L]0 -----------------------------

dir.create("Part1.4"); setwd("Part1.4")

l0NoisyStringList <- c("Noisy [P]0 and [L]0", "Noise-free [P]0 and [L]0")

# simulation of experiments
plNoisy_tlNoisy_array1 <- gen_plNoisy_array2(tl_num_large, paraTrue=para1, repN=repN5, expN=expN, sd_level="medium", tp_cv="medium", tl_cv="medium")
plNoisy_tlNoisy_array2 <- gen_plNoisy_array2(tl_num_large, paraTrue=para2, repN=repN5, expN=expN, sd_level="medium", tp_cv="medium", tl_cv="medium")
plNoisy_tlNoisy_array3 <- gen_plNoisy_array2(tl_num_large, paraTrue=para3, repN=repN5, expN=expN, sd_level="medium", tp_cv="medium", tl_cv="medium")

# fitting process
fitResult_tlNoisy1 <- parApply(clust, plNoisy_tlNoisy_array1, 3, fitfunc, paraTrue=para1, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_tlNoisy2 <- parApply(clust, plNoisy_tlNoisy_array2, 3, fitfunc, paraTrue=para2, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)
fitResult_tlNoisy3 <- parApply(clust, plNoisy_tlNoisy_array3, 3, fitfunc, paraTrue=para3, paraRatio=ratioList[2], repN=repN5, l_conc=tl_num_large)

# violin plots
res2_kd_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(l0NoisyStringList)),
  y=c(fitResult_tlNoisy1["f1_kd",], fitResult_tlNoisy2["f1_kd",], fitResult_tlNoisy3["f1_kd",], 
      fitResult_para1["f1_kd",], fitResult_para2["f1_kd",], fitResult_para3["f1_kd",]),
  f=rep(l0NoisyStringList, each=expN*length(p0StringList))
)
res2_kd_violin <- ggplot(res2_kd_df, aes(x, y, fill=f)) + 
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Noisy [P]0 and [L]0", "Noise-free [P]0 and [L]0"),
                    values=c("Noisy [P]0 and [L]0"=c1_s1, "Noise-free [P]0 and [L]0"=c2_s1)) +
  labs(y="Estimated Kd (nM)") +
  theme(axis.title.x=element_blank())+
  legend_theme

ggsave("kd.png", res2_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", res2_kd_violin, width=7, height=7, units='cm')

write_table2(filename="Output", listname=p0StringList, fitResult_tlNoisy1, fitResult_tlNoisy2, fitResult_tlNoisy3)

setwd("..")

stopCluster(clust)
