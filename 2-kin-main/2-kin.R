source("00_global.R")

# Part 2. Binding Kinetics Experiment ------------------------------------
clust <- makeCluster(32)

clusterExport(clust, c("nlsLM","formula_binding_kinetics", "formula_binding_kinetics_group", 
                       "formula_binding_kinetics_obs", "formula_binding_kinetics_obs_group", 
                       "formula_binding_kinetics_disso", "formula_binding_kinetics_disso_group",
                       "formula_second_order_kin", "formula1", "formula2"))

# * 2.1 Estimation of koff and kon ---------
dir.create("Part2.1.1"); setwd("Part2.1.1")
time_points <- 10

# the unit of p0, koff, and kon is nM, s-1, and nM-1 s-1
kinetics_para1 <- list(p0=100, koff=1e-2, kon=1e-4)
kinetics_para2 <- list(p0=50, koff=1e-2, kon=1e-4)
kinetics_para3 <- list(p0=10, koff=1e-2, kon=1e-4)

l0_list1 <- with(kinetics_para1, l0_calc(p0=p0, kd=koff/kon))
l0_list2 <- with(kinetics_para2, l0_calc(p0=p0, kd=koff/kon))
l0_list3 <- with(kinetics_para3, l0_calc(p0=p0, kd=koff/kon))

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

kon_sec_order_nonlinear_koff_corrected1 <- sapply(1:expN, function(x) fitfunc_kin(plNoisy_array_kinetics1$data_array_asso[,,x], paraTrue=kinetics_para1, paraRatio=ratioList[2], type="sec_order_asso", l0_list1, fitResult_kin_disso1[,x]))
kon_sec_order_nonlinear_koff_corrected2 <- sapply(1:expN, function(x) fitfunc_kin(plNoisy_array_kinetics2$data_array_asso[,,x], paraTrue=kinetics_para2, paraRatio=ratioList[2], type="sec_order_asso", l0_list2, fitResult_kin_disso2[,x]))
kon_sec_order_nonlinear_koff_corrected3 <- sapply(1:expN, function(x) fitfunc_kin(plNoisy_array_kinetics3$data_array_asso[,,x], paraTrue=kinetics_para3, paraRatio=ratioList[2], type="sec_order_asso", l0_list3, fitResult_kin_disso3[,x]))

kd_sec_order_kin1 <- fitResult_kin_disso1["koff_est",]/kon_sec_order_nonlinear_koff_corrected1
kd_sec_order_kin2 <- fitResult_kin_disso2["koff_est",]/kon_sec_order_nonlinear_koff_corrected2
kd_sec_order_kin3 <- fitResult_kin_disso3["koff_est",]/kon_sec_order_nonlinear_koff_corrected3

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
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
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
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Asso_obs", "Nonlinear (koff corrected)"),
                    values=c("Asso_obs"=c1_s1, "Nonlinear (koff corrected)"=c2_s1)) +
  labs(x="[P]0 (nM)", y="Estimated kon (nM-1 s-1)") +
  scale_x_discrete(limits=p0StringList) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
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
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("Linear (koff corrected)", "Nonlinear (koff corrected)"),
                    values=c("Linear (koff corrected)"=c1_s1, "Nonlinear (koff corrected)"=c2_s1)) +
  labs(x="[P]0 (nM)", y="Estimated kon (nM-1 s-1)") +
  scale_x_discrete(limits=p0StringList) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  legend_theme

# violin plot of kon
kinStringList <- c("2nd order (koff corrected)", "Nonlinear (koff corrected)")
kin_asso_kon_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(kinStringList)),
  y=c(kon_sec_order_nonlinear_koff_corrected1, kon_sec_order_nonlinear_koff_corrected2, kon_sec_order_nonlinear_koff_corrected3,
      kon_nonlinear_koff_corrected1, kon_nonlinear_koff_corrected2, kon_nonlinear_koff_corrected3),
  f=rep(kinStringList, each=expN*length(p0StringList))
)
asso_kon_violin3 <- ggplot(kin_asso_kon_df, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
  scale_fill_manual(breaks=c("2nd order (koff corrected)", "Nonlinear (koff corrected)"),
                    values=c("2nd order (koff corrected)"=c1_s1, "Nonlinear (koff corrected)"=c2_s1)) +
  labs(x="[P]0 (nM)", y="Estimated kon (nM-1 s-1)") +
  scale_x_discrete(limits=p0StringList) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  legend_theme

# violin plot of kd
kinStringList <- c("Equilibrium", "Kinetics")
kin_kd_df <- data.frame(
  x=rep(rep(p0StringList, each=expN), times=length(kinStringList)),
  y=c(fitResult_kin_cmp1["f1_kd",], fitResult_kin_cmp2["f1_kd",], fitResult_kin_cmp3["f1_kd",], 
      kd_sec_order_kin1, kd_sec_order_kin2, kd_sec_order_kin3),
  f=rep(kinStringList, each=expN*length(p0StringList))
)
kin_kd_violin <- ggplot(kin_kd_df, aes(x, y, fill=f)) +
  geom_split_violin(nudge=0.01, quantiles=c(0.25, 0.5, 0.75), linewidth=0.1, scale="count") +
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
ggsave("kon3.png", asso_kon_violin3, dpi=800, width=7, height=7, units='cm')
ggsave("kon3.svg", asso_kon_violin3, width=7, height=7, units='cm')
ggsave("kd.png", kin_kd_violin, dpi=800, width=7, height=7, units='cm')
ggsave("kd.svg", kin_kd_violin, width=7, height=7, units='cm')

setwd("..")

dir.create("Part2.1.2"); setwd("Part2.1.2")


kinOutString <- c("[P]0=100nM", "[P]0=50nM", "[P]0=10nM")
kinkoffList <- list(fitResult_kin_disso1, fitResult_kin_disso2, fitResult_kin_disso3)
write_table_kin(filename_base="Output.1.asso_obs", stringList=kinOutString, fitResult_kin_asso_obs_conv1, fitResult_kin_asso_obs_conv2, fitResult_kin_asso_obs_conv3)
write_table_kin_disso(filename_base="Output.2.koff", stringList=kinOutString, fitResult_kin_disso1, fitResult_kin_disso2, fitResult_kin_disso3)
write_table_kin_asso(filename_base="Output.3.kon_linear", stringList=kinOutString, kinkoffList, kon_linear_koff_corrected1, kon_linear_koff_corrected2, kon_linear_koff_corrected3)
write_table_kin_asso(filename_base="Output.4.kon_nonlinear", stringList=kinOutString, kinkoffList, kon_nonlinear_koff_corrected1, kon_nonlinear_koff_corrected2, kon_nonlinear_koff_corrected3)
write_table_kin_asso(filename_base="Output.5.kon_nonlinear_2nd_order", stringList=kinOutString, kinkoffList, kon_sec_order_nonlinear_koff_corrected1, kon_sec_order_nonlinear_koff_corrected2, kon_sec_order_nonlinear_koff_corrected3)
write_table_kin_kd(filename_base="Output.6.kd-equi", stringList=kinOutString, type="equi", fitResult_kin_cmp1, fitResult_kin_cmp2, fitResult_kin_cmp3)
write_table_kin_kd(filename_base="Output.7.kd-kin", stringList=kinOutString, type="kin", kd_sec_order_kin1, kd_sec_order_kin2, kd_sec_order_kin3)
setwd("..")

stopCluster(clust)
