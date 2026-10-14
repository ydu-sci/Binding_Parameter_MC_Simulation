source("00_global.R")

# Part 3. Competitive Binding Experiment ----------------------------------
clust <- makeCluster(32)
clusterExport(clust, c("nlsLM", "formula_competitive_pl", "logistic_func", "drm", "LL.4"))

# * 3.1 Estimation of Apparent IC50, Apparent Ki, and Ki ------------------

dir.create("Part3.1-1"); setwd("Part3.1-1")

para_comp_l1 <- list(p0=5, l0=50, kd=10, ki=1)
para_comp_l2 <- list(p0=5, l0=50, kd=10, ki=10)
para_comp_l3 <- list(p0=5, l0=50, kd=10, ki=100)
# theo. data
plTheo_comp_l1 <- with(para_comp_l1, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp_l2 <- with(para_comp_l2, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
plTheo_comp_l3 <- with(para_comp_l3, formula_competitive_pl(p0, l0, ti_num_large, kd, ki))
# exp. data
plNoisy_comp_l1 <- gen_plNoisy_array(plTheo_comp_l1, paraTrue=para_comp_l1, repN=repN5, expN=expN, sd_level="medium")
plNoisy_comp_l2 <- gen_plNoisy_array(plTheo_comp_l2, paraTrue=para_comp_l2, repN=repN5, expN=expN, sd_level="medium")
plNoisy_comp_l3 <- gen_plNoisy_array(plTheo_comp_l3, paraTrue=para_comp_l3, repN=repN5, expN=expN, sd_level="medium")
# fitted ki
fitki_comp_l1 <- parApply(clust, plNoisy_comp_l1, 3, fitfunc_ki, paraTrue=para_comp_l1, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp_l2 <- parApply(clust, plNoisy_comp_l2, 3, fitfunc_ki, paraTrue=para_comp_l2, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
fitki_comp_l3 <- parApply(clust, plNoisy_comp_l3, 3, fitfunc_ki, paraTrue=para_comp_l3, paraRatio=ratioList[2], repN=repN5, i_conc=ti_num_large)
# fitted results
fitres_comp_l1 <- parApply(clust, plNoisy_comp_l1, 3, fitfunc_app_ic50, paraTrue=para_comp_l1, repN=repN5, i_conc=ti_num_large)
fitres_comp_l2 <- parApply(clust, plNoisy_comp_l2, 3, fitfunc_app_ic50, paraTrue=para_comp_l2, repN=repN5, i_conc=ti_num_large)
fitres_comp_l3 <- parApply(clust, plNoisy_comp_l3, 3, fitfunc_app_ic50, paraTrue=para_comp_l3, repN=repN5, i_conc=ti_num_large)
# fitted apparent IC50
fitappic50_comp_l1 <- fitres_comp_l1["c",]
fitappic50_comp_l2 <- fitres_comp_l2["c",]
fitappic50_comp_l3 <- fitres_comp_l3["c",]
# theoretical apparent IC50
theoappic50_comp_l1 <- with(para_comp_l1, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp_l2 <- with(para_comp_l2, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
theoappic50_comp_l3 <- with(para_comp_l3, inhibition_calc(p0=p0, l0=l0, kd=kd, ki=ki))$app_IC50
# ratio of fitted to theoretical apparent IC50
diffappic50_comp_l1 <- fitappic50_comp_l1 / theoappic50_comp_l1
diffappic50_comp_l2 <- fitappic50_comp_l2 / theoappic50_comp_l2
diffappic50_comp_l3 <- fitappic50_comp_l3 / theoappic50_comp_l3
#ic50 ratio
diffappic50_df1 <- data.frame(ratio=diffappic50_comp_l1, ki=para_comp_l1$ki)
diffappic50_df2 <- data.frame(ratio=diffappic50_comp_l2, ki=para_comp_l2$ki)
diffappic50_df3 <- data.frame(ratio=diffappic50_comp_l3, ki=para_comp_l3$ki)

diffappic50_df <- rbind(diffappic50_df1, diffappic50_df2, diffappic50_df3)
diffappic50_df$ki <- factor(diffappic50_df$ki,
                            levels = sort(unique(diffappic50_df$ki)))

diffappic50_violin <- ggplot(diffappic50_df, aes(x=ki, y=ratio)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.2, colour="grey40") +
  geom_violin(quantiles=c(0.25, 0.5, 0.75), linewidth=0.1,
              quantile.linetype = 1, scale="count", width=0.7, fill=c1_s1) +
  labs(x=expression(K[i]~"(nM)"),
       y="Estimated / Theoretical app. IC50") +
  legend_theme

ggsave("appic50.png", diffappic50_violin, dpi=800, width=6.5, height=6.5, units='cm')
ggsave("appic50.svg", diffappic50_violin, width=6.5, height=6.5, units='cm')

# ki calculated by cheng-prusoff equation
ki_cheng_comp_l1 <- calc_ki_app1(l0=para_comp_l1$l0, i0=fitappic50_comp_l1, kd=para_comp_l1$kd)
ki_cheng_comp_l2 <- calc_ki_app1(l0=para_comp_l2$l0, i0=fitappic50_comp_l2, kd=para_comp_l2$kd)
ki_cheng_comp_l3 <- calc_ki_app1(l0=para_comp_l3$l0, i0=fitappic50_comp_l3, kd=para_comp_l3$kd)

diffki_cheng_comp_l1 <- ki_cheng_comp_l1 / para_comp_l1$ki
diffki_cheng_comp_l2 <- ki_cheng_comp_l2 / para_comp_l2$ki
diffki_cheng_comp_l3 <- ki_cheng_comp_l3 / para_comp_l3$ki

diffki_cheng_df1 <- data.frame(ratio=diffki_cheng_comp_l1, ki=para_comp_l1$ki)
diffki_cheng_df2 <- data.frame(ratio=diffki_cheng_comp_l2, ki=para_comp_l2$ki)
diffki_cheng_df3 <- data.frame(ratio=diffki_cheng_comp_l3, ki=para_comp_l3$ki)

diffki_cheng_df <- rbind(diffki_cheng_df1, diffki_cheng_df2, diffki_cheng_df3)
diffki_cheng_df$ki <- factor(diffki_cheng_df$ki,
                             levels=sort(unique(diffki_cheng_df$ki)))

# ki calculated by wang's group equation
ki_wanggroup_comp_l1 <- with(para_comp_l1, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp_l1, kd=kd))
ki_wanggroup_comp_l2 <- with(para_comp_l2, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp_l2, kd=kd))
ki_wanggroup_comp_l3 <- with(para_comp_l3, calc_ki(p0=p0, l0=l0, i0=fitappic50_comp_l3, kd=kd))

diffki_wanggroup_comp_l1 <- ki_wanggroup_comp_l1 / para_comp_l1$ki
diffki_wanggroup_comp_l2 <- ki_wanggroup_comp_l2 / para_comp_l2$ki
diffki_wanggroup_comp_l3 <- ki_wanggroup_comp_l3 / para_comp_l3$ki

diffki_wanggroup_df1 <- data.frame(ratio=diffki_wanggroup_comp_l1, ki=para_comp_l1$ki)
diffki_wanggroup_df2 <- data.frame(ratio=diffki_wanggroup_comp_l2, ki=para_comp_l2$ki)
diffki_wanggroup_df3 <- data.frame(ratio=diffki_wanggroup_comp_l3, ki=para_comp_l3$ki)

diffki_wanggroup_df <- rbind(diffki_wanggroup_df1, diffki_wanggroup_df2, diffki_wanggroup_df3)
diffki_wanggroup_df$ki <- factor(diffki_wanggroup_df$ki,
                                 levels=sort(unique(diffki_wanggroup_df$ki)))

# ki directly fitted by wang equation
ki_direct_fit_comp_l1 <- fitki_comp_l1["ki",]
ki_direct_fit_comp_l2 <- fitki_comp_l2["ki",]
ki_direct_fit_comp_l3 <- fitki_comp_l3["ki",]

diffki_direct_fit_comp_l1 <- unname(fitki_comp_l1["ki", ]) / para_comp_l1$ki
diffki_direct_fit_comp_l2 <- unname(fitki_comp_l2["ki", ]) / para_comp_l2$ki
diffki_direct_fit_comp_l3 <- unname(fitki_comp_l3["ki", ]) / para_comp_l3$ki

diffki_direct_fit_df1 <- data.frame(ratio=diffki_direct_fit_comp_l1, ki=para_comp_l1$ki)
diffki_direct_fit_df2 <- data.frame(ratio=diffki_direct_fit_comp_l2, ki=para_comp_l2$ki)
diffki_direct_fit_df3 <- data.frame(ratio=diffki_direct_fit_comp_l3, ki=para_comp_l3$ki)

diffki_direct_fit_df <- rbind(diffki_direct_fit_df1, diffki_direct_fit_df2, diffki_direct_fit_df3)
diffki_direct_fit_df$ki <- factor(diffki_direct_fit_df$ki,
                                  levels=sort(unique(diffki_direct_fit_df$ki)))

# 3 in 1 
diffki_df <- rbind(
  transform(diffki_cheng_df,      method="Cheng-Prusoff"),
  transform(diffki_wanggroup_df,  method="Wang's group"),
  transform(diffki_direct_fit_df, method="Direct fit"))
diffki_df$method <- factor(diffki_df$method,
                           levels=c("Cheng-Prusoff", "Wang's group", "Direct fit"))

diffki_violin <- ggplot(diffki_df, aes(x=ki, y=ratio, fill=method)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.2, colour="grey40") +
  geom_violin(quantiles=c(0.25, 0.5, 0.75), linewidth=0.1,
              quantile.linetype = 1, scale="count", width=0.7, 
              position=position_dodge(width=0.8)) +
  scale_fill_manual(values=c("Cheng-Prusoff"=plt_color[1],
                             "Wang's group"=plt_color[2],
                             "Direct fit"=plt_color[3])) +
  labs(x=expression("True"~K[i]~"(nM)"),
       y="Estimated / True Ki", fill=NULL) +
  legend_theme

ggsave("diffki_violin.png", diffki_violin, dpi=800, width=13, height=6.5, units='cm')
ggsave("diffki_violin.svg", diffki_violin, width=13, height=6.5, units='cm')


kiStringList <- c("Ki = 1 nM", "Ki = 10 nM", "Ki = 100 nM")

dfcomp_comp_l1 <- data.frame(x=rep(ti_num_large, each=repN5),
                             y=as.vector(plNoisy_comp_l1[,,1]), group=kiStringList[1])
dfcomp_comp_l2 <- data.frame(x=rep(ti_num_large, each=repN5),
                             y=as.vector(plNoisy_comp_l2[,,1]), group=kiStringList[2])
dfcomp_comp_l3 <- data.frame(x=rep(ti_num_large, each=repN5),
                             y=as.vector(plNoisy_comp_l3[,,1]), group=kiStringList[3])

dfcomp_comp <- rbind(dfcomp_comp_l1, dfcomp_comp_l2, dfcomp_comp_l3)
dfcomp_comp$group <- factor(dfcomp_comp$group, levels=kiStringList)

dfcomp_comp_point <- ggplot(dfcomp_comp, aes(x=x, y=y, colour=group)) +
  stat_summary(fun=mean, geom="line", linewidth=0.3) +
  geom_point(size=0.3, alpha=0.7) +
  scale_x_log10(breaks=trans_breaks("log10", function(x) 10^x),
                labels=trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values=c(plt_color[1], plt_color[2], plt_color[3])) +
  labs(x=expression("[I]"[0]~"(nM)"), y=expression("[PL]"[obs]~"(nM)"), colour=NULL) +
  legend_theme

ggsave("dfcomp_comp_point.png", dfcomp_comp_point, dpi=800, width=6.5, height=6.5, units='cm')
ggsave("dfcomp_comp_point.svg", dfcomp_comp_point, width=6.5, height=6.5, units='cm')

setwd("..")

dir.create("Part3.1-2"); setwd("Part3.1-2")

compOutString <- c("Ki = 1 nM", "Ki = 10 nM", "Ki = 100 nM")

write_table_comp_appic50(filename_base="appic50", stringList=compOutString, 
                         fitres_comp_l1, fitres_comp_l2, fitres_comp_l3)

write_table_comp_cheng(filename_base="cheng", stringList=compOutString, 
                       ki_cheng_comp_l1, ki_cheng_comp_l2, ki_cheng_comp_l3)

write_table_comp_ki(filename_base="wanggroup", stringList=compOutString, 
                    ki_wanggroup_comp_l1, ki_wanggroup_comp_l2, ki_wanggroup_comp_l3)

write_table_comp_ki(filename_base="directfit", stringList=compOutString, 
                    ki_direct_fit_comp_l1, ki_direct_fit_comp_l2, ki_direct_fit_comp_l3)
setwd("..")

stopCluster(clust)
