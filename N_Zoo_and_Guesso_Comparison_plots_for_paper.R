list[Infections.Gueasso.13.fixed.seasonality, attack.rate.5.groups, initial.5.groups] =
  simulate.with.reporting.rate.adjusted.cases (sus.dist =Gueasso.fit1,
                                               time.length = 17,
                                               num.sims= 1000,
                                               waifw.init= t(waifw.5.groups), state.5.group.Gueasso)

epi.curve.dat = convert.data.for.ribbon.plot(Infections.Gueasso.13.fixed.seasonality, Infections.Gueasso.bounded.cases.and.obs ,  Gueasso.seasonal.est.13.bounded$Cases_by_age, Gueasso.bounded_cases_and_obs$Cases_by_age)

epi.curve.dat$run = factor(epi.curve.dat$run, levels=c("13", "17"), labels=c("week 13", "week 17"))
gg.curve.Gueasso <- ggplot(epi.curve.dat, aes(x=Week, y=Cases,id = iter, colour=as.factor(run), fill=as.factor(run))) + 
  stat_summary(fun.y=mean, geom="line") + 
  stat_summary(geom="ribbon", fun.ymin = quantile.05, fun.ymax= quantile.95, alpha=0.25, colour=NA) +
  geom_vline(xintercept = 17, color = 5, linetype = "dashed") +
  geom_vline(xintercept = 13, color = 2, linetype = "dashed")  + 
  theme(legend.position="bottom", legend.title=element_blank()) 
gg.curve.Gueasso



##' Comparison of predictions for week 13 and week 17 data where we simulate with a different simulation
list[Infections.Gueasso.13.fitted.seasonality, attack.rate.NZoo.13, initial.sus.NZoo.13] = simulate.with.R0.and.sus.dist.obs.rates(Gueasso.seasonal.est.13.bounded, Gueasso.seasonal.est.13.bounded$R02, 
                                                                                                                                          time.length = 14, num.sims,
                                                                                                                                          cases.by.age.group = cases.by.age.Gueasso.13,
                                                                                                                                          waifw.init = t(waifw.5.groups), state.Gueasso)

epi.curve.dat = convert.data.for.ribbon.plot(Infections.Gueasso.13.fixed.seasonality, 
                                             Infections.Gueasso.13.fitted.seasonality, 
                                             Gueasso.fit1$Cases_by_age, 
                                             Gueasso.fit2$Cases_by_age)

j = which(epi.curve.dat$run == "17")
epi.curve.dat$Week[j] = epi.curve.dat$Week[j] - 4
epi.curve.dat$run = factor(epi.curve.dat$run, levels=c("13", "17"), labels=c("Fixed Seasonality", "Fitted seasonality"))
gg.curve.Gueasso <- ggplot(epi.curve.dat, aes(x=Week, y=Cases,id = iter, colour=as.factor(run), fill=as.factor(run))) + 
  stat_summary(fun.y=mean, geom="line") + stat_summary(fun.y=mean, geom="point") +
  stat_summary(geom="ribbon", fun.ymin = quantile.05, fun.ymax= quantile.95, alpha=0.25, colour=NA) +
 # geom_vline(xintercept = 17, color = 5, linetype = "dashed") + 
  geom_vline(xintercept = 13, color = 1, linetype = "dashed")  +
  theme(legend.position="bottom", legend.title=element_blank()) 


######################################################################
######################################################################
##' Output figure comparing the predictions for Gueasso using 
pdf("Gueasso_comparison.pdf", height = 5, width = 6)
gg.curve.Gueasso
dev.off()




epi.curve.dat = convert.data.for.ribbon.plot.NZoo(Infections.NZoo.13,
                                                  NZoo.sims3, 
                                                  rowSums(NZoo.cases.report1), 
                                                  rowSums(NZoo.cases.report2))
epi.curve.dat$run = factor(epi.curve.dat$run, levels=c("13", "17"), labels=c("N'Zoo week 13", "N'Zoo week 17"))

gg.curve.NZoo <- ggplot(epi.curve.dat, aes(x=Week, y=Cases, id = iter, colour=as.factor(run), fill=as.factor(run))) + 
  stat_summary(fun.y=mean, geom="line") +
  stat_summary(fun.y=mean, geom="point") + 
  stat_summary(geom="ribbon", fun.ymin = quantile.05, fun.ymax= quantile.95, alpha=0.25, colour=NA) +
  theme(legend.position="bottom", legend.title=element_blank())   +
  geom_vline(xintercept = 13, color = 2,  linetype = "dashed")  +
  geom_vline(xintercept = 17, color = 5,  linetype = "dashed") 
gg.curve.NZoo

pdf("NZoo_comparison.pdf", height = 5, width = 7)
gg.curve.NZoo
dev.off()


