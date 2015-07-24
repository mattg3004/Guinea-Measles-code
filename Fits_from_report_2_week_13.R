estimate.susceptibility.and.observation.5.group.seasonality.report2.week.13 <- function(waifw, R0, iters, chains, 
                                                                                        R0_sd, alpha_mean, alpha_s, sus_prior,
                                                                                        S.prefect, population){
  
  list[ N.5.group, M.5.group, reported.cases, state, weeks] = output.data.for.given.sub.prefecture(S.prefect, population)
  ### Using the later report of cases, we only want to include cases up to week 13 in order to compare with what we were given in report 1 up to week 13.
  j = which(weeks == 13)
  
  ### we want to include seasonality. Estimated seasonality has 0.3 * cos shape, with cos = 1 approximately at weeks 13 and 39. 
  ### Need to work out where the first week of observed cases is on this scale in order to scale correctly for seasonality.
  first.week.adjust = weeks[1] - 13
  weeks = weeks[1:j]
  reported.cases = reported.cases[1:j, ]
  N = length(weeks)
  seasonal.multiplier = matrix(0, N, 1)
  for(i in 1 : N){
    # waifw = output.waifw(waifw.init, R_0 * (1 + cos((3/12 - (7 - i)/48)  * 2 *  pi) * 0.3), state)
    seasonal.multiplier[i] = 1 + cos(pi / 2 +  ((first.week.adjust + i - 1 )/52)* 2 *  pi)
  }
  #waifw.input = output.waifw(waifw = t(waifw), R_0 = R0, state = state)
  waifw = output.waifw(waifw.init, R0, state) / sum(state)
  
  ### Guinea.data holds the data to input into the stan code
  Guinea.data <- list("Obs_cases_age" = reported.cases, "N" = N, "M" = 5,  
                      "state" = state, "R0" = R0, "R0_sd" = R0_sd, "waifw" = waifw, 
                      "seasonal_multiplier" = seasonal.multiplier, "alpha_mean" = alpha_mean, 
                      "alpha_s" = alpha_s, "sus_prior" = sus_prior)
  
  ### Here the stan model is run to estimate the susceptibility in the 5 age groups. 
  fit <- stan(model_code = Guinea_model_tsir_fit_observation_rate_seasonal_forcing , data = Guinea.data, iter = iters , chains = chains) 
  
  ### print the output of the stan model for informational purposes
  print(fit)
  
  ### extract the output of the stan model, which is also the output of this function.
  sus.est = extract(fit)
  
  return(list(fit, sus.est))
}


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

###############         DO MODEL FITS FOR REPORT 2 USING DATA UP TO WEEK 13         ###############

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
sus_prior.seasonal.report2.13 = rbind(c(mean(N.Zoo.seasonality.13.report2$s[, 1]), sd(N.Zoo.seasonality.13.report2$s[, 1])), 
                                      c(mean(N.Zoo.seasonality.13.report2$s[, 2]), sd(N.Zoo.seasonality.13.report2$s[, 2])),
                                      c(mean(N.Zoo.seasonality.13.report2$s[, 3]), sd(N.Zoo.seasonality.13.report2$s[, 3])),
                                      c(mean(N.Zoo.seasonality.13.report2$s[, 4]), sd(N.Zoo.seasonality.13.report2$s[, 4])),
                                      c(mean(N.Zoo.seasonality.13.report2$s[, 5]), sd(N.Zoo.seasonality.13.report2$s[, 5])))

list[Foum.fit.report2.13, Foum.est.report2.13] = estimate.susceptibility.and.observation.5.group.seasonality.report2.week.13(waifw = t(waifw.5.groups), R0 = mean(N.Zoo.seasonality$R02), 
                                                                                                                                     R0_sd = sd(N.Zoo.seasonality$R02), 
                                                                                                                                     iters = 2000, chains = 4, 
                                                                                                                                     alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                                                     alpha_s = sqrt(3 * (sd(N.Zoo.seasonality$alpha))^2 /(pi^2)), 
                                                                                                                                     sus_prior.seasonal.report2.13, S.prefect = "FOUMBADOU",
                                                                                                                                     population = 19438)



#output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.Foum.13.report2, time.length = 14, r, g, b, "Foumbadou", 
#                                                                        Foum.est.report2.week13, obs.cases.by.week = rowSums(cases.by.age.Foum[1:which(weeks == 13),]), 
#                                                                        num.sims, first.week = 14)

list[Gueasso.fit.report2.13, Gueasso.est.report2.13] = estimate.susceptibility.and.observation.5.group.seasonality.report2.week.13(waifw =  t(waifw.5.groups), 
                                                                                                                        R0 = mean(N.Zoo.seasonality.13$R02), 
                                                                                                                        R0_sd = sd(N.Zoo.seasonality.13$R02), 
                                                                                                                        iters = 5000, chains = 4, 
                                                                                                                        alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                                        alpha_s = sqrt(3 * (sd(N.Zoo.seasonality.13$alpha))^2 /(pi^2)), 
                                                                                                                        sus_prior.seasonal.report2.13, S.prefect = "GUEASSO", population = 21116)




list[Gama.fit.report2.13, Gama.est.report2.13] = estimate.susceptibility.and.observation.5.group.seasonality.report2.week.13(waifw =  t(waifw.5.groups), 
                                                                                                                  R0 = mean(N.Zoo.seasonality.13$R02), 
                                                                                                                  R0_sd = sd(N.Zoo.seasonality.13$R02), 
                                                                                                                  iters = 5000, chains = 4, 
                                                                                                                  alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                                  alpha_s = sqrt(3 * (sd(N.Zoo.seasonality.13$alpha))^2 /(pi^2)), 
                                                                                                                  sus_prior.seasonal.report2.13, S.prefect = "GAMA BEREMA", population = 20465)



list[Laine.fit.report2.13, Laine.est.report213] = estimate.susceptibility.and.observation.5.group.seasonality.report2.week.13(waifw =  t(waifw.5.groups), 
                                                                                                                              R0 = mean(N.Zoo.seasonality.13$R02), 
                                                                                                                              R0_sd = sd(N.Zoo.seasonality.13$R02), 
                                                                                                                              iters = 5000, chains = 4, 
                                                                                                                              alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                                              alpha_s = sqrt(3 * (sd(N.Zoo.seasonality.13$alpha))^2 /(pi^2)),  
                                                                                                                              sus_prior.seasonal.report2.13, S.prefect = "LAINE", population = 16591)



list[Kokota.fit.report2.13, Kokota.est.report2.13] = estimate.susceptibility.and.observation.5.group.seasonality.report2.week.13(waifw =  t(waifw.5.groups), 
                                                                                                                                 R0 = mean(N.Zoo.seasonality.13$R02), 
                                                                                                                                 R0_sd = sd(N.Zoo.seasonality.13$R02), 
                                                                                                                                 iters = 5000, chains = 4, 
                                                                                                                                 alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                                                 alpha_s = sqrt(3 * (sd(N.Zoo.seasonality.13$alpha))^2 /(pi^2)),  
                                                                                                                                 sus_prior.seasonal.report2.13, S.prefect = "KOKOTA", population = 14732)



list[Cu.Lola.fit.report2.13, Cu.Lola.est.report2.13] = estimate.susceptibility.and.observation.5.group.seasonality.report2.week.13(waifw =  t(waifw.5.groups), 
                                                                                                                                   R0 = mean(N.Zoo.seasonality.13$R02), 
                                                                                                                                   R0_sd = sd(N.Zoo.seasonality.13$R02), 
                                                                                                                                   iters = 5000, chains = 4, 
                                                                                                                                   alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                                                   alpha_s = sqrt(3 * (sd(N.Zoo.seasonality.13$alpha))^2 /(pi^2)), 
                                                                                                                                   sus_prior.seasonal.report2.13, S.prefect = "CU LOLA", population = 49993)


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

###############         DO SIMULATIONS FOR FITS FROM REPORT 2 USING DATA UP TO WEEK 13       ###############

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

list[ N.5.group, M.5.group, cases.by.age.Foum, state.Foum, weeks] = output.data.for.given.sub.prefecture(S.prefect = "FOUMBADOU", total.pop = 19438)
list[Infections.Foum.13.report2, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Foum.est.report2.13, Foum.est.report2.13$R02, 
                                                                                 time.length = 14, num.sims,
                                                                                 cases.by.age.group = cases.by.age.Foum[1:which(weeks == 13),],
                                                                                 waifw.init = t(waifw.5.groups), state.Foum)



list[ N.5.group, M.5.group, cases.by.age.Gueasso, state.Gueasso, weeks] = output.data.for.given.sub.prefecture(S.prefect = "GUEASSO",total.pop =  21116)
list[Infections.Gueasso.13.report2, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Gueasso.est.report2.13, Gueasso.est.report2.13$R02, 
                                                                                    time.length = 14, num.sims,
                                                                                    cases.by.age.group = cases.by.age.Gueasso[1:which(weeks == 13),],
                                                                                    waifw.init = t(waifw.5.groups), state.Gueasso)




list[ N.5.group, M.5.group, cases.by.age.Gama, state.Gama, weeks] = output.data.for.given.sub.prefecture(S.prefect = "GAMA BEREMA",total.pop =  20465)
list[Infections.Gama.13.report2, attack.rate.NZoo.13, initial.sus.NZoo.13] = simulate.with.R0.and.sus.dist.obs.rates(Gama.est.report2.13, Gama.est.report2.13$R02, 
                                                                                                                     time.length = 14, num.sims,
                                                                                                                     cases.by.age.group = cases.by.age.Gama[1:which(weeks == 13),],
                                                                                                                     waifw.init = t(waifw.5.groups), state.Gama)



list[ N.5.group, M.5.group, cases.by.age.Laine, state.Laine, weeks] = output.data.for.given.sub.prefecture(S.prefect = "LAINE",total.pop =  16591)
list[Infections.Laine.13.report2, a, i] = simulate.with.R0.and.sus.dist.obs.rates( Laine.est.report213,  Laine.est.report213$R02, 
                                                                                   time.length = 14, num.sims,
                                                                                   cases.by.age.group = cases.by.age.Laine[1:which(weeks == 13),],
                                                                                   waifw.init = t(waifw.5.groups), state.Laine)



list[ N.5.group, M.5.group, cases.by.age.Kokota, state.Kokota, weeks] = output.data.for.given.sub.prefecture(S.prefect = "KOKOTA",total.pop =  14732)
list[Infections.Kokota.13.report2, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Kokota.est.report2.13, Kokota.est.report2.13$R02, 
                                                                                   time.length = 14, num.sims,
                                                                                   cases.by.age.group = cases.by.age.Kokota.13[1:which(weeks == 13),],
                                                                                   waifw.init = t(waifw.5.groups), state.Kokota)



list[ N.5.group, M.5.group, cases.by.age.Cu.Lola, state.Cu.Lola, weeks] = output.data.for.given.sub.prefecture(S.prefect = "CU LOLA",total.pop =  49993)
list[Infections.Cu.Lola.13.report2, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Cu.Lola.est.report2.13, Cu.Lola.est.report2.13$R02, 
                                                                                    time.length = 14, num.sims,
                                                                                    cases.by.age.group = cases.by.age.Cu.Lola.13[1:which(weeks == 13),],
                                                                                    waifw.init = t(waifw.5.groups), state.Cu.Lola)


Foum.ests2 = Foum.est.report2.13
Foum.sims2 = Infections.Foum.13.report2

Gueasso.ests2 = Gueasso.est.report2.13
Gueasso.sims2 = Infections.Gueasso.13.report2

Gama.ests2 = Gama.est.report2.13
Gama.sims2 = Infections.Gama.13.report2

Laine.ests2 = Laine.est.report213
Laine.sims2 = Infections.Laine.13.report2

Kokota.ests2 = Kokota.est.report2.13
Kokota.sims2 = Infections.Kokota.13.report2

Cu.Lola.ests2 = Cu.Lola.est.report2.13
Cu.Lola.sims2 = Infections.Cu.Lola.13.report2

save.image("Post_Week13_report2_fits.Rdata")
#save.image("Seasonal_June_12th_Report_2.RData")