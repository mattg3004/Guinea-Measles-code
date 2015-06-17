
Guinea_model_tsir_fit_observation_rate_seasonal_forcing <- '
data {
  int <lower = 0> N;         // week sample size
  int <lower = 0> M;         // number of age groups
  real Obs_cases_age[N, M];  // number of cases observed each week by age group
  real state[M];             // population size of each age group
  real <lower = 0> R0;       // average value of R0 for the model
  real R0_sd;                // standard deviation of R_0
  real waifw[M, M];          // waifw by age groups
  real sus_prior[M, 2];      // priors for susceptibility by age group from fits to NZoo
  real seasonal_multiplier[N, 1];  // this vector gives the scaling of R0 by week, i.e. increasing as we go back in time, looks like a cosine function with a certain wavelength and amplitude. This is scaled later on as we select different amplitudes for the seasonal forcing.
  real alpha_mean;           // mean of seasonal amplitude
  real alpha_s;             // s parameter of logistic function for seasonal amplitude
}

parameters {
  // fit the susceptible proportion by age given a value of R_0
  // s1 is proportion of susceptibles in group 1 (Under 1s), s2 for group 2 etc.
  real <lower = 0, upper = 1> s[M];   // susceptibility in number of age groups specified by M. Initially using Under 1s, 1-4, 5-10, 11-15, 16+
  real <lower = 0.05, upper = 1> delta;       // delta is the reporting rate of cases 
  real <lower = 0> Cases_by_age[N, M];     // true cases by age by week 
  real <lower = 0> R02;                    // value of R_0 chosen from normal distribution with mean and st dev given by fit to NZoo data.
  real <lower = 0.2, upper = 0.6> alpha;   // seasonal forcing amplitude
}

model {
  real C[N-1, M];                  // dummy variable holding the mean of expected cases in week n given cases in previous week.
  real Cum_cases[N, M];            // cumulative cases over age by week
  real beta;                       // transmission term dummy variable
  real sus;                        // number of susceptible dummy variable
  real waifw2[M, M]; //dummy waifw
  real R03;           // dummy R0 value, scaled seasonally
  real phi[N, M];                  // force of infection matrix - calculate this at every step given proposed true cases by age group
  
  s[1] ~ normal(sus_prior[1, 1], sus_prior[1, 2]);         // prior distribution for Under 1s from fitting NZoo.
  s[2] ~ normal(sus_prior[2, 1], sus_prior[2, 2]);         // prior distribution for 1-5 from fitting NZoo.
  s[3] ~ normal(sus_prior[3, 1], sus_prior[3, 2]);         // prior distribution for 6-10 from fitting NZoo.
  s[4] ~ normal(sus_prior[4, 1], sus_prior[4, 2]);         // prior distribution for 11-15 from fitting NZoo.
  s[5] ~ normal(sus_prior[5, 1], sus_prior[5, 2]);         // prior distribution for 16+ from fitting NZoo.
  alpha ~ logistic(alpha_mean, alpha_s);       // prior for amplitude of seasonal scaling. Choose from this distribution for each iteration
  R02 ~ normal(R0, R0_sd);          // prior for R0. Choose a value for each iteration
  
  
  for (m in 1 : M){  
    Obs_cases_age[1, m] ~ normal(delta * Cases_by_age[1, m], sqrt(delta * (1- delta) * Cases_by_age[1, m]));    // for first time step assume that the observed cases are distributed normally with mean (delta * true cases), where delta is the observation rate and true cases = Cases_by_age[1, m]
  
    Cum_cases[1, m] <- Cases_by_age[1, m];  // Update cumulative cases by age
  }
  
  
  // Next we loop over the number of weeks that have observed cases given by N.
  
  for(n in 2: N){
  
  
  // recalculate the waifw for sampled R0 and sampled seasonal amplitude alpha
    R03 <-  (R02 * (1 + ((seasonal_multiplier[n, 1] - 1) * alpha))) ;   // Rescale R0 to account for seasonality
    for(m1 in 1 : M){
      for(m2 in 1 : M){
        waifw2[m1, m2] <- -log(1 - (R03/R0)*(1 - exp(-waifw[m1, m2])))* sum(state);  // rescale the waifw to account for seasonality and change in R0
      }
    }
  
    for(m in 1 : M){
      phi[(n-1), m] <- 1 - exp(-(waifw2[m, 1] * Cases_by_age[n-1, 1] + waifw2[m, 2] * Cases_by_age[n-1, 2] + 
      waifw2[m, 3] * Cases_by_age[n-1, 3] + waifw2[m, 4] * Cases_by_age[n-1, 4] + 
      waifw2[m, 5] * Cases_by_age[n-1, 5])/ sum(state));   // Calculate the force of infection by age given the cases by age and the waifw
    }
  
    for(m in 1 : M){
  
      // Assume that the true cases by age at time n, given by Cases_by_age[n, m], are normally distributed dependent on the true cases at time n-1, which is contained in phi and also on the number of susceptible people left in the age group, given by (s1 * state[1] - Cum_cases[n-1, m]).
      // (1 - (1 - phi[n, m])^(R02/R0) term transforms force of infection to the required level for the sampled value of R0.
      
      beta <- phi[n-1, m];      //  set beta to be the force of infection on the given age group
      sus <- (s[m] * state[m] - Cum_cases[(n-1), m]);  // sus is the remaining susceptibles in the age group
      
      Cases_by_age[n, m] ~ normal(beta * sus, sqrt(beta * sus));
      Cum_cases[n, m] <- Cum_cases[n-1, m] + Cases_by_age[n, m];  // Update cumulative cases
    
    }
    
    for (m in 1 : M){
      Obs_cases_age[n, m] ~ normal(delta * Cases_by_age[n, m], sqrt(delta * (1- delta) * Cases_by_age[n, m]));  // Assume observed cases are normally distributed
    }
  
  }
}
'


sus_prior.seasonal = rbind(c(mean(N.Zoo.seasonality$s[, 1]), sd(N.Zoo.seasonality$s[, 1])), 
                           c(mean(N.Zoo.seasonality$s[, 2]), sd(N.Zoo.seasonality$s[, 2])),
                           c(mean(N.Zoo.seasonality$s[, 3]), sd(N.Zoo.seasonality$s[, 3])),
                           c(mean(N.Zoo.seasonality$s[, 4]), sd(N.Zoo.seasonality$s[, 4])),
                           c(mean(N.Zoo.seasonality$s[, 5]), sd(N.Zoo.seasonality$s[, 5])))

estimate.susceptibility.and.observation.5.group.seasonality <- function(waifw, R0, iters, chains, 
                                                                        R0_sd, alpha_mean, alpha_s, sus_prior,
                                                                        S.prefect, population){
  
  list[ N.5.group, M.5.group, cases.by.age.group.5.group, state.5.group, weeks] = output.data.for.given.sub.prefecture(S.prefect, population)
  waifw.input = output.waifw(waifw = t(waifw), R_0 = R0, state = state.5.group)
  
  ### we want to include seasonality. Estimated seasonality has 0.3 * cos shape, with cos = 1 approximately at weeks 13 and 39. Need to work out where the first week of observed cases is on this scale in order to scale correctly for seasonality.
  first.week.adjust = weeks[1] - 13
  N = length(weeks)
  seasonal.multiplier = matrix(0, N, 1)
  for(i in 1 : N){
    # waifw = output.waifw(waifw.init, R_0 * (1 + cos((3/12 - (7 - i)/48)  * 2 *  pi) * 0.3), state)
    seasonal.multiplier[i] = 1 + cos(pi / 2 +  ((first.week.adjust + i - 1 )/52)* 2 *  pi)
  }
  #waifw.input = output.waifw(waifw = t(waifw), R_0 = R0, state = state.5.group)
  waifw = output.waifw(waifw.init, R0, state) / sum(state)
  
  ### Guinea.data holds the data to input into the stan code
  Guinea.data <- list("Obs_cases_age" = cases.by.age.group.5.group, "N" = N, "M" = 5,  
                      "state" = state.5.group, "R0" = R0, "R0_sd" = R0_sd, "waifw" = waifw, 
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


list[Foum.seasonal.fit, Foum.seasonal.est] = estimate.susceptibility.and.observation.5.group.seasonality(waifw = t(waifw.5.groups), 
                                                                                                         R0 = mean(N.Zoo.seasonality$R02), 
                                                                                                         R0_sd = sd(N.Zoo.seasonality$R02), 
                                                                                                         iters = 2000, chains = 4, 
                                                                                                         alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                         alpha_s = sqrt(3 * (sd(N.Zoo.seasonality$alpha))^2 /(pi^2)), 
                                                                                                         sus_prior.seasonal, S.prefect = "FOUMBADOU", population = 19438)

list[Infections.Foum.seasonal, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Foum.seasonal.est, Foum.seasonal.est$R02, time.length = 14, num.sims,
                                                                               cases.by.age.group = cases.by.age.Foum,
                                                                               waifw.init = t(waifw.5.groups), state.Foum)



list[Gueasso.seasonal.fit, Gueasso.seasonal.est] = estimate.susceptibility.and.observation.5.group.seasonality(waifw =  t(waifw.5.groups), 
                                                                                                               R0 = mean(N.Zoo.seasonality$R02), 
                                                                                                               R0_sd = sd(N.Zoo.seasonality$R02), 
                                                                                                               iters = 5000, chains = 4, 
                                                                                                               alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                               alpha_s = sqrt(3 * (sd(N.Zoo.seasonality$alpha))^2 /(pi^2)), 
                                                                                                               sus_prior.seasonal, S.prefect = "GUEASSO", population = 21116)

list[ N.5.group, M.5.group, cases.by.age.Gueasso, state.Gueasso, weeks] = output.data.for.given.sub.prefecture(S.prefect = "GUEASSO",total.pop =  21116)

list[Infections.Gueasso.seasonal, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Gueasso.seasonal.est, Gueasso.seasonal.est$R02, time.length = 14, num.sims,
                                                                                  cases.by.age.group = cases.by.age.Gueasso,
                                                                                  waifw.init = t(waifw.5.groups), state.Gueasso)



list[Gama.seasonal.fit, Gama.seasonal.est] = estimate.susceptibility.and.observation.5.group.seasonality(waifw =  t(waifw.5.groups), 
                                                                                                         R0 = mean(N.Zoo.seasonality$R02), 
                                                                                                         R0_sd = sd(N.Zoo.seasonality$R02), 
                                                                                                         iters = 5000, chains = 4, 
                                                                                                         alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                         alpha_s = sqrt(3 * (sd(N.Zoo.seasonality$alpha))^2 /(pi^2)), 
                                                                                                         sus_prior.seasonal, S.prefect = "GAMA BEREMA", population = 20465)


list[ N.5.group, M.5.group, cases.by.age.Gama, state.Gama, weeks] = output.data.for.given.sub.prefecture(S.prefect = "GAMA BEREMA",total.pop =  20465)
list[Infections.Gama.seasonal, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Gama.seasonal.est, Gama.seasonal.est$R02, 
                                                                               time.length = 14, num.sims,
                                                                               cases.by.age.group = cases.by.age.Gama,
                                                                               waifw.init = t(waifw.5.groups), state.Gama)


list[Laine.seasonal.fit, Laine.seasonal.est] = estimate.susceptibility.and.observation.5.group.seasonality(waifw =  t(waifw.5.groups), 
                                                                                                           R0 = mean(N.Zoo.seasonality$R02), 
                                                                                                           R0_sd = sd(N.Zoo.seasonality$R02), 
                                                                                                           iters = 5000, chains = 4, 
                                                                                                           alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                           alpha_s = sqrt(3 * (sd(N.Zoo.seasonality$alpha))^2 /(pi^2)), 
                                                                                                           sus_prior.seasonal, S.prefect = "LAINE", population = 16591)


list[ N.5.group, M.5.group, cases.by.age.Laine, state.Laine, weeks] = output.data.for.given.sub.prefecture(S.prefect = "LAINE",total.pop =  16591)
list[Infections.Laine.seasonal, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Laine.seasonal.est, Laine.seasonal.est$R02, 
                                                                                time.length = 14, num.sims,
                                                                                cases.by.age.group = cases.by.age.Laine,
                                                                                waifw.init = t(waifw.5.groups), state.Laine)




list[Kokota.seasonal.fit, Kokota.seasonal.est] = estimate.susceptibility.and.observation.5.group.seasonality(waifw =  t(waifw.5.groups), 
                                                                                                             R0 = mean(N.Zoo.seasonality$R02), 
                                                                                                             R0_sd = sd(N.Zoo.seasonality$R02), 
                                                                                                             iters = 5000, chains = 4, 
                                                                                                             alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                             alpha_s = sqrt(3 * (sd(N.Zoo.seasonality$alpha))^2 /(pi^2)), 
                                                                                                             sus_prior.seasonal, S.prefect = "KOKOTA", population = 14732)



list[ N.5.group, M.5.group, cases.by.age.Kokota, state.Kokota, weeks] = output.data.for.given.sub.prefecture(S.prefect = "KOKOTA",total.pop =  14732)
list[Infections.Kokota.seasonal, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Kokota.seasonal.est, Kokota.seasonal.est$R02, 
                                                                                 time.length = 14, num.sims,
                                                                                 cases.by.age.group = cases.by.age.Kokota,
                                                                                 waifw.init = t(waifw.5.groups), state.Kokota)



list[Cu.Lola.seasonal.fit, Cu.Lola.seasonal.est] = estimate.susceptibility.and.observation.5.group.seasonality(waifw =  t(waifw.5.groups), 
                                                                                                               R0 = mean(N.Zoo.seasonality$R02), 
                                                                                                               R0_sd = sd(N.Zoo.seasonality$R02), 
                                                                                                               iters = 5000, chains = 4, 
                                                                                                               alpha_mean = mean(N.Zoo.seasonality$alpha), 
                                                                                                               alpha_s = sqrt(3 * (sd(N.Zoo.seasonality$alpha))^2 /(pi^2)), 
                                                                                                               sus_prior.seasonal, S.prefect = "CU LOLA", population = 49993)


list[ N.5.group, M.5.group, cases.by.age.Cu.Lola, state.Cu.Lola, weeks] = output.data.for.given.sub.prefecture(S.prefect = "CU LOLA",total.pop =  49993)

list[Infections.Cu.Lola.seasonal, a, i] = simulate.with.R0.and.sus.dist.obs.rates(Cu.Lola.seasonal.est, Cu.Lola.seasonal.est$R02, 
                                                                                  time.length = 14, num.sims,
                                                                                  cases.by.age.group = cases.by.age.Cu.Lola,
                                                                                  waifw.init = t(waifw.5.groups), state.Cu.Lola)


########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################

###############                 LABEL  RESULTS FOR EASE OF USE LATER                 ###############

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################

Foum.ests3 = Foum.seasonal.est
Foum.sims3 = Infections.Foum.seasonal

Gueasso.ests3 = Gueasso.seasonal.est
Gueasso.sims3 = Infections.Gueasso.seasonal

Gama.ests3 = Gama.seasonal.est
Gama.sims3 = Infections.Gama.seasonal

Laine.ests3 = Laine.seasonal.est
Laine.sims3 = Infections.Laine.seasonal

Kokota.ests3 = Kokota.seasonal.est
Kokota.sims3 = Infections.Kokota.seasonal

Cu.Lola.ests3 = Cu.Lola.seasonal.est
Cu.Lola.sims3 = Infections.Cu.Lola.seasonal

save.image("Post_Week17_fits.Rdata")