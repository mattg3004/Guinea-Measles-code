library(rstan)
source("Functions.R")
source("Functions_for_tsir_fits.R")
source("Forward_simulations.R")
Guinea_model_tsir_with_susceptible_prior_and_seasonality <- '
data {
  int<lower=0> N;       // week sample size
  int<lower=0> M;       // number of age groups
  int y[N, M];          // infection data by age group
  int x[N, M];          // infection data by age group, 2 weeks previously
  real state[M];        // population size of each age group
  int cum_cases[N, M];  // cumulative cases by week and age group
  real <lower = 0> R0;  // average value of R0 for the model
  real R0_sd;           // standard deviation for R0
  real waifw[M, M];              // waifw by age groups
  real seasonal_multiplier[N, 1];  // this vector gives the scaling of R0 by week, i.e. increasing as we go back in time, looks like a cosine function with a certain wavelength and amplitude. This is scaled later on as we select different amplitudes for the seasonal forcing.

}

parameters {
// fit the susceptible proportion by age given a value of R_0
  real <lower = 0, upper = 1> s[M];   // susceptibility in number of age groups specified by M. Initially using Under 1s, 1-4, 5-10, 11-15, 16+
  real <lower = 10> R02;                         // fit a value of R0 
  real <lower = 0.2, upper = 0.6> alpha;   // seasonal forcing amplitude
}

model {

  // We assume that the number of cases at time t is dependent on number of cases from two weeks previous (this is contained in phi) and the number of susceptibles in the age group. Assuming that the epidemic the reporting rate is 1 the likelihood of observing the cases we did is assumed to be given by a poisson distribution.
  // s[1] * state[1] gives the initial number of susceptibles in group 1 and subtracting cum_cases from this gives the number of susceptibles at time t. Same applies for other groups. Looping over m is for the different age groups, over n gives the different weeks in the dataset.
  
  real waifw2[M, M]; //dummy waifw
  real R03;           // dummy R0 value, scaled seasonally
  real phi[N, M];                  // force of infection matrix - calculate this at every
  
  s[2] ~ normal(0.28, 0.03);         // prior distribution of susceptibles for 1-5s from Science method. Choose s[2] from this for each iteration
  alpha ~ logistic(0.3, 0.05);       // prior for amplitude of seasonal scaling. Choose from this distribution for each iteration
  R02 ~ normal(R0, R0_sd);          // prior for R0. Choose a value for each iteration
  
  for (m in 1 : M){
    for (n in 1 : N){
  
      R03 <- (R02 * (1 + ((seasonal_multiplier[n, 1] - 1) * alpha))) ;   // Rescale R0 to account for seasonality
      for(m1 in 1 : M){
        for(m2 in 1 : M){
          waifw2[m1, m2] <- -log(1 - (R03/R0)*(1 - exp(-waifw[m1, m2])))* sum(state);  // rescale the waifw to account for seasonality and change in R0
        }
      }
  
      phi[n, m] <- 1 - exp(-(waifw2[m, 1] * x[n, 1] ^ 0.97 + waifw2[m, 2] * x[n, 2] ^ 0.97 + 
      waifw2[m, 3] * x[n, 3] ^ 0.97 + waifw2[m, 4] * x[n, 4] ^ 0.97 + waifw2[m, 5] * x[n, 5] ^ 0.97)/ sum(state)); // Calculate the force of infection by age given the cases by age and the waifw
  
      y[n, m] ~ poisson((s[m] * state[m] - cum_cases[n, m]) * (1 - (1 - phi[n, m])));  // likelihood of observed cases follows Poisson(S_t phi(t))
    }
  }
}
'



estimate.susceptibility.5.group.N.Zoo.with.seasonality <- function(waifw.init, x, y, R0, state, iters, chains, R0_sd){
  
  ### N is the number of weeks of cases, M is the number of age groups. Can get both of these from the dimensions of case data (x or y).
  ### M is also the number of rows and columns of the waifw.
  N = length(x[, 1])
  M = length(waifw.init[, 1])
  
  ### denom is the total population size, used to scale the WAIFW
  denom = sum(state)
  
  ### cum_cases holds cumulative cases at each week by age group. Used in the TSIR code to amend the susceptible population in each age group as observed cases occur in the dataset.
  cum_cases = matrix(0, N, M)
  for(i in 1 : M){
    cum_cases[, i] = cumsum(x[, i])
  }
  
  ### waifw is waifw.init scaled according to the given value of R_0. This is calculated for every week, as we include seasonality with an amplitude of 0.3
  ### phi is the attack rate by week and age group, which is dependent on the waifw and the cases from two weeks previous
  phi = matrix(0, N, M)
  seasonal.multiplier = matrix(0, N, 1)
  for(i in 1 : N){
    # waifw = output.waifw(waifw.init, R_0 * (1 + cos((3/12 - (7 - i)/48)  * 2 *  pi) * 0.3), state)
    seasonal.multiplier[i] = 1 + cos((3/12 - (7 - i)/48)  * 2 *  pi)
    #waifw = output.waifw(waifw.init, R_0, state) / sum(state)
    #phi[i, ] = calc.phi.1.week(waifw, x[i, ], denom)
  }
  waifw = output.waifw(waifw.init, R0, state) / sum(state)
  
  ### Guinea.data holds the data to input into the stan code
  Guinea.data <- list("y" = y, "x" = x, "N" = N, "M" = M,  "state" = state, "cum_cases" = cum_cases, "R0" = R0, "R0_sd" = R0_sd, 
                      "waifw" = waifw, "seasonal_multiplier" = seasonal.multiplier)
  
  ### Here the stan model is run to estimate the susceptibility in the 5 age groups. 
  fit <- stan(model_code = Guinea_model_tsir_with_susceptible_prior_and_seasonality , data = Guinea.data, iter = iters , chains = chains) 
  
  ### print the output of the stan model for informational purposes
  print(fit)
  
  ### extract the output of the stan model, which is also the output of this function.
  sus.est = extract(fit)
  
  return(list(fit, sus.est))
}



waifw.5.groups =  rbind(c(6.9 , 3 , 3/2 , 2/3,  0.5 ) * state.5.group[1],
                        c(3 , 6.9 , 3/2 , 2/3,  0.5 ) * state.5.group[2],
                        c(3/2 , 2/3, 3.5 , 3.5, 0.5 ) * state.5.group[3],
                        c(3/2 , 2/3, 3.5 , 3.5, 0.5 ) * state.5.group[4],
                        c(0.5 , 0.5, 0.5 , 0.5, 1.6 ) * state.5.group[5])

time.length = 10
source("Data_NZoo_week13.R")
list[N.Zoo.seasonality.fit.13, N.Zoo.seasonality.13] = estimate.susceptibility.5.group.N.Zoo.with.seasonality(waifw.init = t(waifw.5.groups), x = x.5.group[1:8,],
                                                                                                              y = y.5.group[1:8,], R0 = 18, state = state.5.group, 
                                                                                                              iters = 2000, chains = 4, R0_sd = 4)

list[Infections.NZoo.13, a, i] = simulate.with.given.R0.and.sus.dist.waifw.groups.5.groups(sus.dist = N.Zoo.seasonality.13, 
                                                                                           time.length = 14, 
                                                                                           num.sims = 1000,
                                                                                           cases.by.age.group = cases.by.age.group.5.group,
                                                                                           t(waifw.5.groups), 
                                                                                           state.5.group)





source("Data_NZoo.R")     


list[N.Zoo.seasonality.fit.13.report2, N.Zoo.seasonality.13.report2] = estimate.susceptibility.5.group.N.Zoo.with.seasonality(t(waifw.5.groups), x.5.group[1:8, ], y.5.group[1:8,], R0 = 18,
                                                                                                                              state.5.group, iters = 2000, chains = 4, R0_sd = 4)

list[N.Zoo.seasonality.fit, N.Zoo.seasonality] = estimate.susceptibility.5.group.N.Zoo.with.seasonality(t(waifw.5.groups), x.5.group, y.5.group, R0 = 18,
                                                                                                        state.5.group, iters = 2000, chains = 4, R0_sd = 4)

list[Infections.NZoo, attack.rate.NZoo, initial.NZoo] = simulate.with.given.R0.and.sus.dist.waifw.groups.5.groups(sus.dist = N.Zoo.seasonality, 
                                                                                                                  time.length = 10, 
                                                                                                                  num.sims = 1000,
                                                                                                                  cases.by.age.group = cases.by.age.group.5.group,
                                                                                                                  t(waifw.5.groups), 
                                                                                                                  state.5.group)


list[Infections.NZoo.13.report.2, a, i] = simulate.with.given.R0.and.sus.dist.waifw.groups.5.groups(sus.dist = N.Zoo.seasonality.13.report2, 
                                                                                                    time.length = 14, 
                                                                                                    num.sims = 1000,
                                                                                                    cases.by.age.group = cases.by.age.group.5.group[1:10, ],
                                                                                                    t(waifw.5.groups), 
                                                                                                    state.5.group)

NZoo.ests1 = N.Zoo.seasonality.13
NZoo.ests2 = N.Zoo.seasonality.13.report2
NZoo.ests3 = N.Zoo.seasonality

NZoo.sims1 = Infections.NZoo.13
NZoo.sims2 = Infections.NZoo.13.report.2
NZoo.sims3 = Infections.NZoo
save.image("Post_NZoo_fits.RData")