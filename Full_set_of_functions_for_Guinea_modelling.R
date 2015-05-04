library(rstan)

Guinea_model_tsir_with_susceptible_prior_2 <-{ '
data {
  int<lower=0> N;       // week sample size
  int<lower=0> M;       // number of age groups
  int y[N, M];          // infection data by age group
  real phi[N, M];       // force of infection on each age group by week, dependent on WAIFW and R_0 - Klepac paper method. This has been calculated for R_0 = 18
  real state[M];        // population size of each age group
  int cum_cases[N, M];  // cumulative cases by week and age group
  real <lower = 0> R0;  // average value of R0 for the model
  real R0_sd;           // standard deviation for R0

}

parameters {
// fit the susceptible proportion by age given a value of R_0
  real <lower = 0, upper = 1> s1;          // s1 - susceptibility in Under 1s
  real <lower = 0, upper = 1> s2;          // s2 - susceptibility in 1-5
  real <lower = 0, upper = 1> s3;          // s3 - susceptibility in 6-10
  real <lower = 0, upper = 1> s4;          // s4 - susceptibility in 11-15
  real <lower = 0, upper = 1> s5;          // s5 - susceptibility in 16+
  real R02;                         // fit a value of R0 
}

model {

// We assume that the number of cases at time t is dependent on number of cases from two weeks previous (this is contained in phi) and the number of susceptibles in the age group. Assuming that the epidemic the reporting rate is 1 the likelihood of observing the cases we did is assumed to be given by a poisson distribution.
// s1 * state[1] gives the initial number of susceptibles in group 1 and subtracting cum_cases from this gives the number of susceptibles at time t. Same applies for other groups. Looping over m is for the different age groups, over n gives the different weeks in the dataset.
  
  s2 ~ normal(0.28, 0.03);         // prior distribution of susceptibles for 1-5s from Science method.
  R02 ~ normal(R0, 4);             // Prior distribution for R_0
  
  for (m in 1 : M){
    for (n in 1 : N){
      if(m == 1){
        y[n, m] ~ poisson((s1 * state[1] - cum_cases[n, m]) * phi[n, m]);  
      }
      if(m == 2){
        y[n, m] ~ poisson((s2 * state[2] - cum_cases[n, m]) * phi[n, m]);  
      }
      if(m == 3){
        y[n, m] ~ poisson((s3 * state[3] - cum_cases[n, m]) * phi[n, m]);  
      }
      if(m == 4){
        y[n, m] ~ poisson((s4 * state[4] - cum_cases[n, m]) * phi[n, m]);  
      }
      if(m == 4){
        y[n, m] ~ poisson((s5 * state[5] - cum_cases[n, m]) * phi[n, m]);  
      }
    }
  }
}
'
}


Guinea_model_tsir_with_susceptible_prior <- '
data {
  int<lower=0> N;       // week sample size
  int<lower=0> M;       // number of age groups
  int y[N, M];          // infection data by age group
  real phi[N, M];       // force of infection on each age group by week, dependent on WAIFW and R_0 - Klepac paper method. This has been calculated for R_0 = 18
  real state[M];        // population size of each age group
  int cum_cases[N, M];  // cumulative cases by week and age group
  real <lower = 0> R0;  // average value of R0 for the model
  real R0_sd;           // standard deviation for R0

}

parameters {
// fit the susceptible proportion by age given a value of R_0
  real <lower = 0, upper = 1> s[M];   // susceptibility in number of age groups specified by M. Initially using Under 1s, 1-4, 5-10, 11-15, 16+
  real R02;                         // fit a value of R0 
}

model {

// We assume that the number of cases at time t is dependent on number of cases from two weeks previous (this is contained in phi) and the number of susceptibles in the age group. Assuming that the epidemic the reporting rate is 1 the likelihood of observing the cases we did is assumed to be given by a poisson distribution.
// s1 * state[1] gives the initial number of susceptibles in group 1 and subtracting cum_cases from this gives the number of susceptibles at time t. Same applies for other groups. Looping over m is for the different age groups, over n gives the different weeks in the dataset.

  s[2] ~ normal(0.28, 0.03);         // prior distribution of susceptibles for 1-5s from Science method.
  R02 ~ normal(R0, 4);             // Prior distribution for R_0

  for (m in 1 : M){
    for (n in 1 : N){
      y[n, m] ~ poisson((s[m] * state[m] - cum_cases[n, m]) * (1 - (1 - phi[n, m])^(R02/R0)));  // likelihood of observed cases follows Poisson(S_t phi(t))
    }
  }
}
'


estimate.susceptibility.5.group.N.Zoo <- function(waifw.init, x, y, R0, state, iters, chains, R0_sd){
  
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
  for(i in 1 : N){
    waifw = output.waifw(waifw.init, R_0 * (1 + cos((3/12 - (7 - i)/48)  * 2 *  pi) * 0.3), state)
    phi[i, ] = calc.phi.1.week(waifw, x[i, ], denom)
  }
  
  ### Guinea.data holds the data to input into the stan code
  Guinea.data <- list("y" = y, "N" = N, "M" = M, "phi" = phi, "state" = state, "cum_cases" = cum_cases, "R0" = R0, "R0_sd" = R0_sd)
  
  ### Here the stan model is run to estimate the susceptibility in the 5 age groups. 
  fit <- stan(model_code = Guinea_model_tsir_with_susceptible_prior, data = Guinea.data, iter = iters , chains = chains) 
  
  ### print the output of the stan model for informational purposes
  print(fit)
  
  ### extract the output of the stan model, which is also the output of this function.
  sus.est = extract(fit)
  
  return(sus.est)
}

source("Data_NZoo.R")

waifw.5.groups =  rbind(c(6.9 , 3 , 3/2 , 2/3,  0.5 ) * state.5.group[1],
                        c(3 , 6.9 , 3/2 , 2/3,  0.5 ) * state.5.group[2],
                        c(3/2 , 2/3, 3.5 , 3.5, 0.5 ) * state.5.group[3],
                        c(3/2 , 2/3, 3.5 , 3.5, 0.5 ) * state.5.group[4],
                        c(0.5 , 0.5, 0.5 , 0.5, 1.6 ) * state.5.group[5])

NZoo.estimates = estimate.susceptibility.5.group.N.Zoo(t(waifw.5.groups), x.5.group, y.5.group, R0 = 18, state.5.group, iters = 2000, chains = 4, R0_sd = 4)




Guinea_model_tsir_fit_observation_rate<- '
data {
  int <lower = 0> N;         // week sample size
  int <lower = 0> M;         // number of age groups
  real Obs_cases_age[N, M];  // number of cases observed each week by age group
  real Cases_week_obs[N];    // cases observed each week
  real state[M];             // population size of each age group
  real <lower = 0> R0;       // average value of R0 for the model
  real R0_sd;                // standard deviation of R_0
  real waifw[M, M];          // waifw by age groups
  real sus_prior[M, 2];      // priors for susceptibility by age group from fits to NZoo
  int Late_n;               // point at which week 14 occurs in data, as this is the point at which R effecrtive seems to shift to a lower value.
}

parameters {
// fit the susceptible proportion by age given a value of R_0
// s1 is proportion of susceptibles in group 1 (Under 1s), s2 for group 2 etc.
  real <lower = 0, upper = 1> s[M];   // susceptibility in number of age groups specified by M. Initially using Under 1s, 1-4, 5-10, 11-15, 16+
  real <lower = 0.05, upper = 1> delta;       // delta is the reporting rate of cases 
  real <lower = 0, upper = 100> Cases_by_age[N, M];     // true cases by age by week 
  real <lower = 0> R02;                    // value of R_0 chosen from normal distribution with mean and st dev given by fit to NZoo data.
  real <lower = 0> R0_Late;                    // R0 for the later stage of the epidemic where spread has decreased.
}

model {
  real C[N-1, M];                  // dummy variable holding the mean of expected cases in week n given cases in previous week.
  real Cum_cases[N, M];            // cumulative cases over age by week
  real beta;                       // transmission term dummy variable
  real sus;                        // number of susceptible dummy variable
  real phi[N, M];                  // force of infection matrix - calculate this at every step given proposed true cases by age group
  R02 ~ normal(R0, R0_sd);         // Input posterior from previous fit for R_0
  R0_Late ~ normal(R0, R0_sd);
  s[1] ~ normal(sus_prior[1, 1], sus_prior[1, 2]);         // prior distribution for Under 1s from fitting NZoo.
  s[2] ~ normal(sus_prior[2, 1], sus_prior[2, 2]);         // prior distribution for 1-5 from fitting NZoo.
  s[3] ~ normal(sus_prior[3, 1], sus_prior[3, 2]);         // prior distribution for 6-10 from fitting NZoo.
  s[4] ~ normal(sus_prior[4, 1], sus_prior[4, 2]);         // prior distribution for 11-15 from fitting NZoo.
  s[5] ~ normal(sus_prior[5, 1], sus_prior[5, 2]);         // prior distribution for 16+ from fitting NZoo.

for (m in 1 : M){  
  Obs_cases_age[1, m] ~ normal(delta * Cases_by_age[1, m], sqrt(delta * (1- delta) * Cases_by_age[1, m]));    // for first time step assume that the observed cases are distributed normally with mean (delta * true cases), where delta is the observation rate and true cases = Cases_by_age[1, m]

  Cum_cases[1, m] <- Cases_by_age[1, m];  // Update cumulative cases by age
}


// Next we loop over the number of weeks that have observed cases given by N.

  for(n in 2: (Late_n - 1)){
  
    for(m in 1 : M){
    phi[(n-1), m] <- 1 - exp(-(waifw[m, 1] * Cases_by_age[n-1, 1] + waifw[m, 2] * Cases_by_age[n-1, 2] + waifw[m, 3] * Cases_by_age[n-1, 3] + waifw[m, 4] * Cases_by_age[n-1, 4] + waifw[m, 5] * Cases_by_age[n-1, 5])/ sum(state));   // Calculate the force of infection by age given the cases by age and the waifw
    }
  
  
    for(m in 1 : M){
  
  // Assume that the true cases by age at time n, given by Cases_by_age[n, m], are normally distributed dependent on the true cases at time n-1, which is contained in phi and also on the number of susceptible people left in the age group, given by (s1 * state[1] - Cum_cases[n-1, m]).
  // (1 - (1 - phi[n, m])^(R02/R0) term transforms force of infection to the required level for the sampled value of R0.
  
    beta <- (1 - (1 - phi[n-1, m])^(R02/R0));      //  set beta to be the force of infection on the given age group
    sus <- (s[m] * state[m] - Cum_cases[(n-1), m]);  // sus is the remaining susceptibles in the age group
  
    Cases_by_age[n, m] ~ normal(beta * sus, sqrt(beta * sus));
    Cum_cases[n, m] <- Cum_cases[n-1, m] + Cases_by_age[n, m];  // Update cumulative cases
  
    }
  
    for (m in 1 : M){
  Obs_cases_age[n, m] ~ normal(delta * Cases_by_age[n, m], sqrt(delta * (1- delta) * Cases_by_age[n, m]));  // Assume observed cases are normally distributed
    }
  }


  // Now attempt to fit for the latter portion of the cases, where R effectvie is reduced.
  for(n in Late_n : N){
  
    for(m in 1 : M){
    phi[(n-1), m] <- 1 - exp(-(waifw[m, 1] * Cases_by_age[n-1, 1] + waifw[m, 2] * Cases_by_age[n-1, 2] + waifw[m, 3] * Cases_by_age[n-1, 3] + waifw[m, 4] * Cases_by_age[n-1, 4] + waifw[m, 5] * Cases_by_age[n-1, 5])/ sum(state));   // Calculate the force of infection by age given the cases by age and the waifw
    }
  
  
    for(m in 1 : M){
  
  // Assume that the true cases by age at time n, given by Cases_by_age[n, m], are normally distributed dependent on the true cases at time n-1, which is contained in phi and also on the number of susceptible people left in the age group, given by (s1 * state[1] - Cum_cases[n-1, m]).
  // (1 - (1 - phi[n, m])^(R0_Late/R0) term transforms force of infection to the required level for the sampled value of R0.
  
    beta <- (1 - (1 - phi[n-1, m])^(R0_Late/R0));      //  set beta to be the force of infection on the given age group
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




output.data.for.given.sub.prefecture <- function(S.prefect, total.pop, week.13.only = 0){
  ### read in data on cases
  data = read.csv("Lola_Measles_week_17.csv")
  #data = read.csv("Lola_Measles_week_13_included.csv")
  
  ### subset to cases in given sub prefecture
  S.prefect.data = subset(data, as.character(data$Sous.Prefecture) == S.prefect)
  week = matrix(0, nrow(S.prefect.data), 1)
  for(i in 1 : nrow(S.prefect.data)){
    week[i] = as.numeric(strsplit(as.character(S.prefect.data$Semaine.N.), "-")[[i]][2])
  }
  S.prefect.data$Week = week
  A =  as.character(unique(S.prefect.data$Week)[order((unique(S.prefect.data$Week)))])
  A = A[A != ""]
  A = seq(min(as.numeric(A)), 17)
  if (week.13.only == 1){
    A = seq(min(as.numeric(A)), 13)
  }
  ## split into cases by week
  cases.by.week = matrix(0, length(A), 2)
  cases.by.week[, 1] = as.character(A)
  for(i in 1:length(A)){
    cases.by.week[i, 2] = length(which(S.prefect.data$Week == A[i]))
  }
  
  ### add in age in years variable
  #S.prefect.data$Age.in.years = floor(S.prefect.data$Age.Mois/ 12)
  S.prefect.data$Age.in.years = floor(S.prefect.data$AGE)
  ### calculate cases by age in years for each week where data exists for given sub prefecture
  cases.by.age.in.years = matrix(0, length(A), (max(18, S.prefect.data$Age.in.years) + 1))
  
  for(i in 1 : (max(S.prefect.data$Age.in.years) + 1)){
    for(j in 1 : length(A)){
      A1 = subset(S.prefect.data, ((S.prefect.data$Age.in.years) == (i-1)) & (S.prefect.data$Week == A[j]))
      if(length(A1[, 1]) > 0){
        cases.by.age.in.years[j, i] = length(A1[, 1])  
      }
    }
  }
  
  ### convert cases by year to cases by age group for 4 group version
  cases.by.age.group = matrix(0, length(A), 4)
  cases.by.age.group[, 1] = cases.by.age.in.years[, 1]
  cases.by.age.group[, 2] = rowSums(cases.by.age.in.years[, 2:6])
  cases.by.age.group[, 3] = rowSums(cases.by.age.in.years[, 7:(max(18,(max(S.prefect.data$Age.in.years) + 1)))])
  cases.by.age.group[, 4] = rowSums(cases.by.age.in.years[, 16:(max(18,(max(S.prefect.data$Age.in.years) + 1)))])
  
  ### specify population in each age group.
  state = matrix(0, 4, 1)
  state[1] = round((664 / 15559)* total.pop)
  state[2] = round((3318 - 664 + 600) * total.pop/15559)
  state[3] = round((round(0.425 * 15559 ) - 3318 - 600) * total.pop/15559)
  state[4] = round((total.pop - sum(state[1:3])) )
  state = as.numeric(state)
  
  
  denom = sum(state)
  
  
  #x = cases.by.age.group[1:(length(unique(S.prefect.data$Semaine.N.)) - 2), ]
  # y = cases.by.age.group[3:length(unique(S.prefect.data$Semaine.N.)), ]
  N = (length(A) - 2)
  M = 4
  
  
  ### calculate cases by age group when we consider 5 group version
  cases.by.age.group.5.group = matrix(0, length(A), 5)
  cases.by.age.group.5.group[, 1] = cases.by.age.in.years[, 1]
  cases.by.age.group.5.group[, 2] = rowSums(cases.by.age.in.years[, 2:6])
  cases.by.age.group.5.group[, 3] = rowSums(cases.by.age.in.years[, 7:11])
  cases.by.age.group.5.group[, 4] = rowSums(cases.by.age.in.years[, 12:16])
  cases.by.age.group.5.group[, 5] = rowSums(cases.by.age.in.years[, 17:(max(18,(max(S.prefect.data$Age.in.years) + 1)))])
  
  ### specify population in each age group.
  state.5.group = matrix(0, 5, 1)
  state.5.group[1] = round(664 * total.pop/15559 )
  state.5.group[2] = round((3318 - 664 + 600) * total.pop/15559)
  state.5.group[3] = round(((round(0.425 * 15559 ) - 3318 - 600) / 2)* total.pop/15559)
  state.5.group[4] = round((((round(0.425 * 15559 ) - 3318- 600) / 2) + 300* total.pop/15559)* total.pop/15559)
  state.5.group[5] = round((total.pop - sum(state.5.group[1:4])))
  state.5.group = as.numeric(state.5.group)
  
  
  denom.5.group = sum(state.5.group)
  
  
  # x.5.group = cases.by.age.group.5.group[1:(length(unique(S.prefect.data$Semaine.N.)) - 2), ]
  # y.5.group = cases.by.age.group.5.group[3:length(unique(S.prefect.data$Semaine.N.)), ]
  N.5.group = (length(A) - 2)
  M.5.group = 5
  
  ### output data related to cases and population for 4 and 5 group models.
  return(list(N, M, cases.by.age.group, state, N.5.group, M.5.group, cases.by.age.group.5.group, state.5.group, A))
}

sus_prior = rbind(c(mean(NZoo.estimates$s[, 1]), sd(NZoo.estimates$s[, 1])), 
                  c(mean(NZoo.estimates$s[, 2]), sd(NZoo.estimates$s[, 2])),
                  c(mean(NZoo.estimates$s[, 3]), sd(NZoo.estimates$s[, 3])),
                  c(mean(NZoo.estimates$s[, 4]), sd(NZoo.estimates$s[, 4])),
                  c(mean(NZoo.estimates$s[, 5]), sd(NZoo.estimates$s[, 5])))
                

fit.for.given.sub.prefecture <- function(S.prefect, population, R0, R0_sd, waifw, iters, chains, sus_prior, week.13.only = 0){
  list[N, M, cases.by.age.group, state, N.5.group, M.5.group, cases.by.age.group.5.group, state.5.group, weeks] = output.data.for.given.sub.prefecture(S.prefect, population, week.13.only)
  waifw.input = output.waifw(waifw = t(waifw), R_0 = R0, state = state.5.group)
  
  S.prefect.data <- list("Obs_cases_age" = cases.by.age.group.5.group, "Cases_week_obs" = rowSums(cases.by.age.group.5.group), "N" = length(cases.by.age.group.5.group[, 1]), "M" = 5, "state" = state.5.group, "R0" = R0, "R0_sd" = R0_sd, "waifw" = waifw.input, "sus_prior" = sus_prior, "Late_n" = which(weeks == 14))
  
  fit <- stan(model_code = Guinea_model_tsir_fit_observation_rate, data = S.prefect.data, iter = iters , chains = chains) 
  print(fit)
  ests = extract(fit)
  plot.susceptibility.estimates.5.group.one.R0(ests, S.prefect)
  plot.reporting.rate.estimates(ests, S.prefect)
  return(list(ests, fit))
}




fit.for.given.sub.prefecture.no.prior <- function(S.prefect, population, R0, R0_sd, waifw, iters, chains, sus_prior){
  list[N, M, cases.by.age.group, state, N.5.group, M.5.group, cases.by.age.group.5.group, state.5.group, weeks] = output.data.for.given.sub.prefecture(S.prefect, population)
  waifw.input = output.waifw(waifw = t(waifw), R_0 = R0, state = state.5.group)
  
  S.prefect.data <- list("Obs_cases_age" = cases.by.age.group.5.group, "Cases_week_obs" = rowSums(cases.by.age.group.5.group), "N" = length(cases.by.age.group.5.group[, 1]), "M" = 5, "state" = state.5.group, "R0" = R0, "R0_sd" = R0_sd, "waifw" = waifw.input, "sus_prior" = sus_prior)
  
  fit <- stan(model_code = Guinea_model_tsir_fit_observation_rate_no_sus_prior, data = S.prefect.data, iter = iters , chains = chains) 
  print(fit)
  ests = extract(fit)
  plot.susceptibility.estimates.5.group.one.R0(ests, S.prefect)
  plot.reporting.rate.estimates(ests, S.prefect)
  return(list(ests, fit))
}
list[Kokota.ests, Kokota.fits] = fit.for.given.sub.prefecture.no.prior("KOKOTA", population = 14732, 
                                                              R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02), 
                                                              waifw = waifw.5.groups, iters = 2000, chains = 4, sus_prior)

list[Foum.ests, Foum.fits] = fit.for.given.sub.prefecture("FOUMBADOU", population = 19438, 
                                                          R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02), 
                                                          waifw = waifw.5.groups, iters = 1000, chains = 4, sus_prior)
list[Kokota.ests, Kokota.fits] = fit.for.given.sub.prefecture("KOKOTA", population = 14732, 
                                                          R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02), 
                                                          waifw = waifw.5.groups, iters = 2000, chains = 4, sus_prior)

list[Laine.ests, Laine.fits] = fit.for.given.sub.prefecture("LAINE", population = 16591, 
                                                              R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02), 
                                                              waifw = waifw.5.groups, iters = 2000, chains = 4, sus_prior)
list[Laine.ests.no.prior, Laine.fits.no.prior] = fit.for.given.sub.prefecture.no.prior("LAINE", population = 16591, 
                                                                     R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02), 
                                                                     waifw = waifw.5.groups, iters = 2000, chains = 4, sus_prior)


list[Gama.Berema.ests, Gama.Berema.fits] = fit.for.given.sub.prefecture("GAMA BEREMA", population = 20465, 
                                                                        R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02),
                                                                        waifw = waifw.5.groups, iters = 1000, chains = 4, sus_prior, week.13.only = 0)

list[Gueasso.ests, Gueasso.fits] = fit.for.given.sub.prefecture("GUEASSO", population = 21116, 
                                                                R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02),
                                                                waifw = waifw.5.groups, iters = 2000, chains = 4, sus_prior)

list[Kokota.ests, Kokota.fits] = fit.for.given.sub.prefecture("KOKOTA", population = 14732, 
                                                                R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02),
                                                                waifw = waifw.5.groups, iters = 1000, chains = 4, sus_prior)

list[Cu.Lola.ests, Cu.Lola.fits] = fit.for.given.sub.prefecture("CU LOLA", population = 49933, 
                                                              R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02),
                                                              waifw = waifw.5.groups, iters = 1000, chains = 4, sus_prior)


list[N, M, cases.by.age.group, state, N.5.group, M.5.group, cases.by.age.group.5.group, state.5.group, weeks] = output.data.for.given.sub.prefecture(S.prefect, population)
waifw.input = output.waifw(waifw = t(waifw), R_0 = R0, state = state.5.group)

S.prefect.data <- list("Obs_cases_age" = cases.by.age.group.5.group, "Cases_week_obs" = rowSums(cases.by.age.group.5.group), "N" = length(cases.by.age.group.5.group[, 1]), "M" = 5, "state" = state.5.group, "R0" = R0, "R0_sd" = R0_sd, "waifw" = waifw.input, "sus_prior" = sus_prior)

fit <- stan(model_code = Guinea_model_tsir_fit_observation_rate, data = S.prefect.data, iter = iters , chains = chains) 



list[N, M, cases.by.age.group, state, N.5.group, M.5.group, cases.by.age.group.5.group.Foum, state.5.group.Foumbadou] = output.data.for.given.sub.prefecture("FOUMBADOU", 19438)
list[N, M, cases.by.age.group, state,N.5.group, M.5.group, cases.by.age.group.5.group.Gueasso, state.5.group.Gueasso] = output.data.for.given.sub.prefecture("GUEASSO", 21116)
list[N, M, cases.by.age.group, state, N.5.group, M.5.group, cases.by.age.group.5.group.Gama, state.5.group.Gama.Berema] = output.data.for.given.sub.prefecture("GAMA BEREMA", 20465)
list[N, M, cases.by.age.group, state, N.5.group, M.5.group, cases.by.age.group.5.group.Kokota, state.5.group.Kokota] = output.data.for.given.sub.prefecture("KOKOTA", 14732)
list[ N, M, cases.by.age.group, state, N.5.group, M.5.group, cases.by.age.group.5.group.C.Lola, state.5.group.Cu.Lola] = output.data.for.given.sub.prefecture("CU LOLA", 49933)
list[N, M, cases.by.age.group, state,  N.5.group, M.5.group, cases.by.age.group.5.group.Laine, state.5.group.Laine] = output.data.for.given.sub.prefecture("LAINE", 16591)


reporting.rate.dist = c(Foum.ests$delta, Gama.Berema.ests$delta, Gueasso.ests$delta)
s1.dist = c(Foum.ests$s[, 1], Gama.Berema.ests$s[, 1], Gueasso.ests$s[, 1])
s2.dist = c(Foum.ests$s[, 2], Gama.Berema.ests$s[, 2], Gueasso.ests$s[, 2])
s3.dist = c(Foum.ests$s[, 3], Gama.Berema.ests$s[, 3], Gueasso.ests$s[, 3])
s4.dist = c(Foum.ests$s[, 4], Gama.Berema.ests$s[, 4], Gueasso.ests$s[, 4])
s5.dist = c(Foum.ests$s[, 5], Gama.Berema.ests$s[, 5], Gueasso.ests$s[, 5])
R0.dist = c(Foum.ests$R02, Gama.Berema.ests$R02, Gueasso.ests$R02)

extract.95.percentile.cases <- function(sus.dist, num.draws, obs.cases){
  cases.by.age.group = matrix(0, dim(sus.dist$Cases_by_age)[2], 5) 
  total.cases = matrix(0, num.draws, 1)
  reps = matrix(0, num.draws, 1)
  for (i in 1 : num.draws){
    for(j in 1 : dim(sus.dist$Cases_by_age)[2]){
      for(k in 1 : 5){
        cases.by.age.group[j, k] =  sample(sus.dist$Cases_by_age[, j, k], 1)  
      }
      total.cases[i] = sum(cases.by.age.group)
    }
  }
  # for(i in 1 : num.draws){
  #    reps[i] = sample(sus.dist$delta, 1)
  #    total.cases[i] = obs.cases/reps[i]
  #  }
  A = total.cases[order(total.cases)][round(0.025*num.draws):round(0.975*num.draws)]
  percentiles.95 = A[c(1,round(length(A)/2),length(A))]
  return(percentiles.95)
}


extract.95.percentile.cases.non.converged <- function(reporting.rate.dist, num.draws, obs.cases){
  cases.by.age.group = matrix(0, length(obs.cases), 5) 
  total.cases = matrix(0, num.draws, 1)
  reps = matrix(0, num.draws, 1)
  
  for(i in 1 : num.draws){
    reps[i] = sample(reporting.rate.dist, 1)
    total.cases[i] = obs.cases/reps[i]
  }
  A = total.cases[order(total.cases)][round(0.025*num.draws):round(0.975*num.draws)]
  percentiles.95 = A[c(1,round(length(A)/2),length(A))]
  return(percentiles.95)
}


Foum.95 = extract.95.percentile.cases(Foum.ests, 10000, sum(cases.by.age.group.5.group.Foum))
Gama.95 = extract.95.percentile.cases(Gama.Berema.ests, 10000,sum(cases.by.age.group.5.group.Gama))
Gueasso.95 = extract.95.percentile.cases(Gueasso.ests, 10000, sum(cases.by.age.group.5.group.Gueasso))
C.Lola.95 = extract.95.percentile.cases.non.converged(reporting.rate.dist, 1000, sum(cases.by.age.group.5.group.C.Lola))
Laine.95 = extract.95.percentile.cases.non.converged(reporting.rate.dist, 1000, sum(cases.by.age.group.5.group.Laine))
Kokota.95 = extract.95.percentile.cases.non.converged(reporting.rate.dist, 1000, sum(cases.by.age.group.5.group.Kokota))


Cases = matrix(0,6,5)
Cases[, 1] = as.character(c("Foumbadou","Kokota", "Gama Berema", "Gueasso", "Laine", "C. Lola"))
Cases[, 2] = as.numeric(c(27,11, 38, 35, 22, 4))
Cases[, 3] = round(as.numeric(c(Foum.95[1], Kokota.95[1], Gama.95[1], Gueasso.95[1], Laine.95[1], C.Lola.95[1])))
Cases[, 4] = round(as.numeric(c(Foum.95[2], Kokota.95[2], Gama.95[2], Gueasso.95[2], Laine.95[2], C.Lola.95[2])))
Cases[, 5] = round(as.numeric(c(Foum.95[3], Kokota.95[3], Gama.95[3], Gueasso.95[3], Laine.95[3], C.Lola.95[3])))
Cases = data.frame(Cases)
colnames(Cases) = c("Sub Prefecture", "Reported cases", "2.5%","50%", "97.5%")
View(Cases)






#######################
### Function runs forward forecasts given a set of estimated susceptibility profiles by age group produced by stan along with R0. Similar to previous function, but this is for when we have 5 groups.
### Inputs are R0; time.length - the number of weeks we project forward for; sus.dist - which is the estimated susceptiblity by age from stan; num.sims - number of forward projections performed; cases.by.age.group - the cases by age group that have been observed in the outbreak; waifw.init - the initial waifw that was used for the stan estimation, and will be used for the forward projections; state - which is the population by age group.
#######################

simulate.with.reporting.rate.adjusted.cases <- function(sus.dist, time.length, num.sims,
                                                        waifw.init, state){ 
  
  
  ### Initialize the variables we want to output from the function. These are, for each simulation, the attack rate over each age group, the initial number of susceptibles in each age group, and the number of infections in each age group.
  attack.rate = matrix(0, num.sims, 5)
  initial.sus = matrix(0, num.sims, 5)
  end.sus = matrix(0, num.sims, 5)
  Infections = array(0, c(length(state), time.length, num.sims))
  
  ### Perform a loop for each of the simulations
  for(i in 1 : num.sims){
    
    ### Select a beginning susceptibility for each age group from the stan output
    start.sus.proportion = c(sample(sus.dist$s[, 1], 1, replace = T), sample(sus.dist$s[, 2], 1, replace = T),
                             sample(sus.dist$s[, 3], 1, replace = T), sample(sus.dist$s[, 4], 1, replace = T),
                             sample(sus.dist$s[, 5], 1, replace = T))
    
    cases.by.age.group = matrix(0, dim(sus.dist$Cases_by_age)[2], 5)
    for(qq in 1 : dim(sus.dist$Cases_by_age)[2]){
      for(jj in 1 : 5){
        cases.by.age.group[qq, jj] =  sample(sus.dist$Cases_by_age[, qq, jj], 1)  
      }
    }
    
    ### Select a value of R_0 to use in the simulation
    R0 = sample(sus.dist$R02, 1)
    
    ### This susceptibility is for before the epidemic began, therefore we subtract any susceptibles who were infected in the outbreak up to the end of the data.
    sus = round(start.sus.proportion * state - colSums(cases.by.age.group))
    
    ### Store these numbers of susceptibles in the initial.sus matrix, which will be output from the function
    initial.sus[i, ] = sus
    
    ### Check that none of the numbers for initial susceptibles is <0.
    for(i1 in 1:length(sus)){
      sus[i1] = max(0, sus[i1])
    }
    
    ### Average the last two weeks of cases by age group to give the initial number of susceptibles by age group. This will be used to calculate the force of infection on each age group through the waifw. 
    ### The averaging is done, as if we continually use the number of cases from two weeks prior to the present week for prediction, then we will get a pattern which is too dependent on the number of cases, for example if we had 10 one week and 30 the next, then projecting this forwards, we would see a pattern of low then high incidence in the predictions.
    ### Another option is to look at bi-weeks, though I haven't looked at this yet.
    I.t = colSums(tail(cases.by.age.group, 2))/2
    
    ### We now loop over the number of weeks that we have chosen to project forwards for.
    for(j in 1 : time.length){
      
      ### Calculate the waifw, which takes into account the seasonality of the outbreak and also calculate the force of infection by age group, which is dependent on the waifw and the number of infecteds in each age group.
      waifw = output.waifw(waifw.init, R0 * (1 + cos((3/12 + j/48)  * 2 *  pi) * 0.3), state)
      phi = calc.phi.1.week(waifw, as.numeric(I.t), denom)
      
      ### Loop over each of the age group to generate the forward projections of number of cases by age group
      for(l in 1 : length(sus)){
        ### The number of infections by age group is assumed to be binomially distributed Bin(n,p), with n = number of susceptibles in the age group; p = force of infection for this age group
        Infections[l, j, i] = rbinom(1, sus[l], phi[l])
        
        ### Update the number of susceptibles by subtracting the number of infections from the number of susceptibles in the age group
        sus[l] = sus[l] - Infections[l, j, i]
      }
      ### Update the number of infections for the next time step
      if (j == 1) {I.t = (tail(cases.by.age.group, 1) + (Infections[, j , i]))/2}
      if (j > 1)  { I.t = ((Infections[, j - 1, i]) + (Infections[, j, i]))/ 2}
      #I.t = Infections[, j, i]
    }
    end.sus[i, ] = sus
    for(pp in 1 : 5){
      ### Calculate the attack rate for each of the age groups
      attack.rate[i, pp] = (sum(Infections[pp, , i]))/state[pp]
    }
  }
  
  ### Output the projections along with the attack rate and initial susceptibles for each simulation. 
  return(list(Infections, attack.rate, initial.sus))
}


output.plots.with.estimates.for.observed.data.sample.for.observed.cases <- function(Infections.array, time.length, r, g, b, S.Prefect, sus.dist, obs.cases.by.week, num.sims){
  
  cases.by.age.group = matrix(0, length(obs.cases.by.week), 5)
  cases.by.week = matrix(0, num.sims, length(obs.cases.by.week))
  
  for(i in 1 : num.sims){
    for(qq in 1 : length(obs.cases.by.week)){
      for(jj in 1 : 5){
        cases.by.age.group[qq, jj] =  sample(sus.dist$Cases_by_age[, qq, jj], 1)  
      }
      cases.by.week[i, qq] = sum(cases.by.age.group[qq,])
    }
  }
  
  
  mean.cases.by.week = matrix(0, length(obs.cases.by.week), 1)
  
  
  mean.cumulative.Infections.array = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){
    mean.cumulative.Infections.array[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  
  
  
  
  plot(seq(14 - length(obs.cases.by.week), 14 + time.length - 1), c(cases.by.week[1, ], colSums( Infections.array[, , 1])),  col = rgb(red = r, green = g, blue = b, alpha = 0.1), type = "l", 
       ylim = c(0,max(c(cases.by.week,colSums(Infections.array)))), xaxt = "n",  xlab = "Week", ylab = "Predicted infections")
  axis(side = 1, at = seq(14 - length(obs.cases.by.week), 14 + time.length - 1))
  
  # plot(seq(14, 14 + time.length - 1),colSums( Infections.array[, , 1]),  col = rgb(red = r, green = g, blue = b, alpha = 0.1), type = "l", 
  #     ylim = c(0,max(colSums(Infections.array))), xaxt = "n",  xlab = "Week", ylab = "Predicted new cases")
  #axis(side = 1, at = seq(14, (14 + time.length - 1)))
  for(i in 1 : 1000){
    lines(seq(14 - length(obs.cases.by.week), 14 + time.length - 1), c(cases.by.week[i, ], colSums(Infections.array[, , i])), col = rgb(red = r, green = g, blue = b, alpha = 0.1))
  }
  means = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means[i] = mean(colSums(Infections.array[, i, ]))
  }
  lines(seq(14, 14 + time.length - 1), means, lwd = 3, col = rgb(red = r-0.05, green = g-0.05, blue = b-0.05))  
  legend("topright", legend = S.Prefect, 
         col = rgb(red = r, green = g, blue = b),
         lwd = 2, bty="n", cex = 1)
  
}


list[Infections.5.groups.Foum, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Foum.ests, time.length, num.sims= 1000, waifw.init = t(waifw.5.groups), state.5.group.Foumbadou)


list[Infections.5.groups.Gama, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Gama.Berema.ests, time.length, num.sims= 1000, waifw.init= t(waifw.5.groups), state.5.group.Gama.Berema)


list[Infections.5.groups.Gueasso, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Gueasso.ests, time.length, num.sims= 1000,  waifw.init= t(waifw.5.groups), state.5.group.Gueasso)






########################
### Simulate forward forecasts for the sub prefectures that haven't converged. We use the concatenated reporting rate distribution from the converged sub prefectures and concatenated susceptible proportions by age group from the converged sub prefectures
########################

simulate.with.reporting.rate.non.converged <- function(s1.dist, s2.dist, s3.dist, s4.dist, s5.dist, 
                                                      reporting.rate.dist, R0.dist,
                                                      time.length, num.sims, obs.cases,
                                                      waifw.init, state){ 
  
  
  ### Initialize the variables we want to output from the function. These are, for each simulation, the attack rate over each age group, the initial number of susceptibles in each age group, and the number of infections in each age group.
  attack.rate = matrix(0, num.sims, 5)
  initial.sus = matrix(0, num.sims, 5)
  end.sus = matrix(0, num.sims, 5)
  Infections = array(0, c(length(state), time.length, num.sims))
  
  ### Perform a loop for each of the simulations
  for(i in 1 : num.sims){
    
    ### Select a beginning susceptibility for each age group from the stan output
    start.sus.proportion = c(sample(s1.dist, 1, replace = T), sample(s2.dist, 1, replace = T),
                             sample(s3.dist, 1, replace = T), sample(s4.dist, 1, replace = T),
                             sample(s5.dist, 1, replace = T))
    
    cases.by.age.group = matrix(0, length(obs.cases[, 1]), 5)
    for(qq in 1 : length(obs.cases[, 1])){
      for(jj in 1 : 5){
        cases.by.age.group[qq, ] =  obs.cases[qq,]/ sample(reporting.rate.dist, 1)
      }
    }
    
    ### Select a value of R_0 to use in the simulation
    R0 = sample(R0.dist, 1)
    
    ### This susceptibility is for before the epidemic began, therefore we subtract any susceptibles who were infected in the outbreak up to the end of the data.
    sus = round(start.sus.proportion * state - colSums(cases.by.age.group))
    
    ### Store these numbers of susceptibles in the initial.sus matrix, which will be output from the function
    initial.sus[i, ] = sus
    
    ### Check that none of the numbers for initial susceptibles is <0.
    for(i1 in 1:length(sus)){
      sus[i1] = max(0, sus[i1])
    }
    
    ### Average the last two weeks of cases by age group to give the initial number of susceptibles by age group. This will be used to calculate the force of infection on each age group through the waifw. 
    ### The averaging is done, as if we continually use the number of cases from two weeks prior to the present week for prediction, then we will get a pattern which is too dependent on the number of cases, for example if we had 10 one week and 30 the next, then projecting this forwards, we would see a pattern of low then high incidence in the predictions.
    ### Another option is to look at bi-weeks, though I haven't looked at this yet.
    I.t = colSums(tail(cases.by.age.group, 2))/2
    
    ### We now loop over the number of weeks that we have chosen to project forwards for.
    for(j in 1 : time.length){
      
      ### Calculate the waifw, which takes into account the seasonality of the outbreak and also calculate the force of infection by age group, which is dependent on the waifw and the number of infecteds in each age group.
      waifw = output.waifw(waifw.init, R0 * (1 + cos((3/12 + j/48)  * 2 *  pi) * 0.3), state)
      phi = calc.phi.1.week(waifw, as.numeric(I.t), denom)
      
      ### Loop over each of the age group to generate the forward projections of number of cases by age group
      for(l in 1 : length(sus)){
        ### The number of infections by age group is assumed to be binomially distributed Bin(n,p), with n = number of susceptibles in the age group; p = force of infection for this age group
        Infections[l, j, i] = rbinom(1, sus[l], phi[l])
        
        ### Update the number of susceptibles by subtracting the number of infections from the number of susceptibles in the age group
        sus[l] = sus[l] - Infections[l, j, i]
      }
      ### Update the number of infections for the next time step
      if (j == 1) {I.t = (tail(cases.by.age.group, 1) + (Infections[, j , i]))/2}
      if (j > 1)  { I.t = ((Infections[, j - 1, i]) + (Infections[, j, i]))/ 2}
      #I.t = Infections[, j, i]
    }
    end.sus[i, ] = sus
    for(pp in 1 : 5){
      ### Calculate the attack rate for each of the age groups
      attack.rate[i, pp] = (sum(Infections[pp, , i]))/state[pp]
    }
  }
  
  ### Output the projections along with the attack rate and initial susceptibles for each simulation. 
  return(list(Infections, attack.rate, initial.sus))
}


output.plots.with.estimates.for.observed.data.non.convergence <- function(Infections.array, time.length, r, g, b, S.Prefect, reporting.rate.dist, obs.cases.by.week, num.sims){
  
  #cases.by.age.group = matrix(0, dim(sus.dist$Cases_by_age)[2], 5)
  cases.by.week = matrix(0, length(obs.cases.by.week), num.sims)
  reps = matrix(0,  num.sims, 1)
  for (i in 1 : num.sims){
    reps[i] = sample(reporting.rate.dist, 1)
  }
  cases.by.week = matrix(0, num.sims, length(obs.cases.by.week))
  for(i in 1 : num.sims){
    cases.by.week[i,] = obs.cases.by.week/reps[i]
  }
  
  mean.cases.by.week = matrix(0, length(obs.cases.by.week), 1)
  
  
  mean.cumulative.Infections.array = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){
    mean.cumulative.Infections.array[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  
  
  
  
  plot(seq(14 - length(obs.cases.by.week), 14 + time.length - 1), c(cases.by.week[1, ], colSums( Infections.array[, , 1])),  col = rgb(red = r, green = g, blue = b, alpha = 0.1), type = "l", 
       ylim = c(0,max(colSums(Infections.array))), xaxt = "n",  xlab = "Week", ylab = "Predicted infections")
  axis(side = 1, at = seq(14 - length(obs.cases.by.week), 14 + time.length - 1))
  
  # plot(seq(14, 14 + time.length - 1),colSums( Infections.array[, , 1]),  col = rgb(red = r, green = g, blue = b, alpha = 0.1), type = "l", 
  #     ylim = c(0,max(colSums(Infections.array))), xaxt = "n",  xlab = "Week", ylab = "Predicted new cases")
  #axis(side = 1, at = seq(14, (14 + time.length - 1)))
  for(i in 1 : 1000){
    lines(seq(14 - length(obs.cases.by.week), 14 + time.length - 1), c(cases.by.week[i, ], colSums(Infections.array[, , i])), col = rgb(red = r, green = g, blue = b, alpha = 0.1))
  }
  means = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means[i] = mean(colSums(Infections.array[, i, ]))
  }
  lines(seq(14, 14 + time.length - 1), means, lwd = 3, col = rgb(red = r-0.05, green = g-0.05, blue = b-0.05))  
  legend("topright", legend = S.Prefect, 
         col = rgb(red = r, green = g, blue = b),
         lwd = 2, bty="n", cex = 1)
  
}




list[Infections.5.groups.Kokota, attack.rate.5.groups, initial.5.groups]  = 
  simulate.with.reporting.rate.non.converged(s1.dist, s2.dist, s3.dist, s4.dist, s5.dist,
                                            reporting.rate.dist, R0.dist,
                                            time.length, num.sims, cases.by.age.group.5.group.Kokota,
                                            t(waifw.5.groups), state = state.5.group.Kokota)



list[Infections.5.groups.Laine, attack.rate.5.groups, initial.5.groups]  = 
  simulate.with.reporting.rate.non.converged(s1.dist, s2.dist, s3.dist, s4.dist, s5.dist,
                                            reporting.rate.dist, R0.dist,
                                            time.length, num.sims, cases.by.age.group.5.group.Laine,
                                            t(waifw.5.groups), state = state.5.group.Laine)


list[Infections.5.groups.C.Lola, attack.rate.5.groups, initial.5.groups]  = 
  simulate.with.reporting.rate.non.converged(s1.dist, s2.dist, s3.dist, s4.dist, s5.dist,
                                            reporting.rate.dist, R0.dist,
                                            time.length, num.sims, cases.by.age.group.5.group.C.Lola,
                                            t(waifw.5.groups), state = state.5.group.Cu.Lola)



par(mfrow=c(3, 2))
r = 0.15
g = 0.15
b = 0.7


output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.Foum, time.length, r, g, b, "Foumbadou", Foum.ests,obs.cases.by.week = rowSums(cases.by.age.group.5.group.Foum), num.sims, first.week = 18)


output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.Gama, time.length, r, g, b, "Gama Berema", Gama.Berema.ests,obs.cases.by.week = rowSums(cases.by.age.group.5.group.Gama), num.sims)

output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.Gueasso, time.length, r, g, b, "Gueasso", Gueasso.ests,obs.cases.by.week = rowSums(cases.by.age.group.5.group.Gueasso), num.sims)


output.plots.with.estimates.for.observed.data.non.convergence(Infections.array = Infections.5.groups.Laine, time.length, r, g, b, S.Prefect = "Laine", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Laine), num.sims)

output.plots.with.estimates.for.observed.data.non.convergence(Infections.array = Infections.5.groups.Kokota, time.length, r, g, b, S.Prefect = "Kokota", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Kokota), num.sims)

output.plots.with.estimates.for.observed.data.non.convergence(Infections.array = Infections.5.groups.C.Lola, time.length, r, g, b, S.Prefect = "Central Lola", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.C.Lola), num.sims)
