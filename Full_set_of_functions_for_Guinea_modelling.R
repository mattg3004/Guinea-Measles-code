library(rstan)
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



estimate.susceptibility.5.group.N.Zoo <- function(waifw.init, x, y, R0, state, iters, R0_sd){
  
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
  fit <- stan(model_code = Guinea_model_tsir_with_susceptible_prior, data = Guinea.data, iter = iters , chains = 2) 
  
  ### print the output of the stan model for informational purposes
  print(fit)
  
  ### extract the output of the stan model, which is also the output of this function.
  sus.est = extract(fit)
  
  return(sus.est)
}


estimates.5.groups = estimate.susceptibility.5.group.N.Zoo(t(waifw.5.groups), x.5.group, y.5.group, R0 = 18, state.5.group, iters, R0_sd = 4)
