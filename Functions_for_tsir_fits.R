source("Data_NZoo.R")
source("Functions.R")
library(rstan)


###############################
# Make it so that it is possible to receive the output of a function which is in a list, and automatically split it into components using list[..]
###############################
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

###############
### Initially, we fit the proportion of susceptibles if we split the population into 4 age groups, <1's, 1-4, 5-14, 15+.
### Given an R_0, we fit the susceptilbe proportion of each of these groups by assuming that the epidemic is fully observed in N'Zoo.
###############

Guinea_model_tsir <- '
data {
  int<lower=0> N;       // week sample size
  int<lower=0> M;       // number of age groups
  int y[N, M];          // infection data by age group
  real phi[N, M];       // force of infection on each age group by week, dependent on WAIFW and R_0 - Klepac paper method
  real state[M];        // population size of each age group
  int cum_cases[N, M];  // cumulative cases by week and age group
}

parameters {
// fit the susceptible proportion by age given a value of R_0
// s1 is proportion of susceptibles in group 1 (Under 1s), s2 for group 2 etc.

  real <upper = 1> s1;    // s1 - susceptibility in Under 1s
  real <upper = 1> s2;    // s2 - susceptibility in 1-4
  real <upper = 1> s3;    // s3 - susceptibility in 5-14
  real <upper = 1> s4;    // s4 - susceptibility in 15+
}
transformed parameters {
}
model {

// We assume that the number of cases at time t is dependent on number of cases from two weeks previous (this is contained in phi) and
// the number of susceptibles in the age group. Assuming that the epidemic is fully observed (as is thought to be true for NZoo),
// the likelihood of observing the cases we did is assumed to be given by a poisson distribution.
// s1 * state[1] gives the initial number of susceptibles in group 1 and subtracting cum_cases from this gives the number of susceptibles at time t.
// Same applies for other groups. Looping over m is for the different age groups, over n gives the different weeks in the dataset.

  
  for (m in 1 : M){
    for (n in 1 : N){
      if(m == 1){
        y[n, m] ~ poisson((s1 * state[1] - cum_cases[n, m])  * phi[n, m]);  
      }
      if(m == 2){
        y[n, m] ~ poisson((s2 * state[2] - cum_cases[n, m])  * phi[n, m]);  
      }
      if(m == 3){
        y[n, m] ~ poisson((s3 * state[3] - cum_cases[n, m] ) * phi[n, m]);  
      }
      if(m == 4){
        y[n, m] ~ poisson((s4 * state[4] - cum_cases[n, m]) * phi[n, m]);  
      }
    }
  }
}
'



#######################
### We split the population into 5 age groups here, <1's, 1-4, 5-10, 11-15, 16+, and use the same method
### as for the 4 age groups to fit susceptible proportions in each group.
#######################

Guinea_model_tsir_5_groups <- '
data {
  int<lower=0> N;       // week sample size
  int<lower=0> M;       // number of age groups
  int y[N, M];          // infection data by age group
  real phi[N, M];       // force of infection on each age group by week, dependent on WAIFW and R_0 - Klepac paper method
  real state[M];        // population size of each age group
  int cum_cases[N, M];  // cumulative cases by week and age group
}

parameters {
// fit the susceptible proportion by age given a value of R_0
// s1 is proportion of susceptibles in group 1 (Under 1s), s2 for group 2 etc.
  real <upper = 1> s1;    // s1 - susceptibility in Under 1s
  real <upper = 1> s2;    // s2 - susceptibility in 1-4
  real <upper = 1> s3;    // s3 - susceptibility in 5-10
  real <upper = 1> s4;    // s4 - susceptibility in 11-15
  real <upper = 1> s5;    // s5 - susceptibility in 16+
}
transformed parameters {
}
model {

// We assume that the number of cases at time t is dependent on number of cases from two weeks previous (this is contained in phi) and
// the number of susceptibles in the age group. Assuming that the epidemic is fully observed (as is thought to be true for NZoo),
// the likelihood of observing the cases we did is assumed to be given by a poisson distribution.
// s1 * state[1] gives the initial number of susceptibles in group 1 and subtracting cum_cases from this gives the number of susceptibles at time t.
// Same applies for other groups. Looping over m is for the different age groups, over n gives the different weeks in the dataset.

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

#######################
### This function takes in a WAIFW and returns a WAIFW that corresponds to a specified R_0
#######################

output.waifw<- function(waifw, R_0, state){ 
  denom <- sum(state)
  next.gen <- as.numeric(state)*(1-exp(-waifw/denom))
  
  #get the first eigen value
  cur.R0 <- Re(eigen(next.gen)$value[1])
  
  #More correct transform
  R.ratio <- R_0/cur.R0 #print(R0); #print(cur.R0); #print(R.ratio)
  waifw <- -log(1-R.ratio*(1-exp(-waifw/denom)))*denom
  return(waifw)
}



#######################
### This function calculates the force of infection on each age group - method as described in Petra Klepac paper.
### Input the waifw, the cases at each time and the size of the population.
#######################

calc.phi <- function(waifw, x, denom){
  
  phi <- (waifw %*% t(x[, ] ) ^ 0.97) /denom
  
  phi <- 1 - exp(-phi)
  
  phi <- t(phi)
  return(phi)
}





#######################
### This function does necessary calculations to input into rstan code for susceptibility estimation.
### The output is the estimates of the susceptibility in each of the 4 groups.
### Inputs are desired R_0, the initial waifw, which will get scaled according to R_0.
### x and y are the cases observed, with the cases in x being two weeks previous to the cases in y, which is used in the stan model to calculate likelihoods for susceptibility
### state is the population of each age group
#######################

estimate.susceptibility <- function(waifw.init, x, y, R_0, state){
  
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
  Guinea.data <- list("y" = y, "N" = N, "M" = M, "phi" = phi, "state" = state, "cum_cases" = cum_cases)
  
  ### Here the stan model is run to estimate the susceptibility in the age groups.
  fit <- stan(model_code = Guinea_model_tsir, data = Guinea.data, iter = 1000 , chains = 2)
  
  ### print the output of the stan model for informational purposes
  print(fit)
  
  ### extract the output of the stan model, which is also the output of this function.
  sus.est = extract(fit)
  
  return(sus.est)
}




#######################
### This function does necessary calculations to input into rstan code for susceptibility estimation. This function is for the 5 groups split of the population.
### The output is the estimates of the susceptibility in each of the 5 groups.
### Inputs are desired R_0, the initial waifw, which will get scaled according to R_0.
### x and y are the cases observed, with the cases in x being two weeks previous to the cases in y, which is used in the stan model to calculate likelihoods for susceptibility
### state is the population of each age group
#######################

estimate.susceptibility.5.group <- function(waifw.init, x, y, R_0, state, iters){
  
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
  Guinea.data <- list("y" = y, "N" = N, "M" = M, "phi" = phi, "state" = state, "cum_cases" = cum_cases)
  
  ### Here the stan model is run to estimate the susceptibility in the 5 age groups. 
  fit <- stan(model_code = Guinea_model_tsir_5_groups, data = Guinea.data, iter = iters , chains = 2) 
  
  ### print the output of the stan model for informational purposes
  print(fit)
  
  ### extract the output of the stan model, which is also the output of this function.
  sus.est = extract(fit)
  
  return(sus.est)
}


#######################
### This function plots boxplots of the susceptibility estimates for the 4 age groups for 3 different values of R_0.
### Inputs are the fits generated by the stan code for each of the 3 different R_0s, along with the values of these R_0s
### This was initially developed for comparing the fits when we used different R_0s, at some point it was decided that only considering high values of R_0s was consistent with what we expected in terms of susceptibility, so this function is semi-defunct.
#######################

plot.susceptibility.estimates <- function(A, B, C, R0.1, R0.2, R0.3){
  par(mfrow = c(3,1))
  par(mar=c(1,1,1,1))
  par(oma = c(3.5,4,2,2))
  boxplot(A$s1, A$s2, A$s3, A$s4,  
          xaxt = "n", cex.axis = 1.2)
  mtext(side = 2, text = "Susceptibility", line = 1.6, outer= T)
  mtext(side = 1, text = "Age", line = 1.6, outer= T)
  axis(side = 1, at = c(1,2,3,4), labels = c("","","",""))
  legend("topright", legend = paste("R0 =", R0.1), border = F)
  boxplot(B$s1, B$s2, B$s3, B$s4, 
          xaxt = "n", cex.axis = 1.2)
  axis(side = 1, at = c(1,2,3,4), labels = c("","","",""))
  legend("topright", legend = paste("R0 =", R0.2))
  boxplot(C$s1, C$s2, C$s3, C$s4, 
          xaxt = "n", cex.axis = 1.2)
  axis(side = 1, at = c(1,2,3,4), labels = c("<1", "1-4", "5-14", "15+"), cex.axis = 1.5)
  legend("topright", legend = paste("R0 =", R0.3))
  
}


#######################
### This function plots boxplots of the susceptibility estimates for the 5 age groups for 3 different values of R_0.
### Inputs are the fits generated by the stan code for each of the 3 different R_0s, along with the values of these R_0s
### This was initially developed for comparing the fits when we used different R_0s, at some point it was decided that only considering high values of R_0s was consistent with what we expected in terms of susceptibility, so this function is semi-defunct.
#######################

plot.susceptibility.estimates.5.group <- function(A, B, C, R0.1, R0.2, R0.3){
  par(mfrow = c(3,1))
  par(mar=c(1,1,1,1))
  par(oma = c(3.5,4,2,2))
  boxplot(A$s1, A$s2, A$s3, A$s4, A$s5,  
          xaxt = "n", cex.axis = 1.2)
  mtext(side = 2, text = "Susceptibility", line = 1.6, outer= T)
  mtext(side = 1, text = "Age", line = 1.6, outer= T)
  axis(side = 1, at = c(1,2,3,4,5), labels = c("","","","",""))
  legend("topright", legend = paste("R0 =", R0.1), border = F)
  boxplot(B$s1, B$s2, B$s3, B$s4, B$s5,
          xaxt = "n", cex.axis = 1.2)
  axis(side = 1, at = c(1,2,3,4,5), labels = c("","","","",""))
  legend("topright", legend = paste("R0 =", R0.2))
  boxplot(C$s1, C$s2, C$s3, C$s4, C$s5,
          xaxt = "n", cex.axis = 1.2)
  axis(side = 1, at = c(1,2,3,4,5), labels = c("<1", "1-4", "5-10","11-15", "16+"), cex.axis = 1.5)
  legend("topright", legend = paste("R0 =", R0.3))
  
}





#######################
### Function runs forward forecasts given a set of estimated susceptibility profiles by age group produced by stan along with R0. 
### Inputs are R0; time.length - the number of weeks we project forward for; sus.dist - which is the estimated susceptiblity by age from stan; num.sims - number of forward projections performed; cases.by.age.group - the cases by age group that have been observed in the outbreak; waifw.init - the initial waifw that was used for the stan estimation, and will be used for the forward projections; state - which is the population by age group.
#######################

simulate.with.given.R0.and.sus.dist.waifw.groups <- function(R0, sus.dist, time.length, num.sims,
                                                             cases.by.age.group = cases.by.age.group,
                                                             waifw.init, state){ 
  
  ### Initialize the variables we want to output from the function. These are, for each simulation, the attack rate over each age group, the initial number of susceptibles in each age group, and the number of infections in each age group.
  attack.rate = matrix(0, num.sims, 4)
  initial.sus = matrix(0, num.sims, 4)
  end.sus = matrix(0, num.sims, 4)
  Infections = array(0, c(length(state), time.length, num.sims))
  
  ### Perform a loop for each of the simulations
  for(i in 1 : num.sims){
    
    ### Select a beginning susceptibility for each age group from the stan output
    start.sus.proportion = c(sample(sus.dist$s1, 1, replace = T), sample(sus.dist$s2, 1, replace = T),
                             sample(sus.dist$s3, 1, replace = T), sample(sus.dist$s4, 1, replace = T))
    
    ### This susceptibility is for before the epidemic began, therefore we subtract any susceptibles who were infected in the outbreak up to the end of the data.
    sus = round(start.sus.proportion * state) - colSums(cases.by.age.group)
    
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
    }
    
    end.sus[i, ] = sus
    for(pp in 1 : 4){
      ### Calculate the attack rate for each of the age groups
      attack.rate[i, pp] = (sum(Infections[pp, , i]))/state[pp]
    }
  }
  
  ### Output the projections along with the attack rate and initial susceptibles for each simulation. 
  return(list(Infections, attack.rate, initial.sus))
}




#######################
### Function runs forward forecasts given a set of estimated susceptibility profiles by age group produced by stan along with R0. Similar to previous function, but this is for when we have 5 groups.
### Inputs are R0; time.length - the number of weeks we project forward for; sus.dist - which is the estimated susceptiblity by age from stan; num.sims - number of forward projections performed; cases.by.age.group - the cases by age group that have been observed in the outbreak; waifw.init - the initial waifw that was used for the stan estimation, and will be used for the forward projections; state - which is the population by age group.
#######################

simulate.with.given.R0.and.sus.dist.waifw.groups.5.groups <- function(sus.dist, time.length, num.sims,
                                                                      cases.by.age.group = cases.by.age.group,
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
    
    ### Select a value of R_0 to use in the simulation
    R0 = sample(sus.dist$R02, 1)
    alpha = sample(sus.dist$alpha, 1)
    ### This susceptibility is for before the epidemic began, therefore we subtract any susceptibles who were infected in the outbreak up to the end of the data.
    sus = round(start.sus.proportion * state) - colSums(cases.by.age.group)
    
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
      waifw = output.waifw(waifw.init, R0 * (1 + cos((4/12 + j/48)  * 2 *  pi) * alpha), state)
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




#######################
### This function plots boxplots of the susceptibility estimates for the 4 age groups for a given value of R_0.
### Inputs are the fits generated by the stan code for the value of R_0, along with the values of R_0
#######################

plot.susceptibility.estimates.4.group.one.R0 <- function(A, R0){
  par(mfrow = c(1,1))
  par(mar=c(1,1,1,1))
  par(oma = c(3.5,4,2,2))
  boxplot(A$s1, A$s2, A$s3, A$s4,  
          xaxt = "n", cex.axis = 1.2)
  mtext(side = 2, text = "Susceptibility", line = 1.6, outer= T)
  mtext(side = 1, text = "Age", line = 1.6, outer= T)
  legend("topright", legend = paste("R0 =", R0), border = F)
  axis(side = 1, at = c(1,2,3,4), labels = c("<1", "1-4", "5-14", "15+"), cex.axis = 1.5)
}




#######################
### This function plots boxplots of the susceptibility estimates for the 5 age groups for a given value of R_0.
### Inputs are the fits generated by the stan code for the value of R_0, along with the values of R_0
#######################

plot.susceptibility.estimates.5.group.one.R0 <- function(A, S.Prefect){
  par(mfrow = c(1,1))
  par(mar=c(1,1,1,1))
  par(oma = c(3.5,4,2,2))
  boxplot(A$s[, 1], A$s[, 2], A$s[, 3], A$s[, 4], A$s[, 5],
          xaxt = "n", cex.axis = 1.2)
  mtext(side = 2, text = "Susceptibility", line = 1.6, outer= T)
  mtext(side = 1, text = "Age", line = 1.6, outer= T)
  legend("topright", legend = S.Prefect, border = F)
  axis(side = 1, at = c(1,2,3,4,5), labels = c("<1", "1-5", "6-10", "11-15", "16+"), cex.axis = 1.5)
}






#######################
### This function plots boxplots of the susceptibility estimates for the 5 age groups for a given value of R_0.
### Inputs are the fits generated by the stan code for the value of R_0, along with the values of R_0
#######################

plot.reporting.rate.estimates <- function(A, S.prefect){
  par(mfrow = c(1,1))
  par(mar=c(1,1,1,1))
  par(oma = c(3.5,4,2,2))
  boxplot(A$delta,
          xaxt = "n", cex.axis = 1.2, ylim = c(0, 1))
  #mtext(side = 2, text = "Susceptibility", line = 1.6, outer= T)
  mtext(side = 1, text = "Reporting rate", line = 1.6, outer= T)
  legend("topright", legend = paste(S.prefect), border = F)
}




#######################
### This function plots the forward projections for two different simulations
#######################

combine.2.plots <- function(Infections.array1, Infections.array2,
                          R0,
                          r1, g1, b1, r2, g2, b2,
                          time.length){
  
  maximum = max(max(colSums(Infections.array1)), max(colSums(Infections.array2)))
  Infections.array = Infections.array1
  
  
  plot(seq(14, 14 + time.length - 1),colSums( Infections.array[, , 1]),  col = rgb(red = r1, green = g1, blue = b1, alpha = 0.1), type = "l", 
       ylim = c(0, maximum), xaxt = "n",  xlab = "Week", ylab = "Predicted new cases")
  axis(side = 1, at = seq(14, (14 + time.length - 1)))
  for(i in 1 : 1000){
    lines(seq(14, 14 + time.length - 1), colSums(Infections.array[, , i]), col = rgb(red = r1, green = g1, blue = b1, alpha = 0.1))
  }
  means1 = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means1[i] = mean(colSums(Infections.array[, i, ]))
  }
  
  
  
  Infections.array = Infections.array2
  
  for(i in 1 : 1000){
    lines(seq(14, 14 + time.length - 1), colSums(Infections.array[, , i]), col = rgb(red = r2, green = g2, blue = b2, alpha = 0.1))
  }
  means2 = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means2[i] = mean(colSums(Infections.array[, i, ]))
  }
  
  
  
  lines(seq(14, 14 + time.length - 1), means1, lwd = 3, col = rgb(red = r1, green = g1, blue = b1))  
  lines(seq(14, 14 + time.length - 1), means2, lwd = 3, col = rgb(red = r2, green = g2, blue = b2))  
  
  
  legend("topright", legend = c(paste("R0 = ", R0,", 4 groups", sep = ""), paste("R0 = ", R0, ", 5 groups", sep = "")),  
         col = c(rgb(red = r1, green = g1, blue = b1), rgb(red = r2, green = g2, blue = b2)),
         lwd = 2, bty="n", cex = 1)
}





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
    start.sus.proportion = c(sample(sus.dist$s1, 1, replace = T), sample(sus.dist$s2, 1, replace = T),
                             sample(sus.dist$s3, 1, replace = T), sample(sus.dist$s4, 1, replace = T),
                             sample(sus.dist$s5, 1, replace = T))
    
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





simulate.with.reporting.rate.no.converged <- function(s1.dist, s2.dist, s3.dist, s4.dist, s5.dist, 
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


