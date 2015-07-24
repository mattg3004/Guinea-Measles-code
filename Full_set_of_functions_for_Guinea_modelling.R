library(rstan)
source("Functions_for_tsir_fits.R")

############################################################################################################################################################
### Do fits for N'Zoo district, as we assume that there is 100% observation in this sub-prefecture
############################################################################################################################################################
source("NZoo_fits.R")

waifw.5.groups =  rbind(c(6.9 , 3 , 3/2 , 2/3,  0.5 ) * state.5.group[1],
                        c(3 , 6.9 , 3/2 , 2/3,  0.5 ) * state.5.group[2],
                        c(3/2 , 2/3, 3.5 , 3.5, 0.5 ) * state.5.group[3],
                        c(3/2 , 2/3, 3.5 , 3.5, 0.5 ) * state.5.group[4],
                        c(0.5 , 0.5, 0.5 , 0.5, 1.6 ) * state.5.group[5])

############################################################################################################################################################
### This fits to N'Zoo data where we assume that there is the same value of R_0 throughout the epidemic
############################################################################################################################################################

list[NZoo.fit, NZoo.estimates] = estimate.susceptibility.5.group.N.Zoo(waifw.init = t(waifw.5.groups), 
                                                                       R0 = 18, 
                                                                       iters = 2000, 
                                                                       chains = 4, 
                                                                       R0_sd = 4)

############################################################################################################################################################
### This fits to N'Zoo data where we assume that there is one value of R_0 at the beginning of the epidemic and a different one from week 14 onwards.
############################################################################################################################################################

list[NZoo.fit.Late.R0, NZoo.estimates.Late.R0] = estimate.susceptibility.5.group.N.Zoo.Late.R0(waifw.init = t(waifw.5.groups), 
                                                                                               R0 = 18, 
                                                                                               iters = 2000, 
                                                                                               chains = 4, 
                                                                                               R0_sd = 4)


source("Forward_simulations.R")
source("Data_NZoo.R")
time.length = 10
num.sims = 1000
list[Infections.NZoo, attack.rate.NZoo, initial.sus.NZoo] = simulate.with.R0.and.sus.dist(NZoo.estimates, NZoo.estimates$R02, time.length, num.sims,
                                                                                          cases.by.age.group = cases.by.age.group.5.group,
                                                                                          waifw.init = t(waifw.5.groups), state.5.group)


output.plots.with.observed.cases(Infections.NZoo, time.length, r = 0.15, g = 0.15, b = 0.7, "N'Zoo", first.week = 18, obs.cases = rowSums(cases.by.age.group.5.group))

list[Infections.NZoo.R0.Late, attack.rate.NZoo.R0.Late, initial.sus.NZoo.R0.Late] = simulate.with.R0.and.sus.dist(sus.dist = NZoo.estimates.Late.R0,
                                                                                                                  R0_distribution = NZoo.estimates.Late.R0$R0_Late, 
                                                                                                                  time.length, 
                                                                                                                  num.sims,
                                                                                                                  cases.by.age.group = cases.by.age.group.5.group,
                                                                                                                  waifw.init = t(waifw.5.groups), state.5.group)

output.plots.with.observed.cases(Infections.NZoo.R0.Late, time.length, r = 0.8, g = 0.15, b = 0.15, "N'Zoo", first.week = 18, obs.cases = rowSums(cases.by.age.group.5.group))
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################



#########
### Output data for each sub-prefecture that has observed cases. We output the number of cases by each different age group: Under 1's, 1-5, 6-10, 11-15, 16+;
### the number of weeks since the first case to week 17, which is the last week of cases in the data so far; the number of age groups = 5; 
### and the population of each age group for the sub-prefecture
#########


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
  # cases.by.age.group = matrix(0, length(A), 4)
  #  cases.by.age.group[, 1] = cases.by.age.in.years[, 1]
  #  cases.by.age.group[, 2] = rowSums(cases.by.age.in.years[, 2:6])
  #  cases.by.age.group[, 3] = rowSums(cases.by.age.in.years[, 7:(max(18,(max(S.prefect.data$Age.in.years) + 1)))])
  #  cases.by.age.group[, 4] = rowSums(cases.by.age.in.years[, 16:(max(18,(max(S.prefect.data$Age.in.years) + 1)))])
  
  ### specify population in each age group.
  #  state = matrix(0, 4, 1)
  #  state[1] = round((664 / 15559)* total.pop)
  #  state[2] = round((3318 - 664 + 600) * total.pop/15559)
  #  state[3] = round((round(0.425 * 15559 ) - 3318 - 600) * total.pop/15559)
  #  state[4] = round((total.pop - sum(state[1:3])) )
  #  state = as.numeric(state)
  
  
  #  denom = sum(state)
  
  
  #x = cases.by.age.group[1:(length(unique(S.prefect.data$Semaine.N.)) - 2), ]
  # y = cases.by.age.group[3:length(unique(S.prefect.data$Semaine.N.)), ]
  #  N = (length(A) - 2)
  #  M = 4
  
  
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
  return(list(N.5.group, M.5.group, cases.by.age.group.5.group, state.5.group, A))
}

list[N.5.group, M.5.group, cases.by.age.group.5.group.Foum, state.5.group.Foumbadou]    =  output.data.for.given.sub.prefecture("FOUMBADOU", 19438)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Gueasso, state.5.group.Gueasso]   =  output.data.for.given.sub.prefecture("GUEASSO", 21116)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Gama, state.5.group.Gama.Berema]  =  output.data.for.given.sub.prefecture("GAMA BEREMA", 20465)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Kokota, state.5.group.Kokota]     =  output.data.for.given.sub.prefecture("KOKOTA", 14732)
list[N.5.group, M.5.group, cases.by.age.group.5.group.C.Lola, state.5.group.Cu.Lola]    =  output.data.for.given.sub.prefecture("CU LOLA", 49933)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Laine, state.5.group.Laine]       =  output.data.for.given.sub.prefecture("LAINE", 16591)


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

### Non_N_Zoo_fits.R contains stan code to fit for other sub-prefectures. 
source("Non_N_Zoo_fits.R")


### We use the posteriors from the fits to the NZoo case data for susceptibility as the priors for the other sub districts. We also try to fit an observation rate, as observation in these sub-prefectures is assumed to be non comprehensive.

sus_prior = rbind(c(mean(NZoo.estimates$s[, 1]), sd(NZoo.estimates$s[, 1])), 
                  c(mean(NZoo.estimates$s[, 2]), sd(NZoo.estimates$s[, 2])),
                  c(mean(NZoo.estimates$s[, 3]), sd(NZoo.estimates$s[, 3])),
                  c(mean(NZoo.estimates$s[, 4]), sd(NZoo.estimates$s[, 4])),
                  c(mean(NZoo.estimates$s[, 5]), sd(NZoo.estimates$s[, 5])))

###############
### Run code for each of the remaining sub-prefectures to fit. Seems like Foumbadou is the only place where this converges fully. The others converge with varying success, with Laine being the next best and Kokota being the worst.
###############
list[Foum.ests, Foum.fits] = fit.for.given.sub.prefecture("FOUMBADOU", population = 19438, 
                                                          R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02), 
                                                          waifw = waifw.5.groups, iters = 2000, chains = 4, 
                                                          sus_prior, week.13.only = 0)
  
list[Kokota.ests, Kokota.fits] = fit.for.given.sub.prefecture("KOKOTA", population = 14732, 
                                                              R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02), 
                                                              waifw = waifw.5.groups, iters = 10000, chains = 4, 
                                                              sus_prior, week.13.only = 0)
                                                            
list[Laine.ests, Laine.fits] = fit.for.given.sub.prefecture("LAINE", population = 16591, 
                                                            R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02), 
                                                            waifw = waifw.5.groups, iters = 10000, chains = 4, 
                                                            sus_prior, week.13.only = 0)
  
list[Gama.Berema.ests, Gama.Berema.fits] = fit.for.given.sub.prefecture("GAMA BEREMA", population = 20465, 
                                                                        R0 = mean(NZoo.estimates$R02), 
                                                                        R0_sd = sd(NZoo.estimates$R02),
                                                                        waifw = waifw.5.groups, iters = 10000, 
                                                                        chains = 4, sus_prior, week.13.only = 0)
  
list[Gueasso.ests, Gueasso.fits] = fit.for.given.sub.prefecture("GUEASSO", population = 21116, 
                                                                R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02),
                                                                waifw = waifw.5.groups, iters = 10000, 
                                                                chains = 4, sus_prior, week.13.only = 0)
  
list[Cu.Lola.ests, Cu.Lola.fits] = fit.for.given.sub.prefecture("CU LOLA", population = 49933, 
                                                                R0 = mean(NZoo.estimates$R02), R0_sd = sd(NZoo.estimates$R02),
                                                                waifw = waifw.5.groups, iters = 10000, 
                                                                chains = 4, sus_prior, week.13.only = 0)
                                                                






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




#######################
### Function runs forward forecasts given a set of estimated susceptibility profiles by age group produced by stan along with R0. Similar to previous function, but this is for when we have 5 groups.
### Inputs are R0_late, the estimated R0 for after week 14; time.length - the number of weeks we project forward for; sus.dist - which is the estimated susceptiblity by age from stan; num.sims - number of forward projections performed; cases.by.age.group - the cases by age group that have been observed in the outbreak; waifw.init - the initial waifw that was used for the stan estimation, and will be used for the forward projections; state - which is the population by age group.
#######################

simulate.with.reporting.rate.adjusted.cases.R0.late <- function(sus.dist, time.length, num.sims,
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
    R0 = sample(sus.dist$R0_Late, 1)
    
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


output.plots.with.estimates.for.observed.data.sample.for.observed.cases <- function(Infections.array, time.length, r, g, b, S.Prefect, 
                                                                                    sus.dist, obs.cases.by.week, num.sims){
  
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


#list[Infections.5.groups.Foum, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Foum.ests, time.length, num.sims= 1000, waifw.init = t(waifw.5.groups), state.5.group.Foumbadou)


#list[Infections.5.groups.Gama, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Gama.Berema.ests, time.length, num.sims= 1000, waifw.init= t(waifw.5.groups), state.5.group.Gama.Berema)


#list[Infections.5.groups.Gueasso, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Gueasso.ests, time.length, num.sims= 1000,  waifw.init= t(waifw.5.groups), state.5.group.Gueasso)



list[Infections.5.groups.Foum, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases.R0.late (sus.dist = Foum.ests, time.length, num.sims= 1000, waifw.init = t(waifw.5.groups), state.5.group.Foumbadou)


list[Infections.5.groups.Gama, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases.R0.late (sus.dist = Gama.Berema.ests, time.length, num.sims= 1000, waifw.init= t(waifw.5.groups), state.5.group.Gama.Berema)


list[Infections.5.groups.Gueasso, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases.R0.late (sus.dist = Gueasso.ests, time.length, num.sims= 1000,  waifw.init= t(waifw.5.groups), state.5.group.Gueasso)


list[Infections.5.groups.C.Lola, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases.R0.late (sus.dist = Cu.Lola.ests, time.length, num.sims= 1000,  waifw.init= t(waifw.5.groups), state.5.group.Cu.Lola)

list[Infections.5.groups.Laine, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases.R0.late (sus.dist = Laine.ests, time.length, num.sims= 1000,  waifw.init= t(waifw.5.groups), state.5.group.Laine)

list[Infections.5.groups.Kokota, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases.R0.late (sus.dist = Kokota.ests, time.length, num.sims= 1000,  waifw.init= t(waifw.5.groups), state.5.group.Kokota)






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

pdf("Forecast_Cases_10_weeks.pdf", height = 12, width = 12)
par(mfrow=c(4, 2))

output.plots.with.observed.cases(Infections.NZoo, time.length, r = 0.15, g = 0.15, b = 0.7, "N'Zoo same R_0", first.week = 18, obs.cases = rowSums(cases.by.age.group.5.group))

output.plots.with.observed.cases(Infections.NZoo.R0.Late, time.length, r = 0.15, g = 0.15, b = 0.7, "N'Zoo varying R_0", first.week = 18, obs.cases = rowSums(cases.by.age.group.5.group))

output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.Foum, time.length, r, g, b, "Foumbadou", Foum.ests,obs.cases.by.week = rowSums(cases.by.age.group.5.group.Foum), num.sims, first.week = 18)

output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.Gueasso, time.length, r, g, b, "Gueasso", Gueasso.ests,obs.cases.by.week = rowSums(cases.by.age.group.5.group.Gueasso), num.sims, first.week = 18)

output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.Gama, time.length, r, g, b, "Gama Berema", Gama.Berema.ests,obs.cases.by.week = rowSums(cases.by.age.group.5.group.Gama), num.sims, first.week = 18)

output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.Laine, time.length, r, g, b, "Laine", Laine.ests, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Laine), num.sims, first.week = 18)

output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.Kokota, time.length, r, g, b, "Kokota", Kokota.ests, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Kokota), num.sims, first.week = 18)

output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.5.groups.C.Lola, time.length, r, g, b, "Central Lola", Cu.Lola.ests, obs.cases.by.week = rowSums(cases.by.age.group.5.group.C.Lola), num.sims, first.week = 18)

dev.off()

#output.plots.with.estimates.for.observed.data.non.convergence(Infections.array = Infections.5.groups.Laine, time.length, r, g, b, S.Prefect = "Laine", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Laine), num.sims,  first.week = 18)

#output.plots.with.estimates.for.observed.data.non.convergence(Infections.array = Infections.5.groups.Kokota, time.length, r, g, b, S.Prefect = "Kokota", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Kokota), num.sims,  first.week = 18)

#output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.array = Infections.5.groups.C.Lola, time.length, r, g, b, S.Prefect = "Central Lola", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.C.Lola), num.sims,  first.week = 18)
