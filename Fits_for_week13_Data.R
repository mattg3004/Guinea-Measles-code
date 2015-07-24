library(rstan)
source("Functions_for_tsir_fits.R")

############################################################################################################################################################
### Do fits for N'Zoo district, as we assume that there is 100% observation in this sub-prefecture
############################################################################################################################################################
source("NZoo_fits.R")

waifw.5.groups =  rbind(c(6.9 , 3 , 3/2 , 2/3,  0.5 ) * state.5.group[1],
                        c(3 , 6.9 , 3/2 , 2/3,  0.5 ) * state.5.group[2],
                        c(3/2 , 3/2, 3.5 , 3.5, 0.5 ) * state.5.group[3],
                        c(2/3 , 2/3, 3.5 , 3.5, 0.5 ) * state.5.group[4],
                        c(0.5 , 0.5, 0.5 , 0.5, 1.6 ) * state.5.group[5])

############################################################################################################################################################
### This fits to N'Zoo data where we assume that there is the same value of R_0 throughout the epidemic
############################################################################################################################################################
source("Data_NZoo_week13.R")
list[NZoo.fit.13, NZoo.estimates.13] = estimate.susceptibility.5.group.N.Zoo.week.13(waifw.init = t(waifw.5.groups), 
                                                                                     R0 = 18, 
                                                                                     iters = 2000, 
                                                                                     chains = 4, 
                                                                                     R0_sd = 4, alpha = 0.3)
source("Forward_simulations.R")

time.length = 14
num.sims = 1000
list[Infections.NZoo.13, attack.rate.NZoo.13, initial.sus.NZoo.13] = simulate.with.R0.and.sus.dist(NZoo.estimates.13, NZoo.estimates.13$R02, time.length, num.sims,
                                                                                                   cases.by.age.group = cases.by.age.group.5.group,
                                                                                                   waifw.init = t(waifw.5.groups), state.5.group)


output.plots.with.observed.cases(Infections.NZoo.13, time.length, r = 0.15, g = 0.15, b = 0.7, "N'Zoo", first.week = 14, obs.cases = rowSums(cases.by.age.group.5.group))

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################



#########
### Output data for each sub-prefecture that has observed cases. We output the number of cases by each different age group: Under 1's, 1-5, 6-10, 11-15, 16+;
### the number of weeks since the first case to week 17, which is the last week of cases in the data so far; the number of age groups = 5; 
### and the population of each age group for the sub-prefecture
#########


output.data.for.given.sub.prefecture.week.13 <- function(S.prefect, total.pop){
  ### read in data on cases
  data = read.csv("Lola_Measles_week_13_included.csv")
  
  ### subset to cases in given sub prefecture
  S.prefect.data = subset(data, as.character(data$Sous.Prefecture) == S.prefect)
  week = matrix(0, nrow(S.prefect.data), 1)
  for(i in 1 : nrow(S.prefect.data)){
    week[i] = as.numeric(strsplit(as.character(S.prefect.data$Semaine.N.), "-")[[i]][2])
  }
  #if (S.prefect == "FOUMBADOU"){
  #  week = week + 1
  #}
  S.prefect.data$Week = week
  A =  as.character(unique(S.prefect.data$Week)[order((unique(S.prefect.data$Week)))])
  A = A[A != ""]
  A = seq(min(as.numeric(A)), 13)
  
  ## split into cases by week
  cases.by.week = matrix(0, length(A), 2)
  cases.by.week[, 1] = as.character(A)
  for(i in 1:length(A)){
    cases.by.week[i, 2] = length(which(S.prefect.data$Week == A[i]))
  }
  
  ### add in age in years variable
  S.prefect.data$Age.in.years = floor(S.prefect.data$Age.Mois/ 12)
  #S.prefect.data$Age.in.years = floor(S.prefect.data$AGE)
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

list[N.5.group, M.5.group, cases.by.age.group.5.group.NZoo.13, state.5.group.Foumbadou.13]    =  output.data.for.given.sub.prefecture.week.13("NZOO", 1)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Foum.13, state.5.group.Foumbadou.13]    =  output.data.for.given.sub.prefecture.week.13("FOUMBADOU", 19438)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Gueasso.13, state.5.group.Gueasso.13]   =  output.data.for.given.sub.prefecture.week.13("GUEASSO", 21116)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Gama.13, state.5.group.Gama.Berema.13]  =  output.data.for.given.sub.prefecture.week.13("GAMA BEREMA", 20465)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Kokota.13, state.5.group.Kokota]     =  output.data.for.given.sub.prefecture.week.13("KOKOTA", 14732)
list[N.5.group, M.5.group, cases.by.age.group.5.group.C.Lola.13, state.5.group.Cu.Lola]    =  output.data.for.given.sub.prefecture.week.13("CU LOLA", 49933)
list[N.5.group, M.5.group, cases.by.age.group.5.group.Laine.13, state.5.group.Laine]       =  output.data.for.given.sub.prefecture.week.13("LAINE", 16591)

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

### Non_N_Zoo_fits.R contains stan code to fit for other sub-prefectures. 
source("Non_N_Zoo_fits.R")


### We use the posteriors from the fits to the NZoo case data for susceptibility as the priors for the other sub districts. We also try to fit an observation rate, as observation in these sub-prefectures is assumed to be non comprehensive.

sus_prior_13 = rbind(c(mean(NZoo.estimates.13$s[, 1]), sd(NZoo.estimates.13$s[, 1])), 
                  c(mean(NZoo.estimates.13$s[, 2]), sd(NZoo.estimates.13$s[, 2])),
                  c(mean(NZoo.estimates.13$s[, 3]), sd(NZoo.estimates.13$s[, 3])),
                  c(mean(NZoo.estimates.13$s[, 4]), sd(NZoo.estimates.13$s[, 4])),
                  c(mean(NZoo.estimates.13$s[, 5]), sd(NZoo.estimates.13$s[, 5])))

###############
### Run code for each of the remaining sub-prefectures to fit. Seems like Foumbadou is the only place where this converges fully. The others converge with varying success, with Laine being the next best and Kokota being the worst.
###############
source("Fits_non_nZoo_week_13.R")
list[Foum.ests.13.cases.limit, Foum.fits.13.cases.limits] = fit.for.given.sub.prefecture.week.13.upper.cases.limits("FOUMBADOU", population = 19438, 
                                                                                                                    R0 = mean(NZoo.estimates.13$R02), R0_sd = sd(NZoo.estimates.13$R02), 
                                                                                                                    waifw = waifw.5.groups, iters = 1000, chains = 4, 
                                                                                                                    sus_prior_13)
                                                                        

list[Foum.ests.13.no.cases.limit, Foum.fits.13.no.cases.limits] = fit.for.given.sub.prefecture.week.13.no.case.limits("FOUMBADOU", population = 19438, 
                                                                                                                       R0 = mean(NZoo.estimates.13$R02), R0_sd = sd(NZoo.estimates.13$R02), 
                                                                                                                       waifw = waifw.5.groups, iters = 2000, chains = 4, 
                                                                                                                       sus_prior_13)

list[Foum.ests.13.upper.cases.limit.no.obs.lower,
     Foum.fits.13.upper.cases.limit.no.obs.lower] = fit.for.given.sub.prefecture.week.13.no.case.limits.no.obs.lower ("FOUMBADOU", population = 19438, 
                                                                                                                             R0 = mean(NZoo.estimates.13$R02), 
                                                                                                                             R0_sd = sd(NZoo.estimates.13$R02),
                                                                                                                             waifw = waifw.5.groups, iters = 3000, 
                                                                                                                             chains = 4, sus_prior_13)


list[Gama.Berema.ests.13.no.cases.limit, Gama.Berema.fits.13.no.cases.limits] = fit.for.given.sub.prefecture.week.13.no.case.limits("GAMA BEREMA", population = 20465, 
                                                                                                                                    R0 = mean(NZoo.estimates.13$R02), 
                                                                                                                                    R0_sd = sd(NZoo.estimates.13$R02), 
                                                                                                                                    waifw = waifw.5.groups, iters = 1000, chains = 4, 
                                                                                                                                    sus_prior_13)

list[Gama.Berema.ests.13.upper.cases.limit, Gama.Berema.fits.13.upper.cases.limit] = fit.for.given.sub.prefecture.week.13.upper.cases.limits("GAMA BEREMA", population = 20465, 
                                                                                                                                             R0 = mean(NZoo.estimates.13$R02), 
                                                                                                                                             R0_sd = sd(NZoo.estimates.13$R02),
                                                                                                                                             waifw = waifw.5.groups, iters = 1000, 
                                                                                                                                             chains = 4, sus_prior_13)

list[Gama.Berema.ests.13.upper.cases.limit.no.obs.lower,
     Gama.Berema.fits.13.upper.cases.limit.no.obs.lower] = fit.for.given.sub.prefecture.week.13.no.case.limits.no.obs.lower ("GAMA BEREMA", population = 20465, 
                                                                                                                             R0 = mean(NZoo.estimates.13$R02), 
                                                                                                                             R0_sd = sd(NZoo.estimates.13$R02),
                                                                                                                             waifw = waifw.5.groups, iters = 1000, 
                                                                                                                             chains = 4, sus_prior_13)


list[Gueasso.ests.13, Gueasso.fits.13] = fit.for.given.sub.prefecture.week.13("GUEASSO", population = 21116, 
                                                                              R0 = mean(NZoo.estimates.13$R02), R0_sd = sd(NZoo.estimates.13$R02),
                                                                              waifw = waifw.5.groups, iters = 1000, 
                                                                              chains = 4, sus_prior_13)








reporting.rate.dist.13 = c(Foum.ests.13$delta, Gama.Berema.ests.13$delta, Gueasso.ests.13$delta)
s1.dist.13 = c(Foum.ests.13$s[, 1], Gama.Berema.ests.13$s[, 1], Gueasso.ests.13$s[, 1])
s2.dist.13 = c(Foum.ests.13$s[, 2], Gama.Berema.ests.13$s[, 2], Gueasso.ests.13$s[, 2])
s3.dist.13 = c(Foum.ests.13$s[, 3], Gama.Berema.ests.13$s[, 3], Gueasso.ests.13$s[, 3])
s4.dist.13 = c(Foum.ests.13$s[, 4], Gama.Berema.ests.13$s[, 4], Gueasso.ests.13$s[, 4])
s5.dist.13 = c(Foum.ests.13$s[, 5], Gama.Berema.ests.13$s[, 5], Gueasso.ests.13$s[, 5])
R0.dist.13 = c(Foum.ests.13$R02, Gama.Berema.ests.13$R02, Gueasso.ests.13$R02)





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




output.plots.for.sims.week.13 <- function(Infections.array, time.length, r, g, b, S.Prefect, sus.dist, obs.cases.by.week, num.sims){
  
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
  #legend("topright", legend = S.Prefect, 
  #       col = rgb(red = r, green = g, blue = b),
  #       lwd = 2, bty="n", cex = 1)
  
}





list[Infections.5.groups.Foum.13, attack.rate.5.groups.13, initial.5.groups.13] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Foum.ests.13.cases.limit, 
                                                                                                                               time.length, num.sims= 1000, waifw.init = t(waifw.5.groups), state.5.group.Foumbadou)


list[Infections.5.groups.Gama.13, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Gama.Berema.ests.13, time.length, num.sims= 1000, waifw.init= t(waifw.5.groups), state.5.group.Gama.Berema)


list[Infections.5.groups.Gueasso.13, attack.rate.5.groups, initial.5.groups] = simulate.with.reporting.rate.adjusted.cases (sus.dist = Gueasso.ests.13, time.length, num.sims= 1000,  waifw.init= t(waifw.5.groups), state.5.group.Gueasso)





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


output.plots.with.estimates.for.observed.data.non.convergence.week.13 <- function(Infections.array, time.length, r, g, b, S.Prefect, reporting.rate.dist, obs.cases.by.week, num.sims){
  
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
  #legend("topright", legend = S.Prefect, 
  #       col = rgb(red = r, green = g, blue = b),
  #       lwd = 2, bty="n", cex = 1)
  return(cases.by.week)
}




list[Infections.5.groups.Kokota.13, attack.rate.5.groups, initial.5.groups]  = 
  simulate.with.reporting.rate.non.converged(s1.dist.13, s2.dist.13, s3.dist.13, s4.dist.13, s5.dist.13,
                                             reporting.rate.dist.13, R0.dist.13,
                                             time.length, num.sims, cases.by.age.group.5.group.Kokota.13,
                                             t(waifw.5.groups), state = state.5.group.Kokota)



list[Infections.5.groups.Laine.13, attack.rate.5.groups, initial.5.groups]  = 
  simulate.with.reporting.rate.non.converged(s1.dist.13, s2.dist.13, s3.dist.13, s4.dist.13, s5.dist.13,
                                             reporting.rate.dist.13, R0.dist.13,
                                             time.length, num.sims, cases.by.age.group.5.group.Laine.13,
                                             t(waifw.5.groups), state = state.5.group.Laine)


list[Infections.5.groups.C.Lola.13, attack.rate.5.groups, initial.5.groups]  = 
  simulate.with.reporting.rate.non.converged(s1.dist.13, s2.dist.13, s3.dist.13, s4.dist.13, s5.dist.13,
                                             reporting.rate.dist.13, R0.dist.13,
                                             time.length, num.sims, cases.by.age.group.5.group.C.Lola.13,
                                             t(waifw.5.groups), state = state.5.group.Cu.Lola)



par(mfrow=c(3, 2))
r = 0.15
g = 0.15
b = 0.7

pdf("Forecast_Cases_10_weeks_week_13.pdf", height = 12, width = 12)
par(mfrow=c(4, 2))

output.plots.with.observed.cases(Infections.NZoo.13, time.length, r = 0.15, g = 0.15, b = 0.7, "N'Zoo", first.week = 14, obs.cases = rowSums(cases.by.age.group.5.group.13))

output.plots.for.sims.week.13(Infections.5.groups.Foum.13, time.length, r, g, b, "Foumbadou", Foum.ests.13, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Foum.13), num.sims)

output.plots.for.sims.week.13(Infections.5.groups.Gueasso.13, time.length, r, g, b, "Gueasso", Gueasso.ests.13, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Gueasso.13), num.sims)

output.plots.for.sims.week.13(Infections.5.groups.Gama.13, time.length, r, g, b, "Gama Berema", Gama.Berema.ests.13, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Gama.13), num.sims)

output.plots.with.estimates.for.observed.data.non.convergence.week.13(Infections.5.groups.Laine.13, num.sims = 1000,time.length, r, g, b, "Laine", reporting.rate.dist.13, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Laine.13))

output.plots.with.estimates.for.observed.data.non.convergence.week.13(Infections.5.groups.Kokota.13, time.length, r, g, b, "Kokota", reporting.rate.dist.13, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Kokota.13), num.sims)

output.plots.with.estimates.for.observed.data.non.convergence.week.13(Infections.5.groups.C.Lola.13, time.length, r, g, b, "Central Lola", reporting.rate.dist.13, obs.cases.by.week = rowSums(cases.by.age.group.5.group.C.Lola.13), num.sims)

dev.off()

#output.plots.with.estimates.for.observed.data.non.convergence(Infections.array = Infections.5.groups.Laine, time.length, r, g, b, S.Prefect = "Laine", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Laine), num.sims,  first.week = 18)

#output.plots.with.estimates.for.observed.data.non.convergence(Infections.array = Infections.5.groups.Kokota, time.length, r, g, b, S.Prefect = "Kokota", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.Kokota), num.sims,  first.week = 18)

#output.plots.with.estimates.for.observed.data.sample.for.observed.cases(Infections.array = Infections.5.groups.C.Lola, time.length, r, g, b, S.Prefect = "Central Lola", reporting.rate.dist, obs.cases.by.week = rowSums(cases.by.age.group.5.group.C.Lola), num.sims,  first.week = 18)
