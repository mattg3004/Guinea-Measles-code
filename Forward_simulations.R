##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##
##' Simulate forward in time. 
##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##
simulate.with.R0.and.sus.dist<- function(sus.dist, R0_distribution, time.length, num.sims,
                                         cases.by.age.group = cases.by.age.group,
                                         waifw.init, state){ 
  
  
  ##' Initialize the variables we want to output from the function. These are, for each simulation, the attack rate over each age group, the initial number of susceptibles in each age group, and the number of infections in each age group.
  attack.rate = matrix(0, num.sims, 5)
  initial.sus = matrix(0, num.sims, 5)
  end.sus = matrix(0, num.sims, 5)
  Infections = array(0, c(length(state), time.length, num.sims))
  denom = sum(state)
  ##' Perform a loop for each of the simulations
  for(i in 1 : num.sims){
    
    ##' Select a beginning susceptibility for each age group from the stan output
    start.sus.proportion = c(sample(sus.dist$s[, 1], 1, replace = T), sample(sus.dist$s[, 2], 1, replace = T),
                             sample(sus.dist$s[, 3], 1, replace = T), sample(sus.dist$s[, 4], 1, replace = T),
                             sample(sus.dist$s[, 5], 1, replace = T))
    
    ##' Select a value of R_0 to use in the simulation
    
    R0 = sample(R0_distribution, 1)
    
    alpha = 0.3
    
    if( length(sus.dist$alpha > 0)){
      alpha = sample(sus.dist$alpha, 1)
    }
    ##' This susceptibility is for before the epidemic began, therefore we subtract any susceptibles who were infected in the outbreak up to the end of the data.
    sus = ceiling(start.sus.proportion * state) - colSums(cases.by.age.group)
    
    ##' Store these numbers of susceptibles in the initial.sus matrix, which will be output from the function
    initial.sus[i, ] = sus
    
    ##' Check that none of the numbers for initial susceptibles is <0.
    for(i1 in 1:length(sus)){
      sus[i1] = max(0, sus[i1])
    }
    
    ##' Average the last two weeks of cases by age group to give the initial number of susceptibles by age group. This will be used to calculate the force of infection on each age group through the waifw. 
    ##' The averaging is done, as if we continually use the number of cases from two weeks prior to the present week for prediction, then we will get a pattern which is too dependent on the number of cases, for example if we had 10 one week and 30 the next, then projecting this forwards, we would see a pattern of low then high incidence in the predictions.
    ##' Another option is to look at bi-weeks, though I haven't looked at this yet.
    I.t = colSums(tail(cases.by.age.group, 2))/2
    
    ##' We now loop over the number of weeks that we have chosen to project forwards for.
    for(j in 1 : time.length){
      
      ##' Calculate the waifw, which takes into account the seasonality of the outbreak and also calculate the force of infection by age group, which is dependent on the waifw and the number of infecteds in each age group.
      waifw = output.waifw(waifw.init, R0 * (1 + cos((3/12 + j/48)  * 2 *  pi) * alpha), state)
      phi = calc.phi.1.week(waifw, as.numeric(I.t), denom)
      
      ##' Loop over each of the age group to generate the forward projections of number of cases by age group
      for(l in 1 : length(sus)){
        ##' The number of infections by age group is assumed to be binomially distributed Bin(n,p), with n = number of susceptibles in the age group; p = force of infection for this age group
        Infections[l, j, i] = rbinom(1, sus[l], phi[l])
        
        ##' Update the number of susceptibles by subtracting the number of infections from the number of susceptibles in the age group
        sus[l] = sus[l] - Infections[l, j, i]
      }
      ##' Update the number of infections for the next time step
      if (j == 1) {I.t = (tail(cases.by.age.group, 1) + (Infections[, j , i]))/2}
      if (j > 1)  { I.t = ((Infections[, j - 1, i]) + (Infections[, j, i]))/ 2}
    }
    end.sus[i, ] = sus
    for(pp in 1 : 5){
      ##' Calculate the attack rate for each of the age groups
      attack.rate[i, pp] = (sum(Infections[pp, , i]))/state[pp]
    }
  }
  
  ##' Output the projections along with the attack rate and initial susceptibles for each simulation. 
  return(list(Infections, attack.rate, initial.sus))
}




##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##
##' Simulate forward in time. 
##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##
  simulate.with.R0.and.sus.dist.obs.rates<- function(sus.dist, R0_distribution, time.length, num.sims,
                                                   cases.by.age.group = cases.by.age.group,
                                                   waifw.init, state, first.week = 14){ 
  
  
  ##' Initialize the variables we want to output from the function. These are, for each simulation, the attack rate over each age group, the initial number of susceptibles in each age group, and the number of infections in each age group.
  attack.rate = matrix(0, num.sims, 5)
  initial.sus = matrix(0, num.sims, 5)
  end.sus = matrix(0, num.sims, 5)
  sample.cases = array(0, c(dim(sus.dist$Cases_by_age)[2], dim(sus.dist$Cases_by_age)[3], num.sims))
  Infections = array(0, c(length(state), time.length, num.sims))
  denom = sum(state)
  ##' Perform a loop for each of the simulations
  for(ll in 1 : dim(sus.dist$Cases_by_age)[2]){
    for(kk in 1 : dim(sus.dist$Cases_by_age)[3]){
      sample.cases[ll, kk, ] =  sample(sus.dist$Cases_by_age[,ll,kk], num.sims, replace = T)
    }
  }
  for(i in 1 : num.sims){
    ##' Select a beginning susceptibility for each age group from the stan output
    start.sus.proportion = c(sample(sus.dist$s[, 1], 1, replace = T), sample(sus.dist$s[, 2], 1, replace = T),
                             sample(sus.dist$s[, 3], 1, replace = T), sample(sus.dist$s[, 4], 1, replace = T),
                             sample(sus.dist$s[, 5], 1, replace = T))
    ##' Select a value of R_0 to use in the simulation
    R0 = sample(R0_distribution, 1)
    alpha = 0.3
    if( length(sus.dist$alpha > 0)){
      alpha = sample(sus.dist$alpha, 1)
    }
    ##' This susceptibility is for before the epidemic began, therefore we subtract any susceptibles who were infected in the outbreak up to the end of the data.
    sus = ceiling(start.sus.proportion * state) - colSums(sample.cases[, , i])
    ##' Store these numbers of susceptibles in the initial.sus matrix, which will be output from the function
    initial.sus[i, ] = sus
    ##' Check that none of the numbers for initial susceptibles is <0.
    for(i1 in 1:length(sus)){
      sus[i1] = round(max(0, sus[i1]))
    }
    ##' Average the last two weeks of cases by age group to give the initial number of susceptibles by age group. This will be used to calculate the force of infection on each age group through the waifw. 
    ##' The averaging is done, as if we continually use the number of cases from two weeks prior to the present week for prediction, then we will get a pattern which is too dependent on the number of cases, for example if we had 10 one week and 30 the next, then projecting this forwards, we would see a pattern of low then high incidence in the predictions.
    ##' Another option is to look at bi-weeks, though I haven't looked at this yet.
    I.t = colSums(tail(sample.cases[, , i], 2))/2
    ##' We now loop over the number of weeks that we have chosen to project forwards for.
    for(j in 1 : time.length){
      
      ##' Calculate the waifw, which takes into account the seasonality of the outbreak and also calculate the force of infection by age group, which is dependent on the waifw and the number of infecteds in each age group.
      waifw = output.waifw(waifw.init, R0 * (1 + cos((3/12 + j/52 + (first.week -14)/52)  * 2 *  pi) * alpha), state)
      phi = calc.phi.1.week(waifw, as.numeric(I.t), denom)
      
      ##' Loop over each of the age group to generate the forward projections of number of cases by age group
      for(l in 1 : length(sus)){
        ##' The number of infections by age group is assumed to be binomially distributed Bin(n,p), with n = number of susceptibles in the age group; p = force of infection for this age group
        Infections[l, j, i] = rbinom(1, sus[l], phi[l])
        
        ##' Update the number of susceptibles by subtracting the number of infections from the number of susceptibles in the age group
        sus[l] = sus[l] - Infections[l, j, i]
      }
      ##' Update the number of infections for the next time step
      if (j == 1) {I.t = (tail(sample.cases[, , i], 1) + (Infections[, j , i]))/2}
      if (j > 1)  { I.t = ((Infections[, j - 1, i]) + (Infections[, j, i]))/ 2}
    }
    end.sus[i, ] = sus
    for(pp in 1 : 5){
      ##' Calculate the attack rate for each of the age groups
      attack.rate[i, pp] = (sum(Infections[pp, , i]))/state[pp]
    }
  }
  
  ##' Output the projections along with the attack rate and initial susceptibles for each simulation. 
  return(list(Infections, attack.rate, initial.sus))
}



##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##
##' Plot forward projections along with observed cases for N'Zoo
##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##

output.plots.with.observed.cases <- function(Infections.array, time.length, r, g, b, S.Prefect, first.week, obs.cases){
  mean.cumulative.Infections.array = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){
    mean.cumulative.Infections.array[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  
  
  means = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means[i] = mean(colSums(Infections.array[, i, ]))
  }
  a = c(obs.cases, means)
  plot(seq(4, 4 + length(a) - 1), a,  col = rgb(red = r, green = g, blue = b), type = "l", 
       ylim = c(0,max(colSums(Infections.array))), xaxt = "n",  xlab = "Week", ylab = "Predicted new cases")
  axis(side = 1, at = seq(4, 4 + length(a) - 1))
  for(i in 1 : 1000){
    lines(seq(first.week - 1, first.week + time.length - 1), c(obs.cases[length(obs.cases)], colSums(Infections.array[, , i])), col = rgb(red = r, green = g, blue = b, alpha = 0.1))
  }
  
  lines(seq(first.week, first.week + time.length - 1), means, lwd = 3, col = rgb(red = r, green = g, blue = b))  
 # legend("topright", legend = S.Prefect, 
 #        col = rgb(red = r, green = g, blue = b),
 #        lwd = 2, bty="n", cex = 1)
  
}




output.plots.with.observed.cases.over.plot <- function(Infections.array, time.length, r, g, b, S.Prefect, first.week, obs.cases){
  mean.cumulative.Infections.array = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){
    mean.cumulative.Infections.array[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  
  
  means = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means[i] = mean(colSums(Infections.array[, i, ]))
  }
  a = c(obs.cases, means)
  lines(seq(4, 4 + length(a) - 1), a,  col = rgb(red = r, green = g, blue = b))
  axis(side = 1, at = seq(4, 4 + length(a) - 1))
  for(i in 1 : 1000){
    lines(seq(first.week-1, first.week + time.length - 1), c(obs.cases[length(obs.cases)], colSums(Infections.array[, , i])), col = rgb(red = r, green = g, blue = b, alpha = 0.1))
  }
  
  lines(seq(first.week, first.week + time.length - 1), means, lwd = 3, col = rgb(red = r, green = g, blue = b))  
 # legend("topright", legend = S.Prefect, 
 #        col = rgb(red = r, green = g, blue = b),
  #       lwd = 2, bty="n", cex = 1)
  
}





##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##
##' Plot forward projections along with observed cases for N'Zoo
##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##'##

output.plots.with.observed.cases.specify.axis <- function(Infections.array, time.length, r, g, b, S.Prefect, first.week, obs.cases, ylim.max){
  mean.cumulative.Infections.array = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){ 
    mean.cumulative.Infections.array[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  
  
  means = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means[i] = mean(colSums(Infections.array[, i, ]))
  }
  a = c(obs.cases, means)
  plot(seq(4, 4 + length(a) - 1), a,  col = rgb(red = r, green = g, blue = b), type = "l", 
       ylim = c(0,ylim.max), xaxt = "n",  xlab = "Week", ylab = "Predicted new cases")
  axis(side = 1, at = seq(4, 4 + length(a) - 1))
  for(i in 1 : 1000){
    lines(seq(first.week, first.week + time.length - 1), colSums(Infections.array[, , i]), col = rgb(red = r, green = g, blue = b, alpha = 0.1))
  }
  
  lines(seq(first.week, first.week + time.length - 1), means, lwd = 3, col = rgb(red = r, green = g, blue = b))  
  legend("topright", legend = S.Prefect, 
         col = rgb(red = r, green = g, blue = b),
         lwd = 2, bty="n", cex = 1)
  
}




