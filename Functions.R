
forecast.forwards <- function(data, r.effective, time.length, simulations){
  Infections = matrix(0, simulations, time.length)
  for(i in 1 : simulations){
    I.t = data[length(data) - 1] 
    for(j in 1 : time.length){
      k = ceiling(1000*runif(1))
      r = r.effective[k]
      Infections[i, j] = rpois(1, I.t * r)
      if (j == 1) {I.t = tail(data, 1)}
      if (j == 2) {I.t = (tail(data, 1) + Infections[i, j -1])/2} 
      if (j > 2)  { I.t = (Infections[i, j -1] + Infections[i, j ])/ 2}
    }
  }
  return(Infections)
}



forecast.forwards.susceptible <- function(data, r.effective, time.length, simulations, susceptibles, pop){
  Infections = matrix(0, simulations, time.length)
  for(i in 1 : simulations){
    I.t = data[length(data) - 1]
    sus = tail(susceptibles, 1)
    for(j in 1 : time.length){
      k = ceiling(1000*runif(1))
      r = r.effective[k]
      Infections[i, j] = rpois(1, I.t * r * sus/ pop)
      sus = sus - Infections[i, j] + 12
      if (j == 1) {I.t = tail(data, 1)}
      if (j == 2) {I.t = (tail(data, 1) + Infections[i, j - 1])/2} 
      if (j > 2)  { I.t = (Infections[i, j -1] + Infections[i, j ])/ 2}
    }
  }
  return(Infections)
}



forecast.forwards.susceptible.by.age <- function(data, r.effective, time.length, simulations, susceptibles, pop.age){
  Infections = array(0, c(length(susceptibles[1, ]), time.length, simulations))
  for(i in 1 : simulations){
    I.t = sum(tail(data, 2))/2
    sus = susceptibles[length(susceptibles[, 1]), ]
    for(j in 1 : time.length){
      for(l in 1 : length(sus)){
        k = ceiling(length(r.effective) * runif(1))
        r = r.effective[k]
        Infections[l, j, i] = rpois(1, I.t * r * sus[l]/ pop.age[l])
        sus[l] = max(0, sus[l] - Infections[l, j, i])
      }
      if (j == 1) {I.t = (tail(data, 1) + sum(Infections[, j , i]))/2}
      #if (j == 2) {I.t = (tail(data, 1) + sum(Infections[, j - 1, i]))/2} 
      if (j > 1)  { I.t = (sum(Infections[, j - 1, i]) + sum(Infections[, j, i]))/ 2}
    }
  }
  return(Infections)
}



forecast.forwards.susceptible.by.age.no.averaging <- function(data, r.effective, time.length, simulations, susceptibles, pop.age){
  Infections = array(0, c(length(susceptibles[1, ]), time.length, simulations))
  for(i in 1 : simulations){
    I.t = data[length(data) - 1]
    sus = susceptibles[length(susceptibles[, 1]), ]
    for(j in 1 : time.length){
      for(l in 1 : length(sus)){
        k = ceiling(length(r.effective) * runif(1))
        r = r.effective[k]
        Infections[l, j, i] = rpois(1, I.t * r * sus[l]/ pop.age[l])
        sus[l] = max(0, sus[l] - Infections[l, j, i])
      }
      if (j == 1) {I.t = tail(data, 1)}
      if (j > 1)  { I.t = sum(Infections[, j - 1, i]) }
    }
  }
  return(Infections)
}



output.plots <- function(Infections.array, time.length, r, g, b, S.Prefect, first.week){
  mean.cumulative.Infections.array = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){
    mean.cumulative.Infections.array[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  
  
  
  
  plot(seq(first.week, first.week+ time.length - 1),colSums( Infections.array[, , 1]),  col = rgb(red = r, green = g, blue = b, alpha = 0.1), type = "l", 
       ylim = c(0,max(colSums(Infections.array))), xaxt = "n",  xlab = "Week", ylab = "Predicted new cases")
  axis(side = 1, at = seq(first.week, (first.week + time.length - 1)))
  for(i in 1 : 1000){
    lines(seq(first.week, first.week + time.length - 1), colSums(Infections.array[, , i]), col = rgb(red = r, green = g, blue = b, alpha = 0.1))
  }
  means = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means[i] = mean(colSums(Infections.array[, i, ]))
  }
  lines(seq(first.week, first.week + time.length - 1), means, lwd = 3, col = rgb(red = r, green = g, blue = b))  
  legend("topright", legend = S.Prefect, 
         col = rgb(red = r, green = g, blue = b),
         lwd = 2, bty="n", cex = 1)

}


combine.plots <- function(Infections.array1, Infections.array2, Infections.array3,
                          R0.1, R0.2, R0.3,
                          r1, g1, b1, r2, g2, b2, r3, g3, b3,
                          time.length){
  
  maximum = max(max(colSums(Infections.array1)), max(colSums(Infections.array2)), max(colSums(Infections.array3)))
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
  
  
  
  Infections.array = Infections.array3
  
  for(i in 1 : 1000){
    lines(seq(14, 14 + time.length - 1), colSums(Infections.array[, , i]), col = rgb(red = r3, green = g3, blue = b3, alpha = 0.1))
  }
  means3 = matrix(0, time.length, 1)
  for(i in 1 : time.length){
    means3[i] = mean(colSums(Infections.array[, i, ]))
  }
  lines(seq(14, 14 + time.length - 1), means1, lwd = 3, col = rgb(red = r1, green = g1, blue = b1))  
  lines(seq(14, 14 + time.length - 1), means2, lwd = 3, col = rgb(red = r2, green = g2, blue = b2))  
  lines(seq(14, 14 + time.length - 1), means3, lwd = 3, col = rgb(red = r3, green = g3, blue = b3))  
  
  legend("topright", legend = c(paste("R0 =", R0.1), paste("R0 =", R0.2), paste("R0 =", R0.3)), 
         col = c(rgb(red = r1, green = g1, blue = b1), rgb(red = r2, green = g2, blue = b2), rgb(red = r3, green = g3, blue = b3)),
         lwd = 2, bty="n", cex = 1)
}


combine.cases.by.age.plots <- function(Infections.array1, Infections.array2, Infections.array3,
                                       r1, g1, b1, r2, g2, b2, r3, g3, b3,
                                       time.length){
 
 
  Infections.array = Infections.array1
  mean.cumulative.Infections.array1 = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){
    mean.cumulative.Infections.array1[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  

  
  Infections.array = Infections.array2
  mean.cumulative.Infections.array2 = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){
    mean.cumulative.Infections.array2[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  
  
  Infections.array = Infections.array3
  mean.cumulative.Infections.array3 = matrix(0, 23, 1)
  for(i in 1 : length(Infections.array[, 1, 1])){
    mean.cumulative.Infections.array3[i] = mean(Infections.array[i, , ] * length(Infections.array[1, , 1]))
  }
  
  p = max(mean.cumulative.Infections.array1, mean.cumulative.Infections.array2, mean.cumulative.Infections.array3)
  
  plot(0:22, mean.cumulative.Infections.array1, xlab = "Age", col = rgb(red = r1, green = g1, blue = b1),
       ylab = paste("Expected infections over next", time.length,"weeks"), ylim = c(0, p), type = "l", lwd = 2)
  lines(0:22, mean.cumulative.Infections.array2,col = rgb(red = r2, green = g2, blue = b2), lwd = 2)
  lines(0:22, mean.cumulative.Infections.array3,col = rgb(red = r3, green = g3, blue = b3), lwd = 2)
  
}
  

calc.phi.1.week <- function(waifw, x, denom){
  
  phi <- (waifw %*% x ^ 0.97) /denom
  
  
  phi <- 1 - exp(-phi)
  
  phi <- t(phi)
  
  return(phi)
}


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