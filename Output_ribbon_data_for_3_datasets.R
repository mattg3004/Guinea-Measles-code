ribbon.plot.NZoo.3.sims <- function(Sims.A, Cases.A, End.Week.A,
                                    Sims.B, Cases.B, End.Week.B,
                                    Sims.C, Cases.C, End.Week.C,
                                    Week1){
  
  A = matrix(0, 1000 * nrow(Cases.A) + 1000 * nrow(Cases.B) + 1000 * nrow(Cases.C) 
             + 1000 * dim(Sims.A)[2] + 1000 * dim(Sims.B)[2] + 1000 * dim(Sims.C)[2], 4)
  count = 1
  for(i in 1 : 1000){
    for(j in 1 : nrow(Cases.A)){
      A[count, 1] = 1
      A[count, 2] = j + Week1 - 1
      A[count, 3] = sum(Cases.A[j, ])
      A[count, 4] = i 
      count = count + 1
    }
    for(j in 1 : dim(Sims.A)[2]){
      A[count, 1] = 1
      A[count, 2] = j + End.Week.A
      A[count, 3] = colSums(Sims.A[ , , i])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  for(i in 1 : 1000){
    for(j in 1 : nrow(Cases.B)){
      A[count, 1] = 2
      A[count, 2] = j + Week1 - 1
      A[count, 3] = sum(Cases.B[j, ])
      A[count, 4] = i 
      count = count + 1
    }
    for(j in 1 : dim(Sims.B)[2]){
      A[count, 1] = 2
      A[count, 2] = j + End.Week.B
      A[count, 3] = colSums(Sims.B[ , , i])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  
  for(i in 1 : 1000){
    for(j in 1 : nrow(Cases.C)){
      A[count, 1] = 3
      A[count, 2] = j + Week1 - 1
      A[count, 3] = sum(Cases.C[j, ])
      A[count, 4] = i 
      count = count + 1
    }
    for(j in 1 : dim(Sims.C)[2]){
      A[count, 1] = 3
      A[count, 2] = j + End.Week.C
      A[count, 3] = colSums(Sims.C[ , , i])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  
  A = data.frame(A)
  colnames(A) = c("run", "Week", "Cases", "iter")
  A$run = factor(A$run, levels=c(1, 2, 3), labels=c("Report 1", "Report 2 week 13", "Full Report 2"))
  gg.curve <- ggplot(A, aes(x=Week, y=Cases,id = iter, colour=as.factor(run), fill=as.factor(run))) + 
    stat_summary(fun.y=mean, geom="point") + stat_summary(fun.y=mean, geom="line") + 
    stat_summary(geom="ribbon", fun.ymin = quantile.05, fun.ymax= quantile.95, alpha=0.25, colour=NA) +
    geom_vline(xintercept = 17, color = rgb(0,0.5,1), linetype = "dashed") + 
    geom_vline(xintercept = 13, color = rgb(1,0.3,0), linetype = "dashed") + 
    theme(legend.position="bottom", legend.title=element_blank()) 
  gg.curve
  return(gg.curve)
}




ribbon.plot.Non.NZoo.3.sims <- function(Sims.A, Cases.A, End.Week.A,
                                         Sims.B, Cases.B, End.Week.B,
                                         Sims.C, Cases.C, End.Week.C,
                                         Week1){
  
  A = matrix(0, dim(Cases.A)[1]*dim(Cases.A)[2] + dim(Cases.B)[1] * dim(Cases.B)[2] + dim(Cases.C)[1] * dim(Cases.C)[2] +
              1000 * dim(Sims.A)[2] + 1000 * dim(Sims.B)[2] + 1000 * dim(Sims.C)[2], 4)
  count = 1
  for(i in 1 : 1000){
    for(j in 1 : dim(Sims.A)[2] ){
      A[count, 1] = 1
      A[count, 2] = j + End.Week.A
      A[count, 3] = colSums(Sims.A[ , , i])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  
  for(i in 1 : dim(Cases.A)[1]){
    for(j in 1 : dim(Cases.A)[2]){
      A[count, 1] = 1
      A[count, 2] = End.Week.A - dim(Cases.A)[2] + j
      A[count, 3] = rowSums(Cases.A[ i, , ])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  for(i in 1 : 1000){
    for(j in 1 : dim(Sims.B)[2] ){
      A[count, 1] = 2
      A[count, 2] = j + End.Week.B
      A[count, 3] = colSums(Sims.B[ , , i])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  
  
  for(i in 1 : dim(Cases.B)[1]){
    for(j in 1 : dim(Cases.B)[2]){
      A[count, 1] = 2
      A[count, 2] = End.Week.B - dim(Cases.B)[2] + j
      A[count, 3] = rowSums(Cases.B[ i, , ])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  for(i in 1 : 1000){
    for(j in 1 : dim(Sims.C)[2] ){
      A[count, 1] = 3
      A[count, 2] = j + End.Week.C
      A[count, 3] = colSums(Sims.C[ , , i])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  
  for(i in 1 : dim(Cases.C)[1]){
    for(j in 1 : dim(Cases.C)[2]){
      A[count, 1] = 3
      A[count, 2] = End.Week.C - dim(Cases.C)[2] + j
      A[count, 3] = rowSums(Cases.C[ i, , ])[j]
      A[count, 4] = i 
      count = count + 1
    }
  }
  
  A = data.frame(A)
  colnames(A) = c("run", "Week", "Cases", "iter")
  
  
  A$run = factor(A$run, levels=c(1, 2, 3), labels=c("Report 1", "Report 2 week 13", "Full Report 2"))
  gg.curve <- ggplot(A, aes(x=Week, y=Cases,id = iter, colour=as.factor(run), fill=as.factor(run))) + 
    stat_summary(fun.y=mean, geom="point") + stat_summary(fun.y=mean, geom="line") + 
    stat_summary(geom="ribbon", fun.ymin = quantile.05, fun.ymax= quantile.95, alpha=0.25, colour=NA) +
    geom_vline(xintercept = 17, color = rgb(0,0.5,1), linetype = "dashed") + 
    geom_vline(xintercept = 13, color = rgb(1,0.3,0), linetype = "dashed")  +
    theme(legend.position="bottom", legend.title=element_blank()) 
  
  return(gg.curve)
}

