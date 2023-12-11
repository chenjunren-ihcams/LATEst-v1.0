simulation.data.func <- function(eta, N, alpha, x.mult = FALSE) {
  
  t = 0 : 100
  
  Z <- sample(c(rep(0, round(N / 2)), 
                rep(1, N - round(N / 2))), N, replace = FALSE)
  x <- sample(c(rep(0, round(N / 2)),
                rep(1, N - round(N / 2))), N, replace = FALSE)
  
  # Set x as one-dimensional co-variate or multi-dimensional co-variates
  if (x.mult) {
    x.array <- array(NA, dim = c(N, 99))
    for (i in 1 : 99) {
      x.array[, i] <- sample(c(rep(0, round(N / 2)), 
                               rep(1, N - round(N / 2))), N, replace = FALSE)
    }
  } else {
    x.array <- array(1 / 2, dim = c(N, 1))
  }
  
  D <- array(0, N)
  D[Z == 1 & x == 1] <- sample(c(rep(0, round(sum(Z == 1 & x == 1) * 0.2)), 
                                 rep(1, sum(Z == 1 & x == 1) -
                                       round(sum(Z == 1 & x == 1) * 0.2))), 
                               sum(Z == 1 & x == 1), replace = FALSE)
  D[Z == 1 & x == 0] <- sample(c(rep(0, round(sum(Z == 1 & x == 0) * 0.8)), 
                                 rep(1, sum(Z == 1 & x == 0) -
                                       round(sum(Z == 1 & x == 0) * 0.8))), 
                               sum(Z == 1 & x == 0), replace = FALSE)
  
  ## Generate survival probabilities
  temp <- runif(N)
  death <- array(1, dim = c(N, length(t) - 1))
  survival.final <- array(NA, N)
  for (i in 1 : N) {
    death[i, ] <- 
      -diff((1 / (1 + exp((t - 60) / 20)) * 0.5 + 1 - 
               1 / (1 + exp(-60 / 20)) * 0.5) ^ 
              (1.5 ^ (eta * x[i]) * (prod(1.01 ^ (2 * x.array[i, ] - 1))) *
                 (1 + (alpha - 1) * D[i])))
    survival.final[i] <- 
      (1 / (1 + exp((100 - 60) / 20)) * 0.5 + 1 - 
         1 / (1 + exp(-60 / 20)) * 0.5) ^ 
      (1.5 ^ (eta * x[i]) * (prod(1.01 ^ (2 * x.array[i, ] - 1))) *
         (1 + (alpha - 1) * D[i]))
    
  }
  
  ## Generate event samples and time to event
  event <- array(NA, N)
  event.time <- array(NA, N)
  for (i in 1:N) {
    if (runif(1) <= survival.final[i]) {
      event[i] <- 0 
      event.time[i] <- 100
    } else {
      event.time[i] <- sample(1:100, prob = death[i, ])
      event[i] <- 1
    }
  }
  
  event[is.na(event)] <- 0
  event.time[is.na(event.time)] <- Inf
  
  ## Generate censored samples and censoring time
  C <- array(Inf, N)
  for (i in 1 : N) {
    if (runif(1) < 0.1) {
      C[i] <- sample(1 : 100, 1)
      if (C[i] < event.time[i]) {
        event.time[i] <- C[i]
        event[i] <- 0
      }
    }
  }
  
  list(Z = Z, D = D, event = event, event.time = event.time)
  
}