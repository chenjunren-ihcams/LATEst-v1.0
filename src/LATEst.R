LATEst.func <- function(Z, D, event, event.time){
  t <- 0 : ceiling(max(event.time, na.rm = T))
  
  f.g.func <- function(I){
    y <- array(NA, length(t))
    y[t == 0] <- 1
    
    for (j in 2 : length(t)) {
      y[j] <- y[j - 1] *
        sum(event.time[I] > t[j] | (event.time[I] == t[j] & event[I] == 0)) / 
        sum(event.time[I] >= t[j])
    }
    II <- is.na(diff(y))
    
    y <- c(1, 1 + cumsum(predict(smooth.spline(diff(y)[!II], 
                                               x = (1 : (length(y) - 1))[!II],
                                               cv = TRUE),
                                 x = 1 : (length(y) - 1))$y))
    
    y
  }
  
  f <- f.g.func(Z == 0)
  g <- f.g.func(Z == 1 & D == 1)
  
  p.complier <- sum(Z == 1 & D == 1) / sum(Z == 1)
  
  deriv <- function(ff) {
    m <- length(ff)
    ff - c(1, ff[1 : (m - 1)])
  }
  
  # Compute log-likelihood
  log.likelihood <- function(alpha) {
    h <- f + p.complier * (g - g ^ (1 / alpha))
    dhdt <- deriv(h)
    if (max(diff(h) / h[1 : (length(h) - 1)], na.rm = T) <= 0.01) {
      dhdt[dhdt >= 0] <- max(dhdt[dhdt < 0], na.rm = T) / 100
      loglikelihood <-
        sum(log(h[event.time[Z == 1 & event == 0]]), na.rm = T) +
        sum(log(-dhdt[event.time[Z == 1 & event == 1]]), na.rm = T) +
        log(dgamma(alpha, 5, 4))
    } else {
      loglikelihood <- -Inf
    }
    loglikelihood
  }
  
  alphas <- c(1 : 500) / 50
  loglikelihoods <- array(NA, length(alphas))
  for (i in 1 : length(alphas)) {
    loglikelihoods[i] <- log.likelihood(alphas[i])
  }
  loglikelihoods[loglikelihoods == Inf] <- max(loglikelihoods[loglikelihoods != Inf], na.rm = TRUE)
  loglikelihoods <- loglikelihoods - max(loglikelihoods, na.rm = TRUE)
  likelihoods <- exp(loglikelihoods)
  
  alpha.ML <- alphas[which.max(likelihoods)]
  # Compute confidence interval
  for (i in seq(min(likelihoods, na.rm = T), max(likelihoods, na.rm = T), length = 500)) {
    area <- sum(likelihoods[likelihoods > i], na.rm = T) / sum(likelihoods, na.rm = T)
    if (area < 0.95) {
      alpha.min <- alphas[min(which(likelihoods > i)) - 1]
      alpha.max <- alphas[max(which(likelihoods > i)) + 1]
      break
    }
  }
  
  list(alpha.ML = alpha.ML, alpha.min = alpha.min, alpha.max = alpha.max,
       alphas = alphas, loglikelihoods = loglikelihoods)
  
}