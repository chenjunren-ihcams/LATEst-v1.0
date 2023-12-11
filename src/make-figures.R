rm(list = ls())
library(survival)
library(survminer)
library(readxl)
# Set path ######################################################################
path <- 'E:/work/project/LATE/LATEst'
datapath <- paste0(path, '/data')
savepath <- paste0(path, '/output')
if (! dir.exists(savepath)) {
  dir.create(savepath,recursive = TRUE)
}
setwd(savepath)
source('../src/LATEst.R') # import "LATEst.func" function
source('../src/simulation-data-func.R') # import "simulation.data.func" function

# Simulation analyses ############################################################
### Figure 1b ####
set.seed(21)
eta <- 1  ###
N <- 2e3  ###
alpha <- 0.5   ###

### Generate simulated data

simulation.data <- simulation.data.func(eta, N, alpha) ###

Z <- simulation.data$Z
D <- simulation.data$D
event <- simulation.data$event
event.time <- simulation.data$event.time

####
df <- data.frame(cbind(event.time, event, Z, D))
## "treatment vs. no treatment" model
model.D <- coxph(Surv(event.time, event) ~ D, data = df)
summary(model.D)
## "before vs. after" model
model.Z <- coxph(Surv(event.time, event) ~ Z, data = df)
summary(model.Z)

####
### Estimate alpha and confidence interval

LATEst.res <- LATEst.func(Z, D, event, event.time) ###

alphas <- LATEst.res$alphas
loglikelihoods <- LATEst.res$loglikelihoods

plot(alphas, loglikelihoods)
plot(alphas, exp(loglikelihoods), 'l')

alpha.ML <- LATEst.res$alpha.ML
alpha.min <- LATEst.res$alpha.min
alpha.max <- LATEst.res$alpha.max

alpha.ML
alpha.min
alpha.max
####
df <- data.frame(cbind(event.time, event, Z, D))
fit1 <- survfit(Surv(event.time, event) ~ Z + D, data = df)

pdf(file = './Fig1-b-surv-v1.pdf', width = 4.5, height = 5)
ggsurvplot(fit1, data = df, ylim = c(0.5, 1), palette = c('blue', 'cyan', 'red'))
dev.off()

write.csv(df, './Fig1-b-21-Z-D-event-event.time-v1.csv', 
          fileEncoding = 'GB2312')

pdf(file = './Fig1-b-likelihood-v1.pdf', width = 4, height = 5)
plot(alphas, exp(loglikelihoods) / sum(exp(loglikelihoods), na.rm = T) * 50, 'l', 
     xlim = c(0, 1), xlab = 'alphas ^ hat', ylab = 'density')
abline(v = c(alpha.min, alpha.max), col = 2)
dev.off()

write.csv(data.frame(alphas = alphas, 
                     likelihoods = exp(loglikelihoods) / 
                       sum(exp(loglikelihoods), na.rm = T) * 50),
          './Fig1-b-21-alphas-likelihoods-v1.csv', 
          fileEncoding = 'GB2312')

# Simulation analyses ############################################################
### Figure 1c ####
set.seed(1)
num <- 100 # run 100 trials for each combination of N, alpha and eta
## Set different parameter values
N.range <- c(500, 1000, 2000, 5000)
eta.range <- c(1, -1)
alpha.range <- c(1 / 2, 1 / 1.5, 1 / 1.3, 1, 1.3, 1.5, 2)

array.func <- function(x) {
  data <- array(NA, dim = c(length(alpha.range), length(N.range)))
  rownames(data) <- round(alpha.range, 2)
  colnames(data) <- N.range
  data
}

for (kk in 1:length(eta.range)) {
  
  alpha.array <- array.func()
  D.array <- array.func()
  Z.array <- array.func()
  
  for (jj in 1:length(N.range)) {
    for (ii in 1:length(alpha.range)) {
      
      eta <- eta.range[kk]
      N <- N.range[jj]
      alpha <- alpha.range[ii]
      
      alpha.ML.value <- array(NA, num)
      alpha.ML.min <- array(NA, num)
      alpha.ML.max <- array(NA, num)
      
      Z.value <- array(NA, num)
      Z.min <- array(NA, num)
      Z.max <- array(NA, num)
      
      D.value <- array(NA, num)
      D.min <- array(NA, num)
      D.max <- array(NA, num)
      
      for (mm in 1 : num) {
        
        simulation.data <- simulation.data.func(eta, N, alpha) ###
        
        Z <- simulation.data$Z
        D <- simulation.data$D
        event <- simulation.data$event
        event.time <- simulation.data$event.time
        ####
        df <- data.frame(cbind(event.time, event, Z, D))
        
        model.D <- coxph(Surv(event.time, event) ~ D, data = df)
        summary(model.D)
        model.D.summary <- summary(model.D)
        D.value[mm] <- model.D.summary$conf.int['D', 'exp(coef)']
        D.min[mm] <- model.D.summary$conf.int['D', 'lower .95']
        D.max[mm] <- model.D.summary$conf.int['D', 'upper .95']
        
        model.Z <- coxph(Surv(event.time, event) ~ Z, data = df)
        model.Z.summary <- summary(model.Z)
        Z.value[mm] <- model.Z.summary$conf.int['Z', 'exp(coef)']
        Z.min[mm] <- model.Z.summary$conf.int['Z', 'lower .95']
        Z.max[mm] <- model.Z.summary$conf.int['Z', 'upper .95']
        
        ####
        try ({
          
          LATEst.res <- LATEst.func(Z, D, event, event.time) ###
          
          alphas <- LATEst.res$alphas
          loglikelihoods <- LATEst.res$loglikelihoods
          
          alpha.ML.value[mm] <- LATEst.res$alpha.ML
          alpha.ML.min[mm] <- LATEst.res$alpha.min
          alpha.ML.max[mm] <- LATEst.res$alpha.max

          
        }, silent = TRUE)
        
        if (mm %% 10 == 0) cat('.')
        
      }
      
      cat(paste0('\n N = ', N, ' ; eta = ', eta, '; alpha = ', alpha,
                 '; alpha.ML = ', mean(alpha.ML.value, na.rm = TRUE),
                 '; alpha.ML.min = ', mean(alpha.ML.min, na.rm = TRUE),
                 '; alpha.ML.max = ', mean(alpha.ML.max, na.rm = TRUE), '\n'
      ))
      
      alpha.array[ii, jj] <- paste0(round(mean(alpha.ML.value, na.rm = T), 2),
                                    '\n (', round(mean(alpha.ML.min, na.rm = T), 2),
                                    ', ', round(mean(alpha.ML.max, na.rm = T), 2), ')')
      
      D.array[ii, jj] <- paste0(round(mean(D.value, na.rm = T), 2),
                                '\n (', round(mean(D.min, na.rm = T), 2),
                                ', ', round(mean(D.max, na.rm = T), 2), ')')
      
      Z.array[ii, jj] <- paste0(round(mean(Z.value, na.rm = T), 2),
                                '\n (', round(mean(Z.min, na.rm = T), 2),
                                ', ', round(mean(Z.max, na.rm = T), 2), ')')
      
    }
  }
  
  write.csv(alpha.array, paste0('./Fig1-c-alpha-array-eta = ', eta, '-v1.csv'), 
            fileEncoding = 'GB2312')
  write.csv(D.array, paste0('./Fig1-c-D-array-eta = ', eta, '-v1.csv'), 
            fileEncoding = 'GB2312')
  write.csv(Z.array, paste0('./Fig1-c-Z-array-eta = ', eta, '-v1.csv'), 
            fileEncoding = 'GB2312')
  
}

# Simulated data: one-dimensional subject co-variate #############################
###Figure 2 ####
set.seed(1)
num <- 100 
N.range <- c(500, 1000, 2000, 5000, 8000, 10000)
eta.range <- c(1, -1)
alpha.range <- c(1 / 2, 1, 2)

for (kk in 1:length(eta.range)) {
  for (jj in 1:length(alpha.range)) {
    df.plot <- data.frame(value = NA, min = NA, max = NA,
                          N = NA, model = NA)
    
    for (ii in 1:length(N.range)) {
      eta <- eta.range[kk]
      alpha <- alpha.range[jj]
      N <- N.range[ii]
      
      alpha.ML.value <- array(NA, num)
      alpha.ML.min <- array(NA, num)
      alpha.ML.max <- array(NA, num)
      
      Z.value <- array(NA, num)
      Z.min <- array(NA, num)
      Z.max <- array(NA, num)
      
      D.value <- array(NA, num)
      D.min <- array(NA, num)
      D.max <- array(NA, num)
      
      for (mm in 1 : num) {
        
        simulation.data <- simulation.data.func(eta, N, alpha) ###
        
        Z <- simulation.data$Z
        D <- simulation.data$D
        event <- simulation.data$event
        event.time <- simulation.data$event.time
        
        ####
        df <- data.frame(cbind(event.time, event, Z, D))
        
        model.D <- coxph(Surv(event.time, event) ~ D, data = df)
        summary(model.D)
        model.D.summary <- summary(model.D)
        D.value[mm] <- model.D.summary$conf.int['D', 'exp(coef)']
        D.min[mm] <- model.D.summary$conf.int['D', 'lower .95']
        D.max[mm] <- model.D.summary$conf.int['D', 'upper .95']
        
        model.Z <- coxph(Surv(event.time, event) ~ Z, data = df)
        model.Z.summary <- summary(model.Z)
        Z.value[mm] <- model.Z.summary$conf.int['Z', 'exp(coef)']
        Z.min[mm] <- model.Z.summary$conf.int['Z', 'lower .95']
        Z.max[mm] <- model.Z.summary$conf.int['Z', 'upper .95']
        ####
        try ({
          
          LATEst.res <- LATEst.func(Z, D, event, event.time) ###
          
          alphas <- LATEst.res$alphas
          loglikelihoods <- LATEst.res$loglikelihoods
          
          alpha.ML.value[mm] <- LATEst.res$alpha.ML
          alpha.ML.min[mm] <- LATEst.res$alpha.min
          alpha.ML.max[mm] <- LATEst.res$alpha.max
          
        }, silent = TRUE)
        
        if (mm %% 10 == 0) cat('.')
        
      }
      
      cat(paste0('\n N = ', N, ' ; eta = ', eta, '; alpha = ', alpha,
                 '; alpha.ML = ', mean(alpha.ML.value, na.rm = TRUE),
                 '; alpha.ML.min = ', mean(alpha.ML.min, na.rm = TRUE),
                 '; alpha.ML.max = ', mean(alpha.ML.max, na.rm = TRUE), '\n'
      ))
      
      df.plot[(ii - 1) * 3 + 1, ] <- c(round(mean(D.value, na.rm = T), 2),
                                       round(mean(D.min, na.rm = T), 2),
                                       round(mean(D.max, na.rm = T), 2),
                                       N, 'D')
      
      df.plot[(ii - 1) * 3 + 2, ] <- c(round(mean(Z.value, na.rm = T), 2),
                                       round(mean(Z.min, na.rm = T), 2),
                                       round(mean(Z.max, na.rm = T), 2),
                                       N, 'Z')
      
      df.plot[(ii - 1) * 3 + 3, ] <- c(round(mean(alpha.ML.value, na.rm = T), 2),
                                       round(mean(alpha.ML.min, na.rm = T), 2),
                                       round(mean(alpha.ML.max, na.rm = T), 2),
                                       N, 'LATEst')
      
    }
    df.plot$model <- factor(df.plot$model, levels = c('a', 'b', 'D', 'Z', 'LATEst', 'c', 'd'))
    df.plot$N <- factor(df.plot$N)
    df.plot$value <- as.numeric(df.plot$value)
    df.plot$min <- as.numeric(df.plot$min)
    df.plot$max <- as.numeric(df.plot$max)
    
    write.csv(df.plot, paste0('./Fig2-asy-', alpha, '-', eta, '-v1.csv'),
              fileEncoding = 'GB2312')
    
    boxp.plot.func <- function(eta, alpha){
      rb <- plot(-100, -100, xlim = c(1, length(N.range) * 7), xaxt = 'n',
                 xlab = 'N', ylab = 'expected 95-percent CI', 
                 ylim = c(floor(min(df.plot$min) * 2) / 2,
                          ceiling(max(df.plot$max) * 2) / 2),
                 main = paste0('eta = ', eta, '\n alpha = ', alpha))
      xi <- rep(c(3, 4, 5), length(N.range)) + rep((1:length(N.range) - 1) * 7, each = 3)
      points(xi, df.plot$value, col = rep(c("cyan", 'blue', 'red'), 3), pch = 20, cex = 2)
      arrows(xi, df.plot$min, xi, df.plot$max,
             code = 3, col = rep(c("cyan", 'blue', 'red'), 3), angle = 90, 
             length = .1, lwd = 2)
      abline(h = alpha, lty = 2, col = 1)
      axis(1, at = 4 + (1:length(N.range) - 1) * 7, labels = N.range)
      legend("topright", c('D', 'Z', 'LATEst'), col = c("cyan", 'blue', 'red'),
             lty = c(1, 1, 1), pch = c(20, 20, 20))
    }
    
    pdf(paste0('./Fig2-asy-', alpha, '-', eta, '-v1.pdf'), width = 5, height = 5)
    boxp.plot.func(eta, alpha)
    dev.off()
    
  }
}

# Simulation analyses ############################################################
###Figure 3a ####
set.seed(1)
eta <- 1 ###
alpha <-  1 / 2   ###
num <- 10
N.range <- c(100, 200, 300, 500, 1000, 1500, 2000, 2500, 3000, 5000)
run.time <- array(NA, length(N.range))
for (iii in 1:length(N.range)) {
  
  N <- N.range[iii]
  
  run.time.1 <- array(NA, num)
  
  for (jjj in 1:num) {
    
    simulation.data <- simulation.data.func(eta, N, alpha) ###
    
    Z <- simulation.data$Z
    D <- simulation.data$D
    event <- simulation.data$event
    event.time <- simulation.data$event.time
    ####
    t1 <- Sys.time()
    
    LATEst.res <- LATEst.func(Z, D, event, event.time) ###
    
    run.time.1[jjj] <- difftime(Sys.time(), t1, units = 's')
  }
  run.time[iii] <- mean(run.time.1)
  
  print(N)
}

plot(N.range, run.time, ylim = c(0, 0.6), pch = 20,
     xlab = 'N', ylab = 'LATEst run time (sec)')

write.csv(data.frame(N = N.range, run.time = run.time), 
          './Fig3-a-LATEst-run-time-v1.csv', fileEncoding = 'GB2312')

pdf('./Fig3-a-LATEst-run-time-v1.pdf', width = 5, height = 5)
plot(N.range, run.time, ylim = c(0, 0.6), pch = 20,
     xlab = 'N', ylab = 'LATEst run time (sec)')
dev.off()

# Causal effect of intravenous catheter type on line durability #####
### Figure 4b, 4e ####
set.seed(1)
## Import PICC dataset
df <- read_xlsx(paste0(datapath, '/PICC dataset.xlsx')) #####
Z <- df$Zi
D <- df$Di
event <- df$`Line failure`
event.time <- df$`Time (week)`

####
df <- data.frame(cbind(event.time, event, Z, D))
## "treatment vs. no treatment" model
model.D <- coxph(Surv(event.time, event) ~ D, data = df)
summary(model.D)
## "before vs. after" model
model.Z <- coxph(Surv(event.time, event) ~ Z, data = df)
summary(model.Z)
####
LATEst.res <- LATEst.func(Z, D, event, event.time) ###

alphas <- LATEst.res$alphas
loglikelihoods <- LATEst.res$loglikelihoods

plot(alphas, loglikelihoods)
plot(alphas, exp(loglikelihoods), 'l')

alpha.ML <- LATEst.res$alpha.ML
alpha.min <- LATEst.res$alpha.min
alpha.max <- LATEst.res$alpha.max

alpha.ML
alpha.min
alpha.max
####
fit1 <- survfit(Surv(event.time, event) ~ Z + D, data = df)

pdf(file = './Fig4-b-picc-surv-v1.pdf', width = 4.5, height = 5)
ggsurvplot(fit1, data = df, ylim = c(0.7, 1), palette = c('blue', 'cyan', 'red'))
dev.off()

pdf(file = './Fig4-e-picc-likelihood-v1.pdf', width = 4, height = 5)
plot(alphas, exp(loglikelihoods) / sum(exp(loglikelihoods), na.rm = T) * 50, 'l', 
     xlim = c(0, 5), ylim = c(0, 0.7), xlab = 'alphas ^ hat', ylab = 'density')
abline(v = c(alpha.min, alpha.max), col = 2)
dev.off()

write.csv(data.frame(alphas = alphas, 
                     likelihoods = exp(loglikelihoods) / 
                       sum(exp(loglikelihoods), na.rm = T) * 50),
          './Fig4-e-picc-alphas-likelihoods-v1.csv', fileEncoding = 'GB2312')


# Causal effect of a new transplant centre on patient survival#################
### Figure 5b, 5e ####
set.seed(1)
## Import new transplant centre dataset
df <- read_xlsx(paste0(datapath, '/New transplant centre dataset.xlsx')) #####
Z <- df$Zi
D <- df$Di
event <- df$Death
event.time <- df$`Time (week)`
####
df <- data.frame(cbind(event.time, event, Z, D))
## "treatment vs. no treatment" model
model.D <- coxph(Surv(event.time, event) ~ D, data = df)
summary(model.D)
## "before vs. after" model
model.Z <- coxph(Surv(event.time, event) ~ Z, data = df)
summary(model.Z)
####
LATEst.res <- LATEst.func(Z, D, event, event.time) ###

alphas <- LATEst.res$alphas
loglikelihoods <- LATEst.res$loglikelihoods

plot(alphas, loglikelihoods)
plot(alphas, exp(loglikelihoods), 'l')

alpha.ML <- LATEst.res$alpha.ML
alpha.min <- LATEst.res$alpha.min
alpha.max <- LATEst.res$alpha.max

alpha.ML
alpha.min
alpha.max
####
fit1 <- survfit(Surv(event.time, event) ~ Z + D, data = df)

pdf(file = './Fig5-b-new-transplant-centre-surv-v1.pdf', width = 4.5, height = 5)
ggsurvplot(fit1, data = df, ylim = c(0.7, 1), palette = c('blue', 'cyan', 'red'))
dev.off()

pdf(file = './Fig5-e-new-transplant-centre-likelihood-v1.pdf', width = 4, height = 5)
plot(alphas, exp(loglikelihoods) / sum(exp(loglikelihoods), na.rm = T) * 50, 'l', 
     xlim = c(0, 2), ylim = c(0, 1.4), xlab = 'alphas ^ hat', ylab = 'density')
abline(v = c(alpha.min, alpha.max), col = 2)
dev.off()

write.csv(data.frame(alphas = alphas, 
                     likelihoods = exp(loglikelihoods) / 
                       sum(exp(loglikelihoods), na.rm = T) * 50), 
          './Fig5-e-new-transplant-centre-alphas-likelihoods-v1.csv', 
          fileEncoding = 'GB2312')

