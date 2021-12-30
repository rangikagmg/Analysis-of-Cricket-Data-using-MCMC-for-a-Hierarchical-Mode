betaBinomial.mcmc <- function(x, n, mu.start, gamma.start, burn.in = 1000, n.draws = 5000, sigma.mu = 0.005, sigma.gamma = 0.001) {
 
  m  = function(mu, gamma, x, n) {
    N = length(x)
    l = sum(lbeta(mu*(1-gamma)/gamma + x, (1-mu)*(1-gamma)/gamma+n-x)) - N*lbeta(mu*(1-gamma)/gamma, (1-mu)*(1-gamma)/gamma)
    p = -0.5*log(mu) - 0.5*log(1-mu) - 0.5*log(gamma) - 0.5*log(1-gamma)
    return(l + p)
  }
  
#Create vectors/matrices of zeroes to store MCMC draws
  
  gamma = rep(0, burn.in + n.draws)
  mu = rep(0, burn.in + n.draws)
  theta = matrix(rep(0, length(n)*(burn.in + n.draws)), length(n), (burn.in + n.draws))
  
#Create counter variables to store frequency of acceptance in MH steps
  
  acceptance.mu = 0
  acceptance.gamma = 0
  
#Define starting values for mu and gamma
  
  mu[1] = mu.start
  gamma[1] = gamma.start
  
#Simulate starting values for thetas using initial mu and gamma values
  
  for(i in 1:length(x)) {
    theta[i, 1] = rbeta(1, mu[1]*(1-gamma)[1]/gamma[1] + x[i], (1-gamma)[1]/gamma[1]*(1-mu[1]) + n[i] - x[i])
  }
  
  
#Begin MCMC chain. 
  
  for(j in 2:(burn.in + n.draws)) {

  #Set current chain values equal to previous chain values
    gamma[j] = gamma[j-1]
    mu[j] = mu[j-1]
    
    
#Metropolis-Hastings step for mu
    
    
    candidate = rnorm(1, mu[j-1], sigma.mu)
    
    #Check if candidateidate is between 0 and 1. Otherwise, discard it
    
    if((candidate > 0) & (candidate < 1)) {
      m.old = m(mu[j-1],gamma[j-1],x,n)
      m.new = m(candidate,gamma[j-1],x,n)
      
      
    #Draw an observation from a uniform(0,1)
      
      u = runif(1)
      
      if((m.new - m.old) > log(u)) {
        mu[j] = candidate
        acceptance.mu = acceptance.mu+1
      }
    }
    
    
 #Metropolis-Hastings step for gamma
    
   candidate = rnorm(1,gamma[j-1],sigma.gamma)
    if( (candidate > 0) & (candidate < 1)) {
      m.old = m(mu[j-1],gamma[j-1],x,n)
      m.new = m(mu[j-1],candidate,x,n)
      u = runif(1)
      
      if((m.new - m.old) > log(u)) {
        gamma[j] = candidate
        acceptance.gamma = acceptance.gamma + 1
      } 
 
    }
    
    
#Gibbs sampling step for theta
    
    
    for(i in 1:length(n)) {
      theta[i, j] = rbeta(1, (1-gamma[j])/gamma[j]*mu[j] + x[i], (1-gamma[j])/gamma[j]*(1-mu[j]) + n[i] - x[i])
    } 
    
    
  } #End MCMC Loop
  
  
#Throw out the burn-in iterations
  
  mu <- mu[(burn.in + 1):(burn.in + n.draws)]
  gamma <- gamma[(burn.in + 1):(burn.in + n.draws)]
  theta <- theta[,(burn.in + 1):(burn.in + n.draws)]
  
  
#Return chain for each parameter and acceptance rates for MH step parameters
  
  return(list(mu = mu, gamma = gamma, theta = theta, acceptance = c(acceptance.mu/(burn.in + n.draws), acceptance.gamma/(burn.in + n.draws))))
  
}

#Read in data
library(readxl)
BattingRecord <- read_excel("C:/Users/Rangika/Desktop/BattingRecord.xlsx")
attach(BattingRecord)
x <- BattingRecord$Runs[BattingRecord$Balls >= 2000]
n <- BattingRecord$Balls[BattingRecord$Balls >= 2000]

length(n)

#Make a histogram

hist(x/n, main = "Histogram of Batting Averages\n (Min 2000 Balls)", xlab = "Observed Batting Average")


#MCMC Function


betaBinomial.mcmc <- function(x, n, mu.start, gamma.start, burn.in = 1000, n.draws = 5000, sigma.mu = 0.005, sigma.gamma = 0.001) {
  
  m  = function(mu, gamma, x, n) {
    N = length(x)
    l = sum(lbeta(mu*(1-gamma)/gamma + x, (1-mu)*(1-gamma)/gamma+n-x)) - N*lbeta(mu*(1-gamma)/gamma, (1-mu)*(1-gamma)/gamma)
    p = -0.5*log(mu) - 0.5*log(1-mu) - 0.5*log(gamma) - 0.5*log(1-gamma)
    return(l + p)
  }

  
  #Create vectors/matrices of zeroes to store MCMC draws
  
  gamma = rep(0, burn.in + n.draws)
  mu = rep(0, burn.in + n.draws)
  theta = matrix(rep(0, length(n)*(burn.in + n.draws)), length(n), (burn.in + n.draws))
  

  #Create counter variables to store frequency of acceptance in MH steps
  
  acceptance.mu = 0
  acceptance.gamma = 0
  
  
  #Define starting values for mu and gamma

  
  mu[1] = mu.start
  gamma[1] = gamma.start
  
  
  #Simulate starting values for thetas using initial mu and gamma values
  
  for(i in 1:length(x)) {
    theta[i, 1] = rbeta(1, mu[1]*(1-gamma)[1]/gamma[1] + x[i], (1-gamma)[1]/gamma[1]*(1-mu[1]) + n[i] - x[i])
  }
  
  
  #Begin MCMC chain. This chain will run for the length of the burn-in period
  
  for(j in 2:(burn.in + n.draws)) {
    
    #Set current chain values equal to previous chain values
    
    gamma[j] = gamma[j-1]
    mu[j] = mu[j-1]
    
    
    #Metropolis-Hastings step for mu
    
    candidate = rnorm(1, mu[j-1], sigma.mu)
    #Check if candidateidate is between 0 and 1. If not, discard it.
    
    if((candidate > 0) & (candidate < 1)) {
      m.old = m(mu[j-1],gamma[j-1],x,n)
      m.new = m(candidate,gamma[j-1],x,n)
      
      #Draw an observation from a uniform(0,1)
      
      u = runif(1)
      
      if((m.new - m.old) > log(u)) {
        mu[j] = candidate
        acceptance.mu = acceptance.mu+1
      }
      
    }

    
    #Metropolis-Hastings step for gamma
    
    candidate = rnorm(1,gamma[j-1],sigma.gamma)
    
    #Check if candidateidate is between 0 and 1. If not, discard it.
    
    if( (candidate > 0) & (candidate < 1)) {
      m.old = m(mu[j-1],gamma[j-1],x,n)
      m.new = m(mu[j-1],candidate,x,n)
      
      #Draw an observation from a uniform(0,1)
  
      u = runif(1)
      
      if((m.new - m.old) > log(u)) {
        gamma[j] = candidate
        acceptance.gamma = acceptance.gamma + 1
      } 

    }
    
    #Gibbs sampling step for theta

    
    for(i in 1:length(n)) {
      theta[i, j] = rbeta(1, (1-gamma[j])/gamma[j]*mu[j] + x[i], (1-gamma[j])/gamma[j]*(1-mu[j]) + n[i] - x[i])
    } 
    
  
  } #End MCMC Loop
  
  
  
  #Throw out the burn-in iterations

  mu <- mu[(burn.in + 1):(burn.in + n.draws)]
  gamma <- gamma[(burn.in + 1):(burn.in + n.draws)]
  theta <- theta[,(burn.in + 1):(burn.in + n.draws)]
  
  #Return chain for each parameter and acceptance rates for MH step parameters
  
  return(list(mu = mu, gamma = gamma, theta = theta, acceptance = c(acceptance.mu/(burn.in + n.draws), acceptance.gamma/(burn.in + n.draws))))
}

#Fitting the model

set.seed(5)

chain.1 <- betaBinomial.mcmc(x,n, 0.265, 0.002)
chain.2 <- betaBinomial.mcmc(x,n, 0.5, 0.1)
chain.3 <- betaBinomial.mcmc(x,n, 0.100, 0.0001)


chain.1$acceptance
chain.2$acceptance
chain.3$acceptance

par(mfrow=c(1,1))

matplot(data.frame(chain.1$mu, chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")
matplot(data.frame(chain.1$gamma, chain.2$gamma, chain.3$gamma), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "gamma", main = "MCMC Chains for Gamma")
matplot(data.frame(chain.1$theta[1,], chain.2$theta[1,], chain.3$theta[1,]), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Theta 1", main = "MCMC Chains for Theta 1")


#Check acceptance rates for mu and gamma. 

chain.1$acceptance
chain.2$acceptance
chain.3$acceptance


#Inference

mu <- c(chain.1$mu, chain.2$mu, chain.3$mu)
gamma <- c(chain.1$gamma, chain.2$gamma, chain.3$gamma)
theta <- cbind(chain.1$theta, chain.2$theta, chain.3$theta)


#Mean, standard deviation, and 95% interval for mu

hist(mu, xlab = 'Mu', freq=F)
mean(mu)
sd(mu)
quantile(mu,c(.025,.975))



#Mean, standard deviation, and 95% interval for gamma


hist(gamma, xlab = 'Gamma', freq=F, main = 'Histogram of Gamma')
mean(gamma)
sd(gamma)

quantile(gamma,c(.025,.975))


#Means, standard deviations, and 95% intervals for thetas, the

#true battings average of  each player in the sample


names <- BattingRecord$Player[BattingRecord$Balls >= 2000]
data.frame("Player" = names, "Mean" = apply(theta,1,mean), "SD" = apply(theta,1,sd), "Lower 95" = apply(theta,1,quantile, 0.025), "Upper 95" = apply(theta,1,quantile, 0.975))




###########################################

set.seed(5)




#Good mixing


chain.1 <- betaBinomial.mcmc(x,n, 0.265, 0.00157)
chain.2 <- betaBinomial.mcmc(x,n, 0.285, 0.00157)
chain.3 <- betaBinomial.mcmc(x,n, 0.245, 0.00157)

chain.1$acceptance
chain.2$acceptance
chain.3$acceptance

matplot(data.frame(chain.1$mu, chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")

matplot(data.frame(chain.2$mu, chain.3$mu), type = 'l', lty = c(1,2,3), xlab = "Iteration", ylab = "Mu", main = "MCMC Chains for Mu")












 
