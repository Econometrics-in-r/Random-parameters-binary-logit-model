rm(list = ls())
require(stats4)
require(maxLik)
require(randtoolbox)

# reading and storing data in a dataframe
dataset <- read.csv(file.choose(),header=T)
#Sample Size
N <- nrow(dataset) 

#Dependent variable: choice variable
DVar <- dataset$HH_USE
# choice indicator: currently for binary, adjust if necessary
I <- DVar

# Halton Draws 
preparedraws=function()
{
  draws1[,1:dimensions]<<- qnorm(draws1[,1:dimensions])
}

Ndraws=50      # set number of draws 
dimensions=2    # define number of random parameters in the type 1

# generate draws (using Halton)
draws1=as.matrix(halton(Ndraws*N,dimensions))

# preparing draws for estimation - this may take a while
preparedraws()

## data preparation
# separating the variables with fixed parameters 
dataF1 =  as.matrix(data.frame(dataset$AGE,dataset$MOBPERC))
# separating the variables with random parameters 
dataR1 = as.matrix(data.frame(dataset$CAMERAS,dataset$GPS))

dataR1H=NULL
for(i in 1:Ndraws){
  dataR1H=rbind(dataR1H,dataR1)
}

draws1 = draws1[,1:dimensions]

# Likelihood function
LL <- function(params){ 
  
  Fbeta1 <- params[1:2] # Fixed parameters in the mean Function
  MRbeta1 <- params[3:4]  # Mean of Random parameters in the mean function
  SDRbeta1 <- params[5:6]  # Std of Random parameters in the mean function
  
  #vector of indipendent variables with fixed parameters
  offset1 = rep.int(dataF1%*%as.matrix(Fbeta1,ncol=1),Ndraws)
  #simulating random parameters from their means and standard deviation
  beta1 = t( t(draws1)*SDRbeta1 + MRbeta1)
  # constructing the utility function of the chosen alternative
  u1 <- offset1+rowSums(dataR1H*beta1)
  # probability of the chosen alternative
  prob1 <- exp(u1)/(1+exp(u1))
  # simulated loglikelihood 
  PR1 <- (rowMeans(matrix(prob1, ncol = Ndraws)))
  PR0 <- 1-PR1
  loglik <- sum(I*log(PR1)+((1-I)*log(PR0)))
  
  return(loglik)
}

init <- c(0.1,0.1,#fixed parameters1
          0.1,0.1,#mean of random parameters1
          1,1)#standard deviation of random parameters1

# optimization (maximization of likelihood function)
fit <- maxLik(LL,start=init,method="BFGS")

summary(fit)


