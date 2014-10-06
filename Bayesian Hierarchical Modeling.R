setwd("~/Dropbox/School/UC Davis Spring Quarter 2014/STA 145/Project")
library(coda) #For MCMC + effectiveSize
library(plyr) #For making means of lambdas

lambdaSampler <- function(lambdaData, alpha.current, beta.current) {
  lambdaCurrent <- numeric()
  for (i in 1:nrow(lambdaData)) {
    lambdaCurrent[i] <- rgamma(1, shape = (alpha.current + lambdaData[i, 2]), rate = (beta.current + lambdaData[i, 3]))
  }
  lambdaCurrent
}

poisFunc <- function(lambda.current, lambdaData, starData) { #This calculates the poisson density for each data value corresponding to the correct lambda
  #browser()
  counter <- 1
  outputVector <- numeric()
  for (i in 1:nrow(lambdaData)) {
    for (j in 1:lambdaData[i, 3]) {
      outputVector[counter] <- dpois(starData[counter], lambda.current[i], log = TRUE)
      counter <- counter + 1
    }
  }
  outputVector
}

alphaGenerator <- function(alpha.current, a.alpha, b.alpha, beta.current, a.beta, b.beta, lambda.current, lambdaData, starData) { #
  #browser()
  poisSum <- sum(poisFunc(lambda.current, lambdaData, starData))
  returnVal <- dgamma(alpha.current, shape = a.alpha, rate = b.alpha, log = TRUE) + dgamma(beta.current, shape = b.alpha, rate = b.beta, log = TRUE) + sum(dgamma(lambda.current, shape = alpha.current, rate = beta.current, log = TRUE)) + poisSum 
  returnVal
}

projectWorker <- function(lambdaData, starData, nsamples = 10000, burnin = 1000, a.alpha = 0.0, b.alpha = 0.0, a.beta = 0.0, b.beta = 0.0, v.alpha.proposal = 1.0,  verbose = FALSE) {
  #browser()
  
  #Set up storage for the posterior draws:
  alpha.samples <- rep(NA, length = nsamples)
  beta.samples <- rep(NA, length = nsamples)
  lambda.samples <- matrix(NA, nrow = nsamples, ncol = 36)
  
  #Start at non-robust estimates:
  alpha.current = 1
  beta.current = 1
  lambda.current = rep(1, 36)
    
  #Acceptance vectors to store how many proposals are accepted (only store these after the burnin period)
  alpha.accepted = rep(0, nsamples)
  
  # Obtain posterior samples:
  for (i in 1:(burnin + nsamples)) {
    
    #Sample beta
    beta.current <- rgamma(1, shape = (a.beta + alpha.current * 36), rate = b.beta + sum(lambda.current))
    
    #Sample lambdas
    lambda.current <- lambdaSampler(lambdaData, alpha.current, beta.current)
    
    #Sample alpha
    #Update alpha:
    alpha.proposed <- rnorm(1, mean = alpha.current, sd = sqrt(v.alpha.proposal))
    if (alpha.proposed < 0) { #SAME CONCEPT AS HW.  ALPHA MUST BE GREATER THAN 0.
      alpha.current <- alpha.current
    }
    else {
      log.alpha.proposed <- alphaGenerator(alpha.proposed, a.alpha, b.alpha, beta.current, b.alpha, b.beta, lambda.current, lambdaData, starData)
      log.alpha.current <- alphaGenerator(alpha.current, a.alpha, b.alpha, beta.current, b.alpha, b.beta, lambda.current, lambdaData, starData)
      log.U <- log(runif(1))
      if (log.U < log.alpha.proposed - log.alpha.current) {
        alpha.current <- alpha.proposed
        if (i > burnin) {
          alpha.accepted[i - burnin] <- 1
        }
      }
    }
    
    # Store samples:
    if (i > burnin) {
      alpha.samples[i - burnin] <- alpha.current
      beta.samples[i - burnin] <- beta.current
      lambda.samples[i - burnin, ] <- lambda.current
    }
    
    
    if (verbose && i%%1000 == 0) {
      cat(paste0("Finished iteration ", i, "...\n"))
    }
    
  }
  
  # Obtain posterior summary statistics:
  "sum.stats" <- function(x) {
    nstats <- 5
    ret <- rep(NA, nstats)
    names(ret) <- c("mean","sd", "median", "q02.5", "q97.5")
    ret[1:2] <- c(mean(x), sd(x))
    ret[3:5] <- quantile(x, probs = c(0.50, 0.025, 0.975))
    return(ret)
  }
  # Combine the posterior samples Note: bad memory management but lets do it for conceptual simplicity... :)
  posterior.samples <- cbind(alpha.samples, beta.samples, lambda.samples)
  colnames(posterior.samples) <- c("alpha.samples", "beta.samples", paste("lambda_", 1:36, sep = ""))
  if (verbose) {
    cat("First few posterior samples:\n")
    print(head(posterior.samples))
    cat("Last few posterior samples:\n")
    print(tail(posterior.samples))
    cat("Posterior summaries:\n")
    print(apply(posterior.samples, 2, summary))
    cat("Acceptance rates:\n")
    print(c("alpha" = mean(alpha.accepted)))
    cat("Effective Sample Sizes:\n")
    print(effectiveSize(posterior.samples))
  }
  posterior.stats <- apply(posterior.samples, 2, sum.stats)
  colnames(posterior.stats) <- c("alpha.samples", "beta.samples", paste("lambda_", 1:36, sep = ""))
  
  return(list("posterior.stats" = t(posterior.stats),
              "posterior.samples" = mcmc(posterior.samples),
              "acc.rates" = c("alpha" = mean(alpha.accepted))))
}

stars <- read.table("stars.txt", header = TRUE)
lambdaSum <- ddply(stars, .(galaxy), summarize, sum(photon.count)) #Counts the sum Y_{i, j} for each galaxy (lambda)
lambdaCount <- ddply(stars, .(galaxy), summarize, length(photon.count)) #Finds the length (n_{j}) for each galaxy (lambda)
names(lambdaSum) <- c("lambdaJ", "sum")
names(lambdaCount) <- c("lambdaJ", "n")
lambdaData <- cbind(lambdaSum, lambdaCount$n)
names(lambdaData) <- c("lambdaJ", "sum", "n")

output <- projectWorker(lambdaData, stars$photon.count, nsamples = 20000, burnin = 2000, a.alpha = 1, b.alpha = 1, a.beta = 1, b.beta = 1, v.alpha.proposal = 0.5, verbose = TRUE)
posterior.samples <- output[[2]]
mcmc.draws <- mcmc(posterior.samples)
plot(mcmc.draws)
par(mfrow=c(1,1))

posterior.stats <- output[[1]]
getStatistics <- function(photonCounts) {
  outputVec <- numeric(5)
  outputVec[1] <- max(photonCounts)
  outputVec[2] <- min(photonCounts)
  outputVec[3] <- median(photonCounts)
  outputVec[4] <- mean(photonCounts)
  outputVec[5] <- sd(photonCounts)
  names <- c("Max", "Min", "Median", "Mean", "SD")
  outputData <- data.frame(names, outputVec)
  outputData
}

pValueGenerator <- function(lambdaStats, n, lambda) {
  sampleHolder <- matrix(NA, nrow = 100, ncol = n)
  for (i in 1:100) {
    sampleHolder[i, ] <- t(rpois(n = n, lambda = lambda))
  }
  statHolder <- matrix(NA, nrow = 100, ncol = 5) 
  for (i in 1:100) {
    statHolder[i, ] <- t(getStatistics(sampleHolder[i, ])[, 2])
  }
  ppc.pvalues <- numeric()
  for (i in 1:length(lambdaStats)) {
    ppc.pvalues[i] = mean(statHolder[, i] <= lambdaStats[i])
  }
  ppc.pvalues
}

par(mfrow=c(3,2))
###lambdaOne###
lambdaOne <- stars[stars$galaxy == 1, ]
truehist(lambdaOne[, 2], xlab = c("Photon Counts"), main = c("Lambda1 Data"), nbins = 15)
lambdaOnePPD <- rpois(n = nrow(lambdaOne), lambda = posterior.stats[3, 1])
truehist(lambdaOnePPD, xlab = c("Photon Counts"), main = c("PPD Lambda1 Data"), nbins = 15)
lambdaOneStats <- getStatistics(lambdaOne[, 2])
lambdaOnePPDStats <- getStatistics(lambdaOnePPD)
pValueGenerator(lambdaOneStats[, 2], nrow(lambdaOne), posterior.stats[3, 1])

###lambdaNine###
lambdaNine <- stars[stars$galaxy == 9, ]
truehist(lambdaNine[, 2], xlab = c("Photon Counts"), main = c("Lambda9 Data"), nbins = 15)
lambdaNinePPD <- rpois(n = nrow(lambdaNine), lambda = posterior.stats[11, 1])
truehist(lambdaNinePPD, xlab = c("Photon Counts"), main = c("PPD Lambda9 Data"), nbins = 15)
lambdaNineStats <- getStatistics(lambdaNine[, 2])
lambdaNinePPDStats <- getStatistics(lambdaNinePPD)
pValueGenerator(lambdaNineStats[, 2], nrow(lambdaNine), posterior.stats[11, 1])

###lambda15###
lambda15 <- stars[stars$galaxy == 15, ]
truehist(lambda15[, 2], xlab = c("Photon Counts"), main = c("Lambda15 Data"), nbins = 15)
lambda15PPD <- rpois(n = nrow(lambda15), lambda = posterior.stats[17, 1])
truehist(lambda15PPD, xlab = c("Photon Counts"), main = c("PPD Lambda15 Data"), nbins = 15)
lambda15Stats <- getStatistics(lambda15[, 2])
lambda15PPDStats <- getStatistics(lambda15PPD)
pValueGenerator(lambda15Stats[, 2], nrow(lambda15), posterior.stats[17, 1])

par(mfrow=c(2,2))
###lambda18###
lambda18 <- stars[stars$galaxy == 18, ]
truehist(lambda18[, 2], xlab = c("Photon Counts"), main = c("Lambda18 Data"), nbins = 15)
lambda18PPD <- rpois(n = nrow(lambda18), lambda = posterior.stats[20, 1])
truehist(lambda18PPD, xlab = c("Photon Counts"), main = c("PPD Lambda18 Data"), nbins = 15)
lambda18Stats <- getStatistics(lambda18[, 2])
lambda18PPDStats <- getStatistics(lambda18PPD)
pValueGenerator(lambda18Stats[, 2], nrow(lambda18), posterior.stats[20, 1])

###lambda35###
lambda35 <- stars[stars$galaxy == 35, ]
truehist(lambda35[, 2], xlab = c("Photon Counts"), main = c("Lambda35 Data"), nbins = 15)
lambda35PPD <- rpois(n = nrow(lambda35), lambda = posterior.stats[37, 1])
truehist(lambda35PPD, xlab = c("Photon Counts"), main = c("PPD Lambda35 Data"), nbins = 15)
lambda35Stats <- getStatistics(lambda35[, 2])
lambda35PPDStats <- getStatistics(lambda35PPD)
pValueGenerator(lambda35Stats[, 2], nrow(lambda35), posterior.stats[37, 1])

par(mfrow=c(2,1))
plot(density(posterior.stats[3:nrow(posterior.stats), 1]), main = "Posterior Lambda Density")
lambdaMean <- lambdaData[, 2] / lambdaData[, 3]
plot(density(lambdaMean), main = "Data Lambda Density")