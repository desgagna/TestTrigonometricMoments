# ====================================================================================
# Script : Test_Trigo_File2.R
# Purpose: Functions for the article:
#          "An omnibus goodness-of-fit test based on trigonometric moments"
# Note   : This file contains only function definitions.
#          See Test_Trigo_File1.R for usage examples.
# Author : Alain Desgagne 
# Date   : 2005-07-24
# ====================================================================================

## ---- Constants h1 to h37

h1 <-  function(lambda)      adaptIntegrate(f = function(v) cos(pi * (1 + pgamma(v, shape = 1 / lambda, scale = 1))) * dgamma(v, shape = 1 / lambda + 1, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral
h2 <-  function(lambda)      adaptIntegrate(f = function(v) sin(pi * (1 + pgamma(v, shape = 1 / lambda, scale = 1))) * dgamma(v, shape = 1, scale = 1) , lower = 0, upper = Inf, tol = 1e-9)$integral
h3 <-  function(lambda)      adaptIntegrate(f = function(v) cos(pi * (1 + pgamma(v, shape = 1 / lambda, scale = 1))) * log(lambda * v) * dgamma(v, shape = 1 / lambda + 1, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral
h4 <-  function(lambda)      adaptIntegrate(f = function(v) cos(pi * (1 + pgamma(v, shape = 1 / lambda, scale = 1))) * dgamma(v, shape = 3 / lambda, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral
h5 <-  function(lambda)      adaptIntegrate(f = function(v) sin(pi * (1 + pgamma(v, shape = 1 / lambda, scale = 1))) * dgamma(v, shape = 2 / lambda, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral 
h6 <-  function(a,b,c)       adaptIntegrate(f = function(v) cos(2 * pi * pgamma(v, shape = a, scale = 1)) * dgamma(v, shape = b, scale = c), lower = 0, upper = Inf, tol = 1e-9)$integral
h7 <-  function(a,b,c)       adaptIntegrate(f = function(v) sin(2 * pi * pgamma(v, shape = a, scale = 1)) * dgamma(v, shape = b, scale = c), lower = 0, upper = Inf, tol = 1e-9)$integral
h8 <-  function(lambda)      adaptIntegrate(f = function(v) (v - lambda) * log(v) * cos(2 * pi * pgamma(v, shape = lambda, scale = 1)) * dgamma(v, shape = lambda, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral
h9 <-  function(lambda)      adaptIntegrate(f = function(v) (v - lambda) * log(v) * sin(2 * pi * pgamma(v, shape = lambda, scale = 1)) * dgamma(v, shape = lambda, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral
h10 <- function(alpha)       adaptIntegrate(f = function(v) log(v) * cos(2 * pi * pgamma(v, shape = alpha, scale = 1)) * dgamma(v, shape = alpha, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral 
h11 <- function(alpha)       adaptIntegrate(f = function(v) log(v) * sin(2 * pi * pgamma(v, shape = alpha, scale = 1)) * dgamma(v, shape = alpha, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral 
h12 <- function(lambda)      adaptIntegrate(f = function(v) cos(pi * (2 - pbeta(v, lambda / 2, 1 / 2))) * dbeta(v, lambda / 2, 3 / 2), lower = 0, upper = 1, tol = 1e-9)$integral 
h13 <- function(lambda)      adaptIntegrate(f = function(v) sin(pi * (2 - pbeta(v, lambda / 2, 1 / 2))) * dbeta(v, (lambda + 1) / 2, 1), lower = 0, upper = 1, tol = 1e-9)$integral 
h14 <- function(lambda)      adaptIntegrate(f = function(v) cos(pi * (2 - pbeta(v, lambda / 2, 1 / 2))) * (log(v) + (lambda + 1) / lambda * (1 - v)) * dbeta(v, lambda / 2, 1 / 2), lower = 0, upper = 1, tol = 1e-9)$integral 
h15 <- function(lambda)      adaptIntegrate(f = function(v) cos(pi * (2 - pbeta(v, lambda / 2, 1 / 2))) * dbeta(v, (lambda - 1) / 2, 1), lower = 0, upper = 1, tol = 1e-9)$integral
h16 <- function(lambda)      adaptIntegrate(f = function(v) sin(pi * (2 - pbeta(v, lambda / 2, 1 / 2))) * dbeta(v, (lambda - 1) / 2, 1), lower = 0, upper = 1, tol = 1e-9)$integral
h17 <- function(lambda)      adaptIntegrate(f = function(v) cos(2 * pi * pgamma(v, shape = 1 / lambda, scale = 1)) * log(lambda * v) * dgamma(v, shape = 1 / lambda + 1, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral
h18 <- function(lambda)      adaptIntegrate(f = function(v) sin(2 * pi * pgamma(v, shape = 1 / lambda, scale = 1)) * log(lambda * v) * dgamma(v, shape = 1 / lambda + 1, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral
h19 <- function(rho)         adaptIntegrate(f = function(v) (log(v)) ^ 2 * v * exp(- rho * v), lower = 1, upper = Inf, tol = 1e-9)$integral 
h20 <- function(rho)         adaptIntegrate(f = function(v) log(v) * v * exp(- rho * v), lower = 1, upper = Inf, tol = 1e-9)$integral 
h21 <- function(rho)         adaptIntegrate(f = function(v) cos(2 * pi * (1 - exp(- rho * (v - 1)))) * log(v) * (1 - rho * v) * exp(- rho * v), lower = 1, upper = Inf, tol = 1e-9)$integral 
h22 <- function(rho)         adaptIntegrate(f = function(v) sin(2 * pi * (1 - exp(- rho * (v - 1)))) * log(v) * (1 - rho * v) * exp(- rho * v), lower = 1, upper = Inf, tol = 1e-9)$integral 
h23 <- function(rho)         adaptIntegrate(f = function(v) cos(2 * pi * (1 - exp(- rho * (v - 1)))) * v * exp(- rho * v), lower = 1, upper = Inf, tol = 1e-9)$integral 
h24 <- function(rho)         adaptIntegrate(f = function(v) sin(2 * pi * (1 - exp(- rho * (v - 1)))) * v * exp(- rho * v), lower = 1, upper = Inf, tol = 1e-9)$integral 
h25 <- function(alpha, beta) adaptIntegrate(f = function(v) log(v) * cos(2 * pi * pbeta(v, alpha, beta)) * dbeta(v, alpha, beta), lower = 0, upper = 1, tol = 1e-9)$integral 
h26 <- function(alpha, beta) adaptIntegrate(f = function(v) log(v) * sin(2 * pi * pbeta(v, alpha, beta)) * dbeta(v, alpha, beta), lower = 0, upper = 1, tol = 1e-9)$integral 
h27 <- function(alpha, beta) adaptIntegrate(f = function(v) log(1 - v) * cos(2 * pi * pbeta(v, alpha, beta)) * dbeta(v, alpha, beta), lower = 0, upper = 1, tol = 1e-9)$integral 
h28 <- function(alpha, beta) adaptIntegrate(f = function(v) log(1 - v) * sin(2 * pi * pbeta(v, alpha, beta)) * dbeta(v, alpha, beta), lower = 0, upper = 1, tol = 1e-9)$integral 
h29 <- function(mu, lambda)  adaptIntegrate(f = function(v) v  * cos(2 * pi * pinvgauss(v, mean = mu, shape = lambda)) * dinvgauss(v, mean = mu, shape = lambda), lower = 0, upper = Inf, tol = 1e-9)$integral
h30 <- function(mu, lambda)  adaptIntegrate(f = function(v) v  * sin(2 * pi * pinvgauss(v, mean = mu, shape = lambda)) * dinvgauss(v, mean = mu, shape = lambda), lower = 0, upper = Inf, tol = 1e-9)$integral
h31 <- function(mu, lambda)  adaptIntegrate(f = function(v) (v ^ 2 + mu ^ 2) / v * cos(2 * pi * pinvgauss(v, mean = mu, shape = lambda)) * dinvgauss(v, mean = mu, shape = lambda), lower = 0, upper = Inf, tol = 1e-9)$integral
h32 <- function(mu, lambda)  adaptIntegrate(f = function(v) (v ^ 2 + mu ^ 2) / v * sin(2 * pi * pinvgauss(v, mean = mu, shape = lambda)) * dinvgauss(v, mean = mu, shape = lambda), lower = 0, upper = Inf, tol = 1e-9)$integral
h33 <- function(beta)        adaptIntegrate(f = function(v) cos(2 * pi * (1 - (1 - v) ^ beta)) * log(v) * (1 - v) ^ (beta - 2) * (1 - beta * v), lower = 0, upper = 1, tol = 1e-9)$integral 
h34 <- function(beta)        adaptIntegrate(f = function(v) sin(2 * pi * (1 - (1 - v) ^ beta)) * log(v) * (1 - v) ^ (beta - 2) * (1 - beta * v), lower = 0, upper = 1, tol = 1e-9)$integral 
h35 <- function(beta)        adaptIntegrate(f = function(v) cos(2 * pi * (1 - (1 - v) ^ beta)) * log(1 - v) * (1 - v) ^ (beta - 1), lower = 0, upper = 1, tol = 1e-9)$integral 
h36 <- function(beta)        adaptIntegrate(f = function(v) sin(2 * pi * (1 - (1 - v) ^ beta)) * log(1 - v) * (1 - v) ^ (beta - 1), lower = 0, upper = 1, tol = 1e-9)$integral 
h37 <- function(lambda)      adaptIntegrate(f = function(v) sin(pi * (1 + pgamma(v, shape = 1 / lambda, scale = 1)))  * dgamma(v, shape = 1 / lambda + 1, scale = 1), lower = 0, upper = Inf, tol = 1e-9)$integral

## ---- Distributions on R

## ---- APD distribution (see section 5.3 local alternative)

dAPD <- function(x, lambda, alpha = 0.5, rho, mu = 0, sigma = 1) {
  y <- (x - mu) / sigma
  delta <- 2 * alpha ^ rho * (1 - alpha) ^ rho / (alpha ^ rho + (1 - alpha) ^ rho)
  A <- (1/2 + sign(y) * (1/2 - alpha)) ^ rho
  f <- rho * (delta / lambda) ^ (1 / rho) / sigma / gamma(1 / rho) * exp(- delta / A / lambda *  abs(y) ^ rho) 
  return(f)
}
pAPD <- function(x, lambda, alpha = 0.5, rho, mu = 0, sigma = 1){
  y <- (x - mu) / sigma
  delta <- 2 * alpha ^ rho * (1 - alpha) ^ rho / (alpha ^ rho + (1 - alpha) ^ rho)
  F <- alpha * (1 - pgamma(delta / lambda * (-pmin(0, y) / alpha) ^ rho, shape = 1 / rho)) +
    (1 - alpha) * pgamma(delta / lambda * (pmax(0, y) / (1 - alpha)) ^ rho, shape = 1 / rho)
  return(F)
}
qAPD <- function(p, lambda, alpha = 0.5, rho, mu = 0, sigma = 1) {
  delta <- 2 * alpha ^ rho * (1 - alpha) ^ rho / (alpha ^ rho + (1 - alpha) ^ rho)
  si <- sign(p - alpha)
  x <- mu + sigma * (si + 1 - 2 * alpha) / 2 * (lambda / delta) ^ (1 / rho) * 
    qgamma(2 * abs(p - alpha) / (1 + si * (1 - 2 * alpha)), 1 / rho) ^ (1 / rho) 
  return(x)
}
qAPD2 <- function(p, lambda, alpha = 0.5, rho, mu = 0, sigma = 1) {
  delta <- 2 * alpha ^ rho * (1 - alpha) ^ rho / (alpha ^ rho + (1 - alpha) ^ rho)
  ind.pos <- sign(p - alpha)
  alpha2 <- 1/2 + ind.pos * (1/2 - alpha) 
  p2 <- 1/2 + ind.pos * (1/2 - p)
  x <- mu + sigma * ind.pos * alpha2 * (lambda / delta) ^ (1 / rho) * qgamma(1 - p2 / alpha2, 1 / rho) ^ (1 / rho) 
  return(x)
}
rAPD <- function(n, lambda, alpha = 0.5, rho, mu = 0, sigma = 1) {
  delta <- 2 * alpha ^ rho * (1 - alpha) ^ rho / (alpha ^ rho + (1 - alpha) ^ rho)
  W <- rgamma(n, shape = 1 / rho)
  V <- rbinom(n, 1, alpha)
  x <- mu + sigma * (lambda * W / delta) ^ (1 / rho) * ((1 - alpha) * (1 - V) - alpha * V)
  return(x)
}
r.central.moment.APD <- function(r, lambda, alpha, rho, sigma) {
  # function to calculate mean((x - mu) ^ r) if x has a APD(lambda, alpha, rho, mu, sigma) distribution
  delta <- 2 * alpha ^ rho * (1 - alpha) ^ rho / (alpha ^ rho + (1 - alpha) ^ rho)
  return(sigma ^ r * gamma((1 + r) / rho) / gamma(1 / rho) * (lambda / delta) ^ (r / rho) * ((1 - alpha) ^ (1 + r) + (-1) ^ r * alpha ^ (1 + r))) 
}  
mean.APD <- function(lambda, alpha, rho, mu, sigma) {
  # function to calculate mean(x) if x has a APD(lambda, alpha, rho, mu, sigma) distribution
  delta <- 2 * alpha ^ rho * (1 - alpha) ^ rho / (alpha ^ rho + (1 - alpha) ^ rho)
  return(mu + sigma * gamma(2 / rho) / gamma(1 / rho) * (lambda / delta) ^ (1 / rho) * (1 - 2 * alpha))
}
mean2.APD <- function(lambda, alpha, rho, mu, sigma) {
  # function to calculate mean(x ^ 2) if x has a APD(lambda, alpha, rho, mu, sigma) distribution
  delta <- 2 * alpha ^ rho * (1 - alpha) ^ rho / (alpha ^ rho + (1 - alpha) ^ rho)
  mm <- mu + sigma * gamma(2 / rho) / gamma(1 / rho) * (lambda / delta) ^ (1 / rho) * (1 - 2 * alpha)
  return(sigma ^ 2 * gamma(3 / rho) / gamma(1 / rho) * (lambda / delta) ^ (2 / rho) * ((1 - alpha) ^ 3 +  alpha ^ 3) - mu ^ 2 + 2 * mu * mm)      
}

## ---- EPD distribution

dEPD <- function(x, lambda, mu = 0, sigma = 1, log = FALSE) {
  y <- (x - mu) / sigma
  log.f <-  - log(2) - log(sigma) - (1 / lambda - 1) * log(lambda) - lgamma(1 / lambda) - abs(y) ^ lambda / lambda 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
dEPD2 <- function(x, lambda, mu = 0, sigma = 1, log = FALSE) {
  y <- (x - mu) / sigma
  log.f <-  - log(2) - log(sigma) - log(lambda) / lambda - lgamma(1 + 1 / lambda) - abs(y) ^ lambda / lambda 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pEPD <- function(x, lambda, mu = 0, sigma = 1){
  y <- (x - mu) / sigma
  F <- (1/2) * (1 + sign(y) * pgamma(abs(y) ^ lambda / lambda, shape = 1 / lambda, scale = 1))
  return(F)
}
qEPD <- function(p, lambda, mu = 0, sigma = 1) {
  x <- mu + sigma  * sign(p - 1 / 2) * lambda ^ (1 / lambda) * qgamma(abs(2 * p - 1), shape = 1 / lambda, scale = 1) ^ (1 / lambda)  
  return(x)
}
qEPD2 <- function(p, lambda, mu = 0, sigma = 1) {
  W <- qgamma(abs(2 * p - 1), shape = 1 / lambda, scale = 1)
  V <- ifelse(p < 1 / 2, 1, 0)
  x <- mu + sigma * 2 * lambda ^ (1 / lambda) * W ^ (1 / lambda) * (1 / 2 -  V)
  return(x)
}
rEPD <- function(n, lambda, mu = 0, sigma = 1) {
  W <- rgamma(n, shape = 1 / lambda, scale = 1)
  V <- rbinom(n, size = 1, prob = 1/2)
  x <- mu + sigma * (2 * V - 1) * lambda ^ (1 / lambda) * W ^ (1 / lambda) 
  return(x)
}
ML.EPD <- function(x, lambda = "unknown", mu = "unknown", sigma = "unknown"){
  if (lambda == "unknown" & mu == "unknown" & sigma == "unknown")
  {
    log_vraisemblance <- function(theta, y) {
      lambda <- exp(theta[1])
      mu <- theta[2]
      sigma <- exp(theta[3])
      somme <- sum(dEPD(y, lambda, mu, sigma, log = TRUE))
      return(-somme) 
    }
    lambda0 <- c(0.5, 1, 1.5, 2)
    mu0 <- rep(median(x), 4)
    sigma0 <- mapply(function(mu, lambda) mean(abs(x - mu) ^ lambda) ^ (1 / lambda), mu0, lambda0)
    est <- matrix(0, nrow = 8, ncol = 4)
    for (i in 1:4){
      res <- optim(par = c(log(lambda0[i]), mu0[i], log(sigma0[i])), fn = log_vraisemblance, y = x, method = "BFGS")
      par <- c(exp(res$par[1]), res$par[2], exp(res$par[3]))
      est[2 * i - 1,] <- c(par, log_vraisemblance(res$par, x))
      res <- optim(par = c(log(lambda0[i]), mu0[i], log(sigma0[i])), fn = log_vraisemblance, y = x, method = "Nelder-Mead")
      par <- c(exp(res$par[1]), res$par[2], exp(res$par[3]))
      est[2 * i,] <- c(par, log_vraisemblance(res$par, x))
    }
    est[est[,1] <= 0.1 | est[,3] <= 0.01 | est[,4] <= 0.01, 4] <- 999999  ## discard some extreme cases 
    pos <- est[,1] > 20
    if (sum(pos) > 0) {
      est[pos,1] <- 20
      g <- function(mu) sum(abs(x - mu) ^ (20 - 1) * sign(x - mu))
      muu <-  uniroot(g, c(min(x) - 1, max(x) + 1), tol = 1e-09)$root
      est[pos,2] <- muu
      est[pos,3] <- (mean((abs(x - muu)) ^ 20)) ^ (1 / 20) 
      est[pos,4] <- log_vraisemblance(c(log(est[pos,1][1]), est[pos,2][1], log(est[pos,3][1])), x)
    }
    best.est = est[which.min(est[,4]),]
    lambda.hat <- best.est[1]
    mu.hat <- best.est[2]
    sigma.hat <- best.est[3]
  }
  if (lambda == "unknown" & mu == "unknown" & sigma != "unknown")
  {
    mu0 <- median(x)
    g0 <- function(lambda)  digamma(1 / lambda + 1) + log(lambda) - 1 - mean((abs(x - mu0) / sigma) ^ lambda * (lambda * log(abs(x - mu0) / sigma) - 1))
    lambda0 <- uniroot(g0, c(.1, 100), tol = 1e-09)$root
    g <- function(vars) {
      lambda <- vars[1]
      mu <- vars[2]
      g1 <- digamma(1 / lambda + 1) + log(lambda) - mean((abs(x - mu) / sigma) ^ lambda * (lambda * log(abs(x - mu) / sigma) - 1)) - 1
      g2 <- mean(abs(x - mu) ^ (lambda - 1) * sign(x - mu))
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, mu0), g)
    lambda.hat <- res$x[1]
    mu.hat <- res$x[2]
    sigma.hat <- sigma
  }
  if (lambda == "unknown" & mu != "unknown" & sigma == "unknown")
  {
    sigma0 <- mad(x)
    g0 <- function(lambda) digamma(1 / lambda + 1) + log(lambda) - 1 - mean((abs(x - mu) / sigma0) ^ lambda * (lambda * log(abs(x - mu) / sigma0) - 1))
    lambda0 <- uniroot(g0, c(.1, 100), tol = 1e-09)$root
    g <- function(vars) {
      lambda <- vars[1]
      sigma <- vars[2]
      g1 <- digamma(1 / lambda + 1) + log(lambda) - mean((abs(x - mu) / sigma) ^ lambda * lambda * log(abs(x - mu) / sigma))
      g2 <- sigma - (mean((abs(x - mu)) ^ lambda)) ^ (1 / lambda)
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, sigma0), g)
    lambda.hat <- res$x[1]
    mu.hat <- mu
    sigma.hat <- res$x[2]
  }
  if (lambda == "unknown" & mu != "unknown" & sigma != "unknown")
  {
    g <- function(lambda) digamma(1 / lambda + 1) + log(lambda)  - mean((abs(x - mu) / sigma) ^ lambda * (lambda * log(abs(x - mu) / sigma) - 1)) - 1
    lambda.hat <- uniroot(g, c(.1, 100), tol = 1e-09)$root
    mu.hat <- mu
    sigma.hat <- sigma
  }
  if (lambda != "unknown")
  {
    if (mu == "unknown") 
    {
      if (lambda == 1) mu.hat <- median(x) else 
        if (lambda == 2) mu.hat <- mean(x) else 
        {
          g <- function(mu) sum(abs(x - mu) ^ (lambda - 1) * sign(x - mu))
          mu.hat <- uniroot(g, c(min(x) - 1, max(x) + 1), tol = 1e-09)$root
        }
    } else
      mu.hat <- mu
    if (sigma == "unknown") 
      sigma.hat <- (mean((abs(x - mu.hat)) ^ lambda)) ^ (1 / lambda) else
        sigma.hat <- sigma  
      lambda.hat <- lambda 
  } 
  return(c(lambda.hat, mu.hat, sigma.hat))
}
ML.EPD.equations <- function(x){
  mu0 <- mean(x)
  sigma0 <- mad(x)
  g0 <- function(lambda) digamma(1 / lambda + 1) + log(lambda) -  mean((abs(x - mu0) / sigma0) ^ lambda * lambda * log(abs(x - mu0) / sigma0)) +
    mean((abs(x - mu0) / sigma0) ^ lambda) - 1          
  lambda0 <- uniroot(g0, c(.1, 100), tol = 1e-09)$root
  g <- function(vars) {
    lambda <- vars[1]
    mu <- vars[2]
    sigma <- (mean((abs(x - mu)) ^ lambda)) ^ (1 / lambda)
    g1 <- digamma(1 / lambda + 1) + log(lambda)  - mean((abs(x - mu) / sigma) ^ lambda * lambda * log(abs(x - mu) / sigma))
    g2 <- mean(abs(x - mu) ^ (lambda - 1) * sign(x - mu))
    c(g1, g2)
  }
  res <- nleqslv(c(lambda0, mu0), g)
  lambda.hat <- res$x[1]
  mu.hat <- res$x[2]
  sigma.hat <- (mean((abs(x - mu.hat)) ^ lambda.hat)) ^ (1 / lambda.hat)
  return(c(lambda.hat, mu.hat, sigma.hat))
}
ML.EPD.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown") 
    {
    if (lambda == 1) mu.hat <- median(x) else 
      if (lambda == 2) mu.hat <- mean(x) else 
      {
      g <- function(mu) sum(abs(x - mu) ^ (lambda - 1) * sign(x - mu))
      mu.hat <- uniroot(g, c(min(x) - 1, max(x) + 1), tol = 1e-09)$root
      }
    } else
    mu.hat <- mu
  if (sigma == "unknown") 
     sigma.hat <- (mean((abs(x - mu.hat)) ^ lambda)) ^ (1 / lambda) else
     sigma.hat <- sigma  
  return(c(mu.hat, sigma.hat))
}
MM.EPD <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown")  
    mu.hat <- mean(x) else 
      mu.hat <- mu
    if (sigma == "unknown"){C2 <- gamma(1 / lambda) / lambda ^ (2 / lambda) / gamma(3 / lambda);
      sigma.hat <- sqrt(C2 * mean((x - mu.hat) ^ 2))} else
        sigma.hat <- sigma  
      return(c(lambda, mu.hat, sigma.hat))
}
EPD.test.ML <- function(x, lambda = "unknown", mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.lambda.mu.sigma <- ML.EPD(x, lambda, mu, sigma)
  lambda.hat <- ML.lambda.mu.sigma[1]
  mu.hat <- ML.lambda.mu.sigma[2]
  sigma.hat <- ML.lambda.mu.sigma[3]
  Fi <- pEPD(x, lambda = lambda.hat, mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dEPD(x, lambda = lambda.hat, mu.hat, sigma = sigma.hat)))
  G <-  matrix(c(1 / lambda.hat ^ 2  * (h1(lambda.hat) - h3(lambda.hat)), 0, h1(lambda.hat) / sigma.hat,  
                 0, h2(lambda.hat) / sigma.hat / lambda.hat ^ (1 / lambda.hat - 1) / gamma(1 / lambda.hat), 0), byrow = T, nrow = 2)
  C1 <- digamma(1/lambda.hat + 1) + log(lambda.hat)
  R <-  matrix(c((1 / lambda.hat ^ 3) * ((1 / lambda.hat + 1) * trigamma(1/lambda.hat + 1) +  C1 ^ 2  - 1), 0, - C1 / sigma.hat / lambda.hat,
                 0, lambda.hat ^ (2 - 2 / lambda.hat) * gamma(2 - 1 / lambda.hat) / sigma.hat ^ 2 / gamma(1 / lambda.hat), 0,  
               - C1 / sigma.hat / lambda.hat, 0, lambda.hat / sigma.hat ^ 2), byrow = T, nrow = 3)      
  U <- (1:3)[c(lambda == "unknown", mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.lambda.mu.sigma = ML.lambda.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
EPD.test.ML.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.mu.sigma <- ML.EPD.2p(x, lambda, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- pEPD(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dEPD(x, lambda = lambda, mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h1(lambda), 
                            h2(lambda) / lambda ^ (1 / lambda - 1) / gamma(1 / lambda), 0), byrow = T, nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(lambda ^ (2 - 2 / lambda) * gamma(2 - 1 / lambda) / gamma(1 / lambda), 0, 0, lambda), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
EPD.test.MM <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  C2 <- gamma(1 / lambda) / lambda ^ (2 / lambda) / gamma(3 / lambda)
  C3 <- gamma(3 / lambda) ^ 2 / (gamma(1 / lambda) * gamma(5 / lambda) - gamma(3 / lambda) ^ 2)
  MM.lambda.mu.sigma <- MM.EPD(x, lambda, mu, sigma)
  mu.hat <- MM.lambda.mu.sigma[2]
  sigma.hat <- MM.lambda.mu.sigma[3]
  Fi <- pEPD(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dEPD(x, lambda = lambda, mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(lambda) / lambda ^ (1 / lambda - 1) / gamma(1 / lambda), h1(lambda), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * diag(c(C2, 4 * C3)) 
  J <- 1 / sigma.hat * matrix(c(0, h5(lambda) * gamma(2 / lambda) / lambda ^ (1 / lambda) / gamma(3 / lambda), 2 * C3 * h4(lambda), 0), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]; J <- as.matrix(J[,U])
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.lambda.mu.sigma = MM.lambda.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
## ---- Laplace distribution

dLaplace <- function(x, mu = 0, sigma = 1, log = FALSE) {
  y <- (x - mu) / sigma
  log.f <-  - log(2) - log(sigma) - abs(y) 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pLaplace <- function(x, mu = 0, sigma = 1){
  y <- (x - mu) / sigma
  F <- (1/2) * (1 + sign(y) * (1 - exp(-abs(y))))
  return(F)
}
qLaplace <- function(p, mu = 0, sigma = 1) {
  W <- qexp(abs(2 * p - 1))
  V <- ifelse(p < 1 / 2, 1, 0)
  x <- mu + sigma * 2 * W  * (1 / 2 -  V)
  return(x)
}
rLaplace <- function(n, mu = 0, sigma = 1) {
  W <- rexp(n, rate = 1)
  V <- rbinom(n, size = 1, prob = 1/2)
  x <- mu + sigma * 2 * W  * (1 / 2 -  V)
  return(x)
}
ML.Laplace <- function(x, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown")  
    mu.hat <- median(x) else 
    mu.hat <- mu
  if (sigma == "unknown") 
    sigma.hat <- mean(abs(x - mu.hat)) else
    sigma.hat <- sigma  
  return(c(mu.hat, sigma.hat))
}
MM.Laplace <- function(x, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown")  
    mu.hat <- mean(x) else 
      mu.hat <- mu
    if (sigma == "unknown"){
      sigma.hat <- sqrt(mean((x - mu.hat) ^ 2) / 2)} else
        sigma.hat <- sigma  
      return(c(mu.hat, sigma.hat))
}
Laplace.test.ML <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.mu.sigma <- ML.Laplace(x, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- pLaplace(x, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dLaplace(x, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(1), h1(1), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(1, 0, 0, 1), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
Laplace.test.MM <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  MM.mu.sigma <- MM.Laplace(x, mu, sigma)
  mu.hat <- MM.mu.sigma[1]
  sigma.hat <- MM.mu.sigma[2]
  Fi <- pLaplace(x, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dLaplace(x, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(1), h1(1), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * diag(c(1 / 2, 4 / 5)) 
  J <- 1 / sigma.hat * matrix(c(0, h5(1) / 2, 2 * h4(1) / 5, 0), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]; J <- as.matrix(J[,U])
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.mu.sigma = MM.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- normal distribution

dnormal <- function(x, mu = 0, sigma = 1, log = FALSE) dnorm(x, mean = mu, sd = sigma, log = log)
pnormal <- function(x, mu = 0, sigma = 1) pnorm(x, mean = mu, sd = sigma)
rnormal <- function(n, mu = 0, sigma = 1) rnorm(n, mean = mu, sd = sigma)
qnormal <- function(p, mu = 0, sigma = 1) qnorm(p, mean = mu, sd = sigma)
ML.normal <- function(x, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown")  
    mu.hat <- mean(x) else 
    mu.hat <- mu
  if (sigma == "unknown") 
    sigma.hat <- sqrt(mean((x - mu.hat) ^ 2)) else
    sigma.hat <- sigma  
  return(c(mu.hat, sigma.hat))
}
normal.test.ML <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.mu.sigma <- ML.normal(x, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- pnormal(x, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dnormal(x, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(2) * sqrt(2 / pi), h1(2), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(1, 0, 0, 2), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- exp-gamma distribution

dexp.gamma <- function(x, lambda = 1, mu = 0, sigma = 1, log = FALSE){
  y <- (x - mu) / sigma
  log.f <- lambda * y - exp(y) - log(sigma) - lgamma(lambda) 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pexp.gamma <- function(x, lambda = 1, mu = 0, sigma = 1){
  pgamma(exp((x - mu) / sigma), shape = lambda, scale = 1)
}
qexp.gamma <- function(p, lambda = 1, mu = 0, sigma = 1){
  y <- qgamma(p, shape = lambda, scale = 1)
  x <- mu + sigma * log(y)
  return(x)
}
rexp.gamma <- function(n, lambda = 1, mu = 0, sigma = 1){
  y <- rgamma(n, shape = lambda, scale = 1)
  x <- mu + sigma * log(y)
  return(x)
}
ML.exp.gamma <- function(x, lambda = "unknown", mu = "unknown", sigma = "unknown"){
  # initial values  
  exp.x <- exp(x)
  fit <- flexsurvreg(Surv(exp.x) ~ 1, dist = "gengamma", data = data.frame(exp.x))
  lambda0 <- 1 / fit$res[3,1] ^ 2
  sigma0 <- max(1, fit$res[2,1] / fit$res[3,1])
  if (lambda == "unknown" & mu == "unknown" & sigma == "unknown") {
    mx <- mean(x)
    g <- function(vars) {
      lambda <- vars[1]
      sigma <- vars[2]
      g1 <- digamma(lambda) + log(mean(exp(x / sigma)) / lambda) - mean(x / sigma) 
      g2 <- mx + sigma / lambda - sum(x * exp(x / sigma)) / sum(exp(x / sigma)) 
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, sigma0), g)
    lambda.hat <- res$x[1]
    sigma.hat <- res$x[2]
    mu.hat <- sigma.hat * log(mean(exp(x / sigma.hat)) / lambda.hat)
  }
  if (lambda == "unknown" & mu == "unknown" & sigma != "unknown") {
    g <- function(lambda) digamma(lambda) + log(mean(exp(x / sigma)) / lambda) - mean(x / sigma) 
    lambda.hat <- uniroot(g, interval = c(lambda0 / 10, lambda0 * 10), tol = 1e-09)$root
    mu.hat <- sigma * log(mean(exp(x / sigma)) / lambda.hat)
    sigma.hat <- sigma
  }  
  if (lambda != "unknown" & mu == "unknown" & sigma == "unknown") {
    mx <- mean(x)
    sigma0 <- sqrt(6) / pi * sd(x)
    g <- function(sigma) mx + sigma / lambda - sum(x * exp(x / sigma)) / sum(exp(x / sigma)) 
    sigma.hat <- uniroot(g, interval = c(sigma0 / 10, sigma0 * 10), tol = 1e-09)$root
    mu.hat <- sigma.hat * log(mean(exp(x / sigma.hat))/ lambda)
    lambda.hat <- lambda
  }  
  if (lambda != "unknown" & mu == "unknown" & sigma != "unknown") {
    mu.hat <- sigma * log(mean(exp(x / sigma))/ lambda)
    lambda.hat <- lambda
    sigma.hat <- sigma 
  }  
  if (lambda == "unknown" & mu != "unknown" & sigma == "unknown") {
    mx <- mean(x)
    g <- function(vars) {
      lambda <- vars[1]
      sigma <- vars[2]
      g1 <- digamma(lambda) - mean((x - mu) / sigma)
      g2 <- mean((x - mu) * exp((x - mu) / sigma)) - lambda * (mx - mu) - sigma 
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, sigma0), g)
    lambda.hat <- res$x[1]
    sigma.hat <- res$x[2]
    mu.hat <- mu
  }    
  if (lambda != "unknown" & mu != "unknown" & sigma == "unknown") {
    mx <- mean(x)
    g <- function(sigma) mean((x - mu) * exp((x - mu) / sigma)) - lambda * (mx - mu) - sigma   
    sigma.hat <- uniroot(g, interval = c(sigma0 / 10, sigma0 * 10), tol = 1e-09)$root
    lambda.hat <- lambda
    mu.hat <- mu 
  }  
  if (lambda == "unknown" & mu != "unknown" & sigma != "unknown") {
    g <- function(lambda) digamma(lambda) - mean((x - mu) / sigma)
    lambda.hat <- uniroot(g, interval = c(lambda0 / 10, lambda0 * 10), tol = 1e-09)$root
    mu.hat <- mu
    sigma.hat <- sigma 
  } 
  if (lambda != "unknown" & mu != "unknown" & sigma != "unknown") {
    lambda.hat <- lambda
    mu.hat <- mu
    sigma.hat <- sigma 
  } 
  return(c(lambda.hat, mu.hat, sigma.hat))
}
ML.exp.gamma.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown" & sigma == "unknown") {
    mx <- mean(x)
    sigma0 <- sqrt(6) / pi * sd(x)
    g <- function(sigma) mx + sigma / lambda - sum(x * exp(x / sigma)) / sum(exp(x / sigma)) 
    sigma.hat <- uniroot(g, interval = c(sigma0 / 5, sigma0 * 5), tol = 1e-09)$root
    mu.hat <- sigma.hat * log(mean(exp(x / sigma.hat))/ lambda)
  }  
  if (mu == "unknown" & sigma != "unknown") {
    mu.hat <- sigma * log(mean(exp(x / sigma))/ lambda)
    sigma.hat <- sigma 
  }  
  if (mu != "unknown" & sigma == "unknown") {
    mx <- mean(x)
    sigma0 <- sqrt(6) / pi * sd(x)
    g <- function(sigma) mean((x - mu) * exp((x - mu) / sigma)) - lambda * (mx - mu) - sigma   
    sigma.hat <- uniroot(g, interval = c(sigma0 / 5, sigma0 * 5), tol = 1e-09)$root
    mu.hat <- mu 
  }  
  if (mu != "unknown" & sigma != "unknown") {
    mu.hat <- mu 
    sigma.hat <- sigma 
  } 
  return(c(mu.hat, sigma.hat))
}
exp.gamma.test.ML <- function(x, lambda = "unknown", mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.lambda.mu.sigma <-  ML.exp.gamma(x, lambda, mu, sigma)
  lambda.hat <- ML.lambda.mu.sigma[1]
  mu.hat <- ML.lambda.mu.sigma[2]
  sigma.hat <- ML.lambda.mu.sigma[3]
  digamma.lambda <- digamma(lambda.hat) 
  Fi <- pgamma(exp((x - mu.hat) / sigma.hat), shape = lambda.hat, scale = 1)
  neg2L <- - 2 * sum(log(dexp.gamma(x, lambda = lambda.hat, mu = mu.hat, sigma = sigma.hat)))
  G <-  matrix(c(h10(lambda.hat), lambda.hat / sigma.hat * h6(lambda.hat, lambda.hat + 1, 1), h8(lambda.hat) / sigma.hat, 
                 h11(lambda.hat), lambda.hat / sigma.hat * h7(lambda.hat, lambda.hat + 1, 1), h9(lambda.hat) / sigma.hat), 
                 byrow = T, nrow = 2) 
  R <- matrix(c(trigamma(lambda.hat), 1 / sigma.hat, digamma.lambda / sigma.hat,
                1 / sigma.hat, lambda.hat / sigma.hat ^ 2, (lambda.hat * digamma.lambda + 1) / sigma.hat ^ 2,
                digamma.lambda / sigma.hat, (lambda.hat * digamma.lambda + 1) / sigma.hat ^ 2, (lambda.hat * digamma.lambda ^ 2 + 2 * digamma.lambda + lambda.hat * trigamma(lambda.hat) + 1) / sigma.hat ^ 2), nrow = 3)
  U <- (1:3)[c(lambda == "unknown", mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.lambda.mu.sigma = ML.lambda.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
exp.gamma.test.ML.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  digamma.lambda <- digamma(lambda) 
  ML.mu.sigma <- ML.exp.gamma.2p(x, lambda, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- pgamma(exp((x - mu.hat) / sigma.hat), shape = lambda, scale = 1)
  neg2L <- - 2 * sum(log(dexp.gamma(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(lambda * h6(lambda, lambda + 1, 1), lambda * h7(lambda, lambda + 1, 1), h8(lambda), h9(lambda)), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(lambda, lambda * digamma.lambda + 1, lambda * digamma.lambda + 1, lambda * digamma.lambda ^ 2 + 2 * digamma.lambda + lambda * trigamma(lambda) + 1), nrow = 2)    
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- exp-Weibull distribution

dexp.Weibull <- function(x, mu = 0, sigma = 1, log = FALSE) {
  y <- (x - mu) / sigma 
  log.f <- y - exp(y) - log(sigma)
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pexp.Weibull <- function(x, mu = 0, sigma = 1) {
  y <- (x - mu) / sigma 
  F <- 1 - exp(-exp(y)) 
  return(F)
}
qexp.Weibull <- function(p, mu = 0, sigma = 1) {
  y <- log(-log(1 - p)) 
  x <- mu + sigma * y 
  return(x)
}
rexp.Weibull <- function(n, mu = 0, sigma = 1) {
  y <- log(-log(1 - runif(n))) 
  x <- mu + sigma * y 
  return(x)
}
ML.exp.Weibull <- function(x, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown" & sigma == "unknown") {
    mx <- mean(x)
    sigma0 <- sqrt(6) / pi * sd(x)
    g <- function(sigma) mx + sigma - sum(x * exp(x / sigma)) / sum(exp(x / sigma)) 
    sigma.hat <- uniroot(g, interval = c(sigma0 / 5, sigma0 * 5), tol = 1e-09)$root
    mu.hat <- sigma.hat * log(mean(exp(x / sigma.hat)))
  }  
  if (mu == "unknown" & sigma != "unknown") {
    mu.hat <- sigma * log(mean(exp(x / sigma)))
    sigma.hat <- sigma 
  }  
  if (mu != "unknown" & sigma == "unknown") {
    mx <- mean(x)
    sigma0 <- sqrt(6) / pi * sd(x)
    g <- function(sigma) mean((x - mu) * exp((x - mu) / sigma)) - mx + mu - sigma   
    sigma.hat <- uniroot(g, interval = c(sigma0 / 5, sigma0 * 5), tol = 1e-09)$root
    mu.hat <- mu 
  }  
  if (mu != "unknown" & sigma != "unknown") {
    mu.hat <- mu 
    sigma.hat <- sigma 
  } 
  return(c(mu.hat, sigma.hat))
}
exp.Weibull.test.ML <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  euler <- - digamma(1) 
  ML.mu.sigma <- ML.exp.Weibull(x, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- 1 - exp(- exp((x - mu.hat) / sigma.hat))
  neg2L <- - 2 * sum(log(dexp.Weibull(x, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(h6(1,2,1), h7(1,2,1), h8(1), h9(1)), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(1, 1 - euler, 1 - euler, (euler - 1) ^ 2 + pi ^ 2 / 6), nrow = 2)    
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Gumbel distribution

dGumbel <- function(x, mu = 0, sigma = 1, log = FALSE) {
  y <- (x - mu) / sigma 
  log.f <- - y - exp(- y) - log(sigma)
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pGumbel <- function(x, mu = 0, sigma = 1) {
  y <- (x - mu) / sigma 
  F <- exp(-exp(-y)) 
  return(F)
}
qGumbel <- function(p, mu = 0, sigma = 1) {
  y <- -log(-log(p)) 
  x <- mu + sigma * y 
  return(x)
}
rGumbel <- function(n, mu = 0, sigma = 1) {
  y <- -log(-log(runif(n))) 
  x <- mu + sigma * y 
  return(x)
}
ML.Gumbel <- function(x, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown" & sigma == "unknown") {
    mx <- mean(x)
    sigma0 <- sqrt(6) / pi * sd(x)
    g <- function(sigma) -mx + sigma + sum(x * exp(-x / sigma)) / sum(exp(-x / sigma)) 
    sigma.hat <- uniroot(g, interval = c(sigma0 / 5, sigma0 * 5), tol = 1e-09)$root
    mu.hat <- - sigma.hat * log(mean(exp(-x / sigma.hat)))
  }  
  if (mu == "unknown" & sigma != "unknown") {
    mu.hat <- - sigma * log(mean(exp(-x / sigma)))
    sigma.hat <- sigma 
  }  
  if (mu != "unknown" & sigma == "unknown") {
    mx <- mean(x)
    sigma0 <- sqrt(6) / pi * sd(x)
    g <- function(sigma) - mean((x - mu) * exp(-(x - mu) / sigma)) + (mx - mu) - sigma   
    sigma.hat <- uniroot(g, interval = c(sigma0 / 5, sigma0 * 5), tol = 1e-09)$root
    mu.hat <- mu 
  }  
  if (mu != "unknown" & sigma != "unknown") {
    mu.hat <- mu 
    sigma.hat <- sigma 
  } 
  return(c(mu.hat, sigma.hat))
}
Gumbel.test.ML <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  euler <- - digamma(1) 
  ML.mu.sigma <- ML.Gumbel(x, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- exp(- exp(- (x - mu.hat) / sigma.hat))
  neg2L <- - 2 * sum(log(dGumbel(x, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(-h6(1,2,1), h7(1,2,1), h8(1), -h9(1)), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(1, euler - 1, euler - 1, (euler - 1) ^ 2 + pi ^ 2 / 6), nrow = 2)    
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- logistic distribution

dlogistic <- function(x, mu = 0, sigma = 1) dlogis(x, location = mu, scale = sigma)
plogistic <- function(x, mu = 0, sigma = 1) plogis(x, location = mu, scale = sigma)
qlogistic <- function(p, mu = 0, sigma = 1) qlogis(p, location = mu, scale = sigma)
rlogistic <- function(n, mu = 0, sigma = 1) rlogis(n, location = mu, scale = sigma)
ML.logistic <- function(x, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown" & sigma == "unknown") {
    mu.hat <- mean(x)
    sigma.hat <- sqrt(3) * sd(x) / pi
    g1 <- function(mu) {y <- (x - mu) / sigma.hat; mean(2 / (1 + exp(y))) - 1}
    g2 <- function(sigma) {y <- (x - mu.hat) / sigma; mean(y) - mean(2 * y / (1 + exp(y))) - 1}
    err.mu <- err.sigma <- 1000000
    while (err.mu > 1e-10 | err.sigma > 1e-10){
      mu.hat.old <- mu.hat
      sigma.hat.old <- sigma.hat
      mu.hat <- uniroot(g1, interval = c(mu.hat - 10, mu.hat + 10), tol = 1e-09)$root
      sigma.hat <- uniroot(g2, interval = c(sigma.hat / 5, sigma.hat * 5), tol = 1e-09)$root
      err.mu <- abs(mu.hat - mu.hat.old)
      err.sigma <- abs(sigma.hat - sigma.hat.old)
    }
  }  
  if (mu == "unknown" & sigma != "unknown") {
    mu0 <- mean(x)
    g3 <- function(mu) {y <- (x - mu) / sigma; mean(2 / (1 + exp(y))) - 1}
    mu.hat <- uniroot(g3, interval = c(mu0 - 10, mu0 + 10), tol = 1e-09)$root
    sigma.hat <- sigma 
  }  
  if (mu != "unknown" & sigma == "unknown") {
    sigma0 <- sqrt(3) * sd(x) / pi
    g4 <- function(sigma) {y <- (x - mu) / sigma; mean(y) - mean(2 * y / (1 + exp(y))) - 1}
    sigma.hat <- uniroot(g4, interval = c(sigma0 / 5, sigma0 * 5), tol = 1e-09)$root
    mu.hat <- mu 
  }  
  if (mu != "unknown" & sigma != "unknown") {
    mu.hat <- mu 
    sigma.hat <- sigma 
  } 
  return(c(mu.hat, sigma.hat))
}
MM.logistic <- function(x, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown")  
    mu.hat <- mean(x) else 
      mu.hat <- mu
    if (sigma == "unknown") 
      sigma.hat <- sqrt(3) / pi * sqrt(mean((x - mu.hat) ^ 2))  else
        sigma.hat <- sigma  
      return(c(mu.hat, sigma.hat))
}
logistic.test.ML <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.mu.sigma <- ML.logistic(x, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- plogis(x, location = mu.hat, scale = sigma.hat)
  neg2L <- - 2 * sum(log(dlogis(x, location = mu.hat, scale = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, -1 / pi, 0.698397593884459, 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(1 / 3, 0, 0, (3 + pi ^ 2) / 9), nrow = 2)    
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
logistic.test.MM <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  MM.mu.sigma <- MM.logistic(x, mu, sigma)
  mu.hat <- MM.mu.sigma[1]
  sigma.hat <- MM.mu.sigma[2]
  Fi <- plogis(x, location = mu.hat, scale = sigma.hat)
  neg2L <- - 2 * sum(log(dlogis(x, location = mu.hat, scale = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, -1 / pi, 0.698397593884459, 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * diag(c(3 / pi ^ 2, 5 / 4))
  J <- 1 / sigma.hat * matrix(c(0, -0.235854187, 0.4909114316,0), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]; J <- as.matrix(J[,U])
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.mu.sigma = MM.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Student's t-distribution

dStudent <- function(x, lambda, mu = 0, sigma = 1) dt((x - mu) / sigma, df = lambda) / sigma
pStudent <- function(x, lambda, mu = 0, sigma = 1) pt((x - mu) / sigma, df = lambda)
pStudent2 <- function(x, lambda, mu = 0, sigma = 1) {y <- (x - mu) / sigma; 1 / 2 * (1 + sign(y) * (1 - pbeta(1 / (1 + y ^ 2 / lambda), lambda / 2, 1 / 2, lower.tail = TRUE)))}    
qStudent <- function(p, lambda, mu = 0, sigma = 1) mu + qt(p, df = lambda) * sigma 
rStudent <- function(n, lambda, mu = 0, sigma = 1) mu + rt(n, df = lambda) * sigma 
ML.Student <- function(x, lambda = "unknown", mu = "unknown", sigma = "unknown"){
  if (lambda == "unknown" & mu == "unknown" & sigma == "unknown")
  {
    log_vraisemblance <- function(theta, y) {
      lambda <- exp(theta[1])
      mu <- theta[2]
      sigma <- exp(theta[3])
      somme <- sum(log(dStudent(y, lambda, mu, sigma)))
      return(-somme) 
    }
    lambda0 <- c(2, 4, 6, 8)
    mu0 <- rep(median(x), 4)
    sigma0 <- rep(mad(x), 4)
    est <- matrix(0, nrow = 8, ncol = 4)
    for (i in 1:4){
      res <- optim(par = c(log(lambda0[i]), mu0[i], log(sigma0[i])), fn = log_vraisemblance, y = x, method = "BFGS")
      par <- c(exp(res$par[1]), res$par[2], exp(res$par[3]))
      est[2 * i - 1,] <- c(par, log_vraisemblance(res$par, x))
      res <- optim(par = c(log(lambda0[i]), mu0[i], log(sigma0[i])), fn = log_vraisemblance, y = x, method = "Nelder-Mead")
      par <- c(exp(res$par[1]), res$par[2], exp(res$par[3]))
      est[2 * i,] <- c(par, log_vraisemblance(res$par, x))
    }
    best.est = est[which.min(est[,4]),]
    lambda.hat <- best.est[1]
    mu.hat <- best.est[2]
    sigma.hat <- best.est[3]
  }
  if (lambda == "unknown" & mu == "unknown" & sigma != "unknown")
  {
    mu0 <- median(x)    
    g0 <- function(lambda) {y <- (x - mu0) / sigma 
       0.5 * (digamma((lambda + 1) / 2) - digamma(lambda / 2) - 1 / lambda -  mean(log(1 + y^2 / lambda)) + (lambda + 1) / lambda * mean(y^2 / lambda / (1 + y^2 / lambda)))
    }
    lambda0 <- uniroot(g0, c(.1, 1000), tol = 1e-09)$root
    g <- function(vars) {
      lambda <- vars[1]
      mu <- vars[2]
      y <- (x - mu) / sigma
      g1 <- 0.5 * (digamma((lambda + 1) / 2) - digamma(lambda / 2) - 1 / lambda -  mean(log(1 + y^2 / lambda)) + (lambda + 1) / lambda * mean(y^2 / lambda / (1 + y^2 / lambda)))
      g2 <- sum((x - mu) / (sigma ^ 2  + (x - mu) ^ 2 / lambda))
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, mu0), g)
    lambda.hat <- res$x[1]
    mu.hat <- res$x[2]
    sigma.hat <- sigma
  }
  if (lambda == "unknown" & mu != "unknown" & sigma == "unknown")
  {
    sigma0 <- mad(x)
    g0 <- function(lambda) {y <- (x - mu) / sigma0 
    0.5 * (digamma((lambda + 1) / 2) - digamma(lambda / 2) - 1 / lambda -  mean(log(1 + y^2 / lambda)) + (lambda + 1) / lambda * mean(y^2 / lambda / (1 + y^2 / lambda)))
    }
    lambda0 <- uniroot(g0, c(.1, 1000), tol = 1e-09)$root
    g <- function(vars) {
      lambda <- vars[1]
      sigma <- vars[2]
      y <- (x - mu) / sigma
      g1 <- 0.5 * (digamma((lambda + 1) / 2) - digamma(lambda / 2)  -  mean(log(1 + y^2 / lambda)))
      #g1 <- 0.5 * (digamma((lambda + 1) / 2) - digamma(lambda / 2) - 1 / lambda -  mean(log(1 + y^2 / lambda)) + (lambda + 1) / lambda * mean(y^2 / lambda / (1 + y^2 / lambda)))
      g2 <-  - 1  + (lambda + 1) / lambda * mean((x - mu) ^ 2 / (sigma ^ 2 + (x - mu) ^ 2 / lambda))
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, sigma0), g)
    lambda.hat <- res$x[1]
    mu.hat <- mu
    sigma.hat <- res$x[2]
  }
  if (lambda == "unknown" & mu != "unknown" & sigma != "unknown")
  {
    g <- function(lambda) {y <- (x - mu) / sigma
     0.5 * (digamma((lambda + 1) / 2) - digamma(lambda / 2) - 1 / lambda -  mean(log(1 + y^2 / lambda)) + (lambda + 1) / lambda * mean(y^2 / lambda / (1 + y^2 / lambda)))
    }
    lambda.hat <- uniroot(g, c(.1, 1000), tol = 1e-09)$root
    mu.hat <- mu
    sigma.hat <- sigma
  }
  if (lambda != "unknown" & mu == "unknown" & sigma == "unknown") {
    mu0 <- median(x)
    sigma0 <- mad(x)
    g <- function(vars) {
      mu <- vars[1]
      sigma <- vars[2]
      g1 <- sum((x - mu) / (sigma ^ 2  + (x - mu) ^ 2 / lambda))
      g2 <-  - 1  + (lambda + 1) / lambda * mean((x - mu) ^ 2 / (sigma ^ 2 + (x - mu) ^ 2 / lambda))
      c(g1, g2)
    }
    res <- nleqslv(c(mu0, sigma0), g)
    lambda.hat <- lambda
    mu.hat <- res$x[1]
    sigma.hat <- res$x[2]
    }
  if (lambda != "unknown" & mu == "unknown" & sigma != "unknown") {
    g1 <- function(mu)  sum((x - mu) / (sigma ^ 2  + (x - mu) ^ 2 / lambda))
    mu.hat <- uniroot(g1, interval = c(min(x), max(x)), tol = 1e-09)$root
    lambda.hat <- lambda 
    sigma.hat <- sigma 
  }  
  if (lambda != "unknown" & mu != "unknown" & sigma == "unknown") {
    sigma0 <- as.numeric(quantile(x, .65)  - quantile(x, .35)) 
    g2 <- function(sigma) - 1  + (lambda + 1) / lambda * mean((x - mu) ^ 2 / (sigma ^ 2 + (x - mu) ^ 2 / lambda))
    sigma.hat <- uniroot(g2, interval = c(.0001, sigma0 * 100), tol = 1e-09)$root
    lambda.hat <- lambda 
    mu.hat <- mu 
  }  
  if (lambda != "unknown" & mu != "unknown" & sigma != "unknown") {
    lambda.hat <- lambda 
    mu.hat <- mu 
    sigma.hat <- sigma 
  } 
  return(c(lambda.hat, mu.hat, sigma.hat))
}
ML.Student.equations <- function(x){
  mu0 <- median(x)
  sigma0 <- mad(x)
  g0 <- function(lambda) {y <- (x - mu0) / sigma0 
     0.5 * (digamma((lambda + 1) / 2) - digamma(lambda / 2) - 1 / lambda -  mean(log(1 + y^2 / lambda)) + (lambda + 1) / lambda * mean(y^2 / lambda / (1 + y^2 / lambda)))
  }
  lambda0 <- uniroot(g0, c(.1, 1000), tol = 1e-09)$root
  g <- function(vars) {
    lambda <- vars[1]
    mu <- vars[2]
    sigma <- vars[3]
    y <- (x - mu) / sigma
    g1 <-  digamma((lambda + 1) / 2) - digamma(lambda / 2)  -  mean(log(1 + y^2 / lambda))
    #g1 <- digamma((lambda + 1) / 2) - digamma(lambda / 2) -  mean(log(1 + y ^ 2 / lambda)) + 1 / lambda * ((lambda + 1) / lambda * mean(y^2 / (1 + y ^ 2 / lambda)) - 1)
    g2 <- sum((x - mu) / (sigma ^ 2  + (x - mu) ^ 2 / lambda))
    g3 <-  - 1  + (lambda + 1) / lambda * mean((x - mu) ^ 2 / (sigma ^ 2 + (x - mu) ^ 2 / lambda))
    c(g1, g2, g3)
  }
  res <- nleqslv(c(lambda0, mu0, sigma0), g)
  lambda.hat <- res$x[1]
  mu.hat <- res$x[2]
  sigma.hat <- res$x[3]
  return(c(lambda.hat, mu.hat, sigma.hat))
}
ML.Student.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  ## Newton-Raphson version
  if (mu == "unknown" & sigma == "unknown") {
    n <- length(x)
    mu.hat <- median(x)
    sigma.hat <- as.numeric(quantile(x, .65)  - quantile(x, .35)) 
    err.mu <- err.sigma <- 10000
    while (err.mu > 1e-10 | err.sigma > 1e-10){
      mu.hat.old <- mu.hat
      sigma.hat.old <- sigma.hat
      y <- (x - mu.hat) / sigma.hat
      y2l <- y^2/lambda
      mu.hat <- mu.hat - sigma.hat * sum(y / (1 + y2l)) / sum((-1 + y2l) / (1 + y2l) ^ 2)
      y2 <- (x - mu.hat) ^ 2 / sigma.hat ^ 2
      sigma.hat <- sigma.hat +  sigma.hat * 
        (- n * lambda / (lambda + 1)  + sum(y2 / (1 + y2 / lambda))) /
        (2  * sum(y2 / (1+y2/lambda) ^2))
      err.mu <- abs(mu.hat - mu.hat.old)
      err.sigma <- abs(sigma.hat - sigma.hat.old)
    }
  }  
  if (mu == "unknown" & sigma != "unknown") {
    g1 <- function(mu)  sum((x - mu) / (sigma ^ 2  + (x - mu) ^ 2 / lambda))
    mu.hat <- uniroot(g1, interval = c(min(x), max(x)), tol = 1e-09)$root
    sigma.hat <- sigma 
  }  
  if (mu != "unknown" & sigma == "unknown") {
    sigma0 <- as.numeric(quantile(x, .65)  - quantile(x, .35)) 
    g2 <- function(sigma) - 1  + (lambda + 1) / lambda * mean((x - mu) ^ 2 / (sigma ^ 2 + (x - mu) ^ 2 / lambda))
    sigma.hat <- uniroot(g2, interval = c(.0001, sigma0 * 100), tol = 1e-09)$root
    mu.hat <- mu 
  }  
  if (mu != "unknown" & sigma != "unknown") {
    mu.hat <- mu 
    sigma.hat <- sigma 
  } 
  return(c(mu.hat, sigma.hat))
}
ML.Student2.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown" & sigma == "unknown") {
    mu.hat <- median(x)
    sigma.hat <- quantile(x, .65) - quantile(x, .35)
    g1 <- function(mu)  {y <- (x - mu) / sigma.hat; sum(y / (1  + y ^ 2 / lambda))}
    g2 <- function(sigma) {y <- (x - mu.hat) / sigma; - 1  + (lambda + 1) / lambda * mean(y ^ 2 / (1  + y ^ 2 / lambda))}
    err.mu <- err.sigma <- 10000
    interv <- c(min(x), max(x))
    while (err.mu > 1e-10 | err.sigma > 1e-10){
      mu.hat.old <- mu.hat
      sigma.hat.old <- sigma.hat
      mu.hat <- uniroot(g1, interval = interv, tol = 1e-09)$root
      sigma.hat <- uniroot(g2, interval = c(.0001, sigma.hat * 100), tol = 1e-09)$root
      err.mu <- abs(mu.hat - mu.hat.old)
      err.sigma <- abs(sigma.hat - sigma.hat.old)
    }
  }  
  if (mu == "unknown" & sigma != "unknown") {
    mu0 <- median(x)
    g1 <- function(mu)  {y <- (x - mu) / sigma; sum(y / (1  + y ^ 2 / lambda))}
    mu.hat <- uniroot(g1, interval = c(min(x), max(x)), tol = 1e-09)$root
    sigma.hat <- sigma 
  }  
  if (mu != "unknown" & sigma == "unknown") {
    sigma0 <- as.numeric(quantile(x, .65)  - quantile(x, .35)) 
    g2 <- function(sigma) {y <- (x - mu) / sigma; - 1  + (lambda + 1) / lambda * mean(y ^ 2 / (1  + y ^ 2 / lambda))}
    sigma.hat <- uniroot(g2, interval = c(.0001, sigma0 * 100), tol = 1e-09)$root
    mu.hat <- mu 
  }  
  if (mu != "unknown" & sigma != "unknown") {
    mu.hat <- mu 
    sigma.hat <- sigma 
  } 
  return(c(mu.hat, sigma.hat))
}
MM.Student <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  if (mu == "unknown")  
    mu.hat <- mean(x) else 
      mu.hat <- mu
    if (sigma == "unknown") 
    {C2 <-  gamma((lambda - 1) / 2) * sqrt(lambda) / gamma(lambda / 2) / sqrt(pi); sigma.hat <- mean(abs(x - mu.hat)) / C2}  else
      sigma.hat <- sigma  
    return(c(lambda, mu.hat, sigma.hat))
}
Student.test.ML <- function(x, lambda = "unknown", mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.lambda.mu.sigma <- ML.Student(x, lambda, mu, sigma)
  ML.lambda.mu.sigma[1] <- min(ML.lambda.mu.sigma[1], 1000) 
  lambda.hat <- ML.lambda.mu.sigma[1]
  mu.hat <- ML.lambda.mu.sigma[2]
  sigma.hat <- ML.lambda.mu.sigma[3]
  C1 <- 1 / sqrt(lambda.hat) / beta(lambda.hat / 2, 1 / 2)
  Fi <- pt((x - mu.hat) / sigma.hat, lambda.hat)
  neg2L <- - 2 * sum(log(dStudent(x, lambda = lambda.hat, mu = mu.hat, sigma = sigma.hat)))
  G <-  matrix(c(h14(lambda.hat) / 2, 0,  h12(lambda.hat) / sigma.hat, 
                 0, 2 * C1 * h13(lambda.hat) / sigma.hat, 0), byrow = T, nrow = 2) 
  R <-  matrix(c(1 / 4 * (trigamma(lambda.hat / 2) - trigamma((lambda.hat + 1) / 2) -2 * (lambda.hat + 5) / (lambda.hat * (lambda.hat + 1) * (lambda.hat + 3))), 0, -2 / (sigma.hat * (lambda.hat + 1) * (lambda.hat + 3)),
                 0, (lambda.hat + 1) / sigma.hat ^ 2 / (lambda.hat + 3), 0,
                 -2 / (sigma.hat * (lambda.hat + 1) * (lambda.hat + 3)), 0, 2 * lambda.hat / sigma.hat ^ 2 / (lambda.hat + 3)), nrow = 3)   
  U <- (1:3)[c(lambda == "unknown", mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.lambda.mu.sigma = ML.lambda.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))

}
Student.test.ML.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  # lambda > 1
  n <- length(x)
  #C1 <- gamma((lambda + 1) / 2) / sqrt(lambda * pi) / gamma(lambda/2)
  C1 <- 1 / sqrt(lambda) / beta(lambda / 2, 1 / 2)
  ML.mu.sigma <- ML.Student(x, lambda, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- pt((x - mu.hat) / sigma.hat, lambda)
  neg2L <- - 2 * sum(log(dStudent(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)))
  G <-  1 / sigma.hat * matrix(c(0, 2 * C1 * h13(lambda),  h12(lambda), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c((lambda + 1) / (lambda + 3), 0, 0, 2 * lambda / (lambda + 3)), nrow = 2)   
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
Student.test.MM <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  # lambda > 2
  n <- length(x)
  C1 <- 1 / sqrt(lambda) / beta(lambda / 2, 1 / 2)
  #C1 <- gamma((lambda + 1) / 2) / sqrt(lambda * pi) / gamma(lambda/2)
  C2 <- beta((lambda - 1) / 2, 1) / beta(lambda / 2, 1 / 2)  * sqrt(lambda) 
  #C2 <-  gamma((lambda - 1) / 2) * sqrt(lambda) / gamma(lambda / 2) / sqrt(pi)
  C3 <- C2  / (lambda / (lambda - 2) - C2 ^ 2)
  MM.lambda.mu.sigma <- MM.Student(x, lambda, mu, sigma)
  mu.hat <- MM.lambda.mu.sigma[2]
  sigma.hat <- MM.lambda.mu.sigma[3]
  Fi <- pt((x - mu.hat) / sigma.hat, lambda)
  neg2L <- - 2 * sum(log(dStudent(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)))
  G <-  1 / sigma.hat * matrix(c(0, 2 * C1 * h13(lambda),  h12(lambda), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * diag(c((lambda - 2) / lambda, C2 * C3))
  J <- C2 / sigma.hat * matrix(c(0, (lambda - 2) / lambda * h16(lambda), C3 * h15(lambda),0), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]; J <- as.matrix(J[,U])
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.lambda.mu.sigma = MM.lambda.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Distributions on (0,infinity)

## ---- log-EPD distribution

dlogEPD <- function(x, lambda, mu = 0, sigma = 1, log = FALSE) {
  y <- (log(x) - mu) / sigma
  log.f <-  - log(x) - log(2) - log(sigma) - (1 / lambda - 1) * log(lambda)  - lgamma(1 / lambda) - abs(y) ^ lambda / lambda 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
dlogEPD2 <- function(x, lambda, mu = 0, sigma = 1, log = FALSE) {
  y <- (log(x) - mu) / sigma
  log.f <-  - log(x) - log(2) - log(sigma) - log(lambda) / lambda - lgamma(1 + 1 / lambda) - abs(y) ^ lambda / lambda 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
plogEPD <- function(x, lambda, mu = 0, sigma = 1){
  y <- (log(x) - mu) / sigma
  F <- (1/2) * (1 + sign(y) * pgamma(abs(y) ^ lambda / lambda, shape = 1 / lambda, scale = 1))
  return(F)
}
qlogEPD <-function(p, lambda, mu = 0, sigma = 1) {
  W <- qgamma(abs(2 * p - 1), shape = 1 / lambda, scale = 1)
  V <- ifelse(p < 1 / 2, 1, 0)
  x <- mu + sigma * 2 * lambda ^ (1 / lambda) * W ^ (1 / lambda) * (1 / 2 -  V)
  return(exp(x))
}
rlogEPD <- function(n, lambda, mu = 0, sigma = 1) {
  W <- rgamma(n, shape = 1 / lambda, scale = 1)
  V <- rbinom(n, size = 1, prob = 1/2)
  x <- mu + 2 * sigma * lambda ^ (1 / lambda) * W ^ (1 / lambda) * (1 / 2 -  V)
  return(exp(x))
}
ML.logEPD <- function(x, lambda = "unknown", mu = "unknown", sigma = "unknown"){
  logx <- log(x)
  if (lambda == "unknown" & mu == "unknown" & sigma == "unknown")
  {
    log_vraisemblance <- function(theta, y) {
      lambda <- exp(theta[1])
      mu <- theta[2]
      sigma <- exp(theta[3])
      somme <- sum(dEPD(y, lambda, mu, sigma, log = TRUE))
      return(-somme) 
    }
    lambda0 <- c(0.5, 1, 1.5, 2)
    mu0 <- rep(median(logx), 4)
    sigma0 <- mapply(function(mu, lambda) mean(abs(logx - mu) ^ lambda) ^ (1 / lambda), mu0, lambda0)
    est <- matrix(0, nrow = 8, ncol = 4)
    for (i in 1:4){
      res <- optim(par = c(log(lambda0[i]), mu0[i], log(sigma0[i])), fn = log_vraisemblance, y = logx, method = "BFGS")
      par <- c(exp(res$par[1]), res$par[2], exp(res$par[3]))
      est[2 * i - 1,] <- c(par, log_vraisemblance(res$par, logx))
      res <- optim(par = c(log(lambda0[i]), mu0[i], log(sigma0[i])), fn = log_vraisemblance, y = logx, method = "Nelder-Mead")
      par <- c(exp(res$par[1]), res$par[2], exp(res$par[3]))
      est[2 * i,] <- c(par, log_vraisemblance(res$par, logx))
    }
    est[est[,1] <= 0.1 | est[,3] <= 0.01 | est[,4] <= 0.01, 4] <- 999999  ## discard some extreme cases 
    pos <- est[,1] > 20
    if (sum(pos) > 0) {
      est[pos,1] <- 20
      g <- function(mu) sum(abs(logx - mu) ^ (20 - 1) * sign(logx - mu))
      muu <-  uniroot(g, c(min(logx) - 1, max(logx) + 1), tol = 1e-09)$root
      est[pos,2] <- muu
      est[pos,3] <- (mean((abs(logx - muu)) ^ 20)) ^ (1 / 20) 
      est[pos,4] <- log_vraisemblance(c(log(est[pos,1][1]), est[pos,2][1], log(est[pos,3][1])), logx)
    }
    best.est = est[which.min(est[,4]),]
    lambda.hat <- best.est[1]
    mu.hat <- best.est[2]
    sigma.hat <- best.est[3]
  }
  if (lambda == "unknown" & mu == "unknown" & sigma != "unknown")
  {
    mu0 <- median(logx)
    g0 <- function(lambda)  digamma(1 / lambda + 1) + log(lambda) - 1 - mean((abs(logx - mu0) / sigma) ^ lambda * (lambda * log(abs(logx - mu0) / sigma) - 1))
    lambda0 <- uniroot(g0, c(.1, 100), tol = 1e-09)$root
    g <- function(vars) {
      lambda <- vars[1]
      mu <- vars[2]
      g1 <- digamma(1 / lambda + 1) + log(lambda) - 1  - mean((abs(logx - mu) / sigma) ^ lambda * (lambda * log(abs(logx - mu) / sigma) - 1))
      g2 <- mean(abs(logx - mu) ^ (lambda - 1) * sign(logx - mu))
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, mu0), g)
    lambda.hat <- res$x[1]
    mu.hat <- res$x[2]
    sigma.hat <- sigma
  }
  if (lambda == "unknown" & mu != "unknown" & sigma == "unknown")
  {
    sigma0 <- mad(logx)
    g0 <- function(lambda) digamma(1 / lambda + 1) + log(lambda) - 1 - mean((abs(logx - mu) / sigma0) ^ lambda * (lambda * log(abs(logx - mu) / sigma0) - 1))
    lambda0 <- uniroot(g0, c(.1, 100), tol = 1e-09)$root
    g <- function(vars) {
      lambda <- vars[1]
      sigma <- vars[2]
      g1 <- digamma(1 / lambda + 1) + log(lambda) - 1 - mean((abs(logx - mu) / sigma) ^ lambda * (lambda * log(abs(logx - mu) / sigma) - 1))
      g2 <- sigma - (mean((abs(logx - mu)) ^ lambda)) ^ (1 / lambda)
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, sigma0), g)
    lambda.hat <- res$x[1]
    mu.hat <- mu
    sigma.hat <- res$x[2]
  }
  if (lambda == "unknown" & mu != "unknown" & sigma != "unknown")
  {
    g <- function(lambda) digamma(1 / lambda + 1) + log(lambda) - 1 - mean((abs(logx - mu) / sigma) ^ lambda * (lambda * log(abs(logx - mu) / sigma) - 1))
    lambda.hat <- uniroot(g, c(.1, 100), tol = 1e-09)$root
    mu.hat <- mu
    sigma.hat <- sigma
  }
  if (lambda != "unknown")
  {
    if (mu == "unknown") 
    {
      if (lambda == 1) mu.hat <- median(logx) else 
        if (lambda == 2) mu.hat <- mean(logx) else 
        {
          g <- function(mu) sum(abs(logx - mu) ^ (lambda - 1) * sign(logx - mu))
          mu.hat <- uniroot(g, c(min(logx) - 1, max(logx) + 1), tol = 1e-09)$root
        }
    } else
      mu.hat <- mu
    if (sigma == "unknown") 
      sigma.hat <- (mean((abs(logx - mu.hat)) ^ lambda)) ^ (1 / lambda) else
        sigma.hat <- sigma  
      lambda.hat <- lambda 
  } 
  return(c(lambda.hat, mu.hat, sigma.hat))
}
ML.logEPD.equations <- function(x){
  logx <- log(x)
  mu0 <- mean(logx)
  sigma0 <- mad(logx)
  g0 <- function(lambda) digamma(1 / lambda + 1) + log(lambda) -  mean((abs(logx - mu0) / sigma0) ^ lambda * lambda * log(abs(logx - mu0) / sigma0)) +
    mean((abs(logx - mu0) / sigma0) ^ lambda) - 1          
  #g0 <- function(lambda) digamma(1 / lambda + 1) + log(lambda)  - mean((abs(logx - mu0) / sigma0) ^ lambda * lambda * log(abs(logx - mu0) / sigma0))
  lambda0 <- uniroot(g0, c(.1, 100), tol = 1e-09)$root
  g <- function(vars) {
    lambda <- vars[1]
    mu <- vars[2]
    sigma <- vars[3]
    g1 <- digamma(1 / lambda + 1) + log(lambda)  - mean((abs(logx - mu) / sigma) ^ lambda * lambda * log(abs(logx - mu) / sigma))
    #g1 <- digamma(1 / lambda + 1) + log(lambda) - mean((abs(logx - mu) / sigma) ^ lambda * lambda * log(abs(logx - mu) / sigma)) +
    #  mean((abs(logx - mu) / sigma) ^ lambda) - 1          
    g2 <- mean(abs(logx - mu) ^ (lambda - 1) * sign(logx - mu))
    g3 <- sigma - (mean((abs(logx - mu)) ^ lambda)) ^ (1 / lambda)
    c(g1, g2, g3)
  }
  res <- nleqslv(c(lambda0, mu0, sigma0), g)
  lambda.hat <- res$x[1]
  mu.hat <- res$x[2]
  sigma.hat <- res$x[3]
  return(c(lambda.hat, mu.hat, sigma.hat))
}
ML.logEPD.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  logx <- log(x)
  if (mu == "unknown") 
  {
    if (lambda == 1) mu.hat <- median(logx) else 
      if (lambda == 2) mu.hat <- mean(logx) else 
      {
        g <- function(mu) sum(abs(logx - mu) ^ (lambda - 1) * sign(logx - mu))
        mu.hat <- uniroot(g, c(min(logx) - 1, max(logx) + 1), tol = 1e-09)$root
      }
  } else
    mu.hat <- mu
  if (sigma == "unknown") 
    sigma.hat <- (mean((abs(logx - mu.hat)) ^ lambda)) ^ (1 / lambda) else
    sigma.hat <- sigma  
  return(c(mu.hat, sigma.hat))
}
MM.logEPD <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  logx <- log(x)
  if (mu == "unknown")  
    mu.hat <- mean(logx) else 
      mu.hat <- mu
    if (sigma == "unknown"){C1 <- gamma(1 / lambda) / lambda ^ (2 / lambda) / gamma(3 / lambda);
      sigma.hat <- sqrt(C1 * mean((logx - mu.hat) ^ 2))} else
        sigma.hat <- sigma  
      return(c(lambda, mu.hat, sigma.hat))
}
logEPD.test.ML <- function(x, lambda = "unknown", mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.lambda.mu.sigma <- ML.logEPD(x, lambda, mu, sigma)
  lambda.hat <- ML.lambda.mu.sigma[1]
  mu.hat <- ML.lambda.mu.sigma[2]
  sigma.hat <- ML.lambda.mu.sigma[3]
  Fi <- plogEPD(x, lambda = lambda.hat, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dlogEPD(x, lambda = lambda.hat, mu = mu.hat, sigma = sigma.hat)))
  G <-  matrix(c(1 / lambda.hat ^ 2  * (h1(lambda.hat) - h3(lambda.hat)), 0, h1(lambda.hat) / sigma.hat,  
                 0, h2(lambda.hat) / sigma.hat / lambda.hat ^ (1 / lambda.hat - 1) / gamma(1 / lambda.hat), 0), byrow = T, nrow = 2)
  C1 <- digamma(1/lambda.hat + 1) + log(lambda.hat)
  R <-  matrix(c((1 / lambda.hat ^ 3) * ((1 / lambda.hat + 1) * trigamma(1/lambda.hat + 1) +  C1 ^ 2  - 1), 0, - C1 / sigma.hat / lambda.hat,
                 0, lambda.hat ^ (2 - 2 / lambda.hat) * gamma(2 - 1 / lambda.hat) / sigma.hat ^ 2 / gamma(1 / lambda.hat), 0,  
                 - C1 / sigma.hat / lambda.hat, 0, lambda.hat / sigma.hat ^ 2), byrow = T, nrow = 3)      
  U <- (1:3)[c(lambda == "unknown", mu == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.lambda.mu.sigma = ML.lambda.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
logEPD.test.ML.2p <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.mu.sigma <- ML.logEPD.2p(x, lambda, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- plogEPD(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dlogEPD(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(lambda) / lambda ^ (1 / lambda - 1) / gamma(1 / lambda), h1(lambda), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(lambda ^ (2 - 2 / lambda) * gamma(2 - 1 / lambda) / gamma(1 / lambda), 0, 0, lambda), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
logEPD.test.MM <- function(x, lambda, mu = "unknown", sigma = "unknown"){
  # lambda must be specified
  n <- length(x)
  C1 <- gamma(1 / lambda) / lambda ^ (2 / lambda) / gamma(3 / lambda)
  C2 <- gamma(3 / lambda) ^ 2 / (gamma(1 / lambda) * gamma(5 / lambda) - gamma(3 / lambda) ^ 2)
  MM.lambda.mu.sigma <- MM.logEPD(x, lambda, mu, sigma)
  mu.hat <- MM.lambda.mu.sigma[2]
  sigma.hat <- MM.lambda.mu.sigma[3]
  Fi <- plogEPD(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dlogEPD(x, lambda = lambda, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(lambda) / lambda ^ (1 / lambda - 1) / gamma(1 / lambda), h1(lambda), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * diag(c(C1, 4 * C2)) 
  J <- 1 / sigma.hat * matrix(c(0, h5(lambda) / lambda ^ (1 / lambda) / gamma(3 / lambda), 2 * C2 * h4(lambda), 0), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]; J <- as.matrix(J[,U])
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.lambda.mu.sigma = MM.lambda.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Log-Laplace distribution

dlogLaplace <- function(x, mu = 0, sigma = 1, log = FALSE) {
  y <- (log(x) - mu) / sigma
  log.f <-  - log(2) - log(sigma) - log(x) - abs(y) 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
plogLaplace <- function(x, mu = 0, sigma = 1){
  y <- (log(x) - mu) / sigma
  F <- (1/2) * (1 + sign(y) * (1 - exp(-abs(y))))
  return(F)
}
qlogLaplace <- function(p, mu = 0, sigma = 1) {
  W <- qexp(abs(2 * p - 1))
  V <- ifelse(p < 1 / 2, 1, 0)
  x <- exp(mu + sigma * 2 * W  * (1 / 2 -  V))
  return(x)
}
rlogLaplace <- function(n, mu = 0, sigma = 1) {
  W <- rexp(n, rate = 1)
  V <- rbinom(n, size = 1, prob = 1/2)
  x <- exp(mu + 2 * sigma * W  * (1 / 2 -  V))
  return(x)
}
ML.logLaplace <- function(x, mu = "unknown", sigma = "unknown"){
  logx <- log(x)
  if (mu == "unknown")  
    mu.hat <- median(logx) else 
    mu.hat <- mu
  if (sigma == "unknown") 
    sigma.hat <- mean(abs(logx - mu.hat)) else
    sigma.hat <- sigma  
  return(c(mu.hat, sigma.hat))
}
MM.logLaplace <- function(x, mu = "unknown", sigma = "unknown"){
  logx <- log(x)
  if (mu == "unknown")  
    mu.hat <- mean(logx) else 
      mu.hat <- mu
    if (sigma == "unknown"){
      sigma.hat <- sqrt(mean((logx - mu.hat) ^ 2) / 2)} else
        sigma.hat <- sigma  
      return(c(mu.hat, sigma.hat))
}
logLaplace.test.ML <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.mu.sigma <- ML.logLaplace(x, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- plogLaplace(x, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dlogLaplace(x, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(1), h1(1), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(1, 0, 0, 1), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
logLaplace.test.MM <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  MM.mu.sigma <- MM.logLaplace(x, mu, sigma)
  mu.hat <- MM.mu.sigma[1]
  sigma.hat <- MM.mu.sigma[2]
  Fi <- plogLaplace(x, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dlogLaplace(x, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(1), h1(1), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * diag(c(1 / 2, 4 / 5)) 
  J <- 1 / sigma.hat * matrix(c(0, h5(1) / 2, 2 * h4(1) / 5, 0), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]; J <- as.matrix(J[,U])
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.mu.sigma = MM.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- lognormal distribution

dlognormal <- function(x, mu = 0, sigma = 1) dlnorm(x, meanlog = mu, sdlog = sigma)
plognormal <- function(x, mu = 0, sigma = 1) plnorm(x, meanlog = mu, sdlog = sigma)
qlognormal <- function(p, mu = 0, sigma = 1) qlnorm(p, meanlog = mu, sdlog = sigma)
rlognormal <- function(n, mu = 0, sigma = 1) rlnorm(n, meanlog = mu, sdlog = sigma) 
ML.lognormal <- function(x, mu = "unknown", sigma = "unknown"){
  lx <- log(x)
  if (mu == "unknown")  
    mu.hat <- mean(lx) else 
    mu.hat <- mu
  if (sigma == "unknown") 
    sigma.hat <- sqrt(mean((lx - mu.hat) ^ 2)) else
    sigma.hat <- sigma  
  return(c(mu.hat, sigma.hat))
}
lognormal.test.ML <- function(x, mu = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.mu.sigma <- ML.lognormal(x, mu, sigma)
  mu.hat <- ML.mu.sigma[1]
  sigma.hat <- ML.mu.sigma[2]
  Fi <- plognormal(x, mu = mu.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dlognormal(x, mu = mu.hat, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(0, h2(2) * sqrt(2 / pi), h1(2), 0), nrow = 2) 
  R <- 1 / sigma.hat ^ 2 * matrix(c(1, 0, 0, 2), nrow = 2) 
  U <- (1:2)[c(mu == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.sigma = ML.mu.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- half-EPD distribution

dhalfEPD <- function(x, lambda, sigma = 1, log = FALSE) {
  log.f <-  - log(sigma) - (1 / lambda - 1) * log(lambda) - lgamma(1 / lambda) - (x / sigma) ^ lambda / lambda 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
dhalfEPD2 <- function(x, lambda, sigma = 1, log = FALSE) {
  log.f <-  - log(sigma) - log(lambda) / lambda - lgamma(1 + 1 / lambda) - (x / sigma) ^ lambda / lambda 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
phalfEPD <- function(x, lambda, sigma = 1){
  F <-  pgamma((x / sigma) ^ lambda / lambda, shape = 1 / lambda, scale = 1)
  return(F)
}
qhalfEPD <- function(p, lambda, sigma = 1) {
  W <- qgamma(p, shape = 1 / lambda, scale = 1)
  x <- sigma * lambda ^ (1 / lambda) * W ^ (1 / lambda) 
  return(x)
}
rhalfEPD <- function(n, lambda, sigma = 1) {
  W <- rgamma(n, shape = 1 / lambda, scale = 1)
  x <- sigma * lambda ^ (1 / lambda) * W ^ (1 / lambda) 
  return(x)
}
ML.halfEPD <- function(x, lambda = "unknown", sigma = "unknown"){
  if (lambda == "unknown" & sigma == "unknown")
  {
    sigma0 <- mad(x)
    g0 <- function(lambda) digamma(1 / lambda + 1) + log(lambda) - 1 - mean((x / sigma0) ^ lambda * (lambda * log(x / sigma0) - 1))
    lambda0 <- uniroot(g0, c(.1, 100), tol = 1e-09)$root
    g <- function(vars) {
      lambda <- vars[1]
      sigma <- vars[2]
      #g1 <- digamma(1 / lambda + 1) + log(lambda) +  mean((x / sigma) ^ lambda) - 1 - mean((x / sigma) ^ lambda * lambda * log(x / sigma))
      g1 <- digamma(1 / lambda + 1) + log(lambda)  - mean((x / sigma) ^ lambda * lambda * log(x / sigma))
      g2 <- sigma - (mean(x ^ lambda)) ^ (1 / lambda)
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, sigma0), g)
    lambda.hat <- res$x[1]
    sigma.hat <- res$x[2]
  }
  if (lambda == "unknown" & sigma != "unknown")
  {
    g <- function(lambda) digamma(1 / lambda + 1) + log(lambda) - 1 - mean((x / sigma) ^ lambda * (lambda * log(x / sigma) - 1))
    lambda.hat <- uniroot(g, c(.1, 100), tol = 1e-09)$root
    sigma.hat <- sigma
  }
  if (lambda != "unknown" & sigma == "unknown")
  {  
    lambda.hat <- lambda 
    sigma.hat <- (mean(x ^ lambda)) ^ (1 / lambda)
  } 
  if (lambda != "unknown" & sigma != "unknown")
  {
    lambda.hat <- lambda 
    sigma.hat <- sigma 
  }
  return(c(lambda.hat, sigma.hat))
}
ML.halfEPD.1p <- function(x, lambda, sigma = "unknown"){
  if (sigma == "unknown") sigma.hat <- (mean(x ^ lambda)) ^ (1 / lambda) else sigma.hat <- sigma
  return(sigma.hat)
}
MM.halfEPD <- function(x, lambda, sigma = "unknown"){
  if (sigma == "unknown") {C1 <- gamma(1 / lambda) / lambda ^ (1 / lambda) / gamma(2 / lambda); sigma.hat <- C1 * mean(x)} else sigma.hat <- sigma
  return(c(lambda, sigma.hat))
}
halfEPD.test.ML <- function(x, lambda = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.lambda.sigma <- ML.halfEPD(x, lambda, sigma)
  lambda.hat <- ML.lambda.sigma[1]
  sigma.hat <- ML.lambda.sigma[2]
  Fi <- phalfEPD(x, lambda = lambda.hat, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dhalfEPD(x, lambda = lambda.hat, sigma = sigma.hat)))
  G <-  matrix(c((h6(1 / lambda.hat, 1 / lambda.hat + 1, 1) - h17(lambda.hat)) / lambda.hat ^ 2, h6(1 / lambda.hat, 1 / lambda.hat + 1, 1) / sigma.hat,
                 (h7(1 / lambda.hat, 1 / lambda.hat + 1, 1) - h18(lambda.hat)) / lambda.hat ^ 2, h7(1 / lambda.hat, 1 / lambda.hat + 1, 1) / sigma.hat), byrow = T, nrow = 2) 
  C1 <- digamma(1/lambda.hat + 1) + log(lambda.hat)
  R <-  matrix(c((1 / lambda.hat ^ 3) * ((1 / lambda.hat + 1) * trigamma(1/lambda.hat + 1) +  C1 ^ 2  - 1), - C1 / sigma.hat / lambda.hat,
                 - C1 / sigma.hat / lambda.hat, lambda.hat / sigma.hat ^ 2), byrow = T, nrow = 2)   
  U <- (1:2)[c(lambda == "unknown", sigma == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.lambda.sigma = ML.lambda.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
halfEPD.test.ML.1p <- function(x, lambda, sigma = "unknown"){
  n <- length(x)
  sigma.hat <- ML.halfEPD.1p(x, lambda, sigma) 
  Fi <- phalfEPD(x, lambda = lambda, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dhalfEPD(x, lambda = lambda, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(h6(1 / lambda, 1 / lambda + 1, 1), h7(1 / lambda, 1 / lambda + 1, 1)), nrow = 2) 
  R <- lambda / sigma.hat ^ 2
  U <- c(1)[sigma == "unknown"] 
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.sigma = sigma.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
halfEPD.test.MM <- function(x, lambda, sigma = "unknown"){
  # lambda must be specified
  n <- length(x)
  C1 <- gamma(1 / lambda) / lambda ^ (1 / lambda) / gamma(2 / lambda)
  C2 <- gamma(2 / lambda) ^ 2 / (gamma(1 / lambda) * gamma(3 / lambda) - gamma(2 / lambda) ^ 2)
  sigma.hat <- MM.halfEPD(x, lambda, sigma)[2]
  Fi <- phalfEPD(x, lambda = lambda, sigma = sigma.hat)
  neg2L <- - 2 * sum(log(dhalfEPD(x, lambda = lambda, sigma = sigma.hat)))
  G <- 1 / sigma.hat * matrix(c(h6(1 / lambda, 1 / lambda + 1, 1), h7(1 / lambda, 1 / lambda + 1, 1)), nrow = 2) 
  R <- C2 / sigma.hat ^ 2
  J <- C2 / sigma.hat * matrix(c(h6(1 / lambda, 2 / lambda, 1), h7(1 / lambda, 2 / lambda, 1)), nrow = 2) 
  U <- c(1)[sigma == "unknown"] 
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.lambda.sigma = c(lambda, sigma.hat), Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- generalized gamma or GG distribution

dGG <- function(x, lambda = 1, beta = 1, rho = 1, log = FALSE){
  w <- x / beta
  log.f <- log(rho) - log(x) - lgamma(lambda) + lambda * rho * log(w) - w ^ rho
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pGG <- function(x, lambda = 1, beta = 1, rho = 1){
  pgamma((x / beta) ^ rho, shape = lambda, scale = 1)
}
qGG <- function(p, lambda = 1, beta = 1, rho = 1){
  y <- qgamma(p, shape = lambda, scale = 1)
  x <- beta * y ^ (1 / rho)
  return(x)
}
rGG <- function(n, lambda = 1, beta = 1, rho = 1){
  y <- rgamma(n, shape = lambda, scale = 1)
  x <- beta * y ^ (1 / rho)
  return(x)
}
ML.GG <- function(x, lambda = "unknown", beta = "unknown", rho = "unknown"){
  # initial values  
  fit <- flexsurvreg(Surv(x) ~ 1, dist = "gengamma", data = data.frame(x))
  sigma <- fit$res[2,1]
  Q <- fit$res[3,1]
  lambda0 <- 1 / Q ^ 2
  rho0 <- Q / sigma
  if (lambda == "unknown" & beta == "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    g <- function(vars) {
      lambda <- vars[1]
      rho <- vars[2]
      g1 <- digamma(lambda) + log(mean(x ^ rho) / lambda) - rho * mean(log(x)) 
      g2 <- mlx + 1 / lambda / rho - sum(x ^ rho * lx) / sum(x ^ rho) 
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, rho0), g)
    lambda.hat <- res$x[1]
    rho.hat <- res$x[2]
    beta.hat <- (mean(x ^ rho.hat) / lambda.hat) ^ (1 / rho.hat)
  }
  if (lambda == "unknown" & beta == "unknown" & rho != "unknown") {
    g <- function(lambda) digamma(lambda) + log(mean(x ^ rho) / lambda) - mean(log(x ^ rho)) 
    lambda.hat <- uniroot(g, interval = c(lambda0 / 10, lambda0 * 10), tol = 1e-09)$root
    beta.hat <- (mean(x ^ rho) / lambda.hat) ^ (1 / rho)
    rho.hat <- rho
  }  
  if (lambda != "unknown" & beta == "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    g <- function(rho)  mlx + 1 / lambda / rho - sum(x ^ rho * lx) / sum(x ^ rho) 
    rho.hat <- uniroot(g, interval = c(rho0 / 10, rho0 * 10), tol = 1e-09)$root
    beta.hat <- (mean(x ^ rho.hat) / lambda) ^ (1 / rho.hat)
    lambda.hat <- lambda
  }  
  if (lambda != "unknown" & beta == "unknown" & rho != "unknown") {
    beta.hat <- (mean(x ^ rho) / lambda) ^ (1 / rho)
    lambda.hat <- lambda
    rho.hat <- rho 
  }  
  if (lambda == "unknown" & beta != "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    g <- function(vars) {
      lambda <- vars[1]
      rho <- vars[2]
      g1 <- digamma(lambda) - rho * mean(log(x / beta))
      g2 <- mean((x / beta) ^ rho * log(x / beta)) - lambda * mean(log(x / beta)) - 1 / rho
      c(g1, g2)
    }
    res <- nleqslv(c(lambda0, rho0), g)
    lambda.hat <- res$x[1]
    rho.hat <- res$x[2]
    beta.hat <- beta
  }    
  if (lambda != "unknown" & beta != "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    g <- function(rho)  mean((x / beta) ^ rho * log(x / beta)) - lambda * mean(log(x / beta)) - 1 / rho
    rho.hat <- uniroot(g, interval = c(rho0 / 10, rho0 * 10), tol = 1e-09)$root
    lambda.hat <- lambda
    beta.hat <- beta 
  }  
  if (lambda == "unknown" & beta != "unknown" & rho != "unknown") {
    g <- function(lambda) digamma(lambda) - rho * mean(log(x / beta))
    lambda.hat <- uniroot(g, interval = c(lambda0 / 10, lambda0 * 10), tol = 1e-09)$root
    beta.hat <- beta
    rho.hat <- rho 
  } 
  if (lambda != "unknown" & beta != "unknown" & rho != "unknown") {
    lambda.hat <- lambda
    beta.hat <- beta
    rho.hat <- rho 
  } 
  return(c(lambda.hat, beta.hat, rho.hat))
}
ML.GG.likelihood <- function(x){
  log_vraisemblance <- function(theta, y) {
    lambda <- exp(theta[1])
    beta <- exp(theta[2])
    rho <- exp(theta[3])
    somme <- sum(log(dGG(y, lambda, beta, rho)))
    return(-somme) 
  }
  lambda0 <- c(2, 4, 6, 8)
  rho0 <- sqrt(trigamma(lambda0)) / sd(log(x)) 
  beta0 <- exp(mean(log(x)) - digamma(lambda0) / rho0)
  est <- matrix(0, nrow = 8, ncol = 4)
  for (i in 1:4){
    res <- optim(par = c(log(lambda0[i]), log(beta0[i]), log(rho0[i])), fn = log_vraisemblance, y = x, method = "BFGS")
    par <- c(exp(res$par[1]), exp(res$par[2]), exp(res$par[3]))
    est[2 * i - 1,] <- c(par, log_vraisemblance(res$par, x))
    res <- optim(par = c(log(lambda0[i]), log(beta0[i]), log(rho0[i])), fn = log_vraisemblance, y = x, method = "Nelder-Mead")
    par <- c(exp(res$par[1]), exp(res$par[2]), exp(res$par[3]))
    est[2 * i,] <- c(par, log_vraisemblance(res$par, x))
  }
  best.est = est[which.min(est[,4]),]
  lambda.hat <- best.est[1]
  beta.hat <- best.est[2]
  rho.hat <- best.est[3]
  return(c(lambda.hat, beta.hat, rho.hat))
}
ML.GG.2p <- function(x, lambda, beta = "unknown", rho = "unknown"){
  if (beta == "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    rho0 <- pi / sqrt(6) / sd(lx)
    g <- function(rho)  mlx + 1 / lambda / rho - sum(x ^ rho * lx) / sum(x ^ rho) 
    rho.hat <- uniroot(g, interval = c(rho0 / 5, rho0 * 5), tol = 1e-09)$root
    beta.hat <- (mean(x ^ rho.hat) / lambda) ^ (1 / rho.hat)
  }  
  if (beta == "unknown" & rho != "unknown") {
    beta.hat <- (mean(x ^ rho) / lambda) ^ (1 / rho)
    rho.hat <- rho 
  }  
  if (beta != "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    rho0 <- pi / sqrt(6) / sd(lx)
    g <- function(rho)  mean((x / beta) ^ rho * log(x / beta)) - lambda * mean(log(x / beta)) - 1 / rho
    rho.hat <- uniroot(g, interval = c(rho0 / 5, rho0 * 5), tol = 1e-09)$root
    beta.hat <- beta 
  }  
  if (beta != "unknown" & rho != "unknown") {
    beta.hat <- beta
    rho.hat <- rho 
  } 
  return(c(beta.hat, rho.hat))
}
GG.test.ML <- function(x, lambda = "unknown", beta = "unknown", rho = "unknown"){
  n <- length(x)
  ML.lambda.beta.rho <-  ML.GG(x, lambda, beta, rho)
  lambda.hat <- ML.lambda.beta.rho[1]
  beta.hat <- ML.lambda.beta.rho[2]
  rho.hat <- ML.lambda.beta.rho[3]
  digamma.lambda <- digamma(lambda.hat) 
  Fi <-  pgamma((x / beta.hat) ^ rho.hat, shape = lambda.hat, scale = 1)
  neg2L <- - 2 * sum(log(dGG(x, lambda = lambda.hat, beta = beta.hat, rho = rho.hat)))
  G <- matrix(c(h10(lambda.hat), rho.hat / beta.hat * lambda.hat * h6(lambda.hat, lambda.hat + 1, 1), - 1 / rho.hat * h8(lambda.hat), 
                h11(lambda.hat), rho.hat / beta.hat * lambda.hat * h7(lambda.hat, lambda.hat + 1, 1), - 1 / rho.hat * h9(lambda.hat)), 
              nrow = 2, byrow = TRUE) 
  R <- matrix(c(trigamma(lambda.hat), rho.hat / beta.hat, -1 / rho.hat * digamma.lambda, 
                rho.hat / beta.hat, rho.hat ^ 2 / beta.hat ^ 2 * lambda.hat, - 1 / beta.hat * (lambda.hat * digamma.lambda + 1),
                -1 / rho.hat * digamma.lambda, - 1 / beta.hat * (lambda.hat * digamma.lambda + 1), 1 / rho.hat ^ 2 * (lambda.hat * digamma.lambda ^ 2 + 2 * digamma.lambda + lambda.hat * trigamma(lambda.hat) + 1)), nrow = 3)    
  U <- (1:3)[c(lambda == "unknown", beta == "unknown", rho == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.lambda.beta.rho = ML.lambda.beta.rho, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
GG.test.2p.ML <- function(x, lambda, beta = "unknown", rho = "unknown"){
  n <- length(x)
  digamma.lambda <- digamma(lambda) 
  ML.beta.rho <-  ML.GG.2p(x, lambda, beta, rho)
  beta.hat <- ML.beta.rho[1]
  rho.hat <- ML.beta.rho[2]
  Fi <-  pgamma((x / beta.hat) ^ rho.hat, shape = lambda, scale = 1)
  neg2L <- - 2 * sum(log(dGG(x, lambda = lambda, beta = beta.hat, rho = rho.hat)))
  G <- matrix(c(rho.hat / beta.hat * lambda * h6(lambda, lambda + 1, 1), rho.hat / beta.hat * lambda * h7(lambda, lambda + 1, 1), - 1 / rho.hat * h8(lambda), - 1 / rho.hat * h9(lambda)), nrow = 2) 
  R <- matrix(c(rho.hat ^ 2 / beta.hat ^ 2 * lambda, - 1 / beta.hat * (lambda * digamma.lambda + 1), - 1 / beta.hat * (lambda * digamma.lambda + 1), 1 / rho.hat ^ 2 * (lambda * digamma.lambda ^ 2 + 2 * digamma.lambda + lambda * trigamma(lambda) + 1)), nrow = 2)    
  U <- (1:2)[c(beta == "unknown", rho == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.beta.rho = ML.beta.rho, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Weibull distribution 

dWeibull <- function(x, beta = 1, rho = 1, log = FALSE) dweibull(x, shape = rho, scale = beta, log = log)
pWeibull <- function(x, beta = 1, rho = 1) pweibull(x, shape = rho, scale = beta)
qWeibull <- function(p, beta = 1, rho = 1) qweibull(p, shape = rho, scale = beta) 
rWeibull <- function(n, beta = 1, rho = 1) rweibull(n, shape = rho, scale = beta)
ML.Weibull <- function(x, beta = "unknown", rho = "unknown"){
  if (beta == "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    rho0 <- pi / sqrt(6) / sd(lx)
    g = function(rho)  mlx + 1 / rho - sum(x ^ rho * lx) / sum(x ^ rho) 
    rho.hat <- uniroot(g, interval = c(rho0 / 5, rho0 * 5), tol = 1e-09)$root
    beta.hat <- (mean(x ^ rho.hat)) ^ (1 / rho.hat)
  }  
  if (beta == "unknown" & rho != "unknown") {
    beta.hat <- (mean(x ^ rho)) ^ (1 / rho)
    rho.hat <- rho 
  }  
  if (beta != "unknown" & rho == "unknown") {
    lx <- log(x)
    rho0 <- pi / sqrt(6) / sd(lx)
    g = function(rho)  mean((x / beta) ^ rho * log(x / beta)) - mean(log(x / beta)) - 1 / rho
    rho.hat <- uniroot(g, interval = c(rho0 / 5, rho0 * 5), tol = 1e-09)$root
    beta.hat <- beta 
  }  
  if (beta != "unknown" & rho != "unknown") {
    beta.hat <- beta
    rho.hat <- rho 
  } 
  return(c(beta.hat, rho.hat))
}
Weibull.test.ML <- function(x, beta = "unknown", rho = "unknown"){
  n <- length(x)
  euler <- - digamma(1) 
  ML.beta.rho <-  ML.Weibull(x, beta, rho)
  beta.hat <- ML.beta.rho[1]
  rho.hat <- ML.beta.rho[2]
  Fi <-  1 - exp(- (x / beta.hat) ^ rho.hat)  
  neg2L <- - 2 * sum(log(dWeibull(x, beta = beta.hat, rho = rho.hat)))
  G <- matrix(c(rho.hat / beta.hat * h6(1,2,1), rho.hat / beta.hat * h7(1,2,1), - 1 / rho.hat * h8(1), - 1 / rho.hat * h9(1)), nrow = 2) 
  R <- matrix(c(rho.hat ^ 2 / beta.hat ^ 2, 1 / beta.hat * (euler - 1),  1 / beta.hat * (euler - 1), 1 / rho.hat ^ 2 * ((euler - 1) ^ 2 + pi ^ 2 / 6 )), nrow = 2)    
  U <- (1:2)[c(beta == "unknown", rho == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.beta.rho = ML.beta.rho, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Frechet distribution 

dFrechet <- function(x, beta, rho) dfrechet(x, scale = beta, shape = rho)
pFrechet <- function(x, beta, rho) pfrechet(x, scale = beta, shape = rho)
qFrechet <- function(p, beta, rho) qfrechet(p, scale = beta, shape = rho)
rFrechet <- function(n, beta, rho) rfrechet(n, scale = beta, shape = rho)
ML.Frechet <- function(x, beta = "unknown", rho = "unknown"){
  if (beta == "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    rho0 <- pi / sqrt(6) / sd(lx)
    g = function(rho) - mlx + 1 / rho + sum(x ^ (- rho) * lx) / sum(x ^ (- rho)) 
    rho.hat <- uniroot(g, interval = c(rho0 / 5, rho0 * 5), tol = 1e-09)$root
    beta.hat <- (mean(x ^ (- rho.hat))) ^ (- 1 / rho.hat)
  }  
  if (beta == "unknown" & rho != "unknown") {
    beta.hat <- (mean(x ^ (- rho))) ^ (- 1 / rho)
    rho.hat <- rho 
  }  
  if (beta != "unknown" & rho == "unknown") {
    lx <- log(x)
    mlx <- mean(lx)
    rho0 <- pi / sqrt(6) / sd(lx)
    g = function(rho)  mean((x / beta) ^ (-rho) * log(x / beta)) - mean(log(x / beta)) + 1 / rho
    rho.hat <- uniroot(g, interval = c(rho0 / 5, rho0 * 5), tol = 1e-09)$root
    beta.hat <- beta 
  }  
  if (beta != "unknown" & rho != "unknown") {
    beta.hat <- beta
    rho.hat <- rho 
  } 
  return(c(beta.hat, rho.hat))
}
Frechet.test.ML <- function(x, beta = "unknown", rho = "unknown"){
  n <- length(x)
  euler <- - digamma(1) 
  ML.beta.rho <-  ML.Frechet(x, beta, rho)
  beta.hat <- ML.beta.rho[1]
  rho.hat <- ML.beta.rho[2]
  Fi <-  exp(- (x / beta.hat) ^ (- rho.hat))  
  neg2L <- - 2 * sum(log(dFrechet(x, beta = beta.hat, rho = rho.hat)))
  G <- matrix(c(- rho.hat / beta.hat * h6(1,2,1), rho.hat / beta.hat * h7(1,2,1), - 1 / rho.hat * h8(1), 1 / rho.hat * h9(1)), nrow = 2) 
  R <- matrix(c(rho.hat ^ 2 / beta.hat ^ 2, 1 / beta.hat * (1 - euler),  1 / beta.hat * (1 - euler), 1 / rho.hat ^ 2 * ((euler - 1) ^ 2 + pi ^ 2 / 6)), nrow = 2)    
  U <- (1:2)[c(beta == "unknown", rho == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.beta.rho = ML.beta.rho, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Gompertz distribution

dGompertz <- function(x, beta, rho) beta * rho * exp(rho + beta * x - rho * exp(beta * x))
pGompertz <- function(x, beta, rho) 1 - exp(-rho * (exp(beta * x) - 1)) 
qGompertz <- function(p, beta, rho) log(1 - log(1 - p) / rho) / beta
rGompertz <- function(n, beta, rho) log(1 - log(1 - runif(n)) / rho) / beta
ML.Gompertz <- function(x, beta = "unknown", rho = "unknown") {
  if (beta == "unknown" & rho == "unknown") {
    xb <- mean(x)
    g.beta <- function(beta)  1 / beta + xb - mean(x * exp(beta * x)) / (mean(exp(beta * x)) -  1)
    beta.ML <- uniroot(g.beta, interval = c(0.0001, 100), tol = 1e-7)$root
    rho.ML <- 1 / (mean(exp(beta.ML * x)) - 1)
  }  
  if (beta == "unknown" & rho != "unknown") {
    xb <- mean(x)
    g.beta <- function(beta)  1 / beta + xb  - rho * mean(x * exp(beta * x))
    beta.ML <- uniroot(g.beta, interval = c(0.0001, 100), tol = 1e-7)$root
    rho.ML <- rho
  }
  if (beta != "unknown" & rho == "unknown") {
    rho.ML <- 1 / (mean(exp(beta * x)) - 1)
    beta.ML <- beta
  }
  if (beta != "unknown" & rho != "unknown") {
    beta.ML <- beta
    rho.ML <- rho
  }
  return(c(beta.ML, rho.ML))
}
Gompertz.test.ML <- function(x, beta = "unknown", rho = "unknown"){
  n <- length(x)
  ML.beta.rho <-  ML.Gompertz(x, beta, rho)
  beta.hat <- ML.beta.rho[1]
  rho.hat <- ML.beta.rho[2]
  h20.rho.hat <- h20(rho.hat)
  Fi <- 1 - exp(-rho.hat * (exp(beta.hat * x) - 1))  
  neg2L <- - 2 * sum(log(dGompertz(x, beta = beta.hat, rho = rho.hat)))
  G <- rho.hat * exp(rho.hat) * matrix(c(h21(rho.hat) / beta.hat, h22(rho.hat) / beta.hat, -h23(rho.hat), -h24(rho.hat)), ncol = 2)
  R <- matrix(c((1 + rho.hat ^ 2 * exp(rho.hat) * h19(rho.hat)) / beta.hat ^ 2, rho.hat * exp(rho.hat) * h20.rho.hat / beta.hat, rho.hat * exp(rho.hat) * h20.rho.hat / beta.hat, 1 / rho.hat ^ 2), ncol = 2)
  U <- (1:2)[c(beta == "unknown", rho == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.beta.rho = ML.beta.rho, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- log-logistic distribution

dloglogistic <- function(x, beta, rho) dllogis(x, rho, scale = beta)
ploglogistic <- function(x, beta, rho) pllogis(x, rho, scale = beta)
qloglogistic <- function(p, beta, rho) qllogis(p, rho, scale = beta)
rloglogistic <- function(n, beta, rho) rllogis(n, rho, scale = beta)
ML.loglogistic <- function(x, beta = "unknown", rho = "unknown"){
  lx <- log(x)
  if (beta == "unknown" & rho == "unknown") {
    beta.hat <- exp(mean(lx))
    rho.hat <- 1 / (sqrt(3) * sd(lx) / pi)
    g1 <- function(beta) {w <- (x / beta) ^ rho.hat;  mean(2 / (1 + w)) - 1}
    g2 <- function(rho) {w <- (x / beta.hat) ^ rho; mean(log(w)) - mean(2 * log(w) / (1 + w)) - 1}
    err.beta <- err.rho <- 1000000
    while (err.beta > 1e-10 | err.rho > 1e-10){
      beta.hat.old <- beta.hat
      rho.hat.old <- rho.hat
      beta.hat <- uniroot(g1, interval = c(beta.hat / 5, beta.hat * 5), tol = 1e-09)$root
      rho.hat <- uniroot(g2, interval = c(rho.hat / 5, rho.hat * 5), tol = 1e-09)$root
      err.beta <- abs(beta.hat - beta.hat.old)
      err.rho <- abs(rho.hat - rho.hat.old)
    }
  }  
  if (beta == "unknown" & rho != "unknown") {
    beta0 <- exp(mean(lx))
    g3 <- function(beta) {w <- (x / beta) ^ rho;  mean(2 / (1 + w)) - 1}
    beta.hat <- uniroot(g3, interval = c(beta0 / 5, beta0 * 5), tol = 1e-09)$root
    rho.hat <- rho 
  }  
  if (beta != "unknown" & rho == "unknown") {
    rho0 <- 1 / (sqrt(3) * sd(lx) / pi)
    g4 <- function(rho) {w <- (x / beta) ^ rho; mean(log(w)) - mean(2 * log(w) / (1 + w)) - 1}
    rho.hat <- uniroot(g4, interval = c(rho0 / 5, rho0 * 5), tol = 1e-09)$root
    beta.hat <- beta 
  }  
  if (beta != "unknown" & rho != "unknown") {
    beta.hat <- beta 
    rho.hat <- rho 
  } 
  return(c(beta.hat, rho.hat))
}
MM.loglogistic <- function(x, beta = "unknown", rho = "unknown"){
  logx <- log(x)
  if (beta == "unknown") 
    beta.hat <- exp(mean(logx)) else 
      beta.hat <- beta
    if (rho == "unknown"){
      rho.hat <- 1 / sqrt(3 / pi^2 * mean((logx - log(beta.hat)) ^ 2))} else
        rho.hat <- rho  
      return(c(beta.hat, rho.hat))
}
loglogistic.test.ML <- function(x, beta = "unknown", rho = "unknown"){
  n <- length(x)
  ML.beta.rho <- ML.loglogistic(x, beta, rho)
  beta.hat <- ML.beta.rho[1]
  rho.hat <- ML.beta.rho[2]
  Fi <- pllogis(x, shape = rho.hat, scale = beta.hat)
  neg2L <- - 2 * sum(log(dllogis(x, shape = rho.hat, scale = beta.hat)))
  G <- matrix(c(0, -rho.hat / beta.hat / pi, - 0.698397593884459 / rho.hat, 0), nrow = 2)
  R <- matrix(c(rho.hat ^ 2 / 3 / beta.hat ^ 2, 0, 0, (3 + pi ^ 2) / 9 / rho.hat ^ 2), nrow = 2)    
  U <- (1:2)[c(beta == "unknown", rho == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.beta.rho = ML.beta.rho, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
loglogistic.test.MM <- function(x, beta = "unknown", rho = "unknown"){
  n <- length(x)
  MM.beta.rho <- MM.loglogistic(x, beta, rho)
  beta.hat <- MM.beta.rho[1]
  rho.hat <- MM.beta.rho[2]
  Fi <- pllogis(x, shape = rho.hat, scale = beta.hat)
  neg2L <- - 2 * sum(log(dllogis(x, shape = rho.hat, scale = beta.hat)))
  G <- matrix(c(0, -rho.hat / beta.hat / pi, - 0.698397593884459 / rho.hat, 0), nrow = 2)
  R <- diag(c(3 * rho.hat ^ 2 / beta.hat ^ 2 / pi ^ 2, 5 / 4 / rho.hat ^ 2))
  J <- matrix(c(0, - rho.hat / beta.hat * 0.235854187, - 1 / rho.hat * 0.4909114316,0), nrow = 2) 
  U <- (1:2)[c(beta == "unknown", rho == "unknown")]
  G <- as.matrix(G[,U]); R <- R[U,U]; J <- as.matrix(J[,U])
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.beta.rho = MM.beta.rho, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Gamma distribution

dGamma <- function(x, alpha, beta) dgamma(x, shape = alpha, scale = beta)
pGamma <- function(x, alpha, beta) pgamma(x, shape = alpha, scale = beta)
qGamma <- function(p, alpha, beta) qgamma(p, shape = alpha, scale = beta)

ML.gamma <- function(x, alpha = "unknown", beta = "unknown") {
  if (alpha == "unknown" & beta == "unknown") {
    xbar <- mean(x)
    sn2 <- mean((x - xbar) ^ 2)
    alpha.MM <- xbar ^ 2 / sn2 
    g.alpha <- function(alpha)  digamma(alpha) + log(xbar / alpha) - mean(log(x))
    alpha.ML <- uniroot(g.alpha, interval = c(alpha.MM / 10, alpha.MM * 10), tol = 1e-7)$root
    beta.ML <- xbar / alpha.ML
  }
  if (alpha == "unknown" & beta != "unknown") {
    xbar <- mean(x)
    sn2 <- mean((x - xbar) ^ 2)
    alpha.MM <- xbar ^ 2 / sn2 
    g.alpha <- function(alpha)  digamma(alpha) + log(beta) - mean(log(x))
    alpha.ML <- uniroot(g.alpha, interval = c(alpha.MM / 10, alpha.MM * 10), tol = 1e-7)$root
    beta.ML <- beta
  }
  if (alpha != "unknown" & beta == "unknown") {
    beta.ML <- mean(x) / alpha
    alpha.ML <- alpha
  }
  if (alpha != "unknown" & beta != "unknown") {
    alpha.ML <- alpha
    beta.ML <- beta
  }
  return(c(alpha.ML, beta.ML))
}
gamma.test.ML <- function(x, alpha = "unknown", beta = "unknown"){
  n <- length(x)
  ML.alpha.beta <-  ML.gamma(x, alpha, beta)
  alpha.hat <- ML.alpha.beta[1]
  beta.hat <- ML.alpha.beta[2]
  Fi <- pgamma(x / beta.hat, shape = alpha.hat, scale = 1)
  neg2L <- - 2 * sum(log(dGamma(x, alpha = alpha.hat, beta = beta.hat)))
  G <- matrix(c(h10(alpha.hat), h11(alpha.hat), alpha.hat / beta.hat * h6(alpha.hat, alpha.hat + 1, 1), alpha.hat / beta.hat * h7(alpha.hat, alpha.hat + 1, 1)), ncol = 2)
  R <- matrix(c(trigamma(alpha.hat), 1 / beta.hat, 1 / beta.hat, alpha.hat / beta.hat ^ 2), ncol = 2)
  U <- (1:2)[c(alpha == "unknown", beta == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.alpha.beta = ML.alpha.beta, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- inverse-gamma distribution

dinversegamma <- function(x, alpha, beta) beta ^ alpha / gamma(alpha) / x ^ (alpha + 1) * exp(- beta / x)
pinversegamma <- function(x, alpha, beta) 1 - pgamma(beta / x, shape = alpha)
qinversegamma <- function(p, alpha, beta) 1 / qgamma(1 - p, shape = alpha) * beta
rinversegamma <- function(n, alpha, beta) 1 / rgamma(n, shape = alpha) * beta
ML.inversegamma <- function(x, alpha = "unknown", beta = "unknown") {
  if (alpha == "unknown" & beta == "unknown") {
    mean.1x <- mean(1 / x)
    mlx <- mean(log(x))
    alpha.MM <- mean.1x ^ 2 / mean((1 / x - mean.1x) ^ 2)
    g.alpha <- function(alpha)  digamma(alpha) - log(alpha) + log(mean.1x)  + mlx
    alpha.ML <- uniroot(g.alpha, interval = c(alpha.MM / 5, alpha.MM * 5), tol = 1e-7)$root
    beta.ML <-  alpha.ML / mean.1x
  }  
  if (alpha == "unknown" & beta != "unknown") {
    mean.1x <- mean(1 / x)
    mlx <- mean(log(x))
    alpha.MM <- mean.1x ^ 2 / mean((1 / x - mean.1x) ^ 2)
    g.alpha <- function(alpha)  digamma(alpha) - log(beta) + mlx
    alpha.ML <- uniroot(g.alpha, interval = c(alpha.MM / 5, alpha.MM * 5), tol = 1e-7)$root
    beta.ML <-  beta
  }
  if (alpha != "unknown" & beta == "unknown") {
    beta.ML <-  alpha / mean(1 / x)
    alpha.ML <- alpha
  }
  if (alpha != "unknown" & beta != "unknown") {
    alpha.ML <- alpha
    beta.ML <-  beta
  }
  return(c(alpha.ML, beta.ML))
}
inversegamma.test.ML <- function(x, alpha = "unknown", beta = "unknown"){
  n <- length(x)
  ML.alpha.beta <-  ML.inversegamma(x, alpha, beta)
  alpha.hat <- ML.alpha.beta[1]
  beta.hat <- ML.alpha.beta[2]
  Fi <- 1 - pgamma(beta.hat / x, shape = alpha.hat) 
  neg2L <- - 2 * sum(log(dinversegamma(x, alpha = alpha.hat, beta = beta.hat)))
  G <- matrix(c(h10(alpha.hat), -h11(alpha.hat), -alpha.hat / beta.hat * h6(alpha.hat, alpha.hat + 1, 1), alpha.hat / beta.hat * h7(alpha.hat, alpha.hat + 1, 1)), ncol = 2)
  R <- matrix(c(trigamma(alpha.hat), -1 / beta.hat,- 1 / beta.hat, alpha.hat / beta.hat ^ 2), ncol = 2)
  U <- (1:2)[c(alpha == "unknown", beta == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.alpha.beta = ML.alpha.beta, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Beta prime distribution

dbetaprime <- function(x, alpha, beta) gamma(alpha + beta) / gamma(alpha) / gamma(beta) * x ^ (alpha - 1) / (1 + x) ^ (alpha + beta)
pbetaprime <- function(x, alpha, beta) pbeta(x / (x + 1), alpha, beta)
qbetaprime <- function(p, alpha, beta) {y <- qbeta(p, alpha, beta); return(y / (1 - y))}
rbetaprime <- function(n, alpha, beta) {y <- rbeta(n, alpha, beta); return(y / (1 - y))}
ML.betaprime <- function(x, alpha = "unknown", beta = "unknown") {
  if (alpha == "unknown" & beta == "unknown") {
    y <- x / (x + 1)
    yb <- mean(y)
    vy <- var(y)
    mly <- mean(log(y))
    mly1 <- mean(log(x+1))
    alpha.ML <- yb * (yb * (1 - yb) / vy - 1)  # Method of moments, initial value
    beta.ML <- (1 - yb) * (yb * (1 - yb) / vy - 1) # Method of moments, initial value
    g.alpha <- function(alpha)  digamma(alpha + beta.ML) - digamma(alpha) + mly
    g.beta <- function(beta)  digamma(alpha.ML + beta) - digamma(beta) - mly1
    err.alpha <- err.beta <- 1000000
    while (err.alpha > 1e-10 | err.beta > 1e-10){
      alpha.ML.old <- alpha.ML
      beta.ML.old  <- beta.ML
      alpha.ML <- uniroot(g.alpha, interval = c(alpha.ML / 5, alpha.ML * 5), tol = 1e-09)$root
      beta.ML <- uniroot(g.beta, interval = c(beta.ML / 5, beta.ML * 5), tol = 1e-09)$root
      err.alpha <- abs(alpha.ML - alpha.ML.old)
      err.beta <- abs(beta.ML - beta.ML.old)
    }
  }  
  if (alpha == "unknown" & beta != "unknown") {
    y <- x / (x + 1)
    yb <- mean(y)
    vy <- var(y)
    mly <- mean(log(y))
    alpha0 <- yb * (yb * (1 - yb) / vy - 1)  # Method of moments, initial value
    g.alpha <- function(alpha)  digamma(alpha + beta) - digamma(alpha) + mly
    alpha.ML <- uniroot(g.alpha, interval = c(alpha0 / 5, alpha0 * 5), tol = 1e-7)$root
    beta.ML <- beta    
  }
  if (alpha != "unknown" & beta == "unknown") {
    y <- x / (x + 1)
    yb <- mean(y)
    vy <- var(y)
    mly1 <- mean(log(x+1))
    beta0 <- (1 - yb) * (yb * (1 - yb) / vy - 1) # Method of moments, initial value
    g.beta <- function(beta)  digamma(alpha + beta) - digamma(beta) - mly1
    beta.ML <- uniroot(g.beta, interval = c(beta0 / 5, beta0 * 5), tol = 1e-7)$root
    alpha.ML <- alpha    
  }
  if (alpha != "unknown" & beta != "unknown") {
    alpha.ML <- alpha    
    beta.ML <- beta    
  }
  return(c(alpha.ML, beta.ML))
}
betaprime.test.ML <- function(x, alpha = "unknown", beta = "unknown"){
  n <- length(x)
  ML.alpha.beta <-  ML.betaprime(x, alpha, beta)
  alpha.hat <- ML.alpha.beta[1]
  beta.hat <- ML.alpha.beta[2]
  Fi <- pbetaprime(x, alpha.hat, beta.hat)
  neg2L <- - 2 * sum(log(dbetaprime(x, alpha.hat, beta.hat)))
  G <- matrix(c(h25(alpha.hat, beta.hat), h26(alpha.hat, beta.hat), h27(alpha.hat, beta.hat), h28(alpha.hat, beta.hat)), ncol = 2)
  R <- matrix(c(trigamma(alpha.hat) - trigamma(alpha.hat + beta.hat), - trigamma(alpha.hat + beta.hat), - trigamma(alpha.hat + beta.hat), trigamma(beta.hat) - trigamma(alpha.hat + beta.hat)), ncol = 2)
  U <- (1:2)[c(alpha == "unknown", beta == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.alpha.beta = ML.alpha.beta, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Lomax distribution, (which is a special case of the Pareto distribution of type II with a location fixed to 0)
dLomax <- function(x, alpha, sigma) alpha / sigma / (1 + x / sigma) ^ (1 + alpha)
pLomax <- function(x, alpha, sigma) 1 - (1 + x / sigma) ^ (- alpha)
qLomax <- function(p, alpha, sigma) sigma * ((1 - p) ^ (- 1 / alpha) - 1)
rLomax <- function(n, alpha, sigma) sigma * ((1 - runif(n)) ^ (- 1 / alpha) - 1)
ML.Lomax <- function(x, alpha = "unknown", sigma = "unknown") {
  if (alpha == "unknown" & sigma == "unknown") {
    vx <- var(x)
    mx <- mean(x)
    alpha.MM <- 2 * vx / (vx - mx ^ 2)
    sigma.MM <-  mx * (alpha.MM - 1)
    g.sigma <- function(sigma)  (1 + (mean(log(1 + x / sigma))) ^ (- 1)) * mean(x / (1 + x / sigma)) - sigma
    sigma.ML <- uniroot(g.sigma, interval = c(sigma.MM / 100, sigma.MM * 100), tol = 1e-7)$root
    alpha.ML <- 1 /  mean(log(1 + x / sigma.ML))
  }
  if (alpha == "unknown" & sigma != "unknown") {
    alpha.ML <- 1 /  mean(log(1 + x / sigma))
    sigma.ML <- sigma
  }    
  if (alpha != "unknown" & sigma == "unknown") {
    sigma.MM <- mean(x) * (alpha - 1)
    g.sigma <- function(sigma)  (alpha + 1) * mean(x / (1 + x / sigma)) - sigma
    sigma.ML <- uniroot(g.sigma, interval = c(sigma.MM / 100, sigma.MM * 100), tol = 1e-7)$root
    alpha.ML <- alpha
  }    
  if (alpha != "unknown" & sigma != "unknown") {
    alpha.ML <- alpha
    sigma.ML <- sigma
  }    
  return(c(alpha.ML, sigma.ML))
}
Lomax.test.ML <- function(x, alpha = "unknown", sigma = "unknown"){
  n <- length(x)
  ML.alpha.sigma <-  ML.Lomax(x, alpha, sigma)
  alpha.hat <- ML.alpha.sigma[1]
  sigma.hat <- ML.alpha.sigma[2]
  Fi <- pLomax(x, alpha.hat, sigma.hat)
  neg2L <- - 2 * sum(log(dLomax(x, alpha.hat, sigma.hat)))
  G <- matrix(c(-h6(1,2,1) / alpha.hat, -h7(1,2,1) / alpha.hat,
                -alpha.hat / sigma.hat * h6(1, 1, alpha.hat / (alpha.hat + 1)), 
                -alpha.hat / sigma.hat * h7(1, 1, alpha.hat / (alpha.hat + 1))), ncol = 2)
  R <- matrix(c(1 / alpha.hat ^ 2, rep(-1 / (alpha.hat + 1) / sigma.hat, 2), alpha.hat / (alpha.hat + 2) / sigma.hat ^ 2), ncol = 2)
  U <- (1:2)[c(alpha == "unknown", sigma == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.alpha.sigma = ML.alpha.sigma, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Nakagami distribution

dNakagami <- function(x, lambda, omega) 2 / gamma(lambda) * (lambda / omega) ^ lambda * x ^ (2 * lambda - 1) * exp(- lambda / omega * x ^ 2)
pNakagami <- function(x, lambda, omega) pgamma(lambda / omega * x ^ 2 , shape = lambda, scale = 1)
qNakagami <- function(p, lambda, omega) sqrt(omega / lambda * qgamma(p, shape = lambda, scale = 1)) 
rNakagami <- function(n, lambda, omega) sqrt(omega / lambda * rgamma(n, shape = lambda, scale = 1)) 
ML.Nakagami <- function(x, lambda = "unknown", omega = "unknown") {
  if (lambda == "unknown" & omega == "unknown") {
    omega.ML <- mean(x ^ 2)
    mlx2 <- mean(log(x ^ 2))
    lambda.MM <- omega.ML ^ 2 / mean((x ^ 2 - omega.ML) ^ 2)
    g.lambda <- function(lambda)  digamma(lambda) - log(lambda) + log(omega.ML)  - mlx2
    lambda.ML <- uniroot(g.lambda, interval = c(lambda.MM / 5, lambda.MM * 5), tol = 1e-09)$root
  }
  if (lambda == "unknown" & omega != "unknown") {
    mx2 <- mean(x ^ 2)
    mlx2 <- mean(log(x ^ 2))
    lambda0 <- omega ^ 2 / mean((x ^ 2 - omega) ^ 2)
    g.lambda <- function(lambda)   digamma(lambda) - log(lambda)  + log(omega) - mlx2 + mx2 / omega  - 1
    lambda.ML <- uniroot(g.lambda, interval = c(lambda0 / 5, lambda0 * 5), tol = 1e-09)$root
    omega.ML <- omega
  }
  if (lambda != "unknown" & omega == "unknown") {
    omega.ML <- mean(x ^ 2)
    lambda.ML <- lambda
  }
  if (lambda != "unknown" & omega != "unknown") {
    lambda.ML <- lambda
    omega.ML <- omega
  }
  return(c(lambda.ML, omega.ML))
}
Nakagami.test.ML <- function(x, lambda = "unknown", omega = "unknown"){
  n <- length(x)
  ML.lambda.omega <-  ML.Nakagami(x, lambda, omega)
  lambda.hat <- ML.lambda.omega[1]
  omega.hat <- ML.lambda.omega[2]
  Fi <- pgamma(lambda.hat / omega.hat * x ^ 2 , shape = lambda.hat, scale = 1)
  neg2L <- - 2 * sum(log(dNakagami(x, lambda = lambda.hat, omega = omega.hat)))
  h4.lambda.hat <- h6(lambda.hat, lambda.hat + 1, 1)
  h5.lambda.hat <- h7(lambda.hat, lambda.hat + 1, 1)
  G <- matrix(c(h10(lambda.hat) - h4.lambda.hat, h11(lambda.hat) - h5.lambda.hat, lambda.hat / omega.hat * h4.lambda.hat, lambda.hat / omega.hat * h5.lambda.hat), ncol = 2)
  R <- matrix(c(trigamma(lambda.hat) - 1 / lambda.hat, 0, 0, lambda.hat / omega.hat ^ 2), ncol = 2)
  U <- (1:2)[c(lambda == "unknown", omega == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.lambda.omega = ML.lambda.omega, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- inverse-Gaussian distribution

dinversegaussian <- function(x, mu = 1, lambda = 1, log = FALSE) {
  log.f <-  log(lambda)/2  - log(2 * pi)/2 - 3 / 2 * log(x) - lambda * (x - mu) ^ 2 / 2 / mu ^ 2 / x
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pinversegaussian <- function(x, mu = 1, lambda = 1) {
  z <- sqrt(lambda / x)     
  F <- pnorm(z * (x / mu - 1)) + exp(2 * lambda / mu) * pnorm(- z * (x / mu + 1)) 
  return(F)
}
qinversegaussian <- function(p, mu = 1, lambda = 1) qinvgauss(p, mean = mu, shape = lambda)
rinversegaussian <- function(n, mu = 1, lambda = 1) {
  y1 <- rnorm(n) ^ 2
  x <- mu  + mu ^ 2 * y1 / 2 / lambda - mu / 2 / lambda * sqrt(4 * mu * lambda * y1 + mu ^ 2 * y1 ^ 2)
  x2 <- mu ^ 2 / x
  pos <- runif(n) > mu / (mu + x)
  x[pos] <- x2[pos] 
  return(x)
}
ML.inversegaussian <- function(x, mu = "unknown", lambda = "unknown") {
  if (mu == "unknown" & lambda == "unknown") {
    mu.ML <- mean(x)
    lambda.ML <- 1 / (mean(1 / x) - 1 / mu.ML)
  }
  if (mu == "unknown" & lambda != "unknown") {
    mu.ML <- mean(x)
    lambda.ML <- lambda
  }
  if (mu != "unknown" & lambda == "unknown") {
    lambda.ML <- 1 / (mean(1 / x) + mean(x) / mu ^ 2 - 2 / mu)
    mu.ML <- mu
  }
  if (mu != "unknown" & lambda != "unknown") {
    mu.ML <- mu
    lambda.ML <- lambda
  }
  return(c(mu.ML, lambda.ML))
}
inversegaussian.test.ML <- function(x, mu = "unknown", lambda = "unknown"){
  n <- length(x)
  ML.mu.lambda <-  ML.inversegaussian(x, mu, lambda)
  mu.hat <- ML.mu.lambda[1]
  lambda.hat <- ML.mu.lambda[2]
  Fi <- pinvgauss(x, mean = mu.hat, shape = lambda.hat)
  neg2L <- - 2 * sum(log(dinvgauss(x, mean = mu.hat, shape = lambda.hat)))
  G <- matrix(c(lambda.hat / mu.hat ^ 3 * h29(mu.hat, lambda.hat), lambda.hat / mu.hat ^ 3 * h30(mu.hat, lambda.hat), - h31(mu.hat, lambda.hat) / 2 / mu.hat ^ 2, - h32(mu.hat, lambda.hat) / 2 / mu.hat ^ 2), ncol = 2)
  R <- matrix(c(lambda.hat / mu.hat ^ 3, 0, 0, 1 / 2 / lambda.hat ^ 2), ncol = 2)
  U <- (1:2)[c(mu == "unknown", lambda == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.mu.lambda = ML.mu.lambda, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- exponential distribution

dExp <- function(x, beta) dexp(x, rate = 1 / beta)
pExp <- function(x, beta) pexp(x, rate = 1 / beta)
qExp <- function(p, beta) qexp(p, rate = 1 / beta)
rExp <- function(n, beta) rexp(n, rate = 1 / beta)

ML.exp <- function(x, beta = "unknown"){
  if (beta == "unknown") beta.hat <- mean(x) else beta.hat <- beta
  return(beta.hat)
}
exp.test.ML <- function(x, beta = "unknown"){
  n <- length(x)
  beta.hat <- ML.exp(x, beta)
  Fi <-  1 - exp(- x / beta.hat)
  neg2L <- - 2 * sum(log(dExp(x, beta = beta.hat)))
  G <- 1 / beta.hat * matrix(c(h6(1,2,1), h7(1,2,1)), nrow = 2) 
  R <- 1 / beta.hat ^ 2
  U <- c(1)[beta == "unknown"] 
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.beta = beta.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- half-normal distribution

dhalfnormal <- function(x, delta, log = FALSE){
  log.f <- log(2 / pi) / 2 - log(delta) - (x / delta) ^ 2 / 2 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
phalfnormal <- function(x, delta) 2 * pnorm(x / delta) - 1
qhalfnormal <- function(p, delta) delta * qnorm((1 + p) / 2)
rhalfnormal <- function(n, delta) abs(rnorm(n, 0, delta))
ML.halfnormal <- function(x, delta = "unknown"){
  if (delta == "unknown") delta.hat <- sqrt(mean(x ^ 2)) else delta.hat <- delta
  return(delta.hat)
}
MM.halfnormal <- function(x, delta = "unknown"){
  if (delta == "unknown") delta.hat <- sqrt(pi / 2) * mean(x) else delta.hat <- delta
  return(delta.hat)
}
halfnormal.test.ML <- function(x, delta = "unknown"){
  n <- length(x)
  delta.hat <- ML.halfnormal(x, delta) 
  Fi <- 2 * pnorm(x / delta.hat) -  1 
  neg2L <- - 2 * sum(log(dhalfnormal(x, delta = delta.hat)))
  G <- 1 / delta.hat * matrix(c(h6(1 / 2, 3 / 2, 1), h7(1 / 2, 3 / 2, 1)), nrow = 2) 
  R <- 2 / delta.hat ^ 2
  U <- c(1)[delta == "unknown"] 
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.delta = delta.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
halfnormal.test.MM <- function(x, delta = "unknown"){
  n <- length(x)
  delta.hat <- MM.halfnormal(x, delta) 
  Fi <- 2 * pnorm(x / delta.hat) -  1 
  neg2L <- - 2 * sum(log(dhalfnormal(x, delta = delta.hat)))
  G <- 1 / delta.hat * matrix(c(h6(1 / 2, 3 / 2, 1), h7(1 / 2, 3 / 2, 1)), nrow = 2) 
  R <- 1 / (pi / 2 - 1) / delta.hat ^ 2
  J <- 1 / (pi / 2 - 1)  / delta.hat * matrix(c(h6(1/2,1,1), h7(1/2,1,1)), nrow = 2) 
  U <- c(1)[delta == "unknown"] 
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.delta = delta.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Rayleigh distribution

dRayleigh <- function(x, delta, log = FALSE){
  log.f <- - 2 * log(delta) + log(x) - (x / delta) ^ 2 / 2 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pRayleigh <- function(x, delta) 1 - exp(- (x / delta) ^ 2 / 2)
qRayleigh <- function(p, delta){
  y <- qexp(p, rate = 1 / (2 * delta ^ 2))
  x <- sqrt(y)
  return(x)
}
rRayleigh <- function(n, delta){
  y <- rexp(n, rate = 1 / (2 * delta ^ 2))
  x <- sqrt(y)
  return(x)
}
ML.Rayleigh <- function(x, delta = "unknown"){
  if (delta == "unknown") delta.hat <- (mean(x ^ 2) / 2) ^ (1 / 2) else delta.hat <- delta
  return(delta.hat)
}
MM.Rayleigh <- function(x, delta = "unknown"){
  if (delta == "unknown") delta.hat <- sqrt(2 / pi) * mean(x) else delta.hat <- delta
  return(delta.hat)
}
Rayleigh.test.ML <- function(x, delta = "unknown"){
  n <- length(x)
  delta.hat <- ML.Rayleigh(x, delta) 
  Fi <- 1 - exp(- (x / delta.hat) ^  2 / 2) 
  neg2L <- - 2 * sum(log(dRayleigh(x, delta = delta.hat)))
  G <- 2 / delta.hat * matrix(c(h6(1,2,1), h7(1,2,1)), nrow = 2) 
  R <- 4 / delta.hat ^ 2
  U <- c(1)[delta == "unknown"] 
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.delta = delta.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
Rayleigh.test.MM <- function(x, delta = "unknown"){
  n <- length(x)
  delta.hat <- MM.Rayleigh(x, delta) 
  Fi <- 1 - exp(- (x / delta.hat) ^  2 / 2) 
  neg2L <- - 2 * sum(log(dRayleigh(x, delta = delta.hat)))
  G <- 2 / delta.hat * matrix(c(h6(1,2,1), h7(1,2,1)), nrow = 2) 
  R <- 1 / (4 / pi - 1) / delta.hat ^ 2
  J <- 1 / (4 / pi - 1)  / delta.hat * matrix(c(h6(1,3/2,1), h7(1,3/2,1)), nrow = 2) 
  U <- c(1)[delta == "unknown"] 
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.delta = delta.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Maxwell-Boltzmann distribution

dMaxwell <- function(x, delta, log = FALSE){
  log.f <- log(2 / pi) / 2 - 3 * log(delta) + 2 * log(x) - (x / delta) ^ 2 / 2 
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
pMaxwell <- function(x, delta) pgamma((x / delta) ^ 2 / 2, shape = 3 / 2, scale = 1)
qMaxwell <- function(p, delta){
  y <- qgamma(p, shape = 3 / 2, scale = 2 * delta ^ 2)
  x <- sqrt(y)
  return(x)
}
rMaxwell <- function(n, delta){
  y <- rgamma(n, shape = 3 / 2, scale = 2 * delta ^ 2)
  x <- sqrt(y)
  return(x)
}
ML.Maxwell <- function(x, delta = "unknown"){
  if (delta == "unknown") delta.hat <- (mean(x ^ 2) / 3) ^ (1 / 2) else delta.hat <- delta
  return(delta.hat)
}
MM.Maxwell <- function(x, delta = "unknown"){
  if (delta == "unknown") delta.hat <- sqrt(pi/8) * mean(x) else delta.hat <- delta
  return(delta.hat)
}
Maxwell.test.ML <- function(x, delta = "unknown"){
  n <- length(x)
  delta.hat <- ML.Maxwell(x, delta) 
  Fi <- pgamma((x / delta.hat) ^  2 / 2, shape = 3 / 2, scale = 1)
  neg2L <- - 2 * sum(log(dMaxwell(x, delta.hat)))
  G <- 3  / delta.hat * matrix(c(h6(3 / 2, 5 / 2, 1), h7(3 / 2, 5 / 2, 1)), nrow = 2) 
  R <- 6 / delta.hat ^ 2
  U <- c(1)[delta == "unknown"] 
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.delta = delta.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
Maxwell.test.MM <- function(x, delta = "unknown"){
  n <- length(x)
  delta.hat <- MM.Maxwell(x, delta) 
  Fi <- pgamma((x / delta.hat) ^  2 / 2, shape = 3 / 2, scale = 1)
  neg2L <- - 2 * sum(log(dMaxwell(x, delta.hat)))
  G <- 3  / delta.hat * matrix(c(h6(3 / 2, 5 / 2, 1), h7(3 / 2, 5 / 2, 1)), nrow = 2) 
  R <- 1 / (3 * pi / 8 - 1) / delta.hat ^ 2
  J <- 1 / (3 * pi / 8 - 1)  / delta.hat * matrix(c(h6(3/2,2,1), h7(3/2,2,1)), nrow = 2) 
  U <- c(1)[delta == "unknown"] 
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.delta = delta.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Chi-squared distribution

ML.chisquared <- function(x, k = "unknown") {
  if (k == "unknown") {
    xbar <- mean(x)
    mlx <- mean(log(x))
    g.k <- function(k)  digamma(k / 2) + log(2) - mlx
    k.ML <- uniroot(g.k, interval = c(xbar / 5, xbar * 5), tol = 1e-7)$root
  } else 
    k.ML <- k
  return(k.ML)
}
MM.chisquared <- function(x, k = "unknown") {
  if (k == "unknown") {
    k.MM <- mean(x)
  } else 
    k.MM <- k
  return(k.MM)
}
chisquared.test.ML <- function(x, k = "unknown"){
  n <- length(x)
  k.hat <-  ML.chisquared(x, k)
  Fi <-   pchisq(x, k.hat)
  neg2L <- - 2 * sum(log(dchisq(x, k.hat)))
  G <- 1 / 2 * matrix(c(h10(k.hat / 2), h11(k.hat / 2)), ncol = 1)
  R <- trigamma(k.hat / 2) / 4
  U <- c(1)[k == "unknown"] 
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.k = k.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
chisquared.test.MM <- function(x, k = "unknown"){
  n <- length(x)
  k.hat <-  MM.chisquared(x, k)
  Fi <-   pchisq(x, k.hat)
  neg2L <- - 2 * sum(log(dchisq(x, k.hat)))
  G <- 1 / 2 * matrix(c(h10(k.hat / 2), h11(k.hat / 2)), ncol = 1)
  R <- 1 / (2 * k.hat)
  J <- 1 / 2 * matrix(c(h6(k.hat/2, k.hat/2 + 1,1), h7(k.hat/2, k.hat/2 + 1,1)), nrow = 2) 
  U <- c(1)[k == "unknown"] 
  if (length(U) > 0)  {Ri <- solve(R);  Sigma <- diag(2)/2 - G %*% Ri  %*% t(J) - 
    J %*% Ri  %*% t(G) + G %*% Ri  %*% t(G)} else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(MM.k = k.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Distributions on (1, infinity)

## ---- Pareto distribution

dPareto <- function(x, alpha = 1.0, log = FALSE){
  indic.domain <- as.numeric(x >= 1)
  x[indic.domain == 0] <- 1 # to avoid log(x) = NaN if x < 0
  log.f <- log(indic.domain) + log(alpha)  - (1 + alpha) * log(x)
  return(if (log == FALSE) {exp(log.f)} else {log.f})
}
qPareto <- function(p, alpha = 1.0){
  return(1 / (1 - p) ^ (1 / alpha))
}
rPareto <- function(n, alpha = 1.0){
  return(1 / (1 - runif(n)) ^ (1 / alpha))
}
pPareto <- function(x, alpha = 1.0, lower.tail = TRUE){
  indic.domain <- as.numeric(x >= 1)
  x[indic.domain == 0] <- 1 # to avoid (1 / x) ^ alpha = NaN if x <= 0
  cdf <- (1 - (1 / x) ^ alpha) * indic.domain
  return(if (lower.tail == TRUE) {cdf} else {1 - cdf})
}
ML.Pareto <- function(x, alpha = "unknown"){
  if (alpha == "unknown") alpha.hat <- 1 / mean(log(x)) else alpha.hat <- alpha
  return(alpha.hat)
}  
Pareto.test.ML <- function(x, alpha = "unknown"){
  n <- length(x)
  alpha.hat <- ML.Pareto(x, alpha)
  Fi <-  1 - (1 / x) ^ alpha.hat
  neg2L <- - 2 * sum(log(dPareto(x, alpha = alpha.hat)))
  G <- -1 / alpha.hat * matrix(c(h6(1,2,1), h7(1,2,1)), nrow = 2) 
  R <- 1 / alpha.hat ^ 2
  U <- c(1)[alpha == "unknown"] 
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(alpha.hat = alpha.hat, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Distributions on (0,1)

## ---- beta distribution

ML.beta <- function(x, alpha = "unknown", beta = "unknown") {
  if (alpha == "unknown" & beta == "unknown") {
    xb <- mean(x)
    vx <- var(x)
    mlx <- mean(log(x))
    mlx1 <- mean(log(1 - x))
    alpha.ML <- xb * (xb * (1 - xb) / vx - 1)  # Method of moments, initial value
    beta.ML <- (1 - xb) * (xb * (1 - xb) / vx - 1) # Method of moments, initial value
    g.alpha <- function(alpha)  mlx  - digamma(alpha) + digamma(alpha + beta.ML)
    g.beta <- function(beta)  mlx1  - digamma(beta) + digamma(alpha.ML + beta)
    err.alpha <- err.beta <- 1000000
    while (err.alpha > 1e-10 | err.beta > 1e-10){
      alpha.ML.old <- alpha.ML
      beta.ML.old  <- beta.ML
      alpha.ML <- uniroot(g.alpha, interval = c(alpha.ML / 5, alpha.ML * 5), tol = 1e-7)$root
      beta.ML <- uniroot(g.beta, interval = c(beta.ML / 5, beta.ML * 5), tol = 1e-7)$root
      err.alpha <- abs(alpha.ML - alpha.ML.old)
      err.beta <- abs(beta.ML - beta.ML.old)
    }
  }
  if (alpha == "unknown" & beta != "unknown") {
      xb <- mean(x)
      vx <- var(x)
      mlx <- mean(log(x))
      alpha0 <- xb * (xb * (1 - xb) / vx - 1)  # Method of moments, initial value
      g.alpha <- function(alpha)  mlx  - digamma(alpha) + digamma(alpha + beta)
      alpha.ML <- uniroot(g.alpha, interval = c(alpha0 / 5, alpha0 * 5), tol = 1e-7)$root
      beta.ML <- beta
  }      
  if (alpha != "unknown" & beta == "unknown") {
      xb <- mean(x)
      vx <- var(x)
      mlx1 <- mean(log(1 - x))
      beta0 <- (1 - xb) * (xb * (1 - xb) / vx - 1) # Method of moments, initial value
      g.beta <- function(beta)  mlx1  - digamma(beta) + digamma(alpha + beta)
      beta.ML <- uniroot(g.beta, interval = c(beta0 / 5, beta0 * 5), tol = 1e-7)$root
      alpha.ML <- alpha
  }      
  if (alpha != "unknown" & beta != "unknown") {
      alpha.ML <- alpha
      beta.ML <- beta
  }      
  return(c(alpha.ML, beta.ML))
}
beta.test.ML <- function(x, alpha = "unknown", beta = "unknown"){
  n <- length(x)
  ML.alpha.beta <-  ML.beta(x, alpha, beta)
  alpha.hat <- ML.alpha.beta[1]
  beta.hat <- ML.alpha.beta[2]
  Fi <- pbeta(x, alpha.hat, beta.hat)
  neg2L <- - 2 * sum(log(dbeta(x, alpha.hat, beta.hat)))
  G <- matrix(c(h25(alpha.hat, beta.hat), h26(alpha.hat, beta.hat), h27(alpha.hat, beta.hat), h28(alpha.hat, beta.hat)), ncol = 2)
  R <- matrix(c(trigamma(alpha.hat) - trigamma(alpha.hat + beta.hat), - trigamma(alpha.hat + beta.hat), - trigamma(alpha.hat + beta.hat), trigamma(beta.hat) - trigamma(alpha.hat + beta.hat)), ncol = 2)
  U <- (1:2)[c(alpha == "unknown", beta == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.alpha.beta = ML.alpha.beta, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- Kumaraswamy distribution

dKumaraswamy <- function(x, alpha, beta) alpha * beta * x ^ (alpha - 1) * (1 - x ^ alpha) ^ (beta - 1) 
pKumaraswamy <- function(x, alpha, beta) 1 - (1 - x ^ alpha) ^ beta 
qKumaraswamy <- function(p, alpha, beta)  (1 - (1 - p) ^ (1 / beta)) ^ (1 / alpha)
rKumaraswamy <- function(n, alpha, beta)  (1 - (1 - runif(n)) ^ (1 / beta)) ^ (1 / alpha)
ML.Kumaraswamy <- function(x, alpha = "unknown", beta = "unknown") {
  if (alpha == "unknown" & beta == "unknown") {
    g.alpha <- function(alpha)  1 + 1 / mean(log(1 - x ^ alpha)) + (1 / alpha + mean(log(x))) / mean(x ^ alpha *log(x) / (1 - x ^ alpha))
    alpha.ML <- uniroot(g.alpha, interval = c(0.0001, 100), tol = 1e-7)$root
    beta.ML <- -1 /  mean(log(1 - x ^ alpha.ML)) 
  }
  if (alpha == "unknown" & beta != "unknown") {
    g.alpha <- function(alpha)  1 - beta + (1 / alpha + mean(log(x))) / mean(x ^ alpha *log(x) / (1 - x ^ alpha))
    alpha.ML <- uniroot(g.alpha, interval = c(0.0001, 100), tol = 1e-7)$root
    beta.ML <- beta
  }
  if (alpha != "unknown" & beta == "unknown") {
    beta.ML <- -1 /  mean(log(1 - x ^ alpha)) 
    alpha.ML <- alpha
  }
  if (alpha != "unknown" & beta != "unknown") {
    alpha.ML <- alpha
    beta.ML <- beta
  }
  return(c(alpha.ML, beta.ML))
}
Kumaraswamy.test.ML <- function(x, alpha = "unknown", beta = "unknown"){
  n <- length(x)
  ML.alpha.beta <-  ML.Kumaraswamy(x, alpha, beta)
  alpha.hat <- ML.alpha.beta[1]
  beta.hat <- ML.alpha.beta[2]
  euler <- -digamma(1)
  Fi <- pKumaraswamy(x, alpha.hat, beta.hat)
  neg2L <- - 2 * sum(log(dKumaraswamy(x, alpha.hat, beta.hat)))
  G <- beta.hat * matrix(c(h33(beta.hat) / alpha.hat, h34(beta.hat) / alpha.hat, h35(beta.hat), h36(beta.hat)), ncol = 2)
  R11 <- ifelse(beta.hat != 2, 1 / alpha.hat ^ 2 + beta.hat / alpha.hat ^ 2 / (beta.hat - 2) * 
                  ((digamma(beta.hat) + euler - 1) ^ 2 - trigamma(beta.hat) + 
                     pi ^ 2 / 6 - 1), (4 * 1.20205690315959 - 3) / alpha.hat ^ 2)
  R12 <- ifelse(beta.hat != 1, (digamma(beta.hat)  + euler - 1 + 1 / beta.hat) / alpha.hat / (1 - beta.hat), (1 - pi ^ 2 / 6) / alpha.hat)
  R <- matrix(c(R11, R12, R12, 1 / beta.hat ^ 2), ncol = 2)
  U <- (1:2)[c(alpha == "unknown", beta == "unknown")] 
  G <- as.matrix(G[,U]); R <- R[U,U]
  if (length(U) > 0)  Sigma <- diag(2)/2 - G %*% solve(R) %*% t(G) else Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.alpha.beta = ML.alpha.beta, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}

## ---- uniform distribution on (a,b)

duniform <- function(x, a = 0, b = 1) dunif(x, a, b)
puniform <- function(x, a = 0, b = 1) punif(x, a, b)
quniform <- function(p, a = 0, b = 1) qunif(p, a, b)
runiform <- function(n, a = 0, b = 1) runif(n, a, b)
ML.uniform <- function(x, a = "unknown", b = "unknown") {
  if (a == "unknown" & b == "unknown") {
    a.ML <- min(x)
    b.ML <- max(x)
  }
  if (a == "unknown" &  b != "unknown") {
    a.ML <- min(x)
    b.ML <- b
  }
  if (a != "unknown" &  b == "unknown") {
    a.ML <- a
    b.ML <- max(x)
  }
  if (a != "unknown" &  b != "unknown") {
    a.ML <- a
    b.ML <- b
  }
  return(c(a.ML, b.ML))
}
uniform.test.ML <- function(x, a = "unknown", b = "unknown"){
  n <- length(x)
  ML.a.b <-  ML.uniform(x, a, b)
  a.hat <- ML.a.b[1]
  b.hat <- ML.a.b[2]
  Fi <- puniform(x, a.hat, b.hat)
  neg2L <- - 2 * sum(log(duniform(x, a.hat, b.hat)))
  Sigma <- diag(2)/2
  Sigma.inv <- solve(Sigma)  
  tau.bar.sqrtn <- sqrt(n) * matrix(c(mean(cos(2 * pi * Fi)), mean(sin(2 * pi * Fi))), ncol = 1)
  rownames(tau.bar.sqrtn) = c("sqrt(n) Cn", "sqrt(n) Sn")
  Tn <- as.vector(t(tau.bar.sqrtn) %*% Sigma.inv %*% tau.bar.sqrtn)
  Z.tau.bar <- tau.bar.sqrtn / sqrt(diag(Sigma)) 
  rownames(Z.tau.bar) = c("Z(Cn)", "Z(Sn)")
  p.value <- pchisq(Tn, df = 2, lower.tail = FALSE)
  return(list(ML.a.b = ML.a.b, Sigma.inv = Sigma.inv, neg2L = neg2L, tau.bar.sqrtn = tau.bar.sqrtn, Z.tau.bar = Z.tau.bar, Tn = Tn, p.value = p.value))
}
