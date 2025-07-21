# ====================================================================================
# Script : Test_Trigo_File1.R
# Purpose: Main code for the article:
#          "An omnibus goodness-of-fit test based on trigonometric moments"
# Note   : This file demonstrates how to use the functions defined in Test_Trigo_File2.R.
# Author : Alain Desgagné
# Date   : 2005-07-21
# ====================================================================================

# ====================================================================================
# Required R packages and function definitions
# ====================================================================================

rm(list=ls())
setwd("C:/Users/desgagne_a/Documents/4Recherche/Articles/Tests LK")
source(file = 'Test_Trigo_File2.R', echo = TRUE)
options(digits = 10)   
library(cubature) #adaptIntegrate function
library(ellipse)
library(flexsurv) # GG, log-logistic
library(nleqslv) # nleqslv
library(statmod) # Inverse Gaussian distribution
library(VGAM) # laplace, Frechet and Gompertz

# ====================================================================================
# Figure 1: The functions cos(2 pi u) and sin(2 pi u) plotted on the interval [0,1]
# ====================================================================================

# pdf("plot_cos_sin.pdf")

cos.f <- function(u) cos(2 * pi * u)
sin.f <- function(u) sin(2 * pi * u)
u <- seq(from = 0, to = 1, length = 100)
u1 <- seq(from = 0, to = 0.25, length = 100)
u2 <- seq(from = 0.25, to = 0.50, length = 100)
u3 <- seq(from = 0.50, to = 0.75, length = 100)
u4 <- seq(from = 0.75, to = 1, length = 100)

par(mfrow = c(1,1))
par(cex = 1.7)
par(cex.axis = 1, cex.lab = 1, cex.main = 1)
par(mgp = c(1.5, 0.5, 0))

plot(u, cos.f(u), type = 'l',  col = "white", lwd = 3, lty = 1,
     ylab = '', xlab = '', xaxt = 'n', yaxt = 'n')
title(xlab= 'u', line=1.3, cex.lab=1.2)
axis(2, at = c(-1, 0, 1), labels = c('-1','0','1'))
axis(1, at = c(0.25, 0.5, 0.75), labels = c('0.25','0.50','0.75'))
axis(1, at = c(0,1), labels = c('0','1'))

lines(u1, cos.f(u1),  col = "black", lwd = 3, lty = 1)
lines(u2, cos.f(u2),  col = "grey65", lwd = 3, lty = 1)
lines(u3, cos.f(u3),  col = "grey65", lwd = 3, lty = 1)
lines(u4, cos.f(u4),  col = "black", lwd = 3, lty = 1)
lines(u1, sin.f(u1),  col = "black", lwd = 3, lty = 2)
lines(u2, sin.f(u2),  col = "black", lwd = 3, lty = 2)
lines(u3, sin.f(u3),  col = "grey65", lwd = 3, lty = 2)
lines(u4, sin.f(u4),  col = "grey65", lwd = 3, lty = 2)
abline(v = c(0.25, 0.5, 0.75), h = c(-1/sqrt(2), 0, 1/sqrt(2)), lty = 2, col = "grey65")
legend('bottomleft', legend = c(expression(paste("cos(2",pi,"u)")), expression(paste("sin(2",pi,"u)"))), 
       bg = 'white', ncol = 1, col = c('black','black'), lty = c(1, 2), cex = .7, lwd = c(3, 3))

#dev.off()

# ====================================================================================
# Asymptotics of the test under local alternatives
# 5.1. First example: the gamma distribution
# ====================================================================================

## ---- Different functions and constants

dens <- function(x) dgamma(x, shape = lambda, scale = beta)
s1 <- function(x, lambda, beta) {- digamma(lambda) - log(beta) + log(x)}
s2 <- function(x, lambda, beta) {- lambda / beta + x / beta ^ 2}
R <- function(lambda, beta) matrix(c(trigamma(lambda), 1 / beta, 1 / beta, lambda / beta ^ 2), ncol = 2)
tau1 <- function(x, lambda, beta) {v <-  (x / beta); cos(2 * pi * pgamma(v, shape = lambda))}
tau2 <- function(x, lambda, beta) {v <-  (x / beta); sin(2 * pi * pgamma(v, shape = lambda))}
G <- function(lambda, beta) matrix(c(h10(lambda), h11(lambda), lambda / beta * h6(lambda, lambda + 1, 1), lambda / beta * h7(lambda, lambda + 1, 1)), ncol = 2)
sk <- function(x, lambda, beta) {v <- x / beta; 1 + lambda * log(v) - v * log(v)}
Gk <- function(lambda) matrix(c(- h8(lambda), - h9(lambda)), nrow = 2) 
S <- function(lambda, beta) matrix(c(-digamma(lambda), - 1 / beta * (lambda * digamma(lambda) + 1)), ncol = 2) 
C1 <- function(lambda) digamma(lambda) - trigamma(lambda) * (lambda * digamma(lambda) + 1)
M <- function(lambda) cbind(c(
    - h8(lambda) - (h10(lambda)  + lambda * C1(lambda) * h6(lambda, lambda + 1, 1)) / (lambda * trigamma(lambda) - 1), 
    - h9(lambda) - (h11(lambda)  + lambda * C1(lambda) * h7(lambda, lambda + 1, 1)) / (lambda * trigamma(lambda) - 1)
   ))
Sigma.inv <- function(lambda, beta) solve(diag(2)/2 - G(lambda, beta) %*% solve(R(lambda, beta)) %*% t(G(lambda, beta)))
V <- function(lambda, beta) t(M(lambda)) %*%  Sigma.inv(lambda, beta) %*% M(lambda)
level <- 0.05 # the nominal significance level is set to 0.05 
# We calculate the critical value for a test with level of 0.05, defined as the 5% quantile (in the right tail) of a chi-square distribution with 2 degrees of freedom
critical.value <- qchisq(level, df = 2, lower = FALSE) 
delta <- matrix(seq(from = 0, to = 30, length.out = 1000), byrow = T, nrow = 1)
lambda <- 1; beta <- 1 # We fix values for validation

## ---- Validation of sk
eps <- 1e-07
x <- rgamma(1, shape = lambda, scale = beta)
(dGG(x, lambda, beta, rho = 1 + eps, log = TRUE) - dGG(x, lambda, beta, rho = 1, log = TRUE)) / eps
sk(x, lambda, beta)

## ---- Validation of Gk
matrix(c(
  adaptIntegrate(f = function(x) tau1(x, lambda, beta) * sk(x, lambda, beta)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral, 
  adaptIntegrate(f = function(x) tau2(x, lambda, beta) * sk(x, lambda, beta)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral),
  nrow = 2
)    
Gk(lambda)

## ---- Validation of S
matrix(c(adaptIntegrate(f = function(x) sk(x, lambda, beta)  * s1(x, lambda, beta) * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral,
         adaptIntegrate(f = function(x) sk(x, lambda, beta)  * s2(x, lambda, beta) * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral), 
       ncol = 2)     
S(lambda, beta)

## ---- Validation of M
Gk(lambda) - G(lambda, beta) %*% solve(R(lambda, beta)) %*% t(S(lambda, beta))
M(lambda)

## ---- Validation of Sigma.inv
gamma.test.ML(rgamma(1e6, lambda, scale = beta))$Sigma.inv
Sigma.inv(lambda, beta)

## ---- power curve for lambda = 0.5
lambda <- .5; beta <- 1
ncp <- diag(t(delta) %*% V(lambda, beta) %*% delta) # non-centrality parameter (ncp) of the noncentral chi-square distribution 
power <- pchisq(critical.value, df = 2, ncp = ncp, lower.tail = FALSE)
power.lambda0.5 <- power

## ---- power curve for lambda = 1
lambda <- 1; beta <- 1
ncp <- diag(t(delta) %*% V(lambda, beta) %*% delta) # non-centrality parameter (ncp) of the noncentral chi-square distribution 
power <- pchisq(critical.value, df = 2, ncp = ncp, lower.tail = FALSE)
power.lambda.1 <- power

## ---- power curve for lambda = 1.5
lambda <- 1.5; beta <- 1
ncp <- diag(t(delta) %*% V(lambda, beta) %*% delta) # non-centrality parameter (ncp) of the noncentral chi-square distribution 
power <- pchisq(critical.value, df = 2, ncp = ncp, lower.tail = FALSE)
power.lambda1.5 <- power

#pdf(file = "Asymptotic_Power_Curves_gamma.pdf", width = 8, height = 8)

par(mfrow = c(1, 1), mar = c(4.5, 5, 1, 1) + 0.1)
color <- c('black', 'grey65', 'black')
lty <- c(1, 1, 2) 
lwd <- c(3, 3, 3)

plot(delta, power.lambda0.5, ylim = c(0, 1), type = 'l',
     col = color[1], lwd = lwd[1], lty = lty[1],
     xlab = expression(delta), ylab = 'Asymptotic Power', cex.lab = 2, cex.axis = 2)

lines(delta, power.lambda.1, ylim = c(0, 1), type = 'l',
      col = color[2], lwd = lwd[2], lty = lty[2], cex.lab = 2, cex.axis = 2)

lines(delta, power.lambda1.5, ylim = c(0, 1), type = 'l',
      col = color[3], lwd = lwd[3], lty = lty[3], cex.lab = 2, cex.axis = 2)

ytick <- seq(0.05, 1, by = 0.05)
yticklab <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
axis(side = 2, at = ytick, labels = yticklab)

legend("bottomright", legend = c(expression(paste(lambda[0], " = 0.5 ")),
                                 expression(paste(lambda[0], " = 1 ")),expression(paste(lambda[0], " = 1.5 "))), 
       bg = "white", ncol = 1, col = color, 
       lty = lty, cex = 2, lwd = lwd)

#dev.off()

## ---- Validation of the curves for each values of lambda = 0.5, 1, 1.5

n <- 100000
delta.test <- 4.0
delta.test <- 10.0
delta.test <- 15.0
delta.test <- 20.0
delta.test <- 30.0
res <- NULL
for (i in 1:1000){
  print(i)
  x <- rGG(n, lambda = lambda, beta = 1, rho = 1 + delta.test / sqrt(n))
  res <- c(res, gamma.test.ML(x)$p.value)
}
abline(v = delta.test)
power.test <- mean(res < .05)
power.test
abline(h =power.test, col = "red")

# ====================================================================================
# Asymptotics of the test under local alternatives
# 5.2. Second example: the Weibull distribution
# ====================================================================================

## ---- Different functions and constants

dens <- function(x) dWeibull(x, beta, rho)
euler <- - digamma(1) 
s1 <- function(x, beta, rho) {y <- rho * log(x / beta); rho / beta * (exp(y) - 1)}
s2 <- function(x, beta, rho) {y <- rho * log(x / beta); -1 / rho * (y * exp(y) - y - 1)}
R <- function(beta, rho) matrix(c(rho ^ 2 / beta ^ 2, 1 / beta * (euler - 1),  1 / beta * (euler - 1), 1 / rho ^ 2 * ((euler - 1) ^ 2 + pi ^ 2 / 6 )), nrow = 2)    
tau1 <- function(x, beta, rho) {v <-  (x / beta) ^ rho; cos(2 * pi * (1 - exp(- v)))}
tau2 <- function(x, beta, rho) {v <-  (x / beta) ^ rho; sin(2 * pi * (1 - exp(- v)))}
G <- function(beta, rho) matrix(c(rho / beta * h6(1,2,1), rho / beta * h7(1,2,1), - 1 / rho * h8(1), - 1 / rho * h9(1)), nrow = 2) 
sk <- function(x, beta, rho) {y <- rho * log(x / beta); y + euler}
Gk <- matrix(c(h10(1), h11(1)), nrow = 2) 
S <- function(beta, rho) matrix(c(rho / beta, euler / rho), ncol = 2) 
M <- cbind(c(
  h10(1) - 6  / pi ^ 2 * ((1 - euler + pi ^ 2 / 6) * h6(1,2,1)  - h8(1)),
  h11(1) - 6  / pi ^ 2 * ((1 - euler + pi ^ 2 / 6) * h7(1,2,1)  - h9(1))
))
Sigma.inv <- function(beta, rho) solve(diag(2)/2 - G(beta, rho) %*% solve(R(beta, rho)) %*% t(G(beta, rho)))
V <- function(beta, rho) t(M) %*%  Sigma.inv(beta, rho) %*% M
level <- 0.05 # the nominal significance level is set to 0.05 
# We calculate the critical value for a test with level of 0.05, defined as the 5% quantile (in the right tail) of a chi-square distribution with 2 degrees of freedom
critical.value <- qchisq(level, df = 2, lower = FALSE) 
delta <- matrix(seq(from = 0, to = 40, length.out = 1000), byrow = T, nrow = 1)
beta <- 1; rho <- 1 # We fix arbitrary values with no impact 

## ---- Validation of sk
eps <- 1e-07
x <- rWeibull(1, beta, rho)
(dGG(x, 1 + eps, beta, rho, log = TRUE) - dGG(x, 1, beta, rho, log = TRUE)) / eps
sk(x, beta, rho)

## ---- Validation of Gk
matrix(c(
  adaptIntegrate(f = function(x) tau1(x, beta, rho) * sk(x, beta, rho)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral, 
  adaptIntegrate(f = function(x) tau2(x, beta, rho) * sk(x, beta, rho)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral),
  nrow = 2)    
Gk

## ---- Validation of S
matrix(c(adaptIntegrate(f = function(x) sk(x, beta, rho)  * s1(x, beta, rho) * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral,
         adaptIntegrate(f = function(x) sk(x, beta, rho)  * s2(x, beta, rho) * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral), 
       ncol = 2)     
S(beta, rho)

## ---- Validation of M
Gk - G(beta, rho) %*% solve(R(beta, rho)) %*% t(S(beta, rho))
M

## ---- Validation of Sigma.inv
Weibull.test.ML(rexp(10))$Sigma.inv
Sigma.inv(beta, rho)

## ---- power curve 
ncp <- diag(t(delta) %*% V(beta, rho) %*% delta) # non-centrality parameter (ncp) of the noncentral chi-square distribution 
power <- pchisq(critical.value, df = 2, ncp = ncp, lower.tail = FALSE)

#pdf(file = "Asymptotic_Power_Curves_Weibull.pdf", width = 8, height = 8)

## ---- we define some parameters of the graph
par(mfrow = c(1, 1), mar = c(4.5, 5, 1, 1) + 0.1)
color <- c('black', 'grey65')
lty <- c(1, 1) 
lwd <- c(3, 3)

plot(delta, power, ylim = c(0, 1), type = 'l',
     col = color[1], lwd = lwd[1], lty = lty[1],
     xlab = expression(delta), ylab = 'Asymptotic Power', cex.lab = 2, cex.axis = 2)

ytick <- seq(0.05, 1, by = 0.05)
yticklab <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
axis(side = 2, at = ytick, labels = yticklab)

#dev.off()

## ---- Validation of the power curve

n <- 100000
delta.test <- 10.0
delta.test <- 15.0
delta.test <- 20.0
delta.test <- 30.0
res  <- NULL
for (i in 1:1000){
  print(i)
  x <- rGG(n, lambda = 1 + delta.test / sqrt(n), beta = 1, rho = 1)
  res <- c(res, Weibull.test.ML(x)$p.value)
}

abline(v = delta.test)
power.test <- mean(res < .05)
power.test
abline(h =power.test, col = "red")

# ====================================================================================
# Asymptotics of the test under local alternatives
# 5.3. Third example: the exponential power distribution (EPD)
# ====================================================================================

## ---- Different functions and constants

dens <- function(x) dEPD(x, lambda, mu, sigma)
C1 <- function(lambda) digamma(1/lambda + 1) + log(lambda)
C2 <- function(lambda) gamma(1 / lambda) / lambda ^ (2 / lambda) / gamma(3 / lambda)
C3 <- function(lambda) gamma(3 / lambda) ^ 2 / (gamma(1 / lambda) * gamma(5 / lambda) - gamma(3 / lambda) ^ 2)
s1 <- function(x, lambda, mu, sigma) {y <- (x - mu) / sigma; (1 / sigma) * (abs(y)) ^ (lambda - 1) * sign(y)}
s2 <- function(x, lambda, mu, sigma) {y <- (x - mu) / sigma; (1 / sigma) * ((abs(y)) ^ lambda - 1)}
R.ML <- function(lambda, sigma) 1 / sigma ^ 2 * matrix(c(lambda ^ (2 - 2 / lambda) * gamma(2 - 1 / lambda) / gamma(1 / lambda), 0, 0, lambda), nrow = 2) 
R.MM <- function(lambda, sigma) 1 / sigma ^ 2 * diag(c(C2(lambda), 4 * C3(lambda))) 
tau1 <- function(x, lambda, mu, sigma) cos(2 * pi * pEPD(x, lambda, mu, sigma))
tau2 <- function(x, lambda, mu, sigma) sin(2 * pi * pEPD(x, lambda, mu, sigma))
G <- function(lambda, sigma) 1 / sigma * matrix(c(0, h2(lambda) / lambda ^ (1 / lambda - 1) / gamma(1 / lambda), h1(lambda), 0), nrow = 2) 
J <- function(lambda, sigma) 1 / sigma * matrix(c(0, h5(lambda) * gamma(2 / lambda) / lambda ^ (1 / lambda) / gamma(3 / lambda), 2 * C3(lambda) * h4(lambda), 0), nrow = 2) 
r1.MM <- function(x, lambda, mu, sigma) {y <- (x - mu) / sigma; 1 / sigma * C2(lambda) * y}
r2.MM <- function(x, lambda, mu, sigma) {y <- (x - mu) / sigma; 1 / sigma * 2 * C3(lambda) * (C2(lambda) * y ^ 2 - 1)}
sk1 <- function(x, lambda, mu, sigma) {y <- (x - mu) / sigma; -2 * (abs(y)) ^ lambda * sign(y)}
sk2 <- function(x, lambda, mu, sigma) {y <- (x - mu) / sigma; -1 / lambda * ((abs(y)) ^ lambda * log(abs(y)) - C1(lambda) / lambda)}
Gk <- function(lambda) matrix(c(0, - 2  * h37(lambda), - 1 / lambda ^ 2 * h3(lambda), 0), nrow = 2) 
S.ML <- function(lambda, sigma) 1 / sigma * matrix(c(- 2 * lambda ^ (2 - 1 / lambda) / gamma(1 / lambda), 0, 0,
                             - 1 / lambda * (C1(lambda) + 1)), nrow = 2) 
S.MM <- function(lambda, sigma) 1 / sigma * matrix(c(- 4 * gamma(2 / lambda) / lambda ^ (1 / lambda) / gamma(3 / lambda), 0, 0,
                             - 2 * C3(lambda) / lambda ^ 2  * (2 * log(lambda) + 3 * digamma(3 / lambda) - digamma(1 / lambda))), nrow = 2) 
M.ML <- function(lambda) matrix(c(0, - 2  * h37(lambda)  + 2 * lambda * h2(lambda) / (gamma(1 / lambda) * gamma(2 - 1 / lambda)), 
                 - 1 / lambda ^ 2 * h3(lambda) +  1 / lambda ^ 2 * h1(lambda) * (C1(lambda) + 1), 0), nrow = 2) 
M.MM <- function(lambda) matrix(c(0, - 2 * h37(lambda) + 4 * lambda * h2(lambda) * gamma(2/lambda) / gamma(1/lambda)^2, 
                 - 1 / lambda ^ 2 * h3(lambda) +  1 / 2 / lambda ^ 2 * h1(lambda)  * (2 * log(lambda) + 3 * digamma(3 / lambda) - digamma(1 / lambda)), 0), nrow = 2) 
Sigma.inv.ML <- function(lambda, sigma) solve(diag(2)/2 - G(lambda, sigma) %*% solve(R.ML(lambda, sigma)) %*% t(G(lambda, sigma)))
Sigma.inv.MM <- function(lambda, sigma) solve(diag(2)/2 - G(lambda, sigma) %*% solve(R.MM(lambda, sigma))  %*% t(J(lambda, sigma)) -  J(lambda, sigma) %*% solve(R.MM(lambda, sigma))  %*% t(G(lambda, sigma)) + G(lambda, sigma) %*% solve(R.MM(lambda, sigma))  %*% t(G(lambda, sigma)))

V.ML <- function(lambda, sigma) t(M.ML(lambda)) %*%  Sigma.inv.ML(lambda, sigma) %*% M.ML(lambda)
V.MM <- function(lambda, sigma) t(M.MM(lambda)) %*%  Sigma.inv.MM(lambda, sigma) %*% M.MM(lambda)
level <- 0.05 # the nominal significance level is set to 0.05 
# We calculate the critical value for a test with level of 0.05, defined as the 5% quantile (in the right tail) of a chi-square distribution with 2 degrees of freedom
critical.value <- qchisq(level, df = 2, lower = FALSE) 
lambda <- 1.5
mu <- 0.2; sigma <- 1.3 

## ---- Validation of sk1
eps <- 1e-07
x <- rEPD(1, lambda, mu, sigma)
(log(dAPD(x, lambda, alpha = .5 + eps, rho = lambda, mu, sigma)) -
    log(dAPD(x, lambda, alpha = .5, rho = lambda, mu, sigma))) / eps
sk1(x, lambda, mu, sigma)

## ---- Validation of sk2
(log(dAPD(x, lambda, alpha = .5, rho = lambda + eps, mu, sigma)) -
    log(dAPD(x, lambda, alpha = .5, rho = lambda, mu, sigma))) / eps
sk2(x, lambda, mu, sigma)

## ---- Validation of Gk
matrix(c(
  round(adaptIntegrate(f = function(x) tau1(x, lambda, mu, sigma) * sk1(x, lambda, mu, sigma)  * dens(x), lower = -Inf, upper = 0, tol = 1e-9)$integral +
          adaptIntegrate(f = function(x) tau1(x, lambda, mu, sigma) * sk1(x, lambda, mu, sigma)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral, 6),
  adaptIntegrate(f = function(x) tau2(x, lambda, mu, sigma) * sk1(x, lambda, mu, sigma) * dens(x), lower = -Inf, upper = Inf, tol = 1e-9)$integral,
  adaptIntegrate(f = function(x) tau1(x, lambda, mu, sigma) * sk2(x, lambda, mu, sigma) * dens(x), lower = -Inf, upper = Inf, tol = 1e-9)$integral,
  round(adaptIntegrate(f = function(x) tau2(x, lambda, mu, sigma) * sk2(x, lambda, mu, sigma)  * dens(x), lower = -Inf, upper = 0, tol = 1e-9)$integral +
          adaptIntegrate(f = function(x) tau2(x, lambda, mu, sigma) * sk2(x, lambda, mu, sigma)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral, 6)), 
  nrow = 2)    
Gk(lambda)

## ---- Validation of S.ML
matrix(c(adaptIntegrate(f = function(x) sk1(x, lambda, mu, sigma) * s1(x, lambda, mu, sigma) * dens(x), lower = -Inf, upper = Inf, tol = 1e-9)$integral,
         round(adaptIntegrate(f = function(x) sk1(x, lambda, mu, sigma) * s2(x, lambda, mu, sigma)  * dens(x), lower = -Inf, upper = 0, tol = 1e-9)$integral +
                 adaptIntegrate(f = function(x) sk1(x, lambda, mu, sigma) * s2(x, lambda, mu, sigma)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral, 6),
         round(adaptIntegrate(f = function(x) sk2(x, lambda, mu, sigma) * s1(x, lambda, mu, sigma)  * dens(x), lower = -Inf, upper = 0, tol = 1e-9)$integral +
                 adaptIntegrate(f = function(x) sk2(x, lambda, mu, sigma) * s1(x, lambda, mu, sigma)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral, 6),
         adaptIntegrate(f = function(x) sk2(x, lambda, mu, sigma) * s2(x, lambda, mu, sigma) * dens(x), lower = -Inf, upper = Inf, tol = 1e-9)$integral), 
       nrow = 2)     
S.ML(lambda, sigma)

## ---- Validation of S.MM
matrix(c(adaptIntegrate(f = function(x) sk1(x, lambda, mu, sigma) * r1.MM(x, lambda, mu, sigma) * dens(x), lower = -Inf, upper = Inf, tol = 1e-9)$integral,
         round(adaptIntegrate(f = function(x) sk1(x, lambda, mu, sigma) * r2.MM(x, lambda, mu, sigma)  * dens(x), lower = -Inf, upper = 0, tol = 1e-9)$integral +
                 adaptIntegrate(f = function(x) sk1(x, lambda, mu, sigma) * r2.MM(x, lambda, mu, sigma)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral, 6),
         round(adaptIntegrate(f = function(x) sk2(x, lambda, mu, sigma) * r1.MM(x, lambda, mu, sigma)  * dens(x), lower = -Inf, upper = 0, tol = 1e-9)$integral +
                 adaptIntegrate(f = function(x) sk2(x, lambda, mu, sigma) * r1.MM(x, lambda, mu, sigma)  * dens(x), lower = 0, upper = Inf, tol = 1e-9)$integral, 6),
         adaptIntegrate(f = function(x) sk2(x, lambda, mu, sigma) * r2.MM(x, lambda, mu, sigma) * dens(x), lower = -Inf, upper = Inf, tol = 1e-9)$integral), 
       nrow = 2)     
S.MM(lambda, sigma)

## ---- Validation of M.ML
Gk(lambda) - G(lambda, sigma) %*% solve(R.ML(lambda, sigma)) %*% t(S.ML(lambda, sigma))
M.ML(lambda)

## ---- Validation of M.MM
Gk(lambda) - G(lambda, sigma) %*% solve(R.MM(lambda, sigma)) %*% t(S.MM(lambda, sigma))
M.MM(lambda)

## ---- Validation of Sigma.inv.ML
EPD.test.ML(rnorm(10), lambda)$Sigma.inv
Sigma.inv.ML(lambda, sigma)

## ---- Validation of Sigma.inv.MM
EPD.test.MM(rnorm(10), lambda)$Sigma.inv
Sigma.inv.MM(lambda, sigma)

## ---- 1st scenario for local alternative, where delta1 varies from 0 to 3.5, delta2 is fixed at 0, with delta = (delta1, delta2)

delta <- matrix(c(seq(from = 0, to = 3.5, length.out = 1000), rep(0, 1000)), byrow = T, nrow = 2)
#delta <- matrix(c(seq(from = 0, to = -3.5, length.out = 1000), rep(0, 1000)), byrow = T, nrow = 2) # same results

## ---- power curves 
ncp.ML <- diag(t(delta) %*% V.ML(lambda, sigma) %*% delta) # non-centrality parameter (ncp) of the noncentral chi-square distribution 
ncp.MM <- diag(t(delta) %*% V.MM(lambda, sigma) %*% delta) # non-centrality parameter (ncp) of the noncentral chi-square distribution 
power.ML <- pchisq(critical.value, df = 2, ncp = ncp.ML, lower.tail = FALSE)
power.MM <- pchisq(critical.value, df = 2, ncp = ncp.MM, lower.tail = FALSE)

#pdf(file = "Asymptotic_Power_Curves_Delta1.pdf", width = 8, height = 8)

## ---- we define some parameters of the graph
par(mfrow = c(1, 1), mar = c(4.5, 5, 1, 1) + 0.1)
color <- c('black', 'grey65')
lty <- c(2, 1) 
lwd <- c(3, 3)

plot(delta[1,], power.ML, ylim = c(0, 1), type = 'l',
     col = color[2], lwd = lwd[2], lty = lty[2],
     xlab = expression(delta[1]), ylab = 'Asymptotic Power', cex.lab = 2, cex.axis = 2)

ytick <- seq(0.05, 1, by = 0.05)
yticklab <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
axis(side = 2, at = ytick, labels = yticklab)

lines(delta[1,], power.MM,  col = color[1], lwd = lwd[1], lty = lty[1])
legend("bottomright", legend = c("MM","ML"), 
       bg = "white", ncol = 1, col = color, 
       lty = lty, cex = 2, lwd = lwd)

#dev.off()

## ---- Validation of the power curve

n <- 1000
delta1.test <- 1.0
delta1.test <- 1.5
delta1.test <- 2.0
res.ML <- res.MM <- NULL
for (i in 1:1000){
  print(i)
  x <- rAPD(n, lambda = lambda, alpha = 0.5 + delta1.test / sqrt(n), rho = lambda, mu = 0, sigma = 1)
  res.ML <- c(res.ML, EPD.test.ML(x, lambda)$p.value)
  res.MM <- c(res.MM, EPD.test.MM(x, lambda)$p.value)
}

abline(v = delta1.test)
power.ML.test <- mean(res.ML < .05)
power.ML.test
abline(h =power.ML.test, col = color[2], lwd = lwd[2], lty = lty[2])
power.MM.test <- mean(res.MM < .05)
power.MM.test
abline(h =power.MM.test, col = color[1], lwd = lwd[1], lty = lty[1])

## ---- 2nd scenario for local alternative, where delta2 varies from 0 to 18, delta1 is fixed at 0, with delta = (delta1, delta2)

delta <- matrix(c(rep(0, 1000), seq(from = 0, to = 18, length.out = 1000)), byrow = T, nrow = 2)
#delta <- matrix(c(rep(0, 1000), seq(from = 0, to = -18, length.out = 1000)), byrow = T, nrow = 2) # same results

## ---- power curves 
ncp.ML <- diag(t(delta) %*% V.ML(lambda, sigma) %*% delta) # non-centrality parameter (ncp) of the noncentral chi-square distribution 
ncp.MM <- diag(t(delta) %*% V.MM(lambda, sigma) %*% delta) # non-centrality parameter (ncp) of the noncentral chi-square distribution 
power.ML <- pchisq(critical.value, df = 2, ncp = ncp.ML, lower.tail = FALSE)
power.MM <- pchisq(critical.value, df = 2, ncp = ncp.MM, lower.tail = FALSE)

#pdf(file = "Asymptotic_Power_Curves_Delta2.pdf", width = 8, height = 8)

## ---- we define some parameters of the graph
par(mfrow = c(1, 1), mar = c(4.5, 5, 1, 1) + 0.1)
color <- c('black', 'grey65')
lty <- c(2, 1) 
lwd <- c(3, 3)

plot(delta[2,], power.ML, ylim = c(0, 1), type = 'l',
     col = color[2], lwd = lwd[2], lty = lty[2],
     xlab = expression(delta[2]), ylab = 'Asymptotic Power', cex.lab = 2, cex.axis = 2)

ytick <- seq(0.05, 1, by = 0.05)
yticklab <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
axis(side = 2, at = ytick, labels = yticklab)

lines(delta[2,], power.MM,  col = color[1], lwd = lwd[1], lty = lty[1])
legend("bottomright", legend = c("MM","ML"), 
       bg = "white", ncol = 1, col = color, 
       lty = lty, cex = 2, lwd = lwd)

#dev.off()

## ---- Validation

n <- 10000
delta2.test <- 5.0
delta2.test <- 7.5
delta2.test <- 10.0
res.ML <- res.MM <- NULL
for (i in 1:1000){
  print(i)
  x <- rAPD(n, lambda = lambda, alpha = 0.5, rho = lambda + delta2.test / sqrt(n), mu = 0, sigma = 1)
  res.ML <- c(res.ML, EPD.test.ML(x, lambda)$p.value)
  res.MM <- c(res.MM, EPD.test.MM(x, lambda)$p.value)
}
abline(v =delta2.test)
power.ML.test <- mean(res.ML < .05)
power.ML.test
abline(h =power.ML.test, col = color[2], lwd = lwd[2], lty = lty[2])
power.MM.test <- mean(res.MM < .05)
power.MM.test
abline(h =power.MM.test, col = color[1], lwd = lwd[1], lty = lty[1])

# ====================================================================================
# Real data example
# ====================================================================================

library(lawstat)
data(bias)
dput(bias)
x.data <- c(
  -3.282,  8.376, -1.546, -0.866,  9.627, -2.311,  1.653,  8.893,  1.015,  2.041, -0.138, -1.593,  0.667, -0.188, 
  0.673, -0.772,  2.435, -1.622, -2.378, -7.266, -1.388,  1.675, -1.858, -0.785, -2.958, -3.804,  3.646, -0.817,
  2.813,  0.071, -0.451,  1.950,  0.663,  5.038,  0.965,  5.398,  1.631,  1.844,  5.154, -2.988,  2.858,  1.215,
  2.495, -4.744,  0.205, -0.875,  0.661, -0.587, -0.409, -2.210, -0.532,  5.888,  4.147,  2.767,  1.382,  2.871,
  0.117,  2.177, -5.753, -1.408, -2.821, -0.238,  2.767, -1.816, -2.062,  7.445,  1.060, -0.570, -1.517,  0.092,
  -1.250, -5.891, -2.888,  0.592,  2.006,  2.276, -3.294, -0.521,  0.846, -8.551,  5.158,  4.849,  0.123,  1.078,
  -1.622, -2.683, -2.583, -3.766, -2.550, -1.647, -1.536, -4.448,  0.865, -1.685,  1.177, -0.729)
n <- length(x.data); n

graph_hist <- function(x, dens, est){
  dens.eval <- function(x) do.call(dens, c(list(x), est))
  couleurs <- c(rep("grey50", 7), rep("grey90", 4), rep("grey70", 8)) # define the color of each of the length(h$counts)=19 rectangles of the histogram
  hist(x.data, xlim = c(-11, 11), ylim = c(0, 0.20), prob = TRUE, breaks = 25, col = couleurs,
       border = "black", main = "Histogram of Temperature Forecast Errors", xlab = "Forecast Error")
  curve(dens.eval, add = TRUE, col = "black", n = 1000, lty = 1, lwd = 2)
}

graph.ellipse <- function(Sigma.inv, x, y, title = "95% Confidence Ellipse"){
  par(cex.axis = 1.6, cex.lab = 1.6, cex.main = 1.6)
  Sigma <- solve(Sigma.inv)
  ell <- ellipse(Sigma, level = 0.95)
  xlimm <- c(min(c(ell[,1], x)), max(c(ell[,1], x)))
  ylimm <- c(min(c(ell[,2], y)), max(c(ell[,2], y)))
  par(mgp = c(2.7, .7, 0))  # Valeurs par défaut sont c(3, 1, 0)
  plot(ell, xlim = xlimm, ylim = ylimm, main = title,
       col = 'black', lwd = 2, type = 'l',
       xlab = expression(sqrt(n) * " " * C[n](hat(bold(theta))[n])),
       ylab = "")
  mtext(expression(sqrt(n) * S[n](hat(bold(theta))[n])), side = 2, line = 1.9, cex = 1.6)
  points(x, y, pch = 16)
  vv <- qnorm(.975) * sqrt(Sigma[1,1]); abline(v = c(-vv, 0, vv), lty = 2)
  vv <- qnorm(.975) * sqrt(Sigma[2,2]); abline(h = c(-vv, 0, vv), lty = 2)
}  

## ---- Histogram

#pdf("plot_hist_example1.pdf")

par(mgp = c(1.8, .7, 0))  # Valeurs par défaut sont c(3, 1, 0)
h <- hist(x.data, breaks = 25, plot = FALSE)
couleurs <- c(rep("grey50", 7), rep("grey90", 4), rep("grey70", 8)) # define the color of each of the length(h$counts)=19 rectangles of the histogram
hist(x.data, xlim = c(-11, 11), ylim = c(0, 0.20), prob = TRUE, breaks = 25, col = couleurs,
     border = "black", main = "Histogram of Temperature Forecast Errors", xlab = "Forecast Error")

dens.eval <- function(x) do.call(dens, c(list(x), ana[[1]]))

ana <- Laplace.test.ML(x = x.data); dens <- dLaplace
curve(dens.eval, add = TRUE, col = "grey60", n = 1000, lty = 4, lwd = 2)

ana <- EPD.test.ML(x = x.data); dens <- dEPD
curve(dens.eval, add = TRUE, col = "grey60", n = 1000, lty = 1, lwd = 2)

ana <- Student.test.ML(x = x.data); dens <- dStudent
curve(dens.eval, add = TRUE, col = "black", n = 1000, lty = 1, lwd = 2)

ana <- logistic.test.ML(x = x.data); dens <- dlogistic
curve(dens.eval, add = TRUE, col = "black", n = 1000, lty = 2, lwd = 2)

ana <- normal.test.ML(x = x.data); dens <- dnorm
curve(dens.eval, add = TRUE, col = "black", n = 1000, lty = 3, lwd = 2)

ana <- exp.Weibull.test.ML(x = x.data); dens <- dexp.Weibull
curve(dens.eval, add = TRUE, col = "grey60", n = 1000, lty = 2, lwd = 2)

ana <- Gumbel.test.ML(x = x.data); dens <- dGumbel
curve(dens.eval, add = TRUE, col = "black", n = 1000, lty = 4, lwd = 2)

legend('topright', legend = c("Laplace","EPD","Student","logistic","normal", "exp-Weibull", "Gumbel"), 
       bg = 'white', ncol = 1, col = c("grey60","grey60", "black","black",  "black","grey60", "black"), 
       lty = c(4, 1, 1, 2, 3, 2, 4), cex = .9, lwd = rep(3, 7), seg.len = 4)

#dev.off()

## ---- End of histogram

## ---- Figure ellipse for the EPD model
#pdf("plot_ellipse_EPD.pdf")
ana <- EPD.test.ML(x.data); ana
graph.ellipse(ana$Sigma.inv, ana$tau.bar.sqrtn[1], ana$tau.bar.sqrtn[2], title = "95% Confidence Ellipse for the EPD Model")
#dev.off()

## ---- Figure ellipse for the Laplace model
ana <- Laplace.test.ML(x.data); ana
graph.ellipse(ana$Sigma.inv, ana$tau.bar.sqrtn[1], ana$tau.bar.sqrtn[2], title = "95% Confidence Ellipse for the Laplace Model")

## ---- Figure ellipse for the normal model
#pdf("plot_ellipse_norm.pdf")
ana <- normal.test.ML(x.data); ana
graph.ellipse(ana$Sigma.inv, ana$tau.bar.sqrtn[1], ana$tau.bar.sqrtn[2], title = "95% Confidence Ellipse for the Normal Model")
#dev.off()

## ---- Figure ellipse for the exp-Weibull model
ana <- exp.Weibull.test.ML(x.data); ana
graph.ellipse(ana$Sigma.inv, ana$tau.bar.sqrtn[1], ana$tau.bar.sqrtn[2], title = "95% Confidence Ellipse for the exp-Weibull Model")

## ---- Figure ellipse for the gumbel model
#pdf("plot_ellipse_Gumbel.pdf")
ana <- Gumbel.test.ML(x.data); ana
graph.ellipse(ana$Sigma.inv, ana$tau.bar.sqrtn[1], ana$tau.bar.sqrtn[2], title = "95% Confidence Ellipse for the Gumbel Model")
#dev.off()

## ---- Figure ellipse for the logistic model
ana <- logistic.test.ML(x.data); ana
graph.ellipse(ana$Sigma.inv, ana$tau.bar.sqrtn[1], ana$tau.bar.sqrtn[2], title = "95% Confidence Ellipse for the logistic Model")

## ---- Figure ellipse for the Student model
ana <- Student.test.ML(x.data); ana
graph.ellipse(ana$Sigma.inv, ana$tau.bar.sqrtn[1], ana$tau.bar.sqrtn[2], title = "95% Confidence Ellipse for the Student Model")

## ---- Table of results

ana1 <- EPD.test.ML(x.data)
ana2 <- Laplace.test.ML(x.data)
ana3 <- normal.test.ML(x.data)
ana4 <- exp.Weibull.test.ML(x.data)
ana5 <- Gumbel.test.ML(x.data)
ana6 <- logistic.test.ML(x.data)
ana7 <- Student.test.ML(x.data)

ML <- list(ML1 = ana1[[1]], ML2 = ana2[[1]], ML3 = ana3[[1]], ML4 = ana4[[1]], ML5 = ana5[[1]], ML6 = ana6[[1]], ML7 = ana7[[1]]); ML
lapply(ML, round, digits = 3)

cbind(
  round(c(ana1$neg2L, ana2$neg2L, ana3$neg2L, ana4$neg2L, ana5$neg2L, ana6$neg2L, ana7$neg2L), 1),
  t(round(cbind(ana1$Z.tau.bar, ana2$Z.tau.bar, ana3$Z.tau.bar, ana4$Z.tau.bar, ana5$Z.tau.bar, ana6$Z.tau.bar, ana7$Z.tau.bar), 2)),
  round(c(ana1$Tn, ana2$Tn, ana3$Tn, ana4$Tn, ana5$Tn, ana6$Tn, ana7$Tn), 2),
  round(c(ana1$p.value, ana2$p.value, ana3$p.value, ana4$p.value, ana5$p.value, ana6$p.value, ana7$p.value), 4)
)

## ---- exp.gamma test (not presented in the paper since similar to the normal model as lamda.hat goes to infinity)

## ---- Figure ellipse for the exp-gamma model
ana <- exp.gamma.test.ML(x.data); ana # computation error since lambda.hat tend to infinity 
ana <- exp.gamma.test.ML(x.data, lambda = 59); ana # we fix lambda as large as possible
graph.ellipse(ana$Sigma.inv, ana$tau.bar.sqrtn[1], ana$tau.bar.sqrtn[2], title = "95% Confidence Ellipse for the exp-Weibull Model")

# similar to the normal ellipse
hist(x.data, xlim = c(-11, 11), ylim = c(0, 0.20), prob = TRUE, breaks = 25, col = couleurs,
     border = "black", main = "Histogram of Temperature Forecast Errors", xlab = "Forecast Error")

## ---- normal curve
ana <- normal.test.ML(x = x.data); dens <- dnorm
curve(dens.eval, add = TRUE, col = "red", n = 1000, lty = 1, lwd = 2)

# exp.gamma curve with lambda.hat <- 2e10, mu.hat and sigma.hat estimated in the limit
# if lambda.hat tend to infinity, mu.hat approx - sd(x.data) * sqrt(lambda.hat) * log(lambda.hat) 
#  and sigma.hat approx sd(x.data) * sqrt(lambda.hat), then the neg2L tend to that of the Normal
lambda.hat <- 2e10 
mu.hat <- - sd(x.data) * sqrt(lambda.hat) * log(lambda.hat)
sigma.hat <-  sd(x.data) * sqrt(lambda.hat)
ana[[1]] <- c(lambda.hat, mu.hat, sigma.hat); ana[[1]]
Fi <- pgamma(exp((x.data - mu.hat) / sigma.hat), shape = lambda.hat, scale = 1)
neg2L <- - 2 * sum(log(dexp.gamma(x.data, lambda = lambda.hat, mu = mu.hat, sigma = sigma.hat)))
neg2L
#496.5056152
dens <- dexp.gamma
curve(dens.eval, add = TRUE, col = "blue", n = 1000, lty = 1, lwd = 2)

## ---- exp.gamma curve with lambda = 59
ana <- exp.gamma.test.ML(x.data, lambda = 59); dens <- dexp.gamma
curve(dens.eval, add = TRUE, col = "black", n = 1000, lty = 1, lwd = 2)

## ---- Simulated p-values (vs chi-squared p-values)

test <- EPD.test.ML
ana <- test(x = x.data)
gen <- function(n) rEPD(n, lambda = ana[[1]][1])

test <- Laplace.test.ML
ana <- test(x = x.data)
gen <- rLaplace

test <- normal.test.ML
ana <- test(x = x.data)
gen <- rnorm

test <- exp.Weibull.test.ML
ana <- test(x = x.data)
gen <- rexp.Weibull

test <- Gumbel.test.ML
ana <- test(x = x.data)
gen <- rGumbel

test <- logistic.test.ML
ana <- test(x = x.data)
gen <- rlogistic

test <- Student.test.ML
ana <- test(x = x.data)
gen <- function(n) rt(n, df = ana[[1]][1])

Tn.xdata <- ana$Tn
pvalue.chi.squared <- ana$p.value
n <- length(x.data); n
Tn.vec <- NULL

## ---- simulation error using 100,000 simulations
nb.sim <- 100000
1.96 * sqrt(pvalue.chi.squared * (1 - pvalue.chi.squared) / 100000) # error using 100,000 simulations
for (i in 1:nb.sim){
  print(i)
  x <- gen(n) 
  res <- test(x)
  Tn.vec <- c(Tn.vec, res$Tn)
}
length(Tn.vec)
hist(Tn.vec)
p.value.sim <- mean(Tn.vec > Tn.xdata)
p.value.sim
pvalue.chi.squared
error.sim <- 1.96 * sqrt(pvalue.chi.squared * (1 - pvalue.chi.squared) / length(Tn.vec))
res <- round(c(length(Tn.vec), error.sim, pvalue.chi.squared, p.value.sim, p.value.sim - error.sim, p.value.sim + error.sim), 4)
names(res) = c("nb.sim", "error", "ch-squared p-value", "simulated p-value", "CI -Left", "CI-Right")
cbind(res)

# EPD 0.3607385722
# 100754.0000  0.0030    0.3853    0.3607    0.3577    0.3637

# Laplace 0.3946612949
# 102759.0000      0.0030      0.4038      0.3947      0.3917      0.3977

# normal 0.02582255503
# 110175.0000      0.0010      0.0271      0.0258      0.0249      0.0268

# exp.Weibull 0
# 94336     0     0     0     0     0

# Gumbel 0.0003826117668
# 101931.0000      0.0001      0.0005      0.0004      0.0002      0.0005

# Logistic 0.3660168462
# 100438.000      0.003      0.362      0.366      0.363      0.369

# Student 0.49982
# 100000         0.0031    0.5090      0.4998     0.4967     0.5029

## ---- End of Real Data Example

# ====================================================================================
# Validation of the distributions
# ====================================================================================

## ---- Validation functions

validation1 <- function(){
  n <- 1e6
  ui <- seq(0 + 1 / n / 2, 1 -  1 / n / 2, length = n)
  xi <- qf(ui) # generates perfectly G distributed sample
  hist(xi, prob = TRUE, breaks = 50)
  curve(dens(x), from = max(left.bound, -1000), add = TRUE, col = "blue", n = 1000, lw = 2)
}
validation2 <- function(est, test){  
  x <- gen(10000)
  val <- matrix(c(
    adaptIntegrate(f = dens, lower = left.bound, upper = right.bound, tol = 1e-9)$integral, 1, 
    qf(cdf(quant)), quant, # validation of qf and cdf
    round(cdf(quant), 6), round(adaptIntegrate(f = dens, lower = left.bound, upper = quant, tol = 1e-9)$integral, 6),
    round(cdf(quant), 6), round(mean(x < quant), 6),
    rbind(round(est(x), 6), theta),
    rbind(round(test(x)[[1]], 6), theta)
  ), byrow = T, ncol = 2) 
  val <- cbind(val, round(abs(val[,1] - val[,2]), 6))
  rownames(val) <- c("density integrates to 1", "CDF and inv-CDF", "CDF", "gen", rep("estimation", nrow(val) - 4))
  colnames(val) <- c("Component 1", "Component2", "Absolute difference")
  return(val)
}
validation3 <- function(test, nb.sim = 5000, n = 1000){
  Tn <- tau.bar.sqrtn <- NULL
  for (i in 1:nb.sim){
    print(i)
    x <- gen(n) 
    res <- test(x)
    Tn <- c(Tn, res$Tn)
    tau.bar.sqrtn <- cbind(tau.bar.sqrtn, res$tau.bar.sqrtn)
  }
  val <- matrix(c(
    cbind(round(apply(tau.bar.sqrtn, 1, mean), 4), c(0, 0)),
    round(rbind(var(t(tau.bar.sqrtn)),  solve(test(x)$Sigma.inv)), 4),
    round(rbind(solve(var(t(tau.bar.sqrtn))),  test(x)$Sigma.inv), 4)
  ), byrow = T, ncol = 4)  
  rownames(val) <- c("mean : sqrt(n) Cn and sqrt(n) Sn", "Sigma[1,]", "Sigma[2,]", 
                     "Sigma.inv[1,]", "Sigma.inv[2,]")
  colnames(val) <- c("Simulated", "Simulated", "Exact", "Exact")
  print(val)
  return(list(Tn = Tn, tau.bar.sqrtn = tau.bar.sqrtn, Sigma.inv = res$Sigma.inv))
}  
plot.ellipse.tau.bar.sqrtn <- function(tau.bar.sqrtn1, tau.bar.sqrtn2, Sigma){
  plot(tau.bar.sqrtn1, tau.bar.sqrtn2)  
  lines(ellipse(Sigma), col = 'red', lwd = 2)
}  
hist.tau.bar.sqrtn <- function(tau.bar.sqrtn){
  hist(tau.bar.sqrtn, freq = F)
  curve(dnorm(x, mean = mean(tau.bar.sqrtn), sd = sd(tau.bar.sqrtn)), add = TRUE, col = "grey0", n = 1000, lw = 2)
  p.value <- normal.test.ML(tau.bar.sqrtn)$p.value
  names(p.value) = "normality test p.value"
  return(p.value)
}
hist.Tn <- function(Tn){
  hist(Tn, prob = TRUE)
  curve(dchisq(x, df = 2), add = TRUE, col = "grey0", n = 1000, lw = 2)
  p.value <- chisquared.test.ML(Tn, k = 2)$p.value
  names(p.value) = "ch-square test p.value"
  return(p.value)
}  
validation4 <- function(Tn){ #Validation of the chi-square distribution
  quant.nominal <- seq(1, 10, 1)
  g <- function(q) mean(Tn> qchisq(q, 2, lower.tail = FALSE))
  quant.sim <- apply(cbind(quant.nominal)/100, 1, g) * 100
  plot(quant.nominal, quant.sim, type = 'l', lwd = 2); abline(a=0, b=1, col = 'red', lwd = 2, lty = 2)
  ress <- round(cbind(quant.nominal, quant.sim, abs(quant.nominal - quant.sim)), 1)
  colnames(ress) <- c("nominal quantiles", "simulated quantiles", "Absolute difference")
  return(ress)
}

## ---- Distributions on R

## ---- EPD
lambda <- 1.5; mu <- -4.3; sigma <- 2.8
theta <- c(lambda, mu, sigma)
dens <- function(x) dEPD(x, lambda, mu, sigma)
cdf  <- function(x) pEPD(x, lambda, mu, sigma)
qf  <- function(p) qEPD(p, lambda, mu, sigma)
gen <- function(n) rEPD(n, lambda, mu, sigma) 
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.EPD(x)
test <- function(x)  EPD.test.ML(x)

est <- function(x)  ML.EPD(x, lambda = lambda)
test <- function(x)  EPD.test.ML(x, lambda = lambda)

est <- function(x)  ML.EPD(x, mu = mu)
test <- function(x)  EPD.test.ML(x, mu = mu)

est <- function(x)  ML.EPD(x, sigma = sigma)
test <- function(x)  EPD.test.ML(x, sigma = sigma)

est <- function(x)  ML.EPD(x, lambda = lambda, mu = mu)
test <- function(x)  EPD.test.ML(x, lambda = lambda, mu = mu)

est <- function(x)  ML.EPD(x, lambda = lambda, sigma = sigma)
test <- function(x)  EPD.test.ML(x, lambda = lambda, sigma = sigma)

est <- function(x)  ML.EPD(x, mu = mu, sigma = sigma)
test <- function(x)  EPD.test.ML(x, mu = mu, sigma = sigma)

est <- function(x)  ML.EPD(x, lambda = lambda, mu = mu, sigma = sigma)
test <- function(x)  EPD.test.ML(x, lambda = lambda, mu = mu, sigma = sigma)

est <- function(x)  MM.EPD(x, lambda)
test <- function(x)  EPD.test.MM(x, lambda)

est <- function(x)  MM.EPD(x, lambda, mu = mu)
test <- function(x)  EPD.test.MM(x, lambda, mu = mu)

est <- function(x)  MM.EPD(x, lambda, sigma = sigma)
test <- function(x)  EPD.test.MM(x, lambda, sigma = sigma)

est <- function(x)  MM.EPD(x, lambda, mu = mu, sigma = sigma)
test <- function(x)  EPD.test.MM(x, lambda, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Laplace
mu <- 2.3; sigma <- 1.8
theta <- c(mu, sigma)
dens <- function(x) dLaplace(x, mu, sigma)
cdf  <- function(x) pLaplace(x, mu, sigma)
qf  <- function(p) qLaplace(p, mu, sigma)
gen <- function(n) rLaplace(n, mu, sigma) 
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.Laplace(x)
test <- function(x)  Laplace.test.ML(x)

est <- function(x)  ML.Laplace(x, mu = mu)
test <- function(x)  Laplace.test.ML(x, mu = mu)

est <- function(x)  ML.Laplace(x, sigma = sigma)
test <- function(x)  Laplace.test.ML(x, sigma = sigma)

est <- function(x)  ML.Laplace(x, mu = mu, sigma = sigma)
test <- function(x)  Laplace.test.ML(x, mu = mu, sigma = sigma)

est <- function(x)  MM.Laplace(x)
test <- function(x)  Laplace.test.MM(x)

est <- function(x)  MM.Laplace(x, mu = mu)
test <- function(x)  Laplace.test.MM(x, mu = mu)

est <- function(x)  MM.Laplace(x, sigma = sigma)
test <- function(x)  Laplace.test.MM(x, sigma = sigma)

est <- function(x)  MM.Laplace(x, mu = mu, sigma = sigma)
test <- function(x)  Laplace.test.MM(x, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Normal
mu <- 2.3; sigma <- 1.8
theta <- c(mu, sigma)
dens <- function(x) dnormal(x, mu, sigma)
cdf  <- function(x) pnormal(x, mu, sigma)
qf  <- function(p) qnormal(p, mu, sigma)
gen <- function(n)  rnormal(n, mu, sigma) 
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.normal(x)
test <- function(x)  normal.test.ML(x)

est <- function(x)  ML.normal(x, mu = mu)
test <- function(x)  normal.test.ML(x, mu = mu)

est <- function(x)  ML.normal(x, sigma = sigma)
test <- function(x)  normal.test.ML(x, sigma = sigma)

est <- function(x)  ML.normal(x, mu = mu, sigma = sigma)
test <- function(x)  normal.test.ML(x, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- exp-gamma
lambda <- 1.5; mu <- 2.3; sigma <- 2.8; 
theta <- c(lambda, mu, sigma)
dens <- function(x) dexp.gamma(x, lambda, mu, sigma)
cdf  <- function(x) pexp.gamma(x, lambda, mu, sigma)
qf  <- function(p) qexp.gamma(p, lambda, mu, sigma)
gen <- function(n) rexp.gamma(n, lambda, mu, sigma) 
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.exp.gamma(x)
test <- function(x)  exp.gamma.test.ML(x)

est <- function(x)  ML.exp.gamma(x, lambda = lambda)
test <- function(x)  exp.gamma.test.ML(x, lambda = lambda)

est <- function(x)  ML.exp.gamma(x, mu = mu)
test <- function(x)  exp.gamma.test.ML(x, mu = mu)

est <- function(x)  ML.exp.gamma(x, sigma = sigma)
test <- function(x)  exp.gamma.test.ML(x, sigma = sigma)

est <- function(x)  ML.exp.gamma(x, lambda = lambda, mu = mu)
test <- function(x)  exp.gamma.test.ML(x, lambda = lambda, mu = mu)

est <- function(x)  ML.exp.gamma(x, lambda = lambda, sigma = sigma)
test <- function(x)  exp.gamma.test.ML(x, lambda = lambda, sigma = sigma)

est <- function(x)  ML.exp.gamma(x, mu = mu, sigma = sigma)
test <- function(x)  exp.gamma.test.ML(x, mu = mu, sigma = sigma)

est <- function(x)  ML.exp.gamma(x, lambda = lambda, mu = mu, sigma = sigma)
test <- function(x)  exp.gamma.test.ML(x, lambda = lambda, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- exp-Weibull
mu <- 2.3; sigma <- 1.8
theta <- c(mu, sigma)
dens <- function(x) dexp.Weibull(x, mu, sigma)
cdf  <- function(x) pexp.Weibull(x, mu, sigma)
qf  <- function(p) qexp.Weibull(p, mu, sigma)
gen <- function(n) rexp.Weibull(n, mu, sigma) 
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.exp.Weibull(x)
test <- function(x)  exp.Weibull.test.ML(x)

est <- function(x)  ML.exp.Weibull(x, mu = mu)
test <- function(x)  exp.Weibull.test.ML(x, mu = mu)

est <- function(x)  ML.exp.Weibull(x, sigma = sigma)
test <- function(x)  exp.Weibull.test.ML(x, sigma = sigma)

est <- function(x)  ML.exp.Weibull(x, mu = mu, sigma = sigma)
test <- function(x)  exp.Weibull.test.ML(x, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Gumbel
mu <- 2.3; sigma <- 1.8
theta <- c(mu, sigma)
dens <- function(x) dGumbel(x, mu, sigma)
cdf  <- function(x) pGumbel(x, mu, sigma)
qf  <- function(p) qGumbel(p, mu, sigma)
gen <- function(n) rGumbel(n, mu, sigma) 
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.Gumbel(x)
test <- function(x)  Gumbel.test.ML(x)

est <- function(x)  ML.Gumbel(x, mu = mu)
test <- function(x)  Gumbel.test.ML(x, mu = mu)

est <- function(x)  ML.Gumbel(x, sigma = sigma)
test <- function(x)  Gumbel.test.ML(x, sigma = sigma)

est <- function(x)  ML.Gumbel(x, mu = mu, sigma = sigma)
test <- function(x)  Gumbel.test.ML(x, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- logistic
mu <- 2.3; sigma <- 1.8
theta <- c(mu, sigma)
dens <- function(x) dlogistic(x, mu, sigma)
cdf  <- function(x) plogistic(x, mu, sigma)
qf  <- function(p) qlogistic(p, mu, sigma)
gen <- function(n) rlogistic(n, mu, sigma)
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.logistic(x)
test <- function(x) logistic.test.ML(x)

est <- function(x)  ML.logistic(x, mu = mu)
test <- function(x) logistic.test.ML(x, mu = mu)

est <- function(x)  ML.logistic(x, sigma = sigma)
test <- function(x) logistic.test.ML(x, sigma = sigma)

est <- function(x)  ML.logistic(x, mu = mu, sigma = sigma)
test <- function(x) logistic.test.ML(x, mu = mu, sigma = sigma)

est <- function(x)  MM.logistic(x)
test <- function(x) logistic.test.MM(x)

est <- function(x)  MM.logistic(x, mu = mu)
test <- function(x) logistic.test.MM(x, mu = mu)

est <- function(x)  MM.logistic(x, sigma = sigma)
test <- function(x) logistic.test.MM(x, sigma = sigma)

est <- function(x)  MM.logistic(x, mu = mu, sigma = sigma)
test <- function(x) logistic.test.MM(x, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Student's-t 
lambda <- 5.3; mu <- 2.3; sigma <- 1.8
theta <- c(lambda, mu, sigma) 
dens <- function(x) dStudent(x, lambda, mu, sigma)
cdf  <- function(x) pStudent(x, lambda, mu, sigma)
qf  <- function(p) qStudent(p, lambda, mu, sigma)
gen <- function(n) rStudent(n, lambda, mu, sigma)
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

est  <- function(x) ML.Student(x)
test <- function(x)  Student.test.ML(x)

est  <- function(x) ML.Student(x, lambda = lambda)
test <- function(x)  Student.test.ML(x, lambda = lambda)

est  <- function(x) ML.Student(x, mu = mu)
test <- function(x)  Student.test.ML(x, mu = mu)

est  <- function(x) ML.Student(x, sigma = sigma)
test <- function(x)  Student.test.ML(x, sigma = sigma)

est  <- function(x) ML.Student(x, lambda = lambda, mu = mu)
test <- function(x)  Student.test.ML(x, lambda = lambda, mu = mu)

est  <- function(x) ML.Student(x, lambda = lambda, sigma = sigma)
test <- function(x)  Student.test.ML(x, lambda = lambda, sigma = sigma)

est  <- function(x) ML.Student(x, mu = mu, sigma = sigma)
test <- function(x)  Student.test.ML(x, mu = mu, sigma = sigma)

est  <- function(x) ML.Student(x, lambda = lambda, mu = mu, sigma = sigma)
test <- function(x)  Student.test.ML(x, lambda = lambda, mu = mu, sigma = sigma)

est  <- function(x) MM.Student(x, lambda)
test <- function(x)  Student.test.MM(x, lambda = lambda)

est  <- function(x) MM.Student(x, lambda, mu = mu)
test <- function(x)  Student.test.MM(x, lambda = lambda, mu = mu)

est  <- function(x) MM.Student(x, lambda, sigma = sigma)
test <- function(x)  Student.test.MM(x, lambda = lambda, sigma = sigma)

est  <- function(x) MM.Student(x, lambda, mu = mu, sigma = sigma)
test <- function(x)  Student.test.MM(x, lambda = lambda, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Distributions on (0, infinity)

## ---- log-EPD
lambda <- 1.5; mu <- 2.3; sigma <- 3.8
theta <- c(lambda, mu, sigma)
dens <- function(x) dlogEPD(x, lambda, mu, sigma)
cdf  <- function(x) plogEPD(x, lambda, mu, sigma)
qf  <- function(p) qlogEPD(p, lambda, mu, sigma)
gen <- function(n) rlogEPD(n, lambda, mu, sigma) 
left.bound <- 0; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.logEPD(x)
test <- function(x)  logEPD.test.ML(x)

est <- function(x)  ML.logEPD(x, lambda = lambda)
test <- function(x)  logEPD.test.ML(x, lambda = lambda)

est <- function(x)  ML.logEPD(x, mu = mu)
test <- function(x)  logEPD.test.ML(x, mu = mu)

est <- function(x)  ML.logEPD(x, sigma = sigma)
test <- function(x)  logEPD.test.ML(x, sigma = sigma)

est <- function(x)  ML.logEPD(x, lambda = lambda, mu = mu)
test <- function(x)  logEPD.test.ML(x, lambda = lambda, mu = mu)

est <- function(x)  ML.logEPD(x, lambda = lambda, sigma = sigma)
test <- function(x)  logEPD.test.ML(x, lambda = lambda, sigma = sigma)

est <- function(x)  ML.logEPD(x, mu = mu, sigma = sigma)
test <- function(x)  logEPD.test.ML(x, mu = mu, sigma = sigma)

est <- function(x)  ML.logEPD(x, lambda = lambda, mu = mu, sigma = sigma)
test <- function(x)  logEPD.test.ML(x, lambda = lambda, mu = mu, sigma = sigma)

est <- function(x)  MM.logEPD(x, lambda)
test <- function(x)  logEPD.test.MM(x, lambda)

est <- function(x)  MM.logEPD(x, lambda, mu = mu)
test <- function(x)  logEPD.test.MM(x, lambda, mu = mu)

est <- function(x)  MM.logEPD(x, lambda, sigma = sigma)
test <- function(x)  logEPD.test.MM(x, lambda, sigma = sigma)

est <- function(x)  MM.logEPD(x, lambda, mu = mu, sigma = sigma)
test <- function(x)  logEPD.test.MM(x, lambda, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- log-Laplace
mu <- 2.3; sigma <- 1.8
theta <- c(mu, sigma)
dens <- function(x) dlogLaplace(x, mu, sigma)
cdf  <- function(x) plogLaplace(x, mu, sigma)
qf  <- function(p) qlogLaplace(p, mu, sigma)
gen <- function(n) rlogLaplace(n, mu, sigma) 
left.bound <- 0; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.logLaplace(x)
test <- function(x)  logLaplace.test.ML(x)

est <- function(x)  ML.logLaplace(x, mu = mu)
test <- function(x)  logLaplace.test.ML(x, mu = mu)

est <- function(x)  ML.logLaplace(x, sigma = sigma)
test <- function(x)  logLaplace.test.ML(x, sigma = sigma)

est <- function(x)  ML.logLaplace(x, mu = mu, sigma = sigma)
test <- function(x)  logLaplace.test.ML(x, mu = mu, sigma = sigma)

est <- function(x)  MM.logLaplace(x)
test <- function(x)  logLaplace.test.MM(x)

est <- function(x)  MM.logLaplace(x, mu = mu)
test <- function(x)  logLaplace.test.MM(x, mu = mu)

est <- function(x)  MM.logLaplace(x, sigma = sigma)
test <- function(x)  logLaplace.test.MM(x, sigma = sigma)

est <- function(x)  MM.logLaplace(x, mu = mu, sigma = sigma)
test <- function(x)  logLaplace.test.MM(x, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- lognormal
mu <- 2.3; sigma <-  1.8
theta <- c(mu, sigma)
dens <- function(x) dlognormal(x, mu, sigma)
cdf  <- function(x) plognormal(x, mu, sigma)
qf  <- function(p) qlognormal(p, mu, sigma)
gen <- function(n) rlognormal(n, mu, sigma) 
left.bound <- 0; right.bound <- Inf
quant <- mu + 1.4

est <- function(x)  ML.lognormal(x)
test <- function(x)  lognormal.test.ML(x)

est <- function(x)  ML.lognormal(x, mu = mu)
test <- function(x)  lognormal.test.ML(x, mu = mu)

est <- function(x)  ML.lognormal(x, sigma = sigma)
test <- function(x)  lognormal.test.ML(x, sigma = sigma)

est <- function(x)  ML.lognormal(x, mu = mu, sigma = sigma)
test <- function(x)  lognormal.test.ML(x, mu = mu, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- half-EPD
lambda <- 1.5; sigma <- 2.3
theta <- c(lambda, sigma)
dens <- function(x) dhalfEPD(x, lambda, sigma)
cdf  <- function(x) phalfEPD(x, lambda, sigma)
qf  <- function(p) qhalfEPD(p, lambda, sigma)
gen <- function(n) rhalfEPD(n, lambda, sigma) 
left.bound <- 0; right.bound <- Inf
quant <- sigma * 1.5

est <- function(x)  ML.halfEPD(x)
test <- function(x)  halfEPD.test.ML(x)

est <- function(x)  ML.halfEPD(x, lambda = lambda)
test <- function(x)  halfEPD.test.ML(x, lambda = lambda)

est <- function(x)  ML.halfEPD(x, sigma = sigma)
test <- function(x)  halfEPD.test.ML(x, sigma = sigma)

est <- function(x)  ML.halfEPD(x, lambda = lambda, sigma = sigma)
test <- function(x)  halfEPD.test.ML(x, lambda = lambda, sigma = sigma)

est <- function(x)  MM.halfEPD(x, lambda = lambda)
test <- function(x)  halfEPD.test.MM(x, lambda = lambda)

est <- function(x)  MM.halfEPD(x, lambda = lambda, sigma = sigma)
test <- function(x)  halfEPD.test.MM(x, lambda = lambda, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- GG 
lambda <- 1.5; beta <- 2.3; rho <- 1.8
theta <- c(lambda, beta, rho)
dens <- function(x) dGG(x, lambda, beta, rho)
cdf  <- function(x) pGG(x, lambda, beta, rho)
qf  <- function(p) qGG(p, lambda, beta, rho)
gen <- function(n) rGG(n, lambda, beta, rho) 
left.bound <- 0; right.bound <- Inf
quant <- beta * 1.5

est <- function(x)  ML.GG(x)
test <- function(x) GG.test.ML(x)

est <- function(x)  ML.GG(x, lambda = lambda)
test <- function(x) GG.test.ML(x, lambda = lambda)

est <- function(x)  ML.GG(x, beta = beta)
test <- function(x) GG.test.ML(x, beta = beta)

est <- function(x)  ML.GG(x, rho = rho)
test <- function(x) GG.test.ML(x, rho = rho)

est <- function(x)  ML.GG(x, lambda = lambda, beta = beta)
test <- function(x) GG.test.ML(x, lambda = lambda, beta = beta)

est <- function(x)  ML.GG(x, lambda = lambda, rho = rho)
test <- function(x) GG.test.ML(x, lambda = lambda, rho = rho)

est <- function(x)  ML.GG(x, beta = beta, rho = rho)
test <- function(x) GG.test.ML(x, beta = beta, rho = rho)

est <- function(x)  ML.GG(x, lambda = lambda, beta = beta, rho = rho)
test <- function(x) GG.test.ML(x, lambda = lambda, beta = beta, rho = rho)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Weibull
beta <- 2.3; rho <- 1.8
theta <- c(beta, rho)
dens <- function(x) dWeibull(x, beta, rho)
cdf  <- function(x) pWeibull(x, beta, rho)
qf  <- function(p) qWeibull(p, beta, rho)
gen <- function(n) rWeibull(n, beta, rho) 
left.bound <- 0; right.bound <- Inf
quant <- beta * 1.5

est <- function(x)  ML.Weibull(x)
test <- function(x) Weibull.test.ML(x)

est <- function(x)  ML.Weibull(x, beta = beta)
test <- function(x) Weibull.test.ML(x, beta = beta)

est <- function(x)  ML.Weibull(x, rho = rho)
test <- function(x) Weibull.test.ML(x, rho = rho)

est <- function(x)  ML.Weibull(x, beta = beta, rho = rho)
test <- function(x) Weibull.test.ML(x, beta = beta, rho = rho)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Frechet
beta <- 2.3; rho <- 1.8
theta <- c(beta, rho)
dens <- function(x) dFrechet(x, beta, rho)
cdf  <- function(x) pFrechet(x, beta, rho)
qf  <- function(p) qFrechet(p, beta, rho)
gen <- function(n) rFrechet(n, beta, rho)
left.bound <- 0; right.bound <- Inf
quant <- beta * 1.5

est <- function(x)  ML.Frechet(x)
test <- function(x) Frechet.test.ML(x)

est <- function(x)  ML.Frechet(x, beta = beta)
test <- function(x) Frechet.test.ML(x, beta = beta)

est <- function(x)  ML.Frechet(x, rho = rho)
test <- function(x) Frechet.test.ML(x, rho = rho)

est <- function(x)  ML.Frechet(x, beta = beta, rho = rho)
test <- function(x) Frechet.test.ML(x, beta = beta, rho = rho)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Gompertz
beta <- 2.45; rho <- 1.45
theta <- c(beta, rho)
dens <- function(x) dGompertz(x, beta, rho)
cdf  <- function(x) pGompertz(x, beta, rho)
qf  <- function(p) qGompertz(p, beta, rho)
gen <- function(n)  rGompertz(n, beta, rho) 
left.bound <- 0; right.bound <- Inf
quant <- .4

est <- function(x)  ML.Gompertz(x) 
test <- function(x)  Gompertz.test.ML(x)

est <- function(x)  ML.Gompertz(x, beta = beta) 
test <- function(x)  Gompertz.test.ML(x, beta = beta)

est <- function(x)  ML.Gompertz(x, rho = rho) 
test <- function(x)  Gompertz.test.ML(x, rho = rho)

est <- function(x)  ML.Gompertz(x, beta = beta, rho = rho) 
test <- function(x)  Gompertz.test.ML(x, beta = beta, rho = rho)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- log-logistic
beta <- 2.3; rho <- 1.8
theta <- c(beta, rho)
dens <- function(x) dloglogistic(x, beta, rho)
cdf  <- function(x) ploglogistic(x, beta, rho)
qf  <- function(p) qloglogistic(p, beta, rho)
gen <- function(n) rloglogistic(n, beta, rho)
left.bound <- 0; right.bound <- Inf
quant <- beta * 1.5

est <- function(x)  ML.loglogistic(x)
test <- function(x) loglogistic.test.ML(x)

est <- function(x)  ML.loglogistic(x, beta = beta)
test <- function(x) loglogistic.test.ML(x, beta = beta)

est <- function(x)  ML.loglogistic(x, rho = rho)
test <- function(x) loglogistic.test.ML(x, rho = rho)

est <- function(x)  ML.loglogistic(x, beta = beta, rho = rho)
test <- function(x) loglogistic.test.ML(x, beta = beta, rho = rho)

est <- function(x)  MM.loglogistic(x)
test <- function(x) loglogistic.test.MM(x)

est <- function(x)  MM.loglogistic(x, beta = beta)
test <- function(x) loglogistic.test.MM(x, beta = beta)

est <- function(x)  MM.loglogistic(x, rho = rho)
test <- function(x) loglogistic.test.MM(x, rho = rho)

est <- function(x)  MM.loglogistic(x, beta = beta, rho = rho)
test <- function(x) loglogistic.test.MM(x, beta = beta, rho = rho)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- gamma
alpha <- 2.3; beta <- 1.8
theta <- c(alpha, beta)
dens <- function(x) dgamma(x, shape = alpha, scale = beta)
cdf  <- function(x) pgamma(x, shape = alpha, scale = beta)
qf  <- function(p) qgamma(p, shape = alpha, scale = beta)
gen <- function(n) rgamma(n, shape = alpha, scale = beta) 
left.bound <- 0; right.bound <- Inf
quant <- 3

est <- function(x)  ML.gamma(x)
test <- function(x)  gamma.test.ML(x)

est <- function(x)  ML.gamma(x, alpha = alpha)
test <- function(x)  gamma.test.ML(x, alpha = alpha)

est <- function(x)  ML.gamma(x, beta = beta)
test <- function(x)  gamma.test.ML(x, beta = beta)

est <- function(x)  ML.gamma(x, alpha = alpha, beta = beta)
test <- function(x)  gamma.test.ML(x, alpha = alpha, beta = beta)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

# inverse-gamma distribution
alpha <- 3.45; beta <- 5.45
theta <- c(alpha,beta)
dens <- function(x) dinversegamma(x, alpha, beta)
cdf <- function(x) pinversegamma(x, alpha, beta)
qf <- function(p) qinversegamma(p, alpha, beta)
gen <- function(n) rinversegamma(n, alpha, beta)
left.bound <- 0; right.bound <- Inf
quant <- 2.3

est <- function(x)  ML.inversegamma(x)
test <- function(x) inversegamma.test.ML(x)

est <- function(x)  ML.inversegamma(x, alpha = alpha)
test <- function(x) inversegamma.test.ML(x, alpha = alpha)

est <- function(x)  ML.inversegamma(x, beta = beta)
test <- function(x) inversegamma.test.ML(x, beta = beta)

est <- function(x)  ML.inversegamma(x, alpha = alpha, beta = beta)
test <- function(x) inversegamma.test.ML(x, alpha = alpha, beta = beta)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- beta-prime distribution
alpha <- 3.45; beta <- 5.45
theta <- c(alpha, beta)
dens <- function(x) dbetaprime(x, alpha, beta)
cdf <- function(x) pbetaprime(x, alpha, beta)
qf <- function(p) qbetaprime(p, alpha, beta)
gen <- function(n) rbetaprime(n, alpha, beta)
left.bound <- 0; right.bound <- Inf
quant <- 1.3

est <- function(x)  ML.betaprime(x) 
test <- function(x)  betaprime.test.ML(x) 

est <- function(x)  ML.betaprime(x, alpha = alpha) 
test <- function(x)  betaprime.test.ML(x, alpha = alpha) 

est <- function(x)  ML.betaprime(x, beta = beta) 
test <- function(x)  betaprime.test.ML(x, beta = beta) 

est <- function(x)  ML.betaprime(x, alpha = alpha, beta = beta) 
test <- function(x)  betaprime.test.ML(x, alpha = alpha, beta = beta) 

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Lomax distribution, (which is a special case of the Pareto distribution of type II with a location fixed to 0)
alpha <- 2.3; sigma <- 1.8
theta <- c(alpha, sigma)
dens <- function(x) dLomax(x, alpha, sigma)
cdf  <- function(x) pLomax(x, alpha, sigma)
qf  <- function(p) qLomax(p, alpha, sigma)
gen <- function(n) rLomax(n, alpha, sigma)
left.bound <- 0; right.bound <- Inf
quant <- .7

est <- function(x)  ML.Lomax(x)
test <- function(x)  Lomax.test.ML(x)

est <- function(x)  ML.Lomax(x, alpha = alpha)
test <- function(x)  Lomax.test.ML(x, alpha = alpha)

est <- function(x)  ML.Lomax(x, sigma = sigma)
test <- function(x)  Lomax.test.ML(x, sigma = sigma)

est <- function(x)  ML.Lomax(x, alpha = alpha, sigma = sigma)
test <- function(x)  Lomax.test.ML(x, alpha = alpha, sigma = sigma)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Nakagami distribution
lambda <- 3.45; omega <- 5.45
theta <- c(lambda,omega)
dens <- function(x) dNakagami(x, lambda, omega)
cdf <- function(x) pNakagami(x, lambda, omega)
qf <- function(p) qNakagami(p, lambda, omega)
gen <- function(n) rNakagami(n, lambda, omega)
left.bound <- 0; right.bound <- Inf 
quant <- 2.3

est <- function(x)  ML.Nakagami(x)
test <- function(x) Nakagami.test.ML(x)

est <- function(x)  ML.Nakagami(x, lambda = lambda)
test <- function(x) Nakagami.test.ML(x, lambda = lambda)

est <- function(x)  ML.Nakagami(x, omega = omega)
test <- function(x) Nakagami.test.ML(x, omega = omega)

est <- function(x)  ML.Nakagami(x, lambda = lambda, omega = omega)
test <- function(x) Nakagami.test.ML(x, lambda = lambda, omega = omega)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- inverse-Gaussian
mu <- 5.45; lambda <- 1.45
theta <- c(mu, lambda)
dens <- function(x) dinversegaussian(x, mu, lambda)
cdf  <- function(x) pinversegaussian(x, mu, lambda)
qf  <- function(p) qinversegaussian(p, mu, lambda)
gen <- function(n) rinversegaussian(n, mu, lambda)
left.bound <- 0; right.bound <- Inf
quant <- 3

est <- function(x)  ML.inversegaussian(x) 
test <- function(x)  inversegaussian.test.ML(x)

est <- function(x)  ML.inversegaussian(x, mu = mu) 
test <- function(x)  inversegaussian.test.ML(x, mu = mu)

est <- function(x)  ML.inversegaussian(x, lambda = lambda) 
test <- function(x)  inversegaussian.test.ML(x, lambda = lambda)

est <- function(x)  ML.inversegaussian(x, mu = mu, lambda = lambda) 
test <- function(x)  inversegaussian.test.ML(x, mu = mu, lambda = lambda)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- exponential
beta <- 2.3
theta <- beta
dens <- function(x) dExp(x, beta)
cdf  <- function(x) pExp(x, beta)
qf  <- function(q) qExp(q, beta)
gen <- function(n) rExp(n, beta) 
left.bound <- 0; right.bound <- Inf
quant <- beta * 1.5

est <- function(x)  ML.exp(x)
test <- function(x) exp.test.ML(x)

est <- function(x)  ML.exp(x, beta = beta)
test <- function(x) exp.test.ML(x, beta = beta)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- half-normal
delta <- 2.3
theta <- delta
dens <- function(x) dhalfnormal(x, delta)
cdf  <- function(x) phalfnormal(x, delta)
qf  <- function(p) qhalfnormal(p, delta)
gen <- function(n) rhalfnormal(n, delta) 
left.bound <- 0; right.bound <- Inf; 
quant <- delta * 1.5

est <- function(x)  ML.halfnormal(x)
test <- function(x) halfnormal.test.ML(x)

est <- function(x)  ML.halfnormal(x, delta = delta)
test <- function(x) halfnormal.test.ML(x, delta = delta)

est <- function(x)  MM.halfnormal(x)
test <- function(x) halfnormal.test.MM(x)

est <- function(x)  MM.halfnormal(x, delta = delta)
test <- function(x) halfnormal.test.MM(x, delta = delta)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Rayleigh
delta <- 2.3
theta <- delta
dens <- function(x) dRayleigh(x, delta)
cdf  <- function(x) pRayleigh(x, delta)
qf  <- function(p) qRayleigh(p, delta)
gen <- function(n) rRayleigh(n, delta) 
left.bound <- 0; right.bound <- Inf; 
quant <- delta * 1.5

est <- function(x)  ML.Rayleigh(x)
test <- function(x) Rayleigh.test.ML(x)

est <- function(x)  ML.Rayleigh(x, delta = delta)
test <- function(x) Rayleigh.test.ML(x, delta = delta)

est <- function(x)  MM.Rayleigh(x)
test <- function(x) Rayleigh.test.MM(x)

est <- function(x)  MM.Rayleigh(x, delta = delta)
test <- function(x) Rayleigh.test.MM(x, delta = delta)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Maxwell-Boltzmann
delta <- 2.3
theta <- delta
dens <- function(x) dMaxwell(x, delta)
cdf  <- function(x) pMaxwell(x, delta)
qf  <- function(p) qMaxwell(p, delta)
gen <- function(n) rMaxwell(n, delta) 
left.bound <- 0; right.bound <- Inf; 
quant <- delta * 1.5

est <- function(x)  ML.Maxwell(x)
test <- function(x) Maxwell.test.ML(x)

est <- function(x)  ML.Maxwell(x, delta = delta)
test <- function(x) Maxwell.test.ML(x, delta = delta)

est <- function(x)  MM.Maxwell(x)
test <- function(x) Maxwell.test.MM(x)

est <- function(x)  MM.Maxwell(x, delta = delta)
test <- function(x) Maxwell.test.MM(x, delta = delta)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- chi-squared distribution 
k <- 5
theta <- k
dens <- function(x) dchisq(x, df = k)
cdf  <- function(x) pchisq(x, df = k)
qf  <- function(p) qchisq(p, df = k)
gen <- function(n) rchisq(n, df = k)
left.bound <- 0; right.bound <- Inf
quant <- 2.3

est <- function(x)  ML.chisquared(x)
test <- function(x)  chisquared.test.ML(x)

est <- function(x)  ML.chisquared(x, k = k)
test <- function(x)  chisquared.test.ML(x, k = k)

est <- function(x)  MM.chisquared(x)
test <- function(x)  chisquared.test.MM(x)

est <- function(x)  MM.chisquared(x, k = k)
test <- function(x)  chisquared.test.MM(x, k = k)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Distribution on (1, infinity)

## ---- Pareto distribution
alpha <- 20.3
theta <- alpha
dens <- function(x) dPareto(x, alpha)
cdf  <- function(x) pPareto(x, alpha)
qf  <- function(p) qPareto(p, alpha)
gen <- function(n) rPareto(n, alpha) 
left.bound <- 1; right.bound <- Inf
quant <- 1 + 1.5 / alpha

est <- function(x)  ML.Pareto(x)
test <- function(x) Pareto.test.ML(x)

est <- function(x)  ML.Pareto(x, alpha = alpha)
test <- function(x) Pareto.test.ML(x, alpha = alpha)

est <- function(x)  ML.Pareto(x / min(x))
test <- function(x) Pareto.test.ML(x / min(x))

est <- function(x)  ML.Pareto(x / min(x), alpha = alpha)
test <- function(x) Pareto.test.ML(x / min(x), alpha = alpha)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Distributions on (0,1)

## ---- beta distribution
alpha <- 2.45; beta <- 1.45
theta <- c(alpha, beta)
dens <- function(x) dbeta(x, alpha, beta)
cdf  <- function(x) pbeta(x, alpha, beta)
qf  <- function(p) qbeta(p, alpha, beta)
gen <- function(n)  rbeta(n, alpha, beta)
left.bound <- 0; right.bound <- 1
quant <- .678

est  <- function(x)  ML.beta(x) 
test <- function(x)  beta.test.ML(x)

est  <- function(x)  ML.beta(x, alpha = alpha) 
test <- function(x)  beta.test.ML(x, alpha = alpha)

est  <- function(x)  ML.beta(x, beta = beta) 
test <- function(x)  beta.test.ML(x, beta = beta)

est  <- function(x)  ML.beta(x, alpha = alpha, beta = beta) 
test <- function(x)  beta.test.ML(x, alpha = alpha, beta = beta)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Kumaraswamy
alpha <- 3.45; beta <- 5.45
theta <- c(alpha, beta)
dens <- function(x) dKumaraswamy(x, alpha, beta)
cdf  <- function(x) pKumaraswamy(x, alpha, beta)
qf <- function(p)   qKumaraswamy(p, alpha, beta)
gen <- function(n)  rKumaraswamy(n, alpha, beta)
left.bound <- 0; right.bound <- 1
quant <- .678

est  <- function(x) ML.Kumaraswamy(x)
test <- function(x)  Kumaraswamy.test.ML(x)

est  <- function(x) ML.Kumaraswamy(x, alpha = alpha)
test <- function(x)  Kumaraswamy.test.ML(x, alpha = alpha)

est  <- function(x) ML.Kumaraswamy(x, beta = beta)
test <- function(x)  Kumaraswamy.test.ML(x, beta = beta)

est  <- function(x) ML.Kumaraswamy(x, alpha = alpha, beta = beta)
test <- function(x)  Kumaraswamy.test.ML(x, alpha = alpha, beta = beta)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Distributions on (a,b)

## ---- uniform
a <- -3.45; b <- 5.45
theta <- c(a, b) 
dens <- function(x) duniform(x, a, b)
cdf  <- function(x) puniform(x, a, b)
qf   <- function(p) quniform(p, a, b)
gen  <- function(n) runiform(n, a, b)
left.bound <- a; right.bound <- b
quant <- .678

est  <- function(x) ML.uniform(x) 
test <- function(x)  uniform.test.ML(x)

est  <- function(x) ML.uniform(x, a = a) 
test <- function(x)  uniform.test.ML(x, a = a)

est  <- function(x) ML.uniform(x, b = b) 
test <- function(x)  uniform.test.ML(x, b = b)

est  <- function(x) ML.uniform(x, a = a, b = b) 
test <- function(x)  uniform.test.ML(x, a = a, b = b)

validation1()
validation2(est, test)
res <- validation3(test, nb.sim = 5000, n = 1000)
#res <- validation3(test, nb.sim = 1000, n = 1000)
plot.ellipse.tau.bar.sqrtn(res$tau.bar.sqrtn[1,], res$tau.bar.sqrtn[2,], solve(res$Sigma.inv))
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[1,])
hist.tau.bar.sqrtn(res$tau.bar.sqrtn[2,])
hist.Tn(res$Tn)
validation4(res$Tn)

## ---- Other densities

## ---- APD used in Section 5.3
lambda <- 1.5; alpha <- 0.4; rho <- 2; mu <- 2.3; sigma <- 1.8
dens <- function(x) dAPD(x, lambda, alpha, rho, mu, sigma)
cdf  <- function(x) pAPD(x, lambda, alpha, rho, mu, sigma)
qf  <- function(p) qAPD(p, lambda, alpha, rho, mu, sigma)
gen <- function(n) rAPD(n, lambda, alpha, rho, mu, sigma) 
left.bound <- -Inf; right.bound <- Inf
quant <- mu + 1.4

validation1()
x <- gen(10000)
val <- matrix(c(
  adaptIntegrate(f = dens, lower = left.bound, upper = right.bound, tol = 1e-9)$integral, 1, 
  qf(cdf(quant)), quant, # validation of qf and cdf
  round(cdf(quant), 6), round(adaptIntegrate(f = dens, lower = left.bound, upper = quant, tol = 1e-9)$integral, 6),
  round(cdf(quant), 6), round(mean(x < quant), 6)
), byrow = T, ncol = 2) 
val <- cbind(val, round(abs(val[,1] - val[,2]), 6))
rownames(val) <- c("density integrates to 1", "CDF and inv-CDF", "CDF", "gen")
colnames(val) <- c("Component 1", "Component2", "Absolute difference")
val

# ====================================================================================
# Remark 2.4: Interpretation of the test statistic as n tends to infinity
# ====================================================================================

## ---- Different functions

## ---- distribution of U = F(X), 0<U<1 (the distribution of X is G, F is the distribution tested in H0)
U.dens <- function(u) {
  x <- F.inv(u)  
  gx.dens(x) / f.dens(x)
} 
## ---- distribution of V = cos(2 * pi * U) = cos(2 * pi * F(X)), -1 < V < 1 (the distribution of X is G)
Cos.dens <- function(v) {
  a <- acos(v) / (2 * pi) # decreasing fonction of -1<v<1 taking values between 1/2 and 0
  res <- (U.dens(a) + U.dens(1 - a)) / (2 *  pi * sqrt(1 - v ^ 2))
  return(res)
}  
## ---- distribution of V = sin(2 * pi * U) = sin(2 * pi * F(X)), -1 < V < 1 (the distribution of X is G)
Sin.dens <- function(v) {
  a <- asin(v) / (2 * pi) # increasing fonction of -1<v<1 taking values between -1/4 and 1/4
  res <- numeric(length(v))   
  pos1 <- (v < 0)
  a1 <- a[pos1] # -1/4 < a1 < 0 
  res[pos1] <- U.dens(1 + a1) + U.dens(1 / 2 - a1)
  pos2 <- (v > 0)
  a2 <- a[pos2] # 0 < a1 < 1/4
  res[pos2] <- U.dens(1 / 2 - a2) + U.dens(a2)
  res <- res / (2 *  pi * sqrt(1 - v ^ 2))
  res[v == 0] <- Inf 
  return(res)
}  
## ---- If F = G, under H0, distribution of V = cos(2 * pi * F(X)) or V = sin(2 * pi * F(X)), -1 < V < 1
darcsine <- function(x, a, b) 1 / pi / sqrt((x - a) * (b - x)) # arcsine distribution on (a,b)
hv.cos.sin.dens.H0 <- function(v) darcsine(v, -1, 1)

## ---- Histogram of the densities f (H0) and g (X)
hist.x <- function() {
  n <- 1e6
  ui <- seq(0 + 1 / n / 2, 1 -  1 / n / 2, length = n)
  xi <- Gx.inv(ui) # generates perfectly G distributed sample
  xlimm <- c(min(Gx.inv(.001), F.inv(.001)), max(Gx.inv(1- .001), F.inv(1 - .001)))
  ylimm <- c(0, max(c(gx.dens(seq(xlimm[1], xlimm[2], .1)), f.dens(seq(xlimm[1], xlimm[2], .1))))) 
  hist(xi, prob = TRUE, breaks = 200, col = 'grey90', border = 'white', ylim = ylimm, xlim = xlimm)
  curve(gx.dens(x), add = TRUE, col = "grey90", n = 1000, lwd = 2)
  curve(f.dens(x), add = TRUE, col = "black", n = 1000, lwd = 2)
  segments(x0 = F.inv(c(.25, 0.5, .75)), y0 = 0, x1 = F.inv(c(.25, 0.5, .75)), y1 = f.dens(F.inv(c(.25, 0.5, .75))), col = "grey0", lwd = 2)
  xx1 <- uniroot(function(x) gx.dens(x) - f.dens(x), c(F.inv(.001), F.inv(.20)), tol = 1e-09)$root
  xx2 <- uniroot(function(x) gx.dens(x) - f.dens(x), c(F.inv(.80), F.inv(.999)), tol = 1e-09)$root
  #gx.dens(xx1); f.dens(xx1)
  #gx.dens(xx2); f.dens(xx2)
  segments(x0 = xx1, y0 = 0, x1 = xx1, y1 = f.dens(xx1), col = "grey0", lwd = 2)
  segments(x0 = xx2, y0 = 0, x1 = xx2, y1 = f.dens(xx2), col = "grey0", lwd = 2)
}
ana.dens <- function() {
  n <- 1e6
  ui <- seq(0 + 1 / n / 2, 1 -  1 / n / 2, length = n)
  xi <- Gx.inv(ui) # generates perfectly G distributed sample
  xx1 <- uniroot(function(x) gx.dens(x) - f.dens(x), c(F.inv(.001), F.inv(.20)), tol = 1e-09)$root
  xx2 <- uniroot(function(x) gx.dens(x) - f.dens(x), c(F.inv(.80), F.inv(.999)), tol = 1e-09)$root
  Z.tau.bar <- round(test(xi)$Z.tau.bar, 3)
  left.mass <- Gx.cdf(F.inv(.5)) # usually > 0.5 if Sn>0 and < 0.5 if Sn<0, or = 0 is Sn = 0
  left.tail.mass <- Gx.cdf(F.inv(.25)) # usually > 0.25 if Cn>0 and < 0.25 if Cn<0
  right.tail.mass <- 1 - Gx.cdf(F.inv(.75)) # usually > 0.25 if Cn>0 and < 0.25 if Cn<0
  center.mass <- Gx.cdf(F.inv(.75)) - Gx.cdf(F.inv(.25)) # usually < 0.25 if Cn>0 and > 0.25 if Cn<0
  extreme.left.tail.mass <- F.cdf(xx1)
  extreme.right.tail.mass <- 1 - F.cdf(xx2)
  return(list(Z.tau.bar = Z.tau.bar, left.mass = left.mass, right.mass = 1 - left.mass, left.tail.mass = left.tail.mass, center.mass = center.mass, right.tail.mass = right.tail.mass, extreme.left.tail.mass = extreme.left.tail.mass,extreme.right.tail.mass =extreme.right.tail.mass))
}
## ---- Histogram of U.dens, the density of U = F(X) , 0<U<1 (the distribution of X is G, F is the distribution tested in H0)
hist.F <- function() {
  n <- 1e6
  ui <- seq(0 + 1 / n / 2, 1 -  1 / n / 2, length = n)
  xi <- Gx.inv(ui) # generates perfectly G distributed sample
  Fi <- F.cdf(xi)
  hist(Fi, prob = T, xlim = c(0,1), breaks = 50, col = 'grey90', border = 'white'); 
  curve(U.dens(x), add = TRUE, col = "grey70", n = 1000, lwd = 2)
  segments(x0 = 0, y0 = 1, x1 = 1, y1 = 1, col = "black", lwd = 2) # uniform distribution under H0
  segments(x0 = c(0, .25, 0.5, .75, 1), y0 = 0, x1 = c(0, .25, 0.5, .75, 1), y1 = 1, col = "grey0", lwd = 2)
}
## ---- Histogram of Cos.dens, the density of V = cos(2 * pi * U) = cos(2 * pi * F(X)), -1 < V < 1
hist.Cos <- function() {
  n <- 1e6
  ui <- seq(0 + 1 / n / 2, 1 -  1 / n / 2, length = n)
  xi <- Gx.inv(ui) # generates perfectly G distributed sample
  Fi <- F.cdf(xi)
  Vi.cos <- cos(2 * pi * Fi)
  hist(Vi.cos, prob = T, xlim = c(-1,1), breaks = 50, col = 'grey90', border = 'white'); 
  curve(Cos.dens(x), add = TRUE, col = "grey70", n = 1000, lwd = 2)
  curve(hv.cos.sin.dens.H0(x), add = TRUE, col = "black", n = 1000, lwd = 2) # under H0
  segments(x0 = 0, y0 = 0, x1 = 0, y1 = hv.cos.sin.dens.H0(0), col = "grey0", lwd = 2)
}
## ---- Histogram of Sin.dens, the density of V = sin(2 * pi * U) = sin(2 * pi * F(X)), -1 < V < 1
hist.Sin <- function() {
  n <- 1e6
  ui <- seq(0 + 1 / n / 2, 1 -  1 / n / 2, length = n)
  xi <- Gx.inv(ui) # generates perfectly G distributed sample
  Fi <- F.cdf(xi)
  Vi.sin <- sin(2 * pi * Fi)
  hist(Vi.sin, prob = T, xlim = c(-1,1), breaks = 50, col = 'grey90', border = 'white'); 
  curve(Sin.dens(x), add = TRUE, col = "grey70", n = 1000, lwd = 2)
  curve(hv.cos.sin.dens.H0(x), add = TRUE, col = "black", n = 1000, lwd = 2) # under H0
  segments(x0 = 0, y0 = 0, x1 = 0, y1 = hv.cos.sin.dens.H0(0), col = "grey0", lwd = 2)
}  
## ---- end of functions

## ---- Validation of U.dens, Cos.dens and Sin.dens 
gx.dens <- function(x) dlaplace(x, 0, 1)
Gx.inv <- function(p) qlaplace(p, 0, 1)
f.dens <- function(x) dnorm(x, 0, sqrt(2))
F.cdf <- function(x) pnorm(x, 0, sqrt(2))
F.inv <- function(p) qnorm(p, 0, sqrt(2))
n <- 1e6; ui <- seq(0 + 1 / n / 2, 1 -  1 / n / 2, length = n); 
xi <- Gx.inv(ui) # generates perfectly G distributed sample
Fi <- F.cdf(xi)
Vi.cos <- cos(2 * pi * Fi)
Vi.sin <- sin(2 * pi * Fi)
a <- .36

## ---- Validation of the U.dens distribution
print(round(matrix(c(
  adaptIntegrate(f = function(u) U.dens(u), lower = 0+1e-10, upper = 1-1e-10, tol = 1e-9)$integral,
  1,
  adaptIntegrate(f = function(u) U.dens(u) * u, lower = 0+1e-10, upper = 0, tol = 1e-9)$integral +
    adaptIntegrate(f = function(u) U.dens(u) * u, lower = 0, upper = 1-1e-10, tol = 1e-9)$integral, 
  mean(Fi),
  adaptIntegrate(f = function(u) U.dens(u) * u ^ 2, lower = 0+1e-10, upper = 1-1e-10, tol = 1e-9)$integral, 
  mean(Fi ^ 2),
  adaptIntegrate(f = function(u) U.dens(u) * u ^ 4, lower = 0+1e-10, upper = 1-1e-10, tol = 1e-9)$integral, 
  mean(Fi ^ 4),
  adaptIntegrate(f = function(u) U.dens(u) * u ^ 6, lower = 0+1e-10, upper = 1-1e-10, tol = 1e-9)$integral, 
  mean(Fi ^ 6),
  round(mean(Fi < a), 5),
  round(adaptIntegrate(f = function(u) U.dens(u), lower = 0+1e-10, upper = a, tol = 1e-9)$integral, 5)
), byrow = T, ncol = 2), 6))

## ---- Validation of the Cos.dens distribution
print(round(matrix(c(
  adaptIntegrate(f = function(v) Cos.dens(v), lower = -1+1e-10, upper = 1-1e-10, tol = 1e-9)$integral, 
  1,
  adaptIntegrate(f = function(v) Cos.dens(v) * v, lower = -1+1e-16, upper = 0, tol = 1e-9)$integral +
    adaptIntegrate(f = function(v) Cos.dens(v) * v, lower = 0, upper = 1-1e-16, tol = 1e-9)$integral ,
  mean(Vi.cos),
  adaptIntegrate(f = function(v) Cos.dens(v) * v ^ 2, lower = -1+1e-10, upper = 1-1e-10, tol = 1e-9)$integral, 
  mean(Vi.cos ^ 2),
  adaptIntegrate(f = function(v) Cos.dens(v) * v ^ 4, lower = -1+1e-10, upper = 1-1e-10, tol = 1e-9)$integral ,
  mean(Vi.cos ^ 4),
  adaptIntegrate(f = function(v) Cos.dens(v) * v ^ 6, lower = -1+1e-10, upper = 1-1e-10, tol = 1e-9)$integral ,
  mean(Vi.cos ^ 6),
  mean(Vi.cos < a),
  adaptIntegrate(f = function(v) Cos.dens(v), lower = -1+1e-10, upper = a, tol = 1e-9)$integral 
), byrow = T, ncol = 2), 6))

## ---- Validation of the Sin.dens distribution
print(round(matrix(c(
  adaptIntegrate(f = function(v) Sin.dens(v), lower = -1+1e-10, upper = 0-1e-10, tol = 1e-9)$integral +
    adaptIntegrate(f = function(v) Sin.dens(v), lower = 0+1e-10, upper = 1-1e-10, tol = 1e-9)$integral ,
  1,
  adaptIntegrate(f = function(v) Sin.dens(v) * v, lower = -1+1e-10, upper = 0 - 1e-10, tol = 1e-9)$integral +
    adaptIntegrate(f = function(v) Sin.dens(v) * v, lower = 0 +1e-10, upper = 1-1e-10, tol = 1e-9)$integral ,
  mean(Vi.sin),
  adaptIntegrate(f = function(v) Sin.dens(v) * v ^ 2, lower = -1+1e-10, upper = 0 - 1e-10, tol = 1e-9)$integral +
    adaptIntegrate(f = function(v) Sin.dens(v) * v ^ 2, lower = 0 +1e-10, upper = 1-1e-10, tol = 1e-9)$integral ,
  mean(Vi.sin ^ 2),
  adaptIntegrate(f = function(v) Sin.dens(v) * v ^ 4, lower = -1+1e-10, upper = 0 - 1e-10, tol = 1e-9)$integral +
    adaptIntegrate(f = function(v) Sin.dens(v) * v ^ 4, lower = 0 +1e-10, upper = 1-1e-10, tol = 1e-9)$integral ,
  mean(Vi.sin ^ 4),
  adaptIntegrate(f = function(v) Sin.dens(v) * v ^ 6, lower = -1+1e-10, upper = 0, tol = 1e-9)$integral +
    adaptIntegrate(f = function(v) Sin.dens(v) * v ^ 6, lower = 0, upper = 1-1e-10, tol = 1e-9)$integral ,
  mean(Vi.sin ^ 6),
  mean(Vi.sin > a),
  adaptIntegrate(f = function(v) Sin.dens(v), lower = a, upper = 1 - 1e-10, tol = 1e-9)$integral 
), byrow = T, ncol = 2), 6))
## ---- End of validations

## ---- Case 1: H0 : F = N(mu, sigma), MLE, X ~  EPD(lambda = 3, mu = 0, sigma = 1)
gx.dens <- function(x) dEPD(x, 3, 0, 1)
Gx.inv <- function(p) qEPD(p, 3, 0, 1)
Gx.cdf <- function(x) pEPD(x, 3, 0, 1)

mu.ML <- mean.APD(mu = 0, sigma = 1, lambda = 3, alpha = 0.5, rho = 3); mu.ML
sigma.ML <- sqrt(mean2.APD(mu = 0, sigma = 1, lambda = 3, alpha = 0.5, rho = 3) - mu.ML ^ 2); sigma.ML
#Validation
round(ML.normal(rEPD(10^6, 3, 0, 1)), 4); round(c(mu.ML, sigma.ML), 4)

f.dens <- function(x) dnorm(x, mu.ML, sigma.ML)
F.cdf <- function(x) pnorm(x, mu.ML, sigma.ML)
F.inv <- function(p) qnorm(p, mu.ML, sigma.ML)
test <- function(x) normal.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 2: H0 : F = N(mu, sigma), MLE, X ~  EPD(lambda = 1.5, mu = 0, sigma = 1)

gx.dens <- function(x) dEPD(x, 1.5, 0, 1)
Gx.inv <- function(p) qEPD(p, 1.5, 0, 1)
Gx.cdf <- function(x) pEPD(x, 1.5, 0, 1)

mu.ML <- mean.APD(mu = 0, sigma = 1, lambda = 1.5, alpha = 0.5, rho = 1.5); mu.ML
sigma.ML <- sqrt(mean2.APD(mu = 0, sigma = 1, lambda = 1.5, alpha = 0.5, rho = 1.5) - mu.ML ^ 2); sigma.ML
#Validation
round(ML.normal(rEPD(10^6, 1.5, 0, 1)), 4); round(c(mu.ML, sigma.ML), 4)

f.dens <- function(x) dnorm(x, mu.ML, sigma.ML)
F.cdf <- function(x) pnorm(x, mu.ML, sigma.ML)
F.inv <- function(p) qnorm(p, mu.ML, sigma.ML)
test <- function(x) normal.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 3: H0 : F = N(mu, sigma), MLE, X ~  APD(lambda = 2, alpha = 0.3, rho = 2, mu = 0, sigma = 1)

gx.dens <- function(x) dAPD(x, 2, 0.3, 2, mu = 0, sigma = 1)
Gx.inv <- function(p) qAPD(p, 2, 0.3, 2, mu = 0, sigma = 1)
Gx.cdf <- function(x) pAPD(x, 2, 0.3, 2, mu = 0, sigma = 1)

mu.ML <- mean.APD(mu = 0, sigma = 1, lambda = 2, alpha = 0.3, rho = 2); mu.ML
sigma.ML <- sqrt(mean2.APD(mu = 0, sigma = 1, lambda = 2, alpha = 0.3, rho = 2) - mean.APD(mu = 0, sigma = 1, lambda = 2, alpha = 0.3, rho = 2) ^ 2); sigma.M
#Validation
round(ML.normal(rAPD(10^6, 2, 0.3, 2, 0, 1)), 4); round(c(mu.ML, sigma.ML), 4)

f.dens <- function(x) dnorm(x, mu.ML, sigma.ML)
F.cdf <- function(x) pnorm(x, mu.ML, sigma.ML)
F.inv <- function(p) qnorm(p, mu.ML, sigma.ML)
test <- function(x) normal.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 4: H0 : F = N(mu, sigma), MLE, X ~  N(0,1)

f.dens <- gx.dens <- function(x) dnorm(x, 0, 1)
F.cdf <- Gx.cdf <-  function(x) pnorm(x, 0, 1)
F.inv <- Gx.inv <- function(p) qnorm(p, 0, 1)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 5: H0 : F = N(mu, sigma), MLE, X ~  Laplace(x, 0, 1)

gx.dens <- function(x) dlaplace(x, 0, 1)
Gx.inv <- function(p) qlaplace(p, 0, 1)
Gx.cdf <- function(x) plaplace(x, 0, 1)

mu.ML <- 0
sigma.ML <- sqrt(2); sigma.ML
#Validation
round(ML.normal(rLaplace(10^6, 0, 1)), 4); round(c(mu.ML, sigma.ML), 4)

f.dens <- function(x) dnorm(x, mu.ML, sigma.ML)
F.cdf <- function(x) pnorm(x, mu.ML, sigma.ML)
F.inv <- function(p) qnorm(p, mu.ML, sigma.ML)
test <- function(x) normal.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 6: H0 : F = Laplace(mu, sigma), MLE, X ~  N(x, 0, 1)

gx.dens <- function(x) dnorm(x, 0, 1)
Gx.inv <- function(p) qnorm(p, 0, 1)
Gx.cdf <- function(x) pnorm(x, 0, 1)

mu.ML <- 0
sigma.ML <- sqrt(2 / pi); sigma.ML
#Validation
round(ML.Laplace(rnorm(10^6, 0, 1)), 4); round(c(mu.ML, sigma.ML), 4)

f.dens <- function(x) dlaplace(x, mu.ML, sigma.ML)
F.cdf <- function(x) plaplace(x, mu.ML, sigma.ML)
F.inv <- function(p) qlaplace(p, mu.ML, sigma.ML)
test <- function(x) Laplace.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 7: H0 : F = N(mu, sigma), MLE, X ~  logistic(x, 0, 1)

gx.dens <- function(x) dlogistic(x, 0, 1)
Gx.inv <- function(p) qlogis(p, 0, 1) 
Gx.cdf <- function(x) plogis(x, 0, 1) 

mu.ML <- 0
sigma.ML <- pi / sqrt(3); sigma.ML
#Validation
round(ML.normal(rlogis(10^6)), 4); round(c(mu.ML, sigma.ML), 4)

f.dens <- function(x) dnorm(x, mu.ML, sigma.ML)
F.cdf <- function(x) pnorm(x, mu.ML, sigma.ML)
F.inv <- function(p) qnorm(p, mu.ML, sigma.ML)
test <- function(x) normal.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 8: H0 : F = logistic(mu, sigma), MLE, X ~  N(0, 1)

gx.dens <- function(x) dnorm(x, 0, 1)
Gx.inv <- function(p) qnorm(p, 0, 1)
Gx.cdf <- function(x) pnorm(x, 0, 1)

mu.ML <- 0
g <- function(sigma) sigma + adaptIntegrate(f = function(x) 2 * x / (1 + exp(x / sigma)) * dnorm(x), lower = -Inf, upper = Inf, tol = 1e-9)$integral
sigma.ML <- uniroot(g, interval = c(.1, 10), tol = 1e-09)$root; sigma.ML
#Validation
round(ML.logistic(rnorm(10^6, 0, 1)), 4); round(c(mu.ML, sigma.ML), 4)

f.dens <- function(x) dlogistic(x, mu.ML, sigma.ML)
F.cdf <- function(x) plogistic(x, mu.ML, sigma.ML)
F.inv <- function(p) qlogis(p, mu.ML, sigma.ML)
test <- function(x) logistic.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 9: H0 : F = N(mu, sigma), MLE, X ~  Gumbel(0, 1)

gx.dens <- function(x) dGumbel(x, 0, 1)
Gx.inv <- function(p) qGumbel(p, 0, 1)
Gx.cdf <- function(x) pGumbel(x, 0, 1)

mu.ML <- -digamma(1); mu.ML
sigma.ML <- pi / sqrt(6); sigma.ML
#Validation
round(ML.normal(rGumbel(10^6)), 4); round(c(mu.ML, sigma.ML), 4)

f.dens <- function(x) dnorm(x, mu.ML, sigma.ML)
F.cdf <- function(x) pnorm(x, mu.ML, sigma.ML)
F.inv <- function(p) qnorm(p, mu.ML, sigma.ML)
test <- function(x) normal.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 10: H0 : F = X2(k), MLE, X ~  Weibull(beta = 7, rho = 2)

gx.dens <- function(x) dWeibull(x, beta = 7, rho = 2)
Gx.inv <- function(p) qWeibull(p, 7, 2)
Gx.cdf <- function(x) pWeibull(x, 7, 2)

g <- function(k) adaptIntegrate(f = function(x) log(x / 2) * dWeibull(x, 7, 2), lower = 0, upper = Inf, tol = 1e-9)$integral - digamma(k / 2)
k.ML <- uniroot(g, interval = c(.1, 100), tol = 1e-09)$root; k.ML
#Validation
ML.chisquared(rWeibull(10^6, 7, 2)) # Validation

f.dens <- function(x) dchisq(x, df = k.ML) 
F.cdf <- function(x) pchisq(x, df = k.ML)
F.inv <- function(p) qchisq(p, df = k.ML) 
test <- function(x) chisquared.test.ML(x)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- Case 11: H0 : F =EPD(lambda = 3, mu, sigma), MLE, X ~  EPD(lambda = 1, mu = 0, sigma = 1)

gx.dens <- function(x) dEPD(x, 1, 0, 1)
Gx.inv <- function(p) qEPD(p, 1, 0, 1)
Gx.cdf <- function(x) pEPD(x, 1, 0, 1)

mu.ML <- 0
sigma.ML <- 1.818 # estimation
round(c(ML.EPD(rEPD(10^6, 1, 0, 1), lambda = 3), 3, mu.ML, sigma.ML), 4) #Validation

f.dens <- function(x) dEPD(x, 3, mu.ML, sigma.ML)
F.cdf <- function(x) pEPD(x, 3, mu.ML, sigma.ML)
F.inv <- function(p) qEPD(p, 3, mu.ML, sigma.ML)
test <- function(x) EPD.test.ML(x, lambda = 3)

hist.x()
ana.dens()
hist.F()    
hist.Cos()
hist.Sin()

## ---- End of the section of intepretation
