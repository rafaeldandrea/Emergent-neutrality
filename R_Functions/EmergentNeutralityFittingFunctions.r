## Functions for fitting James O'Dwyer's species abundance probability distribution

library(dplyr)
library(pracma)
library(dgof)
library(zipfR)
library(dgof)
library(hypergeo)
library(orthopolynom)

## -------------------- Log-series (LS), P(n) = -1/log(1-p)*p^n/n --------------------------------------------

## Probability mass function
LS_pmf = function(p, n)
  - 1 / log(1 - p) * p ^ n / n

## Cumulative distribution function
LS_cdf = function(p, nmax = 5e4) {
  n = seq(nmax)
  pmf = LS_pmf(p, n)
  cdf = 1 + zipfR::Ibeta(p, n + 1, 1e-16) / log(1 - p)
  cum = 1 - cdf + pmf
  return(data.frame(x = n, cdf = cdf, cum = cum))
}

## Maximum likelihood estimation of parameter p
LS_p = function(mu)
  uniroot(
    f = function(p)
      - 1 / log(1 - p) * p / (1 - p) - mu,
    lower = 1e-16,
    upper = 1 - 1e-16
  )$root

## Fit LS, get p-value from discrete Kolmogorov-Smirnov test (implemented via library dgof)
LS_fit = function(data, method) {
  ## transform data into count data frame
  w = if (is.vector(data))
    plyr::count(data)
  else
    data
  
  ## fit parameter p
  p = LS_p(weighted.mean(w$x, w$freq))
  
  ## calculate ecdf given the data
  dtf = LS_cdf(p = p, nmax = 10 * max(w$x))
  myfun = stepfun(w$x, c(0, subset(dtf, x %in% w$x)$cdf))
  
  ## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
  if(method == 'ks') mod = dgof::ks.test(with(w, rep(x, freq)), myfun)
  if(method == 'CVM') mod = dgof::cvm.test(with(w, rep(x, freq)), myfun)
  
  ## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
  return(list(
    dtf = dtf,
    parameter = p,
    statistic = mod$statistic,
    p.value = mod$p.value
  ))
}

## ----------------- Inverse quadratic distribution (IQ), P(n) ~ 1/n/(1+a*n) -----------------------------------------

## Probability mass function
IQ_pmf = function(a, n) {
  k = 1 / (digamma(1 + 1 / a) - digamma(1))
  k * 1 / (n * (1 + a * n))
}

## Cumulative distribution function
IQ_cdf = function(a, nmax = 5e4) {
  n = seq(nmax)
  pmf = IQ_pmf(a, n)
  cdf = cumsum(pmf)
  cum = sum(pmf) - cdf + pmf
  return(data.frame(x = n, cdf = cdf, cum = cum))
}

## Maximum likelihood estimation of parameter a
IQ_a = function(data, amin, amax) {
  w = if (is.vector(data))
    plyr::count(data)
  else
    data
  x = w$x
  y = w$freq
  a = uniroot(
    f = function(a)
      1 / a ^ 2 * pracma::psi(1, 1 + 1 / a) / (digamma(1 + 1 / a) - digamma(1)) *
      sum(y) - sum(x * y / (1 + a * x)),
    lower = amin,
    upper = amax
  )$root
  return(c(a = a))
}

## Fit IQ, get p-value from discrete Kolmogorov-Smirnov test (implemented via library dgof)
IQ_fit = function(data,
                  amin = 1e-16,
                  amax = 1 - 1e-16) {
  ## fit parameter a
  a = IQ_a(data, amin, amax)
  
  ## transform data into count data frame
  w = if (is.vector(data))
    plyr::count(data)
  else
    data
  
  ## calculate ecdf given the data
  dtf = IQ_cdf(a = a, nmax = 10 * max(w$x))
  myfun = stepfun(w$x, c(0, subset(dtf, x %in% w$x)$cdf))
  
  ## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
  mod = dgof::ks.test(with(w, rep(x, freq)), myfun)
  
  ## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
  return(list(
    dtf = dtf,
    parameter = a,
    statistic = mod$statistic,
    p.value = mod$p.value
  ))
}


IQ_KS = function(data,
                 amin = 1e-16,
                 amax = 1 - 1e-16,
                 numnulls = 100) {
  ## fit parameter a
  a = IQ_a(data, amin, amax)
  
  ## transform data into count data frame
  w = if (is.vector(data))
    plyr::count(data)
  else
    data
  
  ## calculate ecdf given the data
  dtf = IQ_cdf(a = a, nmax = 10 * max(w$x))
  myfun = stepfun(w$x, c(0, subset(dtf, x %in% w$x)$cdf))
  
  ## nulls
  n = seq(1e5)
  nulls = sapply(seq(numnulls), function(null) {
    set.seed(null)
    sample(
      n,
      size = sum(w$freq),
      replace = TRUE,
      prob = IQ_pmf(a, n)
    )
  })
  
  
  ## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
  x = dgof::ks.test(rep(w$x, w$freq), myfun)$statistic
  z = apply(nulls, 2, function(null)
    dgof::ks.test(null, myfun)$statistic)
  
  return(list(obs.D = x, exp.D = z))
}


## ------------- Poisson distribution (Pois), P(n) = exp(-lambda)*lambda^n/n! ---------------------------------------

## Probability mass function
Pois_pmf = function(n)
  dpois(n, mean(n))

## Cumulative distribution function
Pois_cdf = function(lambda, nmax = 5e4) {
  n = 0:nmax
  pmf = dpois(n, lambda)
  cdf = zipfR::Rgamma(lambda, n + 1)
  cum = sum(pmf) - cdf + pmf
  return(tibble(x = n, cdf = cdf, cum = cum))
}

## Maximum likelihood estimation of parameter lambda
Pois_lambda <- function(data){
  if(class(data) %in% c('integer', 'numeric')){ 
    return(mean(data))
  } else if(class(data) == 'data.frame'){ 
    return(weighted.mean(data$x, data$freq))
  } else stop('Data is not a vector of integers or a plyr::count data frame')
}

## Fit Poisson, get p-value from discrete Kolmogorov-Smirnov test (implemented via library dgof)
Pois_fit = function(data, max_value = 5e4, method) {
  ## fit parameter lambda
  lambda = Pois_lambda(data)
  
  ## transform data into count data frame
  w = if (is.vector(data)) plyr::count(data) else data
  
  ## calculate ecdf given the data
  dtf = Pois_cdf(lambda = lambda, nmax = max_value)
  myfun = stepfun(w$x, c(0, subset(dtf, x %in% w$x)$cdf))
  
  ## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
  if(method == 'ks') mod = dgof::ks.test(with(w, rep(x, freq)), myfun)
  if(method == 'CVM') mod = cvm_test(with(w, rep(x, freq)), myfun)
  
  ## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
  return(list(
    dtf = dtf,
    parameter = lambda,
    statistic = mod$statistic,
    p.value = mod$p.value
  ))
}


## ----------------- Negative Binomial distribution (NB), P(n) ~ Γ(x+size)/(Γ(size) x!) prob^size (1-prob)^x -------------------------------------
## mu = size * prob/(1 - prob)
NB_fit = function(data, lambda, eta) {
  ## transform data into count data frame
  w = if (is.vector(data))
    plyr::count(data)
  else
    data
  
  ## fit the negative binomial
  nbfit = sads::fitnbinom(data,
                          trunc = NULL,
                          start.value = c(prob = lambda / eta, mu = mean(data)))
  mu = as.numeric(nbfit@details$par['mu'])
  size = as.numeric(nbfit@details$par['size'])
  
  ## calculate ecdf given the data
  dtf = data.frame(x = 0:(2 * max(data)),
                   cdf = pnbinom(0:(2 * max(data)), size = size, mu = mu)) %>% mutate(cum =
                                                                                        1 - cdf)
  myfun = stepfun(w$x, c(0, subset(dtf, x %in% w$x)$cdf))
  
  ## call function ks.test from package dgof, which performs the KS test on discrete data/null distribution
  mod = dgof::ks.test(with(w, rep(x, freq)), myfun)
  
  ## return the fitted parameter, the test statistic (the KS distance) and the p-value of the test
  return(list(
    dtf = dtf,
    parameter = c(size = size, mu = mu),
    statistic = mod$statistic,
    p.value = mod$p.value
  ))
}



## ---------------- Wrapper ----------------
fitSAD = function(data, dbn, method = 'CVM', ...) {
  stopifnot(dbn %in% c('LS', 'IQ', 'Pois', 'NB'))
  if (dbn == 'LS')
    return(LS_fit(data, method = method, ...))
  if (dbn == 'Pois')
    return(Pois_fit(data, method = method, ...))
  if (dbn == 'IQ')
    return(IQ_fit(data, ...))
  if (dbn == 'NB')
    return(NB_fit(data, ...))
}


