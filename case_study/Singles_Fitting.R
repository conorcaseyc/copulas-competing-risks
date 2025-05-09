
------------------------------------------------------
  ####### Fitting the singular distributions #######
------------------------------------------------------
library(fitdistrplus)

# target models
Fit <- c("exp", "weibull", "gamma", "lnorm")

# tolerances 
shape_tol   <- 0.1
aic_pct_tol <- 0.05

# fitting function
fit_dist <- function(x) {
  # fit all models
  fits <- setNames(lapply(fit_candidates, function(d) fitdist(x, d)), 
                   fit_candidates)
  # extract AIC
  aics <- sapply(fits, ⁠ [[ ⁠, "aic")
  best <- names(which.min(aics))
  list(
    best       = best,
    parameters = fits[[best]]$estimate,
    aics       = aics
  )
}



# data
t0 <- na.omit(t_default)

# 2) loglikelihood funciton
ll_weib_mix <- function(par, t, k){
  logits  <- par[1:(k-1)]
  shapes  <- exp(par[ (k):(k + k - 1) ])        # α1…αk
  scales  <- exp(par[(k + k):(k + 2*k - 1)])    #  λ1…λk
  # mixing weights
  pis0<- c(plogis(logits), 1)               # last weight = 1
  pis_hat     <- pis_hat / sum(pis0)
  # component densities
  density <- sapply(1:k, function(j)
    dweibull(t, shape=shapes[j], scale=scales[j])
  )
  # mixture density (floor away zeros)
  f_mix <- pmax(density %*% pis, 10e-6)
  -sum(log(f_mix))
}

 # 3 fit k=1,2,3
fits <- list()
aics <- numeric(3)
for(k in 1:3){
  cat("fitting k =",k,"\n")
 
  if(k==1){
    init <- c(log(1), log(median(t_default)))
    par.names <- c("log(shape1)","log(scale1)")
  } else {
    # for k>1: (k−1) logits + k shapes + k scales
    init <- c(
      rep(0, k-1),                 # logits ~ equal weights
      rep(log(1), k),              # shapes
      log(quantile(t_default, probs= seq(1/(k+1),1-1/(k+1), length.out=k))) # scales
    )
  }
  fit <- optim(
    par    = init,
    fn     = nll_weib_mix,
    t      = t_default,
    k      = k,
    method = "BFGS",
    control= list(maxit=5e3)
  )
  # logLikelihood and AIC
  logLik <- -fit$value
  paramo   <- length(fit$par)
  aics[k] <- 2*paramo - 2*logLik
  fits[[k]] <- list(k=k, fit=fit, logLik=logLik, AIC=aics[k])
}

# 4) compare  models
tab <- data.frame(
  components = 1:3,
  logLik     = sapply(fits, `[[`, "logLik"),
  npar       = sapply(fits, function(x) length(x$fit$par)),
  AIC        = aics
)

print(tab, row.names = FALSE)

# announce best by AIC
cat("\nBest model by AIC:", which.min(aics), "component(s)\n")
