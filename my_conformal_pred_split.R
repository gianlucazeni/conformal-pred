####################################################
############ Conformal prediction ##################
########### author: Gianluca Zeni ##################
####################################################

# This file contains the function my_conformal.pred.split
# it is an extension of the basic function in the conformalInference package

# useful when normalized nonconformity measure is chosen

# indeed now it is able to hande the case where there is no uncertainty in the prediction
# and so the residual is 0 (the scores are not divided by 0, but they are already 0)
### my_conformal.pred.split
# added the check on residuals = 0
### my_conformal.pred.split3
# -> sigma_i = exp(mu_i) + beta (and the check on residuals)

#----------- Auxiliary functions ----------


check.args = function(x=NULL, y=NULL, x0=NULL, alpha=NULL,
                      train.fun=NULL, predict.fun=NULL, mad.train.fun=NULL, mad.predict.fun=NULL,
                      special.fun=NULL) {
  
  if (is.null(x) || !is.numeric(x)) stop("x must be a numeric matrix")
  if (is.null(y) || !is.numeric(y)) stop("y must be a numeric vector")
  if (nrow(x) != length(y)) stop("nrow(x) and length(y) must match")
  if (is.null(x0) || !is.numeric(x0)) stop("x0 must be a numeric matrix")
  if (ncol(x) != ncol(x0)) stop("ncol(x) and ncol(x0) must match")
  check.num.01(alpha)
  if (is.null(train.fun) || !is.function(train.fun))
    stop("train.fun must be a function")
  if (is.null(predict.fun) || !is.function(predict.fun))
    stop("predict.fun must be a function")
  if (!is.null(mad.train.fun) && !is.function(mad.train.fun)) 
    stop("mad.train.fun must be a function")
  if (!is.null(mad.predict.fun) && !is.function(mad.predict.fun)) 
    stop("mad.predict.fun must be a function")
  if ((!is.null(mad.train.fun) && is.null(mad.predict.fun)) ||
      (is.null(mad.train.fun) && !is.null(mad.predict.fun)))
    stop("mad.train.fun and mad.predict.fun must both be provided")
  if (!is.null(special.fun) && !is.function(special.fun)) 
    stop("special.fun must be a function")
}

check.bool = function(b) {
  if (is.null(b) || length(b)!=1 || !is.logical(b))
    stop(paste(deparse(substitute(b)),"must be a Boolean"))
}

check.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a))
    stop(paste(deparse(substitute(a)),"must be a number"))
}

check.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i)
    stop(paste(deparse(substitute(i)),"must be an integer"))
}

check.pos.num = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0)
    stop(paste(deparse(substitute(a)),"must be a positive number"))
}

check.pos.int = function(i) {
  if (is.null(i) || length(i)!= 1 || !is.numeric(i) || round(i) != i || i<1)
    stop(paste(deparse(substitute(i)),"must be a positive integer"))
}

check.num.01 = function(a) {
  if (is.null(a) || length(a)!= 1 || !is.numeric(a) || a<0 || a>1)
    stop(paste(deparse(substitute(a)),"must be a number between 0 and 1"))
}


# Compute threshold for split conformal, at level alpha.
# Residuals must be sorted (and positive)
conformal.quantile = function(res, alpha) {
  n = length(res)
  if (ceiling((n+1)*alpha) <= 1) { return(Inf) }
  return(res[ceiling((n+1)*(1-alpha))])
}

#----------- Main function ----------


my_conformal.pred.split <-
  function (x, y, x0, train.fun, predict.fun, alpha = 0.1, mad.train.fun = NULL, 
            mad.predict.fun = NULL, split = NULL, seed = NULL, verbose = FALSE) 
{
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0, ncol = p)
  n0 = nrow(x0)
  check.args(x = x, y = y, x0 = x0, alpha = alpha, train.fun = train.fun, 
             predict.fun = predict.fun, mad.train.fun = mad.train.fun, 
             mad.predict.fun = mad.predict.fun)
  if (verbose == TRUE) 
    txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  if (!is.null(split)) 
    i1 = split
  else {
    if (!is.null(seed)) 
      set.seed(seed)
    i1 = sample(1:n, floor(n/2))
  }
  i2 = (1:n)[-i1]
  n1 = length(i1)
  n2 = length(i2)
  if (verbose) {
    cat(sprintf("%sSplitting data into parts of size %i and %i ...\n", 
                txt, n1, n2))
    cat(sprintf("%sTraining on first part ...\n", txt))
  }
  out = train.fun(x[i1, , drop = F], y[i1])
  fit = matrix(predict.fun(out, x), nrow = n)
  pred = matrix(predict.fun(out, x0), nrow = n0)
  m = ncol(pred)
  if (verbose) {
    cat(sprintf("%sComputing residuals and quantiles on second part ...\n", 
                txt))
  }
  res = abs(y[i2] - matrix(predict.fun(out, x[i2, , drop = F]), 
                           nrow = n2))
  lo = up = matrix(0, n0, m)
  for (l in 1:m) {
    if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
      res.train = abs(y[i1] - fit[i1, l])
      mad.out = mad.train.fun(x[i1, , drop = F], res.train)
      mad.x2 = mad.predict.fun(mad.out, x[i2, , drop = F])
      mad.x0 = mad.predict.fun(mad.out, x0)
      if(all(mad.x2!=0)){
        res[, l] = res[, l]/mad.x2
      }
    }
    else {
      mad.x0 = rep(1, n0)
    }
    q = conformal.quantile(sort(res[, l]), alpha)
    lo[, l] = pred[, l] - q * mad.x0
    up[, l] = pred[, l] + q * mad.x0
  }
  return(list(pred = pred, lo = lo, up = up, fit = fit, split = i1))
  }


#----------- Sigma correction ----------

# beta version - bad results for small n?
# recommended to test it on large dataset to check if e.thing is ok
my_conformal.pred.split3 <-
  function (x, y, x0, train.fun, predict.fun, alpha = 0.1, beta=0, mad.train.fun = NULL, 
            mad.predict.fun = NULL, split = NULL, seed = NULL, verbose = FALSE) 
  {
    x = as.matrix(x)
    y = as.numeric(y)
    n = nrow(x)
    p = ncol(x)
    x0 = matrix(x0, ncol = p)
    n0 = nrow(x0)
    check.args(x = x, y = y, x0 = x0, alpha = alpha, train.fun = train.fun, 
               predict.fun = predict.fun, mad.train.fun = mad.train.fun, 
               mad.predict.fun = mad.predict.fun)
    if (verbose == TRUE) 
      txt = ""
    if (verbose != TRUE && verbose != FALSE) {
      txt = verbose
      verbose = TRUE
    }
    if (!is.null(split)) 
      i1 = split
    else {
      if (!is.null(seed)) 
        set.seed(seed)
      i1 = sample(1:n, floor(n/2))
    }
    i2 = (1:n)[-i1]
    n1 = length(i1)
    n2 = length(i2)
    if (verbose) {
      cat(sprintf("%sSplitting data into parts of size %i and %i ...\n", 
                  txt, n1, n2))
      cat(sprintf("%sTraining on first part ...\n", txt))
    }
    out = train.fun(x[i1, , drop = F], y[i1])
    fit = matrix(predict.fun(out, x), nrow = n)
    pred = matrix(predict.fun(out, x0), nrow = n0)
    m = ncol(pred)
    if (verbose) {
      cat(sprintf("%sComputing residuals and quantiles on second part ...\n", 
                  txt))
    }
    res = abs(y[i2] - matrix(predict.fun(out, x[i2, , drop = F]), 
                             nrow = n2))
    lo = up = matrix(0, n0, m)
    for (l in 1:m) {
      if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
        res.train = abs(y[i1] - fit[i1, l])
        mad.out = mad.train.fun(x[i1, , drop = F], log(res.train))
        mad.x2 = exp(mad.predict.fun(mad.out, x[i2, , drop = F]))+beta
        mad.x0 = exp(mad.predict.fun(mad.out, x0))+beta
        if(all(!is.nan(mad.x2))){
          res[, l] = res[, l]/mad.x2
        }
        if(is.nan(mad.x0)) mad.x0 <- 0
      }
      else {
        mad.x0 = rep(1, n0)
      }
      q = conformal.quantile(sort(res[, l]), alpha)
      lo[, l] = pred[, l] - q * mad.x0
      up[, l] = pred[, l] + q * mad.x0
    }
    return(list(pred = pred, lo = lo, up = up, fit = fit, split = i1))
  }