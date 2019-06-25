####################################################
############ Conformal prediction ##################
########### author: Gianluca Zeni ##################
####################################################

# This file contains the function my_conformal.pred.jack
# it is an extension of the basic function in the conformalInference package

# useful when normalized nonconformity measure is chosen

# indeed now it is able to hande the case where:
# there is no uncertainty in the prediction
# and so the residual is 0 (the scores are not divided by 0, but they are already 0)
### my_conformal.pred.jack
# added the check on residuals = 0
### my_conformal.pred.jack3
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


#----------- Main function ----------


my_conformal.pred.jack <- 
  function (x, y, x0, train.fun, predict.fun, alpha = 0.1, special.fun = NULL, 
          mad.train.fun = NULL, mad.predict.fun = NULL, verbose = FALSE) 
{
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0, ncol = p)
  n0 = nrow(x0)
  check.args(x = x, y = y, x0 = x0, alpha = alpha, train.fun = train.fun, 
             predict.fun = predict.fun, mad.train.fun = mad.train.fun, 
             mad.predict.fun = mad.predict.fun, special.fun = special.fun)
  if (verbose == TRUE) 
    txt = ""
  if (verbose != TRUE && verbose != FALSE) {
    txt = verbose
    verbose = TRUE
  }
  if (verbose) 
    cat(sprintf("%sInitial training on full data set ...\n", 
                txt))
  out = train.fun(x, y)
  fit = matrix(predict.fun(out, x), nrow = n)
  pred = matrix(predict.fun(out, x0), nrow = n0)
  m = ncol(pred)
  if (!is.null(special.fun) && (is.null(mad.train.fun) && is.null(mad.predict.fun))) {
    res = abs(special.fun(x, y, out))
  }
  else {
    res = matrix(0, n, m)
    for (i in 1:n) {
      if (verbose) {
        cat(sprintf("\r%sProcessing training point %i (of %i) ...", 
                    txt, i, n))
        flush.console()
      }
      xx = x[-i, ]
      yy = y[-i]
      out.i = train.fun(xx, yy)
      res[i, ] = abs(y[i] - predict.fun(out.i, x[i, , drop = F]))
      if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
        r = abs(y[-i] - matrix(predict.fun(out.i, x[-i, 
                                                    , drop = F]), nrow = n - 1))
        for (l in 1:m) {
          out.mad.i = mad.train.fun(x[-i, , drop = F], 
                                    r[, l])
          if(mad.predict.fun(out.mad.i, x[i, , drop = F])){
          res[i, l] = res[i, l]/mad.predict.fun(out.mad.i, 
                                                x[i, , drop = F])}
        }
      }
    }
    if (verbose) 
      cat("\n")
  }
  lo = up = matrix(0, n0, m)
  for (l in 1:m) {
    if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
      res.train = abs(y - fit[, l])
      mad.out = mad.train.fun(x, res.train)
      mad.x0 = mad.predict.fun(mad.out, x0)
    }
    else {
      mad.x0 = rep(1, n0)
    }
    q = quantile(res[, l], 1 - alpha)
    lo[, l] = pred[, l] - q * mad.x0
    up[, l] = pred[, l] + q * mad.x0
  }
  return(list(pred = pred, lo = lo, up = up, fit = fit))
  }


#----------- Sigma correction ----------

# # Test
# x <- x
# y <- loggy
# x0 <- as.matrix(x0)
# alpha=0.1
# beta=mybeta
# train.fun=lm.funs()$train
# predict.fun=lm.funs()$predict
# mad.train.fun=lm.funs()$train
# mad.predict.fun=lm.funs()$predict
# special.fun = NULL
# verbose = FALSE


my_conformal.pred.jack3 <- 
  function (x, y, x0, train.fun, predict.fun, alpha = 0.1, beta=0, special.fun = NULL, 
            mad.train.fun = NULL, mad.predict.fun = NULL, verbose = FALSE) 
  {
    x = as.matrix(x)
    y = as.numeric(y)
    n = nrow(x)
    p = ncol(x)
    x0 = matrix(x0, ncol = p)
    n0 = nrow(x0)
    check.args(x = x, y = y, x0 = x0, alpha = alpha, train.fun = train.fun, 
               predict.fun = predict.fun, mad.train.fun = mad.train.fun, 
               mad.predict.fun = mad.predict.fun, special.fun = special.fun)
    if (verbose == TRUE) 
      txt = ""
    if (verbose != TRUE && verbose != FALSE) {
      txt = verbose
      verbose = TRUE
    }
    if (verbose) 
      cat(sprintf("%sInitial training on full data set ...\n", 
                  txt))
    out = train.fun(x, y)
    fit = matrix(predict.fun(out, x), nrow = n)
    pred = matrix(predict.fun(out, x0), nrow = n0)
    m = ncol(pred)
    if (!is.null(special.fun) && (is.null(mad.train.fun) && is.null(mad.predict.fun))) {
      res = abs(special.fun(x, y, out))
    }
    else {
      res = matrix(0, n, m)
      for (i in 1:n) {
        if (verbose) {
          cat(sprintf("\r%sProcessing training point %i (of %i) ...", 
                      txt, i, n))
          flush.console()
        }
        xx = x[-i, ]
        yy = y[-i]
        out.i = train.fun(xx, yy)
        res[i, ] = abs(y[i] - predict.fun(out.i, x[i, , drop = F]))
        if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
          r = abs(y[-i] - matrix(predict.fun(out.i, x[-i, 
                                                      , drop = F]), nrow = n - 1))
          for (l in 1:m) {
            out.mad.i = mad.train.fun(x[-i, , drop = F], 
                                      log(r[, l]))
            if(!is.nan(mad.predict.fun(out.mad.i, x[i, , drop = F]))){
              res[i, l] = res[i, l]/( exp(mad.predict.fun(out.mad.i, 
                                                    x[i, , drop = F])) + beta)
            }
          }
        }
      }
      if (verbose) 
        cat("\n")
    }
    lo = up = matrix(0, n0, m)
    for (l in 1:m) {
      if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
        res.train = abs(y - fit[, l])
        mad.out = mad.train.fun(x, log(res.train))
        mad.x0 = exp(mad.predict.fun(mad.out, x0))+beta
        if(is.nan(mad.x0)) mad.x0 <- 0
      }
      else {
        mad.x0 = rep(1, n0)
      }
      q = quantile(res[, l], 1 - alpha)
      lo[, l] = pred[, l] - q * mad.x0
      up[, l] = pred[, l] + q * mad.x0
    }
    return(list(pred = pred, lo = lo, up = up, fit = fit))
  }