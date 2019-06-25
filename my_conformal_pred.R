####################################################
############ Conformal prediction ##################
########### author: Gianluca Zeni ##################
####################################################

# This file contains the function my_conformal.pred
# it is an extension of the basic function in the conformalInference package

# useful when normalized nonconformity measure is chosen

# indeed now it is able to hande the case where: 
### my_conformal.pred
# -> sigma_i = something_i + beta 
# (beta defining the relative weigth of normalization)
### my_conformal.pred2
# -> sigma_i = exp(mu_i), mu_i estimate of ln(abs(y_i - yhat_i)) 
### my_conformal.pred3
# -> sigma_i = exp(mu_i) + beta



#----------- Auxiliary functions ----------

# Compute a coverage interval from a grid of values
grid.interval = function(pts, vals, threshold, right.side=TRUE,
                         method=c("linear","floor","ceil")) {
  
  method = match.arg(method)
  n = length(pts)
  
  above = vals >= threshold
  if (right.side) ii = which(above)
  else ii = which(!above)
  
  # If no vals are on the correct side of the threshold, then
  # make up a trivial interval
  if (length(ii)==0) return(list(lo=0,up=0))
  
  i1 = min(ii)
  i2 = max(ii)
  
  if (method=="floor") {
    return(list(lo=pts[i1],up=pts[i2]))
  }
  else if (method=="ceil") {
    return(list(lo=pts[max(i1-1,1)],up=pts[min(i2+1,n)]))
  }
  else {
    if (i1==1) xl = pts[i1]
    else {
      slope = (vals[i1]-vals[i1-1])/(pts[i1]-pts[i1-1])
      xl = pts[i1-1] + (threshold-vals[i1-1])/slope
    }
    
    if (i2==n) xr = pts[i2]
    else {
      slope = (vals[i2+1]-vals[i2])/(pts[i2+1]-pts[i2])
      xr = pts[i2] + (threshold-vals[i2])/slope
    }
    
    return(list(lo=xl,up=xr))
  }
}

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

# # TESTING CODE
# tstep <- 78 # 78 27
# y <- log(1+y[,tstep])
# 
# x_dtemp <- f_delta_temp_mean[,tstep]
# x_rain <- f_rain[,tstep]
# s_we_x_rain <- s_we * x_rain
# x <- cbind(s_we, x_rain, x_dtemp, s_we_x_rain)
# x0 <- as.matrix(data.frame(s_we=my_we, x_rain=my_rain[tstep], x_dtemp=my_dtemp[tstep], s_we_x_rain=my_we_rain[tstep]))
# alpha=0.1
# train.fun=lm.funs()$train
# predict.fun=lm.funs()$predict
# mad.train.fun=lm.funs()$train
# mad.predict.fun=lm.funs()$predict
# num.grid.pts = 100
# grid.factor = 1.25
# grid.method = "linear"
# verbose = FALSE
# # beta = 0.0075
# a <- which(r[, l]/mad.predict.fun(out.mad, xx) >= (r[, l]/mad.predict.fun(out.mad, xx))[43,])
# r[a,]
# mad.predict.fun(out.mad, xx)[a]
# (r[, l]/mad.predict.fun(out.mad, xx))[a]
# (r[, l]/(mad.predict.fun(out.mad, xx)+0.005))[a]

my_conformal.pred <- 
  function (x, y, x0, train.fun, predict.fun, alpha = 0.1, mad.train.fun = NULL, 
            mad.predict.fun = NULL, beta = 0, num.grid.pts = 100, grid.factor = 1.25, 
            grid.method = c("linear", "floor", "ceiling"), verbose = FALSE) 
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
    if (length(num.grid.pts) != 1 || !is.numeric(num.grid.pts) || 
        num.grid.pts <= 1 || num.grid.pts >= 1000 || round(num.grid.pts) != 
        num.grid.pts) {
      stop("num.grid.pts must be an integer between 1 and 1000")
    }
    check.pos.num(grid.factor)
    grid.method = match.arg(grid.method)
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
    ymax = max(abs(y))
    yvals = seq(-grid.factor * ymax, grid.factor * ymax, length = num.grid.pts)
    lo = up = matrix(0, n0, m)
    pvals = matrix(0, num.grid.pts, m)
    xx = rbind(x, rep(0, p))
    for (i in 1:n0) {
      if (verbose) {
        cat(sprintf("\r%sProcessing prediction point %i (of %i) ...", 
                    txt, i, n0))
        flush.console()
      }
      xx[n + 1, ] = x0[i, ]
      for (j in 1:num.grid.pts) {
        yy = c(y, yvals[j])
        if (j == 1) 
          out = train.fun(xx, yy)
        else out = train.fun(xx, yy, out)
        r = abs(yy - matrix(predict.fun(out, xx), nrow = n + 
                              1))
        if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
          for (l in 1:m) {
            if (j == 1 && l == 1) 
              out.mad = mad.train.fun(xx, r[, l])
            else out.mad = mad.train.fun(xx, r[, l], out.mad)
            r[, l] = r[, l]/(mad.predict.fun(out.mad, xx)+beta)
          }
        }
        rr = matrix(rep(r[n + 1, ], each = n + 1), nrow = n + 
                      1)
        pvals[j, ] = colMeans(rr <= r)
      }
      for (l in 1:m) {
        int = grid.interval(yvals, pvals[, l], alpha, right = T, 
                            method = grid.method)
        lo[i, l] = int$lo
        up[i, l] = int$up
      }
    }
    if (verbose) 
      cat("\n")
    return(list(pred = pred, lo = lo, up = up, fit = fit))
  }


#-----------------------------------------

my_conformal.pred2 <- 
  function (x, y, x0, train.fun, predict.fun, alpha = 0.1, mad.train.fun = NULL, 
            mad.predict.fun = NULL, num.grid.pts = 100, grid.factor = 1.25, 
            grid.method = c("linear", "floor", "ceiling"), verbose = FALSE) 
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
    if (length(num.grid.pts) != 1 || !is.numeric(num.grid.pts) || 
        num.grid.pts <= 1 || num.grid.pts >= 1000 || round(num.grid.pts) != 
        num.grid.pts) {
      stop("num.grid.pts must be an integer between 1 and 1000")
    }
    check.pos.num(grid.factor)
    grid.method = match.arg(grid.method)
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
    ymax = max(abs(y))
    yvals = seq(-grid.factor * ymax, grid.factor * ymax, length = num.grid.pts)
    lo = up = matrix(0, n0, m)
    pvals = matrix(0, num.grid.pts, m)
    xx = rbind(x, rep(0, p))
    for (i in 1:n0) {
      if (verbose) {
        cat(sprintf("\r%sProcessing prediction point %i (of %i) ...", 
                    txt, i, n0))
        flush.console()
      }
      xx[n + 1, ] = x0[i, ]
      for (j in 1:num.grid.pts) {
        yy = c(y, yvals[j])
        if (j == 1) 
          out = train.fun(xx, yy)
        else out = train.fun(xx, yy, out)
        r = abs(yy - matrix(predict.fun(out, xx), nrow = n + 
                              1))
        if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
          for (l in 1:m) {
            if (j == 1 && l == 1) 
              out.mad = mad.train.fun(xx, log(r[, l]))
            else out.mad = mad.train.fun(xx, log(r[, l]), out.mad)
            r[, l] = r[, l]/exp(mad.predict.fun(out.mad, xx))
          }
        }
        rr = matrix(rep(r[n + 1, ], each = n + 1), nrow = n + 
                      1)
        pvals[j, ] = colMeans(rr <= r)
      }
      for (l in 1:m) {
        int = grid.interval(yvals, pvals[, l], alpha, right = T, 
                            method = grid.method)
        lo[i, l] = int$lo
        up[i, l] = int$up
      }
    }
    if (verbose) 
      cat("\n")
    return(list(pred = pred, lo = lo, up = up, fit = fit))
  }


#-----------------------------------------

my_conformal.pred3 <- 
  function (x, y, x0, train.fun, predict.fun, alpha = 0.1, mad.train.fun = NULL, 
            mad.predict.fun = NULL, beta = 0, num.grid.pts = 100, grid.factor = 1.25, 
            grid.method = c("linear", "floor", "ceiling"), verbose = FALSE) 
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
    if (length(num.grid.pts) != 1 || !is.numeric(num.grid.pts) || 
        num.grid.pts <= 1 || num.grid.pts >= 1000 || round(num.grid.pts) != 
        num.grid.pts) {
      stop("num.grid.pts must be an integer between 1 and 1000")
    }
    check.pos.num(grid.factor)
    grid.method = match.arg(grid.method)
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
    ymax = max(abs(y))
    yvals = seq(-grid.factor * ymax, grid.factor * ymax, length = num.grid.pts)
    lo = up = matrix(0, n0, m)
    pvals = matrix(0, num.grid.pts, m)
    xx = rbind(x, rep(0, p))
    for (i in 1:n0) {
      if (verbose) {
        cat(sprintf("\r%sProcessing prediction point %i (of %i) ...", 
                    txt, i, n0))
        flush.console()
      }
      xx[n + 1, ] = x0[i, ]
      for (j in 1:num.grid.pts) {
        yy = c(y, yvals[j])
        if (j == 1) 
          out = train.fun(xx, yy)
        else out = train.fun(xx, yy, out)
        r = abs(yy - matrix(predict.fun(out, xx), nrow = n + 
                              1))
        if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
          for (l in 1:m) {
            if (j == 1 && l == 1) 
              out.mad = mad.train.fun(xx, log(r[, l]))
            else out.mad = mad.train.fun(xx, log(r[, l]), out.mad)
            r[, l] = r[, l]/(exp(mad.predict.fun(out.mad, xx))+beta)
          }
        }
        rr = matrix(rep(r[n + 1, ], each = n + 1), nrow = n + 
                      1)
        pvals[j, ] = colMeans(rr <= r)
      }
      for (l in 1:m) {
        int = grid.interval(yvals, pvals[, l], alpha, right = T, 
                            method = grid.method)
        lo[i, l] = int$lo
        up[i, l] = int$up
      }
    }
    if (verbose) 
      cat("\n")
    return(list(pred = pred, lo = lo, up = up, fit = fit))
  }
  
