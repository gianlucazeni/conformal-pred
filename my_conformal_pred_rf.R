####################################################
############ Conformal prediction ##################
########### author: Gianluca Zeni ##################
####################################################

###########################
# CP with random forests reg as underlying algorithm
###########################

#-------------------------------------
# INDEX
### my_conformal.pred.rf1
#       sigma_i = res_i / (exp(mu.i) + beta)
#       mu.i : knn regression
 

# for the notation, refer to the work of
# Johansson et al (2014)
#     Regression conformal prediction with random forests.
#-------------------------------------


library("FNN")  # knn.reg function
# library("DMwR2")  # kNN function
library("robustHD")   # standardize function

#---------auxiliary-----------

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Compute a coverage interval from a grid of values
grid.interval = function(pts, vals, threshold, right.side=TRUE,
                         method=c("linear","floor","ceil")) {
  
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
#-------------------------------------------------

# # TEST
# conformal.pred(x, loggy, , alpha=0.1,
#                train.fun=lm.funs()$train, predict.fun=lm.funs()$predict)

# y <- loggy
# x <- x[,1:3]
# x0 <- as.matrix(x0[1:3], ncol=3)
# alpha = 0.1
# mad.train.fun = NULL
# mad.predict.fun = NULL
# num.grid.pts = 100
# grid.factor = 1.25
# grid.method = "linear"
# verbose = FALSE
# train.fun <- rf.funs(ntree=500,varfrac=1)$train
# predict.fun <- rf.funs(ntree=500)$pred


# function (x, y, x0, train.fun, predict.fun, alpha = 0.1, mad.train.fun = NULL, 
#           mad.predict.fun = NULL, num.grid.pts = 100, grid.factor = 1.25, 
#           grid.method = c("linear", "floor", "ceiling"), verbose = FALSE) 
# {
#   x = as.matrix(x)
#   y = as.numeric(y)
#   n = nrow(x)
#   p = ncol(x)
#   x0 = matrix(x0, ncol = p)
#   n0 = nrow(x0)
#   check.args(x = x, y = y, x0 = x0, alpha = alpha, train.fun = train.fun, 
#              predict.fun = predict.fun, mad.train.fun = mad.train.fun, 
#              mad.predict.fun = mad.predict.fun)
#   if (length(num.grid.pts) != 1 || !is.numeric(num.grid.pts) || 
#       num.grid.pts <= 1 || num.grid.pts >= 1000 || round(num.grid.pts) != 
#       num.grid.pts) {
#     stop("num.grid.pts must be an integer between 1 and 1000")
#   }
#   check.pos.num(grid.factor)
#   grid.method = match.arg(grid.method)
#   if (verbose == TRUE) 
#     txt = ""
#   if (verbose != TRUE && verbose != FALSE) {
#     txt = verbose
#     verbose = TRUE
#   }
#   if (verbose) 
#     cat(sprintf("%sInitial training on full data set ...\n", 
#                 txt))
#   out = train.fun(x, y)
#   fit = matrix(predict.fun(out, x), nrow = n)
#   pred = matrix(predict.fun(out, x0), nrow = n0)
#   m = ncol(pred)
#   ymax = max(abs(y))
#   yvals = seq(-grid.factor * ymax, grid.factor * ymax, length = num.grid.pts)
#   lo = up = matrix(0, n0, m)
#   pvals = matrix(0, num.grid.pts, m)
#   xx = rbind(x, rep(0, p))
#   for (i in 1:n0) {
#     if (verbose) {
#       cat(sprintf("\r%sProcessing prediction point %i (of %i) ...", 
#                   txt, i, n0))
#       flush.console()
#     }
#     xx[n + 1, ] = x0[i, ]
#     for (j in 1:num.grid.pts) {
#       yy = c(y, yvals[j])
#       if (j == 1) 
#         out = train.fun(xx, yy)
#       else out = train.fun(xx, yy, out)
#       r = abs(yy - matrix(predict.fun(out, xx), nrow = n + 
#                             1))
#       if (!is.null(mad.train.fun) && !is.null(mad.predict.fun)) {
#         for (l in 1:m) {
#           if (j == 1 && l == 1) 
#             out.mad = mad.train.fun(xx, r[, l])
#           else out.mad = mad.train.fun(xx, r[, l], out.mad)
#           r[, l] = r[, l]/mad.predict.fun(out.mad, xx)
#         }
#       }
#       rr = matrix(rep(r[n + 1, ], each = n + 1), nrow = n + 
#                     1)
#       pvals[j, ] = colMeans(rr <= r)
#     }
#     for (l in 1:m) {
#       int = grid.interval(yvals, pvals[, l], alpha, right = T, 
#                           method = grid.method)
#       lo[i, l] = int$lo
#       up[i, l] = int$up
#     }
#   }
#   if (verbose) 
#     cat("\n")
#   return(list(pred = pred, lo = lo, up = up, fit = fit))
# }


#--------- NORMALIZED ONLY ------------

# y <- loggy
# x <- x[,1:3]
# x0 <- as.matrix(x0[1:3], ncol=3)
# alpha = 0.1
# beta=0.1
# num.grid.pts = 100
# grid.factor = 1.25
# grid.method = "linear"
# verbose = FALSE
# train.fun <- my.rf.funs$train
# predict.fun <- my.rf.funs$pred
# regressors=2
# k=3



# sigma_i = res_i / (exp(mu.i) + beta)
# mu.i : knn regression
my_conformal.pred.rf1 <- function (x, y, x0, train.fun, predict.fun, k=3, alpha = 0.1, beta=0,
                                    num.grid.pts = 100, grid.factor = 1.25, 
                                    grid.method = "linear", verbose = FALSE, 
                                    regressors = 1) 
  # input parameter REGRESSORS
  # 1: std
  # 2: in [0,1]
{
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0, ncol = p)
  n0 = nrow(x0)
  out = train.fun(x, y)
  m = 1
  ymax = max(abs(y))
  yvals = seq(-grid.factor * ymax, grid.factor * ymax, length = num.grid.pts)
  lo = up = matrix(0, n0, m)
  pvals = matrix(0, num.grid.pts, m)
  xx = rbind(x, rep(0, p))
  for (i in 1:n0) {
    xx[n + 1, ] = x0[i, ]
    # standardizing regressors
    if(regressors==1)
      xx_std <- apply(xx, 2, standardize)
    if(regressors==2)
      xx_std <- apply(xx, 2, range01)
    for (j in 1:num.grid.pts) {
      yy = c(y, yvals[j])
      r = abs(yy - matrix(predict.fun(out, xx), nrow = n + 1))
      # Normalization
      for (l in 1:m) {
        mu.i <- knn.reg(train=xx_std, y=log(r), k=k)$pred # z_i OUT of the bag
        if(all((mu.i)!=-Inf))
          r[, l] = r[, l]/(exp(mu.i)+beta)
      }
      rr = matrix(rep(r[n + 1, ], each = n + 1), nrow = n + 1)
      pvals[j, ] = colMeans(rr <= r)
    }
    for (l in 1:m) {
      int = grid.interval(yvals, pvals[, l], alpha, right = T, 
                          method = grid.method)
      lo[i, l] = int$lo
      up[i, l] = int$up
    }
  }
  
  return(list(lo = lo, up = up))
}

# # test code
# x <- x
# y <- loggy
# x0 <- as.matrix(x0)
# k=5
# alpha = 0.1
# beta=0.1
# num.grid.pts = 100
# grid.factor = 1.25
# grid.method = "linear"
# verbose = FALSE
# regressors = 2

