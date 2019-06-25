####################################################
############ Conformal prediction ##################
########### author: Gianluca Zeni ##################
####################################################

###########################
# CP with knn reg as underlying algorithm
###########################

#-------------------------------------
# INDEX
# my_conformal.pred.knn0 : not normalized residuals
# my_conformal.pred.knn1 : sigma_i = res_i / (beta + lambda_i_k)
# my_conformal.pred.knn2 : sigma_i = res_i / exp(beta * lambda_i_k)
# my_conformal.pred.knn3 : sigma_i = res_i / (beta + xi_i_k)
# my_conformal.pred.knn4 : sigma_i = res_i / exp(beta * xi_i_k)
# my_conformal.pred.knn5 : sigma_i = res_i / (beta + lamba_i_k + xi_i_k)
# my_conformal.pred.knn6 : sigma_i = res_i / exp(beta * lamba_i_k) + exp(rho * xi_i_k)

# for the notation, refer to the work of
# Papadopoulos, Vovk, and Gammerman (2011), 
#    Regression conformal prediction with nearest neighbours.
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
# 
# x <- x
# x0 <- as.matrix(x0)
# alpha = 0.1
# mad.train.fun = NULL
# mad.predict.fun = NULL
# num.grid.pts = 100
# grid.factor = 1.25
# grid.method = "linear"
# verbose = FALSE


#--------- NOT NORMALIZED ------------

my_conformal.pred.knn0 <- function (x, y, x0, k=3, alpha = 0.1, 
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
      xx <- apply(xx, 2, standardize)
    if(regressors==2)
      xx <- apply(xx, 2, range01)
    for (j in 1:num.grid.pts) {
      yy = c(y, yvals[j])
      predic <- knn.reg(train=xx, y=yy, k=k) # z_i OUT of the bag
      # predic <- knn.reg(train=xx, test = xx, y=yy, k=k) # z_i IN the bag
      r = abs(yy - matrix(predic$pred))
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

#--------- NORMALIZED ONLY ------------


# sigma_i = res_i / (beta + lambda_i_k)
my_conformal.pred.knn1 <- function (x, y, x0, k=3, alpha = 0.1, beta=0,
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
      xx <- apply(xx, 2, standardize)
    if(regressors==2)
      xx <- apply(xx, 2, range01)
    for (j in 1:num.grid.pts) {
      yy = c(y, yvals[j])
      # predic <- knn.reg(train=xx, y=yy, k=k) # z_i OUT of the bag
      predic <- knn.reg(train=xx, test = xx, y=yy, k=k) # z_i IN the bag
      r = abs(yy - matrix(predic$pred))
      # Normalization
      for (l in 1:m) {
        dneigh <- get.knn(xx,k)  # z_i OUT of the bag (othwise i have k-1)
        # equivalent to 
        # get.knnx(xx[-43,],matrix(xx[43,], ncol=4),k)
        dik <- rowSums(dneigh$nn.dist)
        if(median(dik)!=0)
          lambdaik <- dik / median(dik)
        if(median(dik)==0)
          lambdaik <- 0
        r[, l] = r[, l]/(lambdaik+ beta)
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

# sigma_i = res_i / exp(beta * lambda_i_k)
my_conformal.pred.knn2 <- function (x, y, x0, k=3, alpha = 0.1, beta=0,
                                    num.grid.pts = 100, grid.factor = 1.25, 
                                    grid.method = "linear", verbose = FALSE, 
                                    regressors = 1) 
{
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x0 = matrix(x0, ncol = p)
  n0 = nrow(x0)
  
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
      xx <- apply(xx, 2, standardize)
    if(regressors==2)
      xx <- apply(xx, 2, range01)
    for (j in 1:num.grid.pts) {
      yy = c(y, yvals[j])
      # predic <- knn.reg(train=xx, y=yy, k=k) # z_i OUT of the bag
      predic <- knn.reg(train=xx, test = xx, y=yy, k=k) # z_i IN the bag
      r = abs(yy - matrix(predic$pred))
      # Normalization
      for (l in 1:m) {
        dneigh <- get.knn(xx,k)  # z_i OUT of the bag (othwise i have k-1)
        # equivalent to 
        # get.knnx(xx[-43,],matrix(xx[43,], ncol=4),k)
        dik <- rowSums(dneigh$nn.dist)
        if(median(dik)!=0)
          lambdaik <- dik / median(dik)
        if(median(dik)==0)
          lambdaik <- 0
        r[, l] = r[, l]/exp(lambdaik * beta)
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


# sigma_i = res_i / (beta + xi_i_k)
my_conformal.pred.knn3 <- function (x, y, x0, k=3, alpha = 0.1, beta=0,
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
      xx <- apply(xx, 2, standardize)
    if(regressors==2)
      xx <- apply(xx, 2, range01)
    for (j in 1:num.grid.pts) {
      yy = c(y, yvals[j])
      # predic <- knn.reg(train=xx, y=yy, k=k) # z_i OUT of the bag
      predic <- knn.reg(train=xx, test = xx, y=yy, k=k) # z_i IN the bag
      r = abs(yy - matrix(predic$pred))
      # Normalization
      for (l in 1:m) {
        dneigh <- get.knn(xx,k)  # z_i OUT of the bag (othwise i have k-1)
        # equivalent to 
        # get.knnx(xx[-43,],matrix(xx[43,], ncol=4),k)
        sik <- sqrt(1-1/k)*apply(dneigh$nn.dist,1,sd)
        if(median(sik)!=0)
          xiik <- sik / median(sik)
        if(median(sik)==0)
          xiik <- 0
        r[, l] = r[, l]/(xiik+ beta)
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

# sigma_i = res_i / exp(beta * xi_i_k)
my_conformal.pred.knn4 <- function (x, y, x0, k=3, alpha = 0.1, beta=0,
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
      xx <- apply(xx, 2, standardize)
    if(regressors==2)
      xx <- apply(xx, 2, range01)
    for (j in 1:num.grid.pts) {
      yy = c(y, yvals[j])
      # predic <- knn.reg(train=xx, y=yy, k=k) # z_i OUT of the bag
      predic <- knn.reg(train=xx, test = xx, y=yy, k=k) # z_i IN the bag
      r = abs(yy - matrix(predic$pred))
      # Normalization
      for (l in 1:m) {
        dneigh <- get.knn(xx,k)  # z_i OUT of the bag (othwise i have k-1)
        # equivalent to 
        # get.knnx(xx[-43,],matrix(xx[43,], ncol=4),k)
        sik <- sqrt(1-1/k)*apply(dneigh$nn.dist,1,sd)
        if(median(sik)!=0)
          xiik <- sik / median(sik)
        if(median(sik)==0)
          xiik <- 0
        r[, l] = r[, l]/exp(xiik * beta)
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

# sigma_i = res_i / (beta + lamba_i_k + xi_i_k)
my_conformal.pred.knn5 <- function (x, y, x0, k=3, alpha = 0.1, beta=0,
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
      xx <- apply(xx, 2, standardize)
    if(regressors==2)
      xx <- apply(xx, 2, range01)
    for (j in 1:num.grid.pts) {
      yy = c(y, yvals[j])
      # predic <- knn.reg(train=xx, y=yy, k=k) # z_i OUT of the bag
      predic <- knn.reg(train=xx, test = xx, y=yy, k=k) # z_i IN the bag
      r = abs(yy - matrix(predic$pred))
      # Normalization
      for (l in 1:m) {
        dneigh <- get.knn(xx,k)  # z_i OUT of the bag (othwise i have k-1)
        # equivalent to 
        # get.knnx(xx[-43,],matrix(xx[43,], ncol=4),k)
        
        dik <- rowSums(dneigh$nn.dist)
        if(median(dik)!=0)
          lambdaik <- dik / median(dik)
        if(median(dik)==0)
          lambdaik <- 0
        
        sik <- sqrt(1-1/k)*apply(dneigh$nn.dist,1,sd)
        if(median(sik)!=0)
          xiik <- sik / median(sik)
        if(median(sik)==0)
          xiik <- 0
        r[, l] = r[, l]/(lambdaik + xiik + beta)
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


# sigma_i = res_i / exp(beta * lamba_i_k) + exp(rho * xi_i_k)
my_conformal.pred.knn6 <- function (x, y, x0, k=3, alpha = 0.1, beta=0, rho=0,
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
      xx <- apply(xx, 2, standardize)
    if(regressors==2)
      xx <- apply(xx, 2, range01)
    for (j in 1:num.grid.pts) {
      yy = c(y, yvals[j])
      # predic <- knn.reg(train=xx, y=yy, k=k) # z_i OUT of the bag
      predic <- knn.reg(train=xx, test = xx, y=yy, k=k) # z_i IN the bag
      r = abs(yy - matrix(predic$pred))
      # Normalization
      for (l in 1:m) {
        dneigh <- get.knn(xx,k)  # z_i OUT of the bag (othwise i have k-1)
        # equivalent to 
        # get.knnx(xx[-43,],matrix(xx[43,], ncol=4),k)
        
        dik <- rowSums(dneigh$nn.dist)
        if(median(dik)!=0)
          lambdaik <- dik / median(dik)
        if(median(dik)==0)
          lambdaik <- 0
        
        sik <- sqrt(1-1/k)*apply(dneigh$nn.dist,1,sd)
        if(median(sik)!=0)
          xiik <- sik / median(sik)
        if(median(sik)==0)
          xiik <- 0
        
        r[, l] = r[, l]/(exp(beta*lambdaik) + exp(rho*xiik))
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
