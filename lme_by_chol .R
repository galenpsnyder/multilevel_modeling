# ## model module -- prepares user input data for optimization-----------------------------------
# lme.model.mod <- function(formula, random, id, dat, nointFE, nointRE, REML, getcorr)
# {
#   
#   ## initializing and formatting function call
#   df <- model.frame(formula, data = dat)
#   N <- nrow(df)
#   
#   ## saving names of fixed and random effects for later
#   namesFE <- names(df)[-1]
#   namesRE <- random
#   df <- as.matrix(df)
#   obs.names <- levels(factor(dat[,id]))
#   
#   ## formatting response vector
#   y <- df[,1]
#   attr(y, "names") <- NULL
#   y <- cbind(y, dat[,id])
#   split.y <- list()
#   for(g in 1:length(unique(y[, 2]))){
#     split.y[[g]] <- y[, 1][which(y[, 2] == g)]
#   }
#   y <- split.y
#   
#   ## getting cluster size, "n", and number of clusters, "j"
#   n <- sapply(y, length)
#   j <- length(y)
#   
#   ## formatting FE model matrix
#   if(ncol(df) == 1){
#     x <- matrix(1, nrow = N, ncol = 1)
#     p <- 1
#     x <- cbind(x, dat[, id])
#     namesFE <- "(intercept)"
#     split.x <- list()
#     for(g in 1:length(unique(x[, ncol(x)]))){
#       x.init <- c()
#       for(m in 1:(ncol(x) - 1)){
#         x.init <- cbind(x.init, x[, m][which(x[, ncol(x)] == g)])
#       }
#       split.x[[g]] <- x.init
#     }
#     x <- split.x
#   } else {
#     x <- df[,-1]
#     attr(x, "names") <- NULL
#     p <- ncol(x)
#     ifelse(is.null(p), p <- 1, p <- p)
#     x <- cbind(x, dat[, id])
#     split.x <- list()
#     for(g in 1:length(unique(x[, ncol(x)]))){
#       x.init <- c()
#       for(m in 1:(ncol(x) - 1)){
#         x.init <- cbind(x.init, x[, m][which(x[, ncol(x)] == g)])
#       }
#       split.x[[g]] <- x.init
#     }
#     x <- split.x
#     
#     # some front end work handling inclusion of intercepts for fixed effects
#     if(nointFE == FALSE){
#       x <- lapply(1:j, function(k) cbind(1, x[[k]]))
#       p <- p + 1
#       namesFE <- c("(intercept)", namesFE)
#     } 
#   }
#   
#   ## formatting RE model matrix
#   if(random == 1){
#     q <- 1
#     z <- matrix(1, nrow = N, ncol = 1)
#     namesRE <- "(intercept)"
#     z <- cbind(z, dat[, id])
#     split.z <- list()
#     for(g in 1:length(unique(z[, ncol(z)]))){
#       z.init <- c()
#       for(m in 1:(ncol(z) - 1)){
#         z.init <- cbind(z.init, z[, m][which(z[, ncol(z)] == g)])
#       }
#       split.z[[g]] <- z.init
#     }
#     z <- split.z
#   } else {
#     q <- length(random)
#     z <- as.matrix(dat[random])
#     z <- cbind(z, dat[, id])
#     split.z <- list()
#     for(g in 1:length(unique(z[, ncol(z)]))){
#       z.init <- c()
#       for(m in 1:(ncol(z) - 1)){
#         z.init <- cbind(z.init, z[, m][which(z[, ncol(z)] == g)])
#       }
#       split.z[[g]] <- z.init
#     }
#     z <- split.z
#     
#     # some front end work handling inclusion of intercepts for random effects
#     if(nointRE == FALSE){
#       z <- lapply(1:j, function(k) cbind(1, z[[k]]))
#       q <- q + 1
#       namesRE <- c("(intercept)", namesRE)
#     } 
#   }
#   
#   # forming necessary crossprods for objective function
#   xtx <- lapply(1:j, function(k) crossprod(x[[k]]))
#   xty <- lapply(1:j, function(k) crossprod(x[[k]], y[[k]]))
#   xtz <- lapply(1:j, function(k) crossprod(x[[k]], z[[k]]))
#   ztx <- lapply(1:j, function(k) crossprod(z[[k]], x[[k]]))
#   zty <- lapply(1:j, function(k) crossprod(z[[k]], y[[k]]))
#   zt  <- lapply(1:j, function(k) t(z[[k]]))
#   
#   
#   # diff <- 1
#   # obsloglik <- -Inf
#   # ll <- NULL
#   # iter <- 1
#   # restarts <- 0
#   
#   # creating objects that will be used in objective function
#   lambdat     <- matrix(0, q, q)
#   lind        <- which(upper.tri(lambdat, diag = T), arr.ind = T)
#   theta.start <- diag(1, q)[lind]
#   upper       <- rep(Inf, q*(q+1)/2)
#   lower       <- ifelse(theta.start, 0, -Inf)
#   
#   mod.out <- list(
#     df = df, N = N, namesFE = namesFE, namesRE = namesRE, obs.names = obs.names,
#     y = y, n = n, j = j, x = x, p = p, z = z, q = q, xtx = xtx, xty = xty, xtz = xtz, ztx = ztx,
#     zty = zty, zt = zt, lambdat = lambdat, lind = lind, theta.start = theta.start, 
#     upper = upper, lower = lower, REML = REML, getcorr = getcorr
#   )
#   mod.out
# }



## formula parsing module -- extracts variables and preps for model matrix creation
parse.lme.formula <- function(formula){
  # coercing user input to formula
  ff <- formula(formula)
  
  # extracting y component
  y <- as.character(ff[2])
  
  # extracting rhs
  rhs <- as.character(ff[3])
  
  # splitting formula by fixed and random components--also detects model type
  split.form <- unlist(strsplit(rhs, split = "[(]"))
  ifelse(length(split.form) == 1, lme <- FALSE, lme <- TRUE)
  
  # fixed effects part
  fe <- split.form[1]
  fe <- unlist(strsplit(fe, "[+]"))
  fe <- fe[fe != " "]
  fe <- gsub("\\s", "", fe)
  ifelse("-1" %in% fe || "0" %in% fe, nointFE <- TRUE, nointFE <- FALSE)
  
  if(lme){
    # random effects part
    re.id <- split.form[2]
    re.id <- unlist(strsplit(re.id, split = "[|]"))
    re <- re.id[1]
    re <- unlist(strsplit(re, split = "[+]"))
    re <- gsub("\\s", "", re)
    ifelse("1" %in% re, nointRE <- FALSE, nointRE <- TRUE)
    
    # id component
    id <- re.id[2]
    id <- gsub(")", "", id)
    id <- gsub(" ", "", id)
    
    out <- list(y = y, fe = fe, re = re, id = id, nointFE = nointFE, nointRE = nointRE, lme = lme, lme.form = ff)
  } else {
    out <- list(y = y, fe = fe, re = NULL, id = NULL, nointFE = nointFE, nointRE = NULL, lme = lme, lme.form = ff)
  }
  out
}



## model module -- prepares user input data for optimization-----------------------------------
lme.model.mod <- function(y, fixed, random, id, dat, nointFE, nointRE, REML, getcorr)
{
  
  # conditional formatting
  # remove variable names that are not present in dataframe
  # if fixed effects contains intercept only, assume empty model
  # if random effects contains intercept only, assume random intercepts only model
  # then remove terms referring to intercept so model frame can be created
  miss <- fixed %in% c(names(dat), "1")
  if(any(miss == F)){
    warning(gettext("in fixed effects specification -- the following variables dropped because they are not present in dat:\n", paste(fixed[!miss], collapse = ", ")), call. = F, immediate. = T)
    fixed <- fixed[miss]
  }
  
  miss <- random %in% c(names(dat), "1")
  if(any(miss == F)){
    warning(gettext("in random effects specification -- the following variables dropped because they are not present in dat:\n", paste(random[!miss], collapse = ", ")), call. = F, immediate. = T)
    random <- random[miss]
  }
  
  ifelse("1" %in% fixed & length(fixed) == 1, empty <- TRUE, empty <- FALSE)
  ifelse("1" %in% random & length(random) == 1, ri <- TRUE, ri <- FALSE)
  fixed <- fixed[fixed != "-1" & fixed != "0" & fixed != "1"]
  random <- random[random != "1"]
  
  if(empty == F & length(fixed) == 0) stop(gettext("model cannot be estimated -- no variables specified")) 
  
  ## saving names of fixed and random effects for later
  namesFE <- fixed
  namesRE <- random
  obs.names <- levels(factor(dat[,id]))
  
  ## formatting response vector
  y <- dat[,y]
  N <- length(y)
  y <- cbind(y, dat[,id])
  split.y <- list()
  for(g in 1:length(unique(y[, 2]))){
    split.y[[g]] <- y[, 1][which(y[, 2] == g)]
  }
  y <- split.y
  
  ## getting cluster size, "n", and number of clusters, "j"
  n <- sapply(y, length)
  j <- length(y)
  
  ## formatting FE model matrix
  if(empty){
    x <- matrix(1, nrow = N, ncol = 1)
    p <- 1
    x <- cbind(x, dat[, id])
    namesFE <- "(intercept)"
    split.x <- list()
    for(g in 1:length(unique(x[, ncol(x)]))){
      x.init <- c()
      for(m in 1:(ncol(x) - 1)){
        x.init <- cbind(x.init, x[, m][which(x[, ncol(x)] == g)])
      }
      split.x[[g]] <- x.init
    }
    x <- split.x
  } else {
    x <- dat[, fixed]
    p <- ncol(x)
    ifelse(is.null(p), p <- 1, p <- p)
    x <- cbind(x, dat[, id])
    split.x <- list()
    for(g in 1:length(unique(x[, ncol(x)]))){
      x.init <- c()
      for(m in 1:(ncol(x) - 1)){
        x.init <- cbind(x.init, x[, m][which(x[, ncol(x)] == g)])
      }
      split.x[[g]] <- x.init
    }
    x <- split.x
    
    # some front end work handling inclusion of intercepts for fixed effects
    if(nointFE == FALSE){
      x <- lapply(1:j, function(k) cbind(1, x[[k]]))
      p <- p + 1
      namesFE <- c("(intercept)", namesFE)
    }
  }
  
  ## formatting RE model matrix
  if(ri){
    q <- 1
    z <- matrix(1, nrow = N, ncol = 1)
    z <- cbind(z, dat[, id])
    namesRE <- "(intercept)"
    split.z <- list()
    for(g in 1:length(unique(z[, ncol(z)]))){
      z.init <- c()
      for(m in 1:(ncol(z) - 1)){
        z.init <- cbind(z.init, z[, m][which(z[, ncol(z)] == g)])
      }
      split.z[[g]] <- z.init
    }
    z <- split.z
  } else {
    q <- length(random)
    z <- as.matrix(dat[, random])
    z <- cbind(z, dat[, id])
    split.z <- list()
    for(g in 1:length(unique(z[, ncol(z)]))){
      z.init <- c()
      for(m in 1:(ncol(z) - 1)){
        z.init <- cbind(z.init, z[, m][which(z[, ncol(z)] == g)])
      }
      split.z[[g]] <- z.init
    }
    z <- split.z
    
    # some front end work handling inclusion of intercepts for random effects
    if(nointRE == FALSE){
      z <- lapply(1:j, function(k) cbind(1, z[[k]]))
      q <- q + 1
      namesRE <- c("(intercept)", namesRE)
    }
  }
  
  # forming necessary crossprods for objective function
  xtx <- lapply(1:j, function(k) crossprod(x[[k]]))
  xty <- lapply(1:j, function(k) crossprod(x[[k]], y[[k]]))
  xtz <- lapply(1:j, function(k) crossprod(x[[k]], z[[k]]))
  ztx <- lapply(1:j, function(k) crossprod(z[[k]], x[[k]]))
  zty <- lapply(1:j, function(k) crossprod(z[[k]], y[[k]]))
  zt  <- lapply(1:j, function(k) t(z[[k]]))
  
  # creating objects that will be used in objective function
  lambdat     <- matrix(0, q, q)
  lind        <- which(upper.tri(lambdat, diag = T), arr.ind = T)
  theta.start <- diag(1, q)[lind]
  upper       <- rep(Inf, q*(q+1)/2)
  lower       <- ifelse(theta.start, 0, -Inf)
  
  mod.out <- list(
    df = df, N = N, namesFE = namesFE, namesRE = namesRE, obs.names = obs.names,
    y = y, n = n, j = j, x = x, p = p, z = z, q = q, xtx = xtx, xty = xty, xtz = xtz, ztx = ztx,
    zty = zty, zt = zt, lambdat = lambdat, lind = lind, theta.start = theta.start,
    upper = upper, lower = lower, REML = REML, getcorr = getcorr
  )
  mod.out
}



## optimization module: deviance function (objective function)---------------------------------
pls <- function(mod.out){
  lambdat <- mod.out$lambdat
  lind    <- mod.out$lind
  x       <- mod.out$x
  y       <- mod.out$y
  N       <- mod.out$N
  q       <- mod.out$q
  p       <- mod.out$p
  j       <- mod.out$j
  xtx     <- mod.out$xtx
  xty     <- mod.out$xty
  xtz     <- mod.out$xtz
  ztx     <- mod.out$ztx
  zty     <- mod.out$zty
  zt      <- mod.out$zt
  REML    <- mod.out$REML
  
  function(theta){
    # updating upper triangular elements of relative covariance factor
    lambdat[lind] <- theta 
    # forming lower cholesky factor of variance components
    L    <- lapply(1:j, function(k) t(chol(tcrossprod(lambdat %*% zt[[k]]) + diag(1, q))))
    # taking 
    rzx  <- lapply(1:j, function(k) solve(L[[k]], lambdat %*% ztx[[k]]))
    
    rx   <- matrix(0, p, p)
    rx.vec <- sapply(1:j, function(i) as.vector(xtx[[i]] - crossprod(rzx[[i]])))
    if(is.null(dim(rx.vec))){
      rx.vec <- sum(rx.vec)
    } else {
      rx.vec <- as.vector(rowSums(rx.vec))
    }
    rx[] <- rx.vec
    rx   <- chol(rx)
    
    cu   <- lapply(1:j, function(k) solve(L[[k]], lambdat %*% zty[[k]]))
    
    cb   <- sapply(1:j, function(i) as.vector(xty[[i]] - crossprod(rzx[[i]], cu[[i]])))
    if(is.null(dim(cb))){
      cb <- sum(cb)
    } else {
      cb   <- as.vector(rowSums(cb))
    }
    cb   <- solve(t(rx), cb)
    
    beta <- solve(rx, cb)
    u    <- lapply(1:j, function(k) solve(t(L[[k]]), cu[[k]] - rzx[[k]] %*% beta))
    mu   <- lapply(1:j, function(k) crossprod(lambdat %*% zt[[k]], u[[k]]) + x[[k]] %*% beta)
    
    prss <- sapply(1:j, function(i) c(y[[i]] - mu[[i]], u[[i]]))
    prss <- sum(prss^2)
    
    ldL  <- sapply(1:j, function(i) diag(L[[i]]))
    ldL  <- sum(log(ldL))
    
    fn   <- N
    if(REML == TRUE){
      ldL <- ldL + sum(log(diag(rx)))
      fn  <- fn - p
    }
    2 * ldL + N * (1 + log(2 * pi * prss) - log(fn))
  }
}


## output module-------------------------------------------------------------------------------
lme.output.mod <- function(mod.out, dev.out, etm){
  setClass("lmeObject", representation(model = "list", deviance = "list", fit = "list"))
  
  lambdat <- mod.out$lambdat
  lind    <- mod.out$lind
  x       <- mod.out$x
  y       <- mod.out$y
  N       <- mod.out$N
  q       <- mod.out$q
  p       <- mod.out$p
  j       <- mod.out$j
  xtx     <- mod.out$xtx
  xty     <- mod.out$xty
  xtz     <- mod.out$xtz
  ztx     <- mod.out$ztx
  zty     <- mod.out$zty
  zt      <- mod.out$zt
  REML    <- mod.out$REML
  getcorr <- mod.out$getcorr
  namesFE <- mod.out$namesFE
  namesRE <- mod.out$namesRE
  
  theta   <- dev.out$par
  dev     <- dev.out$value
  iter    <- dev.out$iterations
  
  lambdat[lind] <- theta
  
  L   <- lapply(1:j, function(k) t(chol(tcrossprod(lambdat %*% zt[[k]]) + diag(1, q))))
  
  rzx <- lapply(1:j, function(k) solve(L[[k]], lambdat %*% ztx[[k]]))
  
  rx   <- matrix(0, p, p)
  rx.vec <- sapply(1:j, function(i) as.vector(xtx[[i]] - crossprod(rzx[[i]])))
  if(is.null(dim(rx.vec))){
    rx.vec <- sum(rx.vec)
  } else {
    rx.vec <- as.vector(rowSums(rx.vec))
  }
  rx[] <- rx.vec
  rx   <- chol(rx)
  
  cu   <- lapply(1:j, function(k) solve(L[[k]], lambdat %*% zty[[k]]))
  
  cb   <- sapply(1:j, function(i) as.vector(xty[[i]] - crossprod(rzx[[i]], cu[[i]])))
  if(is.null(dim(cb))){
    cb <- sum(cb)
  } else {
    cb   <- as.vector(rowSums(cb))
  }
  cb   <- solve(t(rx), cb)
  
  beta <- solve(rx, cb)
  
  u <- lapply(1:j, function(k) solve(t(L[[k]]), cu[[k]] - rzx[[k]] %*% beta))
  
  mu <- lapply(1:j, function(k) crossprod(lambdat %*% zt[[k]], u[[k]]) + x[[k]] %*% beta)
  
  prss <- sapply(1:j, function(i) c(y[[i]] - mu[[i]], u[[i]]))
  prss <- sum(prss^2)
  
  npar     <- p + ((q*(q+1))/2) + 1
  LL       <- dev / -2
  aic      <- dev + 2*npar
  bic      <- dev + npar*log(N)
  s        <- prss / ifelse(REML == F, N, N-p)
  sigma    <- s * crossprod(lambdat)
  var.beta <- s * solve(rx) %*% solve(t(rx))
  std.err  <- sqrt(diag(var.beta))
  t.value  <- (beta / std.err)
  p.value  <- pt(abs(t.value), df = (N-npar), lower.tail = FALSE)
  # p.value  <- sapply(p.value, function(i) ifelse(i < .Machine$double.eps, i <- "< 2e-16", i <- i))
  
  beta <- as.data.frame(matrix(c(beta, std.err, t.value, p.value), nrow = p), row.names = namesFE)
  beta <- format(beta, digits = 4)
  beta[, 4] <- sapply(p.value, function(i) ifelse(i < .Machine$double.eps, i <- paste0("< ", format(.Machine$double.eps, digits = 1)), i <- format(i, digits = 4)))
  names(beta) <- c("Estimate", "Std. error", "t value", "p(>|t|)")
  
  if(getcorr == T){
    vcor  <- matrix(0, q, q)
    dvcor <- sqrt(diag(sigma))
    for(i in 1:q){
      for(k in 1:q){
        vcor[i, k] <- sigma[i, k] / (dvcor[i] * dvcor[k])
      }
    }
    diag(vcor) <- dvcor
    sigma <- vcor
  }
  
  sigma <- data.frame(sigma, row.names = namesRE)
  names(sigma) <- namesRE
  
  fit <- list(
    n.obs  = N,
    npar   = npar,
    LL     = LL,
    AIC    = aic,
    BIC    = bic,
    fixef  = beta,
    vcov   = sigma,
    sigma2 = s
  )
  
  out <- new("lmeObject", model = mod.out, deviance = dev.out, fit = fit)
  
  cat("\n")
  cat("Linear mixed-effects model fit by", if(REML == F) "ML\n" else "REML\n")
  cat("\n")
  cat("Statistical model:\n")
  cat("  Number of groups:", j, "\n")
  cat("  Number of observations:", N, "\n")
  # cat("  Number of components:", k, "\n")
  cat("  Number of parameters:", npar, "\n")
  cat("  Degrees of freedom:", N - npar, "\n")
  cat("\n")
  cat("Iteration process:\n")
  # cat("  Estimation procedure required", restarts, "sets of starting values\n")
  cat("  Model converged normally after", iter, "iterations\n")
  cat("  Time elapsed:", unname(etm[3]), "seconds\n")
  # cat("  Convergence criteria: Log-likelihood difference <", ll.diff.crit, "\n")
  cat("\n")
  cat("Goodness-of-fit statistics\n")
  cat("  Deviance:", dev, "\n")
  cat("  Log-likelihood:", LL, "\n")
  cat("  AIC:", aic, "\n")
  cat("  BIC:", bic, "\n")
  cat("\n")
  cat("Maximum likelihood estimates:")
  cat("\n")
  cat("Fixed effects:\n")
  print(beta)
  cat("\n")
  cat("Residual variance:", s, "\n")
  # print(s)
  cat("\n")
  # cat("Average random effects:\n")
  # print(mu.b)
  # cat("\n")
  # cat("\n")
  cat("Random effects", if(getcorr == F) "(co)variance:\n" else "s.d.-correlation:\n")
  print(sigma)
  cat("\n")
  # cat("\n")
  # cat("Component mixing proportions:\n")
  # print(lambda)
  # cat("\n")
  invisible(out)
}



lm.model.mod <- function(formula, dat, nointFE = FALSE, REML){
  ## initializing and formatting function call
  df <- model.frame(formula, data = dat)
  N <- nrow(df)
  
  namesFE <- names(df)[-1]
  df <- as.matrix(df)
  
  ## formatting response vector
  y <- matrix(df[,1], nrow = N, ncol = 1)
  
  ## formatting predictor matrix
  x <- matrix(df[, -1], nrow = N)
  if(ncol(x) == 0){
    empty <- TRUE
    x <- matrix(1, nrow = N)
    namesFE <- "(intercept)"
  } else {
    empty <- FALSE
    if(nointFE == FALSE){
      x <- cbind(1, x)
      namesFE <- c("(intercept)", namesFE)
    }
  }
  
  qr.x  <- qr(x)
  p     <- qr.x$rank
  r     <- qr.R(qr.x)
  c1    <- qr.qty(qr.x, y)
  c2    <- c1[-(1:p)]
  c1    <- c1[1:p]
  beta  <- solve(r, c1)
  rss   <- sum(c2^2)
  sigma <- rss/(N-p)
  
  out <- list(
    x       = x,
    qr.x    = qr.x,
    y       = y,
    beta    = beta,
    sigma   = sigma,
    rss     = rss,
    empty   = empty,
    N       = N,
    namesFE = namesFE,
    nointFE = nointFE,
    REML    = REML
  )
  out
}



lm.output.mod <- function(x, qr.x, y, beta, sigma, rss, empty, N, namesFE, nointFE, REML){
  r     <- qr.R(qr.x)
  p     <- qr.x$rank
  var.b <- sqrt(diag(sigma * chol2inv(r)))
  t.val <- beta / var.b
  pred  <- x %*% beta
  ifelse(nointFE == F, diff  <- pred- mean(y), diff <- pred)
  ssr   <- sum(diff^2)
  f.val <- (ssr / (p)) / sigma
  r2    <- 1 - (rss / (rss + ssr))
  adj.r <- 1 - (1-r2)*((N-1)/(N-p))
  if(REML == TRUE){
    LL  <- 0.5*(-(N-p)*(log(2*pi) + 1 - log(N-p) + log(rss))) - sum(log(abs(diag(qr.x$qr)[1L:p])))
  } else {
    LL  <- 0.5*(-N*(log(2*pi) + 1 - log(N) + log(rss))) 
  }
  dev   <- -2*LL
  aic   <- dev + 2*(p+1)
  bic   <- dev + (p+1)*log(N)
  
  b <- beta
  beta <- cbind(beta, var.b, t.val)
  beta <- data.frame(beta, row.names = namesFE)
  names(beta) <- c("Estimate", "Std. error", "t value")
  
  out <- list(
    n.obs = N, fix.ef = b, sigma = sigma, LL = LL, aic = aic, bic = bic, pred = pred, res = diff
  )
  
  cat("Linear model fit by QR factorization\n")
  cat("Statistical model\n")
  cat("  Number of observations:", N, "\n")
  cat("  Number of parameters:", p, "\n")
  cat("  Degrees of freedom residual:", N - p, "\n")
  if(empty == FALSE) cat("  Degrees of freedom F:", p, "and", N - p, "\n")
  cat("\n")
  cat("Goodness-of-fit statistics\n")
  cat("  Log-likelihood:", LL, "\n")
  cat("  AIC:", aic, "\n")
  cat("  BIC:", bic, "\n")
  cat("\n")
  cat("Coefficients\n")
  print(beta)
  cat("\n")
  cat("Residual variance:", sigma, "\n")
  if(empty == FALSE) cat("R-squared:", r2, "\n")
  if(empty == FALSE) cat("Adjusted R-squared:", adj.r, "\n")
  if(empty == FALSE) cat("F-statistic:", f.val, "\n")
  invisible(out)
}


lme <- function(formula = NULL, random = NULL, id = NULL, dat = NULL, 
                nointFE = FALSE, nointRE = FALSE, REML = FALSE, getcorr = FALSE)
{
  ## some safety coding
  if(is.null(formula)) stop("No formula specified!")
  form.out <- parse.lme.formula(formula)
  # if(is.null(random)) mod.type <- "lm" else mod.type <- "lmm"
  
  if(form.out$lme == TRUE){
    ## executing each module----------------------------------------------------------------------
    # model module
    
    mod.out <- lme.model.mod(y = form.out$y, fixed = form.out$fe, random = form.out$re, 
                             id = form.out$id, dat = dat,
                             nointFE = form.out$nointFE, nointRE = form.out$nointRE, 
                             REML = REML, getcorr = getcorr)
    
    devfun  <- pls(mod.out)
    
    # timing and executing the optimization module
    stm <- proc.time()
    dev.out <- nelder.mead(parm = mod.out$theta.start, objective = devfun, 
                           lower = mod.out$lower, upper = mod.out$upper)
    etm <- proc.time() - stm
    
    # output module
    lme.output.mod(mod.out, dev.out, etm = etm)
  } else {
    mod.out <- lm.model.mod(formula = formula, dat = dat, nointFE = nointFE, REML = REML)
    lm.output.mod(x = mod.out$x, qr.x = mod.out$qr.x, y = mod.out$y, beta = mod.out$beta, 
                  sigma = mod.out$sigma, rss = mod.out$rss, empty = mod.out$empty, N = mod.out$N, 
                  namesFE = mod.out$namesFE, nointFE = mod.out$nointFE,
                  REML = mod.out$REML)
  }
  
}

lme(y ~ x + (1 | id), dat = test)
# lme(formula = y~1, random = NULL, id = "id", dat = test, nointFE = F, REML = F)
# lme(formula = y~1, random = 1, id = "id", dat = test, nointFE = F, REML = F)


