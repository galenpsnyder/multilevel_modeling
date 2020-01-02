lmeEM <- function(formula = NULL, random = NULL, id = NULL, k = 1, dat = NULL, nointFE = FALSE, nointRE = FALSE, REML = FALSE, getcorr = FALSE){
  # initializing and formatting function call
  df <- model.frame(formula, data = dat)
  namesFE <- names(df)[-1]
  df <- as.matrix(df)
  id <- factor(dat[,id])
  y <- df[,1]
  
  # saving full response vector for reweighting
  yy <- y      
  
  x <- df[,-1]
  z <- as.matrix(dat[random])
  N <- nrow(df)
  
  namesRE <- attributes(z)$dimnames[[-1]]
  
  y <- split(x = y, f = id)
  n <- sapply(y, length)
  j <- length(y)
  
  if(nointFE == FALSE){
    one <- rep(1, N)
    x <- cbind(one, x)
    p <- ncol(x)
    x <- split(x = x, f = id)
    x <- lapply(1:j, function(i) matrix(x[[i]], nrow = n[[i]], ncol = p))
    namesFE <- c("(intercept)", namesFE)
  } else {
    p <- ncol(x)
    x <- split(x = x, f = id)
  }
  
  # saving full FE model matrix for reweighting
  xx <- do.call(rbind, x)                             
  
  if(nointRE == FALSE){
    one <- rep(1, N)
    z <- cbind(one, z)
    q <- ncol(z)
    z <- split(x = z, f = id)
    z <- lapply(1:j, function(i) matrix(x[[i]], nrow = n[[i]], ncol = p))
    namesRE <- c("(intercept)", namesRE)
  } else {
    q <- ncol(as.matrix(z))
    z <- split(x = z, f = id)
  }
  
  zz <- do.call(rbind, z)
  
  ztz <- lapply(1:j, function(k) crossprod(z[[k]]))
  
  ll <- vector(mode = "numeric", length = 1000L)
  a <- 2
  
  # initializing beta
  beta <- solve(crossprod(xx)) %*% crossprod(xx, yy)
  
  # initializing random effects covariance and residual variance
  theta <- runif(n = 1, min = 0.25, max = 0.5) * sqrt(colSums(do.call(rbind, lapply(1:j, function(x) colSums(z[[x]]^2) / j))))
  psi <- diag(theta, q, q)
  s <- runif(1, min = .75 * var(yy), max = var(yy))
  
  # initial residual
  res <- lapply(1:j, function(k) y[[k]] - x[[k]]%*%beta)
  
  sigma <- lapply(1:j, function(k) z[[k]]%*%tcrossprod(psi, z[[k]]) + diag(s, n[[k]]))
  
  if(REML == FALSE){
    ll[a] <- do.call(sum ,lapply(1:j, function(k) -0.5 * n[[k]] * log(2*pi) - 0.5*log(det(sigma[[k]])) - 0.5*crossprod(res[[k]], solve(sigma[[k]]))%*%res[[k]]))
  } else {
    ll[a] <- do.call(sum ,lapply(1:j, function(k) -0.5 * n[[k]]-p * log(2*pi) - 0.5*log(det(sigma[[k]])) - 0.5*crossprod(res[[k]], solve(sigma[[k]]))%*%res[[k]]))
  }
  
  while(abs(ll[a] - ll[a-1]) > 1.0e-6){ # convergence criteria set to be commensurate with lme
    # inverting the RE covariance matrix
    psi <- solve(psi)

    # calculating gamma, a within-cluster RE covariance matrix
    gamma <- lapply(1:j, function(k) solve(ztz[[k]]/s + psi))

    # calculating the expected values of the random effects given y
    b <- lapply(1:j, function(k) (gamma[[k]]  %*% crossprod(z[[k]], res[[k]]))/s)

    # conditional estimate of fixed effects
    bu <- lapply(1:j, function(k) z[[k]]%*%b[[k]])
    beta <- solve(crossprod(xx))%*%crossprod(xx, yy - do.call(c, bu))

    # conditional residual--relative to fixed effects
    res <- lapply(1:j, function(k) y[[k]] - x[[k]]%*%beta)

    # contiditional estimate of RE covariance
    psi <- 0
    for(k in 1:j){
      psi <- psi + (gamma[[k]] + tcrossprod(b[[k]]))
    }
    psi <- psi / j
    
    # # calculating trace--second term in log-likelihood
    # tr <- do.call(sum, lapply(1:j, function(k) sum(diag(ztz[[k]]%*%(gamma[[k]] + tcrossprod(b[[k]]))))))
    # 
    # reszb <- do.call(sum, lapply(1:j, function(k) crossprod(res[[k]], z[[k]])%*%b[[k]]))
    # 
    # # conditional estimate of residual variance
    # if(REML == FALSE){
    #   s <- as.numeric((sum(do.call(c, res)^2) + tr - 2*reszb) / N)
    # } else {
    #   s <- as.numeric((sum(do.call(c, res)^2) + tr - 2*reszb) / (N-p))
    # }
    
    s <- do.call(sum, lapply(1:j, function(k) sum(diag(tcrossprod(y[[k]] - x[[k]] %*% beta - bu[[k]]) + z[[k]] %*% tcrossprod(gamma[[k]], z[[k]]))))) / N

    # conditional estimate of joint variance
    sigma <- lapply(1:j, function(k) z[[k]]%*%tcrossprod(psi, z[[k]]) + diag(s, n[[k]]))

    if(REML == FALSE){
      ll[a + 1] <- do.call(sum ,lapply(1:j, function(k) -0.5 * n[[k]] * log(2*pi) - 0.5*log(det(sigma[[k]])) - 0.5*crossprod(res[[k]], solve(sigma[[k]]))%*%res[[k]]))
    } else {
      ll[a + 1] <- do.call(sum ,lapply(1:j, function(k) -0.5 * n[[k]] * log(2*pi) - 0.5*log(det(sigma[[k]])) - 0.5*crossprod(res[[k]], solve(sigma[[k]]))%*%res[[k]]))
    }

    cat("Iteration",a,": Log-likelihood =",ll[a], "\n")
    cat("Log-likelihood difference = ", ll[a] - ll[a-1], "\n")

    a <- a + 1
  }

  # formatting output
  beta <- data.frame(matrix(beta, nrow = 1), row.names = "Estimates")
  names(beta) <- namesFE

  psi <- data.frame(psi, row.names = namesRE)
  names(psi) <- namesRE

  cat("\n")
  cat("Linear Mixed-Effects Model Fit by ")
  cat(if(REML == FALSE) "ML\n" else "REML\n")
  cat("Model converged normally after", a, "iterations", "\n")
  cat("  Log-",if(REML == FALSE) "likelihood" else "restricted-likelihood" , "at convergence:", ll[a], "\n")
  cat("\n")
  cat("Fixed Effects:", "\n")
  print(beta)
  cat("\n")
  cat("Residual Variance:", s, "\n")
  cat("\n")
  cat("Random Effects Covariance:", "\n")
  print(psi)
  cat("\n")
  cat("Number of Observations:", N, "\n")
  cat("Number of Groups:", j,"\n")
  cat("\n")
}
(lmeEM(formula = y ~ x, random = c("x"), id = "id", k = 1, dat = test))


lme(fixed = y~x, data = test, random = ~x|id, method = "ML")
summary(out)
