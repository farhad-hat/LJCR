# Simple model (no clustering) for high-dimensional LMM regression
# High-dimensional fixed effects are penalized using L2-penalization
# to allow for closed-form updates of the beta parameters

library(expm)
library(Matrix)

# Function to convert from unconstrained parameter space to
# positive definite covariance matrix of the form D*sigma^{-2}
# See Pinheiro and Bates for details on this parameterization.
unconstrained_params_to_cov <- function(param_vec, q) {
  param_matrix <- matrix(0, q, q)
  param_matrix[upper.tri(param_matrix, diag=TRUE)] <- param_vec
  param_matrix[lower.tri(param_matrix, diag=TRUE)] <- param_vec
  
  return(expm(param_matrix, method='R_Eigen'))
}

# Function to calculate the profile log likelihood
# param_vec is unconstrained (parameters can take any values)
profile_log_likelihood <- function(param_vec, y, X, Z, N, m, pen_matrix) {
  
  p = dim(X[[1]])[2]
  q = dim(Z[[1]])[2]
  M = length(y)
  
  # Convert param vector to positive definite matrix D*sigma^{-2}
  D_sigma <- unconstrained_params_to_cov(param_vec, q)
  
  # In case close to not positive-definite, reject. 
   if(any(eigen(D_sigma)$values < 1e-5)) return(list(ll=-Inf))
  
  # Invert to get D^{-1}*sigma^2
  inv_D_sigma <- chol2inv(chol(D_sigma)) 
  
  # Just to be sure, also reject if the inverse is not positive-definite
  if(any(eigen(inv_D_sigma)$values < 1e-5)) return(list(ll=-Inf))
  
  # Cholesky decompose to get Delta
  Delta <- chol(inv_D_sigma)
  
  RC_cp_temp <- matrix()
  RC_cp <- matrix(0, nrow = p+1, ncol = p+1)
  
  log_det_delta = determinant(Delta, logarithm=TRUE)$modulus
  log_det_ratio = 0
 
  # Loop over patients
  for(i in 1:M) {
    # No need to include the penalization matrix here, as it will be included at the end.
    Z_star <- rbind(Z[[i]], Delta)
    y_star <- rbind(y[[i]], matrix(0, q, 1))
    X_star <- rbind(X[[i]], matrix(0, q, p))
    
    Z_qr <- qr(Z_star)
    R_11 <- qr.R(Z_qr)[1:q,]
    Q_t <- t(qr.Q(Z_qr, complete=TRUE))
    R_i <- Q_t %*% X_star
    c_i <- Q_t %*% y_star
    
    RC <- cbind(sqrt(m[i]) * R_i[(q+1):dim(R_i)[1],],
                sqrt(m[i]) * c_i[(q+1):length(c_i)])
    
    RC_cp_temp <- crossprod(RC)
    
    # Add penalization
    RC_cp_temp[1:p,1:p] <- RC_cp_temp[1:p,1:p] + as.matrix(pen_matrix)
    
    log_det_ratio <- log_det_ratio + (m[i] * log_det_delta) -
      (m[i] * determinant(R_11, logarithm=TRUE)$modulus)
    
    RC_cp <- RC_cp + RC_cp_temp
  }
 
  # Sum the t(RC_i) %*% RC_i to get full matrix 
  # t(RC) %*% RC. Then use choleski decomposition
  # of this to get matrix R from the QR decomposition.
  # See: https://www.cs.cornell.edu/~arb/papers/mrtsqr-bigdata2013.pdf
  # if(any(eigen(Reduce('+', RC_cp))$values < 1e-5)) stop('\n', 'Increase sample size M or decrease number of groups K')
  # if(any(eigen(Reduce('+', RC_cp))$values < 1e-5)) browser()
  R_all = chol(RC_cp)
  
  # Estimate of residual sum 
  resids <- R_all[(p+1):dim(R_all)[1],dim(R_all)[2]]
  
  # Log likelihood
  ll <- 0.5*N*(log(N) - log(2*pi) - log(sum(resids^2)) - 1) + log_det_ratio
  
  return(list(ll=ll, R_all=R_all))
}

# Function to calculate the profile log likelihood
# param_vec is unconstrained (parameters can take any values)
# Does not use QR decompositions
profile_log_likelihood_naive <- function(param_vec, y, X, Z, N, pen_matrix, 
                                         beta_hat=NULL, sigma_hat=NULL, D_sigma=NULL) {
  
  p = dim(X[[1]])[2]
  q = dim(Z[[1]])[2]
  M = length(y)
  
  # Convert param vector to positive definite matrix D*sigma^{-2}
  if(is.null(D_sigma)) {
    D_sigma = unconstrained_params_to_cov(param_vec, q)
  }
  
  # In case close to not positive-definite, reject. 
  if(any(eigen(D_sigma)$values < 1e-5)) return(list(ll=-Inf))
  
  XX_term = matrix(0, p, p)
  Xy_term = matrix(0, p, 1)
  
  inv_V = list()
  V = list()
  
  # Loop over patients
  for(i in 1:M) {
    Z_star = rbind(Z[[i]], matrix(0, p, q))
    y_star = rbind(y[[i]], matrix(0, p, 1))
    X_star = rbind(X[[i]], pen_matrix)
    
    V[[i]] = diag(n_i+p) + Z_star %*% tcrossprod(D_sigma, Z_star)
    
    # Invert V_i
    inv_V[[i]] = chol2inv(chol(V[[i]])) 
    
    # Just to be sure, also reject if the inverse is not positive-definite
    if(any(eigen(inv_V[[i]])$values < 0)) return(list(ll=-Inf))
    
    pre_term = crossprod(X_star, inv_V[[i]])
    
    XX_term = XX_term + pre_term %*% X_star
    Xy_term = Xy_term + pre_term %*% y_star
    
  }
  
  if(is.null(beta_hat)) {
    beta_hat = chol2inv(chol(XX_term)) %*% Xy_term
  }
  
  
  sigma_hat_temp = 0
  
  resid_prod = list()
  
  # Loop over patients
  for(i in 1:M) {
    y_star = rbind(y[[i]], matrix(0, p, 1))
    X_star = rbind(X[[i]], pen_matrix)
    
    resid = y_star - X_star %*% beta_hat
    
    resid_prod[[i]] = crossprod(resid, inv_V[[i]]) %*% resid
    
    sigma_hat_temp = sigma_hat_temp + resid_prod[[i]]
  }
  
  if(is.null(sigma_hat)) {
    sigma_hat = sigma_hat_temp / N
  }
  
  ll = 0
  
  resid_sum = 0
  det_sum = 0
  
  # Loop over patients
  for(i in 1:M) {
    
    det_V = determinant(V[[i]], logarithm=TRUE)$modulus
    det_sum = det_sum + det_V
    
    # Log likelihood
    resid_sum = resid_sum + resid_prod[[i]]/sigma_hat
    
    ll = ll - 0.5*(n_i*log(2*pi*sigma_hat) + resid_prod[[i]]/sigma_hat + 
                     det_V)
  }
  
  #print(det_sum)
  #print(resid_sum)
  
  return(list(ll=ll, beta_hat=beta_hat, sigma_hat=c(sigma_hat)))
}

