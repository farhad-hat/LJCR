# Run on the HEC
t = Sys.getenv('SGE_TASK_ID') 
K = as.numeric(t)


set.seed(100)

library(here)
library(Metrics)
library(mvtnorm)
library(glmnet)
library(MASS)
library(foreach)
library(doParallel)
library(nlme)
library(clusterGeneration)
library(MCMCglmm)
library(BBmisc)
library(scales)
library(glasso)

source(here('lmm_high_dim_optimization_more_than_one_group.R'))

ljcr_on_ALS_data <- function(K = 2, N_EM_iteration = 40, width = 1){
  
  # Load ALS dataset and
  df_lme_reduced <- readRDS('df_lme_reduced.RData')
  df_lme_reduced <- df_lme_reduced[, c(1, 3, 2, 4:dim(df_lme_reduced)[2])]
  df_lme_reduced <- data.frame(apply(df_lme_reduced, 2, function(x) as.numeric(as.character(x))))
  df_lme_reduced <- cbind(df_lme_reduced[,c(1,2)], as.data.frame(scale(df_lme_reduced[,-c(1,2)])))
  
  # remove colinear fatures
  df_lme_reduced <- within(df_lme_reduced, rm(GenderF, RaceBlack, Urine_AppearanceCLOUDY, Urine_ColorYELLOW))
  
  df_lme_reduced <- df_lme_reduced[,c(1:56,62)]
  
  SubjectIDs <- sort(unique(df_lme_reduced$SubjectID))
  
  # Normalizing onset-delta feature column
  # df_lme_reduced[,"onset_delta"] <- normalize(-df_lme_reduced[,"onset_delta"], method = "standardize")
  
  # removing those patients who have less than 3 observations, due to the
  # error it produces in RC_cp[[i]][1:p,1:p] (subscript out of bound)
  more_than_3_obs_IDs <- vector()
  for (ID in SubjectIDs) {
    if(dim(df_lme_reduced[df_lme_reduced$SubjectID == ID, ])[1] >= 2){
      more_than_3_obs_IDs <- append(more_than_3_obs_IDs, ID)
    }
  }
  
  # Re-defining df_lme_reduced
  df_lme_reduced <- subset(df_lme_reduced, SubjectID %in% more_than_3_obs_IDs)
  SubjectIDs <- sort(unique(df_lme_reduced$SubjectID))
  
  # Inputs
  M <- length(SubjectIDs) # Number of patients
  K <- K # Number of groups 
  N_EM_iteration <- N_EM_iteration # Number EM iterations
  
  # Count the number of data points for patinet i
  n_i_vec <- vector()
  n_i_vec <- sapply(1:length(SubjectIDs), function(i)
    dim(df_lme_reduced[df_lme_reduced$SubjectID == SubjectIDs[i],])[1])
  
  N <- dim(df_lme_reduced)[1] # Total number of data points
  p <- dim(df_lme_reduced)[2]-2 # Number of covariates
  q <- 2 # Number of random effects
  y <- list(); X <- list()
  Z <- list(); b <- list()
  size = list()
  size = sapply(1:N_EM_iteration, function(i) list())
  beta_list <- list()
  sigma_list <- list()
  D_list <- list()
  beta_inferred_list <- list()
  rmse_sigma_list <- list()
  rmse_D_list <- list()
  m_list <- sapply(1:K, function(z) vector())
  m_list_pre <- list()
  ll_value_list <- list()
  
  # Run joint lme to start with some initial values for beta
  lme_joint <- function(df_lme_reduced){
    xnam <- colnames(subset(df_lme_reduced, select=-c(SubjectID, feature_delta, ALSFRS_Total)))
    fmla <- as.formula(paste("ALSFRS_Total ~ (", paste(xnam, collapse= "+"), ")*feature_delta"))
    
    lme_joint <- lme(fmla, random=~1|SubjectID, data=df_lme_reduced, control=lmeControl(returnObject=TRUE))
    return("lme_joint"=lme_joint)
  }
  
  lme_joint <- lme_joint(df_lme_reduced)
  
  # Initial values for parameters
  for (k in 1:K) { 
    ## Random effect covariance matrix
    D_list[[k]] <- genPositiveDefMat("eigen", dim = q, rangeVar=c(1,10))$Sigma
    
    ## fixed effects
    beta_list[[k]] <- as.vector(coefficients(summary(lme_joint))[c(1:p+1), "Value"])
    beta_list[[k]] <- beta_list[[k]][c(length(beta_list[[k]]),1:(length(beta_list[[k]])-1))]
    
    ## noise
    sigma_list[[k]] <- runif(1)
  }
  
  # Store initial values of the parameters
  D_list_init <- D_list
  beta_list_init <- beta_list
  sigma_list_init <- sigma_list
  
  # ALS data
  for(i in 1:M) {
    Z[[i]] <- cbind(1, df_lme_reduced[df_lme_reduced$SubjectID == SubjectIDs[i],][,3])
    X[[i]] <- as.matrix(df_lme_reduced[df_lme_reduced$SubjectID == SubjectIDs[i],][,-c(1,2)])
    y[[i]] <- as.matrix(df_lme_reduced[df_lme_reduced$SubjectID == SubjectIDs[i],][,2])
  }
  
  # Rescaling the ALSFRS values (responses; y)
  y_scaled <- scale(unlist(y[c(1:M)])) * width
  for(i in 1:M) {
    y[[i]] <- as.matrix(y_scaled[1:n_i_vec[i]])
    y_scaled <- y_scaled[-c(1:n_i_vec[i])]
  }

  # Responsibilities
  responsibilities_initial_iteration <- function(m_list, y, X, beta_list, Z, D_list, sigma_list){
    
    p_y_temp <- list()
    p_y_temp <- sapply(1:K, function(i) vector())
    for (k in 1:K){
      p_y_temp[[k]] <- runif(n=M, min=0, max=1)
    }
    
    for (k in 1:K){
      m_list[[k]] <- p_y_temp[[k]]/(Reduce('+', p_y_temp))
    }
    
    return(m_list=m_list)
    
  }
  
  responsibilities <- function(m_list, y, X, beta_list, Z, D_list, sigma_list){
    
    p_y <- list()
    p_y <- sapply(1:K, function(z) vector())
    
    for (k in 1:K) {
      for (i in 1:M) {
        p_y[[k]][i] <- dmvnorm(c(y[[i]]), 
                               mean = X[[i]] %*% beta_list[[k]], 
                               sigma = (Z[[i]] %*% D_list[[k]] %*% t(Z[[i]])) + (sigma_list[[k]] * diag(n_i_vec[i])))
      }
    }
    log_p_x <- list()
    log_p_x <- sapply(1:K, function(z) vector())
    
    for (k in 1:K) {
      for (i in 1:M) {
        log_p_x[[k]][i] <- 0
        for (t in 1:n_i_vec[i]) {
          log_p_x[[k]][i] <- log_p_x[[k]][i] + 
            log(dmvnorm(c(X[[i]][t,]), mean = mean_prime[[k]], sigma = Sigma_prime[[k]]))
        }
      }
    }
    
    px_zero_ind <- vector()
    for(k in 1:K){
      px_zero_ind <- append(px_zero_ind, which(lapply(log_p_x[[k]], exp)==0))
    }
    px_zero_ind <- unique(px_zero_ind)
    
    denominator_vec <- vector()
    for (i in 1:M) {
      denominator <- 0
      for (k in 1:K) {
        if(i %in% px_zero_ind){
          denominator <- 1
        }else{
          denominator <- denominator + (p_y[[k]][i] * exp(log_p_x[[k]][i]) * sum(m_list[[k]]))
          
          if(denominator==0){
            px_zero_ind <- sort(unique(append(px_zero_ind, i)))
            denominator <- 1
          }
        }
      }
      denominator_vec[i] <- denominator
    }
    
    m_list_temp <- list()
    m_list_temp <- sapply(1:K, function(z) vector())
    
    for (i in 1:M) {
      for (k in 1:K) {
        if(i %in% px_zero_ind){
          m_list_temp[[k]][i] <- 1/K
        }else{
          m_list_temp[[k]][i] <- 
            (p_y[[k]][i] * exp(log_p_x[[k]][i]) * sum(m_list[[k]])) / denominator_vec[i]
        }
      }
    }
    
    m_list <- m_list_temp
    return(m_list=m_list)
  }
  
  # Initialization for optim
  init_vals = c(10,1,1)
  
  # Objective needs to be negative log likelihood
  objective_fun <- function(param_vec, y, X, Z, N, m, pen_matrix) {
    -profile_log_likelihood(param_vec, y, 
                            X, Z, N, m, pen_matrix)$ll
  } 
  
  # Define penalization
  lambda = 0
  pen_matrix = sqrt(lambda/M)*diag(p)
  
  # EM algorithm
  for (iteration in 1:N_EM_iteration) {
    
    beta_inferred_list[[iteration]] <- list()
    rmse_sigma_list[[iteration]] <- list()
    rmse_D_list[[iteration]] <- list()
    
    ## E step (Update responsibilitis)
    
    if(iteration == 1){
      m_list <- responsibilities_initial_iteration(m_list, y, X, beta_list, Z, D_list, sigma_list)
      m_list_pre <- m_list
    }else{
      m_list <- responsibilities(m_list, y, X, beta_list, Z, D_list, sigma_list)
    }
    
    ### Log_likelihood value of the previous iteration step
    if(iteration > 3) ll_value_last_three <- ll_value_last_two
    if(iteration > 2) ll_value_last_two <- ll_value_last
    if(iteration > 1) ll_value_last <- ll_value
    
    ll_value <- 0
    
    # M step (estimate parameters (beta, sigma, D))
    
    # Estimate mean and covariance parameters for X
    ########## updated   
    X_t_m <- sapply(1:K, function(z) list())
    mean_prime <- sapply(1:K, function(z) list(0))
    Sigma_prime <- list()
    
    for (k in 1:K) {
      m_vec <- vector()
      for(i in 1:M) {
        mean_prime[[k]] <- mean_prime[[k]] + as.vector(m_list[[k]][i] * colSums(X[[i]]))
        X_t_m[[k]][[i]] <- t(X[[i]]) * m_list[[k]][i]
        m_vec <- append(m_vec, rep(m_list[[k]][i], n_i_vec[i]))
      }
      mean_prime[[k]] <- mean_prime[[k]] / sum(m_list[[k]] * n_i_vec)
      Sigma_prime[[k]] <- cov.wt(do.call(rbind, X), m_vec)$cov
      Sigma_prime[[k]] <- glasso(Sigma_prime[[k]], sqrt(2*M*log(p))/2, penalize.diagonal=FALSE)$w
    }
    ########## end updated   
    
    ## Run optimization in parallel
    cores=detectCores()
    cl <- makeCluster(cores[1])
    registerDoParallel(cl)
    # optim_results <- list()
    # for (k in 1:K) {
    optim_results_parallel <- foreach(k = 1:K) %dopar% {
      library(expm)
      library(Matrix)
      library(matrixcalc)
      library(here)
      source(here('lmm_high_dim_optimization_more_than_one_group.R'))
      
      # for (k in 1:K) {  
      optim_results = optim(init_vals, objective_fun, y=y, X=X, Z=Z, N=N, 
                            m=m_list[[k]], pen_matrix=pen_matrix^2)
    }
    
    stopCluster(cl) 
    
    for (k in 1:K) {
      
      ll_result = profile_log_likelihood(optim_results_parallel[[k]]$par, y, X, Z, N,
                                         m=m_list[[k]], pen_matrix)
      
      R_all = as.matrix(ll_result$R_all)
      
      ## Recover parameters
      beta_inferred = solve(R_all[1:p,1:p], tol = 1e-20)%*%R_all[1:p,p+1]
      sigma2_inferred = R_all[(p+1):dim(R_all)[1],dim(R_all)[2]]^2/N
      D_sigma = unconstrained_params_to_cov(optim_results_parallel[[k]]$par, q)
      ll_value <- ll_value + as.numeric(ll_result$ll)
      
      # print(D_sigma*sigma2_inferred) # D recovered?
      # plot(beta_list_init[[k]], beta_inferred) # Beta recovered?
      # print(sqrt(sigma2_inferred)) # Sigma recovered?
      
      ## Save results
      beta_inferred_list[[iteration]][k] <- list(beta_inferred)
      rmse_sigma_list[[iteration]][k] <- list(rmse(sqrt(sigma2_inferred), sigma_list[[k]]))
      rmse_D_list[[iteration]][k] <- list(rmse((D_sigma*sigma2_inferred), D_list[[k]]))
      
      ## Update parameters
      sigma_list[[k]] <- sqrt(sigma2_inferred)
      beta_list[[k]] <- beta_inferred
      D_list[[k]] <- D_sigma*sigma2_inferred
      
      # Size of the groups
      size[[iteration]][[k]] <- list(sum(m_list[[k]] >= 1/K))
    }
    
    ## Save results of log-likelihood for different iteration
    ll_value_list[[iteration]] <- ll_value
    
    print(size[[iteration]])
    cat("log_likelihood =           ", as.numeric(ll_value), '\n')
    
    # Stop EM if log_likelihood has almost converged
    if(iteration > 1){
      if(abs(ll_value - ll_value_last) < 1) break()
      # if(ll_value > 54000) break()
    }
    
    # Stop EM if log_likelihood convergence has trapped in a loop
    if(iteration > 3){
      if(round(ll_value, 1) == round(ll_value_last_two, 1) &&
         round(ll_value_last, 1) == round(ll_value_last_three, 1)) break()
    }
    
    # Stop EM if size of one of the groups is small
    stop_loop = FALSE
    if(iteration > 1 ){
      for (k in 1:K) {
        if(sum(m_list[[k]] > (1/K)) <= (M/(10*K))) {
          stop_loop = TRUE
          break()
        }
        if(stop_loop) break()
      }
      if(stop_loop) break()
    }
    if(stop_loop) break()
  }
  
  return(list('beta_list_init'=beta_list_init,
              'beta_list_inferred'=beta_inferred_list,
              'sigma_list_init'=sigma_list_init,
              'sigma_list'=sigma_list,
              'D_list_init'=D_list_init,
              'D_list'=D_list,
              'rmse_sigma_list'=rmse_sigma_list,
              'rmse_D_list'=rmse_D_list,
              'll_value_list'=ll_value_list,
              'm_list_pre'=m_list_pre,
              'm_list'=m_list,
              'size'=size,
              'K'=K,
              'iteration'=iteration))
}

storage = Sys.getenv("global_scratch") 

ljcr_results <- ljcr_on_ALS_data(K=K, N_EM_iteration=300, width=2)
saveRDS(ljcr_results, file.path(storage, paste0("ALS_1_K_",ljcr_results$K,".RData")))
