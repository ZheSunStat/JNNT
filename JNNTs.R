# library(mvtnorm)
# library(invgamma)
# library(truncnorm)
# library(dplyr)
# library(igraph)

# MCMC sampling based on with a combination of Gibbs sampler and MH steps for Joint Node and Network Thresholded selection (JNNTs) model

JNNTs <- function(Y, W, X, Z, d, intercept = TRUE,
                      rho, a_rho, b_rho,
                      lambda, a_lambda, b_lambda,
                      sig2_i = 10,
                      MH_step = 50L, MH_stop = 2000L, MH_acptrate = 0.234,
                      burnin = 5000L, M = 10000L,
                      disp = FALSE, disp_step = 100L,
                      Init_seed, Init_list) {
  ############################################################
  #'@param input
  # Y n x 1 matrix or vector or length n, for response
  # W n x q matrix or missing, for confounder variables
  # X n x p matrix, for region data
  # Z p x p x n array, for connectivity data #! symmetric to save space
  # d integer indicates number of sub-networks involved #! fixed
  # intercept boolean variable, to indicate whether include intercept or not. Default is TRUE
  # rho If rho provided, set it as fixed value; otherwise update rho via MH step with prior distribution U(a_rho, b_rho)
  # a_rho, b_rho, interval range for uniform distribution for rho
  # lambda If lambda provided, set it as fixed value; otherwise update lambda via MH step with prior distribution U(a_lambda, b_lambda)
  # a_lambda, b_lambda, interval range for uniform distribution for lambda
  # sig2_i positive scaler indicates sd for the normal distribution of ita coefficients for confounders #! fixed
  # MH_step, MH_stop, MH_acptrate
  # M, burnin number of MCMC steps after and for burn-in
  # disp, disp_step
  # Init_seed, Init_list

  #'@param output
  # ita_save   q1 x M # FIXME
  # beta_save  p x M
  # gamma_save p x M
  # Alpha_save p x d x M
  # theta_save p x M
  # Theta_save p x d x M
  # FIXME sig2's - 5 x M vectors
  # rho_save if rho provided, rho_save is a scaler with the same value of rho; otherwise rho_save is a vector with length M.
  # lambda_save if lambda provided, lambda_save is a scaler with the same value of lambda; otherwise lambda_save is a vector with length M.
  # rhoss_save, lambdass_save are vectors of length M, respectively, if rho is not provided; otherwise they are not included in the output results.
  # rhoacpt_save, lambdaacpt_save are vectors of length M, respectively, if lambda is not provided; otherwise they are not included in the output results.
  #  M vector
  # elapse_time
  # initial_val Initial values

  ############################################################

  #### >>> ANCHOR convert input into right mode  <<< ####
  ## logical scaler
  ##      Alpha_include, beta_include
  ##      rho_update, lambda_update
  # NOTE feature region or network predictor
  Alpha_include = beta_include = FALSE
  rho_update = lambda_update = FALSE

  Y = matrix(Y, ncol = 1)
  mode(Y) = "double"
  # NOTE feature region or network predictor
  if(missing(W) && missing(X) && missing(Z))
    stop("Provide at least one of W, X and Z")
  # NOTE feature region or network predictor
  if(!missing(X)) {
    mode(X) = "double"
    beta_include = TRUE
    if(dim(X)[1] != nrow(Y)) stop("X is of different subject number from that of Y")
  }
  if(!missing(Z)) {
    mode(Z) = "double"
    Alpha_include = TRUE
    if(dim(Z)[3] != nrow(Y)) stop("Z is of different subject number from that of Y")
  }
  if(!Alpha_include && !beta_include)
    stop("Missing both regional and network predictors") # a trivia case without either regional nor network predictors in the model
  if(Alpha_include && beta_include) {
    if(dim(X)[2] != dim(Z)[2]) stop("X and Z are with different feature dimension")
  }
  if(!missing(W)) {
    if(is.vector(W)) W = matrix(W, ncol = 1)
    if(nrow(W) != nrow(Y)) stop("W is of different subject number from that of Y")
    mode(W) = "double"
  }
  #// sig2V_include = beta_include || Alpha_include
  joint_include = beta_include && Alpha_include
  if(joint_include) { # if(!missing(X) && !missing(Z)) {
    if(dim(X)[2] != dim(Z)[2]) stop("X and Z are with different dimensions")
  }
  if(!Alpha_include) d = 1

  d = as.integer(d)
  intercept = as.logical(intercept)

  if(missing(rho)) {
    rho = double(length = 1)
    if(missing(a_rho) || missing(b_rho)) stop("Please provide lower and upper bound for rho")
    a_rho = as.double(a_rho)
    b_rho = as.double(b_rho)
    rho_update = TRUE
  }
  if(missing(lambda)) {
    lambda = double(length = 1)
    if(missing(a_lambda) || missing(b_lambda)) stop("Please provide lower and upper bound for lambda")
    a_lambda = as.double(a_lambda)
    b_lambda = as.double(b_lambda)
    lambda_update = TRUE
  }
  if(missing(Init_seed)) {
    Init_seed <- as.integer(Sys.times())
  } else {
    Init_seed = as.integer(Init_seed)
  }


  #REVIEW Adjust rho for either only regional predictors nor network predictors
  if(!joint_include){ # if(!beta_include || !Alpha_include) {
    rho_update = FALSE
    rho = a_rho = b_rho = 0
  }
  #? Init_list
  sig2_i = as.double(sig2_i)
  MH_step = as.integer(MH_step)
  MH_stop = as.integer(MH_stop)
  MH_acptrate = as.double(MH_acptrate)
  burnin  = as.integer(burnin)
  M = as.integer(M)
  disp = as.logical(disp)

  #### <<< convert input into right mode  >>> ####


  #### >>> ANCHOR define variables <<< ####

  ## info double scaler
  ##      a_update, b_update;
  ##      diff_norm, double_old,
  ##      mu_update, tau2_update
  ##      n0, p0, v0,
  ##      rho0, sig2_V0, sig2_V1, sig2_e0, sig2_i0
  ##      rho_new, lambda_new, rho_ss, lambda_ss
  ##      ratio, acptrate, lu
  a_update = b_update = double(length = 1)
  diff_norm = double(length = 1) # non-negative
  double_old = double(length = 1)
  mu_update = double(length = 1)
  tau2_update = double(length = 1) # non-negative
  n0 = p0 = v0 = double(length = 1) # positive
  sig2_V0 = sig2_V1 = sig2_e0 = sig2_i0 = double(length = 1)
  # sig2_b0 = double(length = 1)
  rho0_new = rho0 = double(length = 1)
  rho_new = lambda_new = double(length = 1) # positive
  rho_ss = lambda_ss = double(length = 1) # positive
  acptrate = ratio = double(length = 1) # positive
  lu = double(length = 1) # negative


  ## info integer
  ##      n, p, q, q0, q1
  ##      MH_flag
  ##      j_id_n(>=0), gamma_nz(>=0)
  ##      i, j, r, m
  n = as.integer(nrow(Y))
  # NOTE feature region or network predictor
  if(beta_include) {
    p = as.integer(ncol(X))
  } else {
    # p = ifelse(Alpha_include, as.integer(dim(Z)[1]), 0)
    p = as.integer(dim(Z)[1])
  }
  q = as.integer(ifelse(missing(W), 0, ncol(W)))
  q0 = integer(length = 1)
  MH_flag = q1 = as.integer(1)
  gamma_nz = j_id_n = as.integer(0)
  i = j = r = m = integer(length = 1)

  # FIXME
  if(intercept) {
    q0 = as.integer(q + 1)
    if(missing(W)) {
      W = matrix(1, nrow = n, ncol = 1)
    } else {
      W = cbind2(1, W)
    }
  } else {
    q0 = q
  }
  q1 = max(q0, q1)


  ## info double vector
  ##      diff_vec(n), R(n)
  ##      diff2_vec(p), bMu_update(p)
  ##!FIXME      diff3_vec(q1), iMu_update(q1), W_V(q1)
  ##      threshold(2)
  ##      X_V(p), mvtnorm_Z_b(p), mvtnorm_Z_i(q1)
  R = diff_vec = vector(mode = "double", length = n)
  bMu_update = diff2_vec = vector(mode = "double", length = p)
  W_V = iMu_update = diff3_vec = vector(mode = "double", length = q1) #! length max(1, q0)
  threshold = vector(mode = "double", length = 2)
  # info added by Jan 10, 2022
  X_V = mvtnorm_Z_b = double(length = p)
  # FIXME length of mvtnorm_Z_i
  mvtnorm_Z_i = double(length = q1)


  ## info double matrix
  #// ##      bTau_update(p, p),
  #// ##!FIXME     iTau_update(q1, q1),
  ##      diff_mat(n, d)
  ##TODO      V(p, 2)
  ##      Xt(p, n), X0(p, p)
  ##!FIXME     W_Q(q1, n), Wt(q1, n)
  ##      I_p(p, p)
  # ##      Sig0(2,2)
  #// bTau_update = matrix(NA, nrow = p, ncol = p)
  #// mode(bTau_update) <- "double"
  #// iTau_update = matrix(NA, nrow = q1, ncol = q1)
  #// mode(iTau_update) <- "double"
  diff_mat = matrix(NA, nrow = n, ncol = d)
  mode(diff_mat) = "double"
  # TODO feature V
  V = matrix(NA, nrow = p, ncol = 2)
  mode(V) <- "double"
  colnames(V) = c("gamma", "theta")
  Xt = matrix(NA, nrow = p, ncol = n)
  mode(Xt) = "double"
  X0 = matrix(NA, nrow = p, ncol = p)
  mode(X0) = "double"
  W_Q = matrix(NA, nrow = q1, ncol = q1)
  mode(W_Q) = "double"
  Wt = matrix(NA, nrow = q1, ncol = n)
  mode(Wt) = "double"
  I_p = matrix(0, nrow = p, ncol = p)
  diag(I_p) = 1
  mode(I_p) = "double"
  #// Sig0 = matrix(NA, nrow = 2, ncol = 2)
  #// mode(Sig0) = "double"


  ## info double array
  ##      Z0[n, p-1, p]
  Z0 = array(NA, dim = c(n, p-1, p))
  mode(Z0) <- "double"


  ## info logic scalar
  ##      bool_old
  bool_old = logical(length = 1)

  ## info logic vector
  ##      gamma_mask(p), theta_mask(p)
  ##      gamma_mask_new(p), theta_mask_new(p)
  gamma_mask_new = gamma_mask = logical(length = p)
  theta_mask_new = theta_mask = logical(length = p)

  ## logic matrix
  ##      Theta_mask(p,d), Delta(p,d)
  ##      Theta_mask_new(p,d), Delta_new(p,d)
  Theta_mask_new = Theta_mask = matrix(FALSE, nrow = p, ncol = d)
  mode(Theta_mask_new) <- mode(Theta_mask) <- "logical"
  Delta_new = Delta = matrix(FALSE, nrow = p, ncol = d)
  mode(Delta_new) <- mode(Delta) <- "logical"


  ## coefficients
  ita = vector(mode = "double", length = q1) #FIXME
  beta = vector(mode = "double", length = p)
  gamma = vector(mode = "double", length = p)
  theta = vector(mode = "double", length = p)
  Theta = matrix(NA, nrow = p, ncol = d)
  mode(Theta) = "double"
  Alpha = matrix(NA, nrow = p, ncol = d)
  mode(Alpha) = "double"
  #FIXME sig2's - 5 x M vectors
  sig2_e = sig2_b = sig2_A = sig2_T = sig2_V = double(length = 1)
  a_e = b_e = a_b = b_b = a_A = b_A = a_T = b_T = a_V = b_V = double(length = 1)

  ##TODO output <--- pointer
  # ita_save q1 x M
  ita_save = matrix(NA, nrow = q1, ncol = M) #FIXME
  mode(ita_save) = "double"
  # gamma_save p x M
  gamma_save = matrix(NA, nrow = p, ncol = M)
  mode(gamma_save) = "double"
  # beta_save  p x M
  beta_save = matrix(NA, nrow = p, ncol = M)
  mode(beta_save) = "double"
  # Alpha_save p x d x M
  Alpha_save = array(NA, dim = c(p, d, M))
  mode(Alpha_save) = "double"
  # theta_save p x M
  theta_save = matrix(NA, nrow = p, ncol = M)
  mode(theta_save) = "double"
  # Theta_save p x d x M
  Theta_save = array(NA, dim = c(p, d, M))
  mode(Theta_save) = "double"
  #FIXME sig2's - 5 x M vectors
  sig2_e_save = vector(mode = "double", length = M)
  sig2_b_save = vector(mode = "double", length = M)
  sig2_A_save = vector(mode = "double", length = M)
  sig2_T_save = vector(mode = "double", length = M)
  sig2_V_save = vector(mode = "double", length = M)

  # rho and lambda
  # REVIEW length of rho_accept and lambda_accept
  #? M+1 or burnin+M+1
  ##      rho_accept(burnin + M + 1), lambda_accept(burnin + M + 1)
  #! actual value for rho_accept and lambda_rho start from the second index
  if(rho_update) {
    rho_save = vector(mode = "double", length = M)
    rhoss_save = vector(mode = "double", length = MH_stop%/% MH_step + 1)
    rho_accept = rep(0, burnin + M + 1)
  } else {
    rho_save = rho
    rhoss_save = 0
    rho_accept = 1
  }
  if(lambda_update) {
    lambda_save = vector(mode = "double", length = M)
    lambdass_save = vector(mode = "double", length = MH_stop%/% MH_step + 1)
    lambda_accept = rep(0, burnin + M + 1)
  } else {
    lambda_save = lambda
    lambdass_save = 0
    lambda_accept = 1
  }


  # elapse_time
  start_time = double(length = 1)
  end_time = double(length = 1)
  elapse_time = double(length = 1)

  #### <<< define variables >>> ####


  #### >>> ANCHOR Initialize <<< ####
  if(missing(Init_list)) {
    #* Case 3) -- totally random
    ## pre-fixed parameters
    #### Hyperparameters in IG distribution for variance
    a_e = b_e = 1 # epsilon
    a_b = b_b = 1 # beta
    a_A = b_A = 1 # Alpha
    a_T = b_T = 1 # Theta
    a_V = b_V = 1 # gamma and theta ---> #TODO feature V
    rho_ss = lambda_ss = 1

    ## set.seed
    if(disp) cat("seed for random initial values:", Init_seed, "\r\n")
    set.seed(Init_seed)

    sig2_e = rgamma(1, shape = a_e, rate = b_e)
    sig2_b = rgamma(1, shape = a_b, rate = b_b)
    sig2_A = rgamma(1, shape = a_A, rate = b_A)
    sig2_T = rgamma(1, shape = a_T, rate = b_T)
    sig2_V = rgamma(1, shape = a_V, rate = b_V)


    ## generate initial values
    if(rho_update) {
      rho = runif(1, min = a_rho, max = b_rho)
      # rho = a_rho
    }

    if(lambda_update) {
      lambda= runif(1, min = a_lambda, max = b_lambda)
      # lambda = a_lambda
    }

    if(q0 > 0) {
      # ita = rnorm(q0, mean = 0, sd = sig2_i)
      ita = rnorm(q0, mean = 0, sd = 1)
    } else {
      ita = 0
    }

    #// V = mvtnorm::rmvnorm(p, mean = rep(0, 2), sigma = sig2_V * matrix(c(1, rho, rho, 1), nrow = 2))
    #// colnames(V) <- c("gamma", "theta")
    #// beta = mvtnorm::rmvnorm(1, mean = V[, "gamma"], sigma = sig2_b * I_p)
    #// Alpha = t(mvtnorm::rmvnorm(d, mean = V[, "theta"], sigma = sig2_A * I_p))
    #// Theta = t(mvtnorm::rmvnorm(d, mean = V[, "theta"], sigma = sig2_T * I_p))

    V = matrix(rnorm(2*p), nrow = p)
    colnames(V) <- c("gamma", "theta")
    beta =  rnorm(p, mean = 0, sd = 1)
    Alpha = matrix(rnorm(p*d, mean = 0, sd = 1), nrow = p)
    Theta = matrix(rnorm(p*d, mean = 0, sd = 1), nrow = p)
  } else {
    a_e = Init_list$e["a"]
    b_e = Init_list$e["b"]
    sig2_e = Init_list$e["sig2"]
    a_b = Init_list$b["a"]
    b_b = Init_list$b["b"]
    sig2_b = Init_list$b["sig2"]
    a_A = Init_list$A["a"]
    b_A = Init_list$A["b"]
    sig2_A = Init_list$A["sig2"]
    a_T = Init_list$T["a"]
    b_T = Init_list$T["b"]
    sig2_T = Init_list$T["sig2"]
    a_V = Init_list$V["a"]
    b_V = Init_list$V["b"]
    sig2_V   = Init_list$V["sig2"]
    rho_ss = Init_list$rho["ss"]
    lambda_ss = Init_list$lambda["ss"]
    if(rho_update) rho = Init_list$rho["rho"]
    if(lambda_update) lambda = Init_list$lambda["lambda"]
    ita = Init_list$ita
    V = Init_list$V
    beta =  Init_list$beta
    Alpha = Init_list$Alpha
    Theta = Init_list$Theta
  }

  # REVIEW Adjust initial values for only Region or only connective case
  if (!Alpha_include) {
    #// rho_update = FALSE
    #// rho = 0
    #// rho_ss = 0
    Alpha = matrix(rep(0, p*d), nrow = p)
    Theta = matrix(rep(0, p*d), nrow = p)
    V[, "theta"] = rep(0, p)
    sig2_A = 0
    sig2_T = 0
  }

  if (!beta_include) {
    #// rho_update = FALSE
    #// rho = 0
    #// rho_ss = 0
    beta = rep(0, p)
    V[, "gamma"] = rep(0, p)
    sig2_b = 0
  }


  # TODO feature V
  gamma = V[, "gamma", drop = TRUE]
  theta = V[, "theta", drop = TRUE]
  rhoss_save[1] = rho_ss
  lambdass_save[1] = lambda_ss

  Initial_list = list(ita = ita,
                      beta = beta,
                      Alpha = Alpha,
                      Theta = Theta,
                      theta = theta,
                      gamma = gamma,
                      rho = rho,
                      lambda = lambda,
                      sig2_vec = c(e = sig2_e, beta = sig2_b, Alpha = sig2_A, Theta = sig2_T, V = sig2_V, sig2_i = sig2_i),
                      ss = c(rho = rho_ss, lambda = lambda_ss))

  #### <<< initial values >>> ####



  #### >>> ANCHOR helper variables <<< ####
  # variable masks - indicator
  gamma_mask = abs(gamma) > lambda
  theta_mask = abs(theta) > lambda
  Theta_mask = abs(Theta) > lambda
  Delta <- theta_mask * Theta_mask

  #* calculate current residual R(n,1)
  if (q0 > 0) {
    R = W %*% ita
  } else {
    R = matrix(0, nrow = n)
  }

  if (beta_include)
    R = R + X[, gamma_mask, drop = FALSE] %*% beta[gamma_mask]

  if (Alpha_include) {
    for (i in 1:n) {
      for (r in 1:d) {
        diff2_vec = Alpha[, r] * Delta[, r]
        R[i] <- R[i] + diff2_vec %*% Z[, , i] %*% diff2_vec
      }
      #* organize array Z into a n x(p-1) x p array Z0
      #! Z[, , i] is symmetric and with diagonal elements being 0s
      for (j in 1:p) {
        Z0[i, , j] <- Z[j, -j, i]
      }
    }
  }
  R = Y - R

  Delta <- Delta == 1
  gamma_nz <- length(which(gamma_mask))

  # pre-calculate some intermediate variables
  p0 = p/2
  n0 = n/2
  v0 = p*d/2
  sig2_i0 = 1/sig2_i

  # Wt: (q1)-by-n matrix;, i.e., transpose(\W)
  if (q0 > 0) {
    Wt = t(W)
    W_eigen = eigen(Wt %*% W) # FIXME
    W_Q = W_eigen$vectors # W1$vectors # q1-by-q1
    W_V = W_eigen$values  # W1$values  # q1 length
  }

  # Xt: p-by-n matrix; i.e., transpose(\X)
  if (beta_include) {
    Xt = t(X)
    X0 <- Xt %*% X
  }

  #NOTE initialize rho0 and threshold
  ############################################################
  ##* for fixed rho and lambda, put before MCMC run;
  ##* O.W, update rho0 and threshold at the end of MCMC after updating rho and lambda
  ############################################################
  rho0 = 1 - rho^2
  threshold = c(-lambda, lambda)
  #// Sig0 = toeplitz(c(1, -rho))


  #### <<< helper variables >>> ####

  #### >>> ANCHOR MCMC starts <<< ####
  start_time <- Sys.time()
  for (m in 1:(burnin + M)) {
    if (disp && m %% disp_step == 0) {
      cat("iter=", m, "---", sep = "", "\r\n")
    }

    # NOTE updating order
    ############################################################
    #* ita, gamma, beta, theta, Theta, Alpha,
    #* hyperparameters: sig2_e, sig2_b, sig2_A, sig2_T, sig2_A, sig2_V,
    #* rho, lambda, rho0, threshold
    #* If missing regional predictors
    #* SKIP gamma, beta, sig2_b, rho, rho0
    #* IF missing network predictors
    #* SKIP theta, Theta, Alpha, sig2_A, sig2_T, rho, rho0
    ############################################################

    # update some helping intermediate scalars
    sig2_V0 = 1/(sig2_V * rho0)
    sig2_V1 = rho * sig2_V0
    sig2_e0 = 2 * sig2_e
    # sig2_b0 = 1 / sig2_b

    #### >>> ANCHOR Update ita - vector of length q1 <<< ####
    if (disp) cat("ita;")
    if(q0 > 1) {
      # info simplified by eigen-decomposition
      diff3_vec = 1/(W_V/sig2_e + sig2_i0) #* sig2_i0 = 1/sig2_i
      diff_vec = R + W %*% ita
      iMu_update = colSums(W_Q * drop(Wt %*% diff_vec))/sig2_e
      iMu_update = iMu_update * diff3_vec
      iMu_update[] = W_Q %*% iMu_update

      #info update ita
      mvtnorm_Z_i = rnorm(q1)
      ita[] = W_Q %*% (colSums(W_Q * mvtnorm_Z_i) * sqrt(diff3_vec))
      ita = ita + iMu_update

      #info update R
      R = diff_vec - W %*% ita
    } else if (q0 == 1) {
      # FIXME q0 == 1 && !intercept
      #? replace n with drop(Wt %*% W)
      tau2_update = n/sig2_e + sig2_i0 # sig2_i0 = 1/sig2_i
      tau2_update = 1/tau2_update
      double_old = ita[1]
      mu_update = tau2_update * (sum(R) + n * double_old)
      mu_update = mu_update/sig2_e
      ita[1] = rnorm(1, mean = mu_update, sd = sqrt(tau2_update))
      #! update R
      R = R + (double_old - ita[1]) # element-wise
    }
    #### <<< Finish updating for ita >>> ####

    if (beta_include) {
      #### >>> ANCHOR Update gamma - vector of length p <<< ####
      # if (beta_include) {
        if (disp) cat("gamma;")
        # >>> Start loop for gamma <<<
        #info tau2_update is fixed across j=1,...,p
        tau2_update = 1/sig2_b + sig2_V0 #* sig2_V0 = 1/(sig2_V * (1 - rho^2))
        tau2_update = 1/tau2_update
        for (j in 1:p) {
          mu_update = beta[j]/sig2_b + sig2_V1 * theta[j] #* sig2_V1 = rho/(sig2_V * (1 - rho^2))
          mu_update = tau2_update * mu_update
          diff_vec = beta[j] * X[, j] # element-wise ---> vector of length n
          bool_old = gamma_mask[j]
          if (bool_old) {
            # Case I: |beta[j]| > lambda ---> falls in tail area
            a_update = sum(R^2)
            diff_vec = R + diff_vec
            b_update = sum(diff_vec^2)
          } else {
            # Case II: |beta[j]| <= lambda ---> falls in center area
            diff_vec = R - diff_vec
            a_update = sum(diff_vec^2)
            b_update = sum(R^2)
          }
          a_update = a_update/sig2_e0 #* sig2_e0 = 2 * sig2_e
          b_update = b_update/sig2_e0

          #! update gamma[j]
          gamma[j] = ThreeMixNorm(nm = mu_update, ns = sqrt(tau2_update), threshold, a = a_update, b = b_update)
          #! update gamma_mask[j]
          gamma_mask[j] = abs(gamma[j]) > lambda
          #! update R
          if(bool_old != gamma_mask[j]) {
            R = diff_vec
            gamma_nz = gamma_nz + ifelse(bool_old, -1, 1)
          } # else {R = R; gamma_nz = gamma_nz}
        }
        # <<< End loop for gamma >>>
      # }
      #### <<< Finish updating for gamma >>> ####


      #### >>> ANCHOR Update beta - vector of length p; i.e. \tilde(beta) ####
      # if (beta_include) {
        if (disp) cat("beta;")
        # >>> Start loop for beta <<<
        if (gamma_nz > 0) { # if (TRUE %in% gamma_mask) {
          # FIXME case when gamma_nz = 1
          # Case I: exit none-zero beta
          # info edited by Jan 10, 2022
          # REVIEW simplified by eigen decomposition
          # FIXME convert into Rcpp
          #! dimension and declaration of X_eigen and X_Q in Rcpp;
          #// bTau_update = 1/sig2_b * I_p # diag(rep(1/sig2_b, p)); sig2_b0 = 1/sig2_b
          #// bTau_update[gamma_mask, gamma_mask] = bTau_update[gamma_mask, gamma_mask] + X0[gamma_mask, gamma_mask]/sig2_e # X0 = t(X) %*% X
          #// bTau_update = solve(bTau_update)
          #// diff_vec = R + X[, gamma_mask, drop = FALSE] %*% beta[gamma_mask] # dimension -- vector of length n
          #// bMu_update[gamma_mask] = (Xt[gamma_mask, ] %*% diff_vec)/sig2_e + gamma[gamma_mask]/sig2_b
          #// bMu_update[!gamma_mask] = gamma[!gamma_mask]/sig2_b
          #// bMu_update = bTau_update %*% bMu_update
          X_eigen = eigen(X0[gamma_mask, gamma_mask])
          X_Q = X_eigen$vectors
          # X_V = X_eigen$values
          # X_V = 1/(X_V/sig2_e + 1/sig2_b)
          X_V[1:gamma_nz] = 1/(X_eigen$values/sig2_e + 1/sig2_b) # X_V > 0 always true
          diff_vec = R + X[, gamma_mask, drop = FALSE] %*% beta[gamma_mask] # dimension -- vector of length n
          bMu_update[1:gamma_nz] = colSums(X_Q * (drop(Xt[gamma_mask, ] %*% diff_vec)/sig2_e + gamma[gamma_mask]/sig2_b))
          bMu_update[1:gamma_nz] = bMu_update[1:gamma_nz] * X_V[1:gamma_nz]
          bMu_update[gamma_mask] = X_Q %*% bMu_update[1:gamma_nz]
          bMu_update[!gamma_mask] = gamma[!gamma_mask]

          #info update beta
          #// beta =  mvtnorm::rmvnorm(1, mean = bMu_update, sigma = bTau_update)
          mvtnorm_Z_b = rnorm(p)
          beta[gamma_mask] = X_Q %*% (colSums(X_Q * mvtnorm_Z_b[gamma_mask]) * sqrt(X_V[1:gamma_nz]))
          beta[!gamma_mask] = mvtnorm_Z_b[!gamma_mask] * sqrt(sig2_b)
          beta = beta + bMu_update

          #info update R
          R = diff_vec - X[, gamma_mask, drop = FALSE] %*% beta[gamma_mask] # dimension -- vector of length n
        } else {
          # Case II: all beta's are zeros
          #* bTau_update = sig2_b * I_p # diag(rep(sig2_b, p))
          #* bMu_update = gamma
          #info update beta
          #// beta = mvtnorm::rmvnorm(1, mean = gamma, sigma = sig2_b * I_p)
          beta = rnorm(p, mean = 0, sd = sqrt(sig2_b)) + gamma

          #info update R
          # R = R
        }
        # <<< End loop for beta >>>
      # }
      #### <<< Finish updating for beta >>> ####
    }

    if (Alpha_include) {
      #### >>> ANCHOR Update theta - vector of length p; i.e., \theta <<< ####
      # if (Alpha_include) {
        if (disp) cat("theta;")
        # >>> Start loop for theta <<<
        #* Calculate mean and variance of truncated normal distribution
        #* variance keeps the same across j=1,...,p
        tau2_update = d/sig2_A + d/sig2_T + sig2_V0 # sig2_V0 = 1/(sig2_V * (1 - rho^2))
        tau2_update = 1/tau2_update

        for (j in 1:p) {
          mu_update = sum(Alpha[j, ])/sig2_A + sum(Theta[j, ])/sig2_T + sig2_V1 * gamma[j] # sig2_V1 = rho/(sig2_V * (1 - rho^2))
          mu_update = tau2_update * mu_update

          if (TRUE %in% Theta_mask[j, ] && TRUE %in% Delta[-j, ]) { # O.W. a_update = b_update = 1 very race case

            diff_mat <- Z0[, , j] %*% (Delta[-j, ] * Alpha[-j, ]) # element-wise ---> n x d matrix
            diff_vec <- diff_mat[, Theta_mask[j, ], drop = FALSE] %*% Alpha[j, Theta_mask[j, ]]
            diff_vec <- 2 * diff_vec #! symmetric property brings the factor 2

            bool_old <- theta_mask[j]
            if (bool_old) {
              #* Case I: |theta[j]| > lambda
              a_update = sum(R^2) # t(R) %*% R
              diff_vec = R + diff_vec
              b_update = sum(diff_vec^2) # t(diff_vec) %*% diff_vec
            } else {
              #* Case II: |theta[j]| <= lambda
              diff_vec = R - diff_vec
              a_update = sum(diff_vec^2) # t(diff_vec) %*% diff_vec
              b_update = sum(R^2) # t(R) %*% R
            }
            a_update = a_update / sig2_e0
            b_update = b_update / sig2_e0

            #! sample theta[j]
            theta[j] = ThreeMixNorm(nm = mu_update, ns = sqrt(tau2_update), threshold, a = a_update, b = b_update)
            #! update theta_mask[j]
            theta_mask[j] = abs(theta[j]) > lambda
            #! update R
            #* if theta_mask[j](new) != heta_mask[j](old) ---> 1) case I: R = R+__; case II: R = R-__; else ---> R = R
            if (theta_mask[j] != bool_old) {
              R = diff_vec
            }
          } else {
            #* Theta_mask[j, ] == (FALSE, ..., FALSE) or Delta[-j, ] = (FALSE, ..., FALSE) ---> diff_vec = 0 ---> a_update = b_update = t(R) %*% R; R = R
            #! sample theta[j]
            theta[j] = ThreeMixNorm(nm = mu_update, ns = sqrt(tau2_update), threshold)
            #! update theta_mask[j]
            theta_mask[j] = abs(theta[j]) > lambda
            #! update R
            # R = R
          }

          #! update Delta
          Delta[j, ] = as.logical(Theta_mask[j, ] * theta_mask[j])
        }
        # <<< End loop for theta >>>
      # }
      #### <<< Finish updating for theta <<< ####


      #### >>> ANCHOR Update Theta - matrix of size p-by-d; i.e., [\theta^1, \theta^d] <<< ####
      # if (Alpha_include) {
        if (disp) cat("Theta;")
        # >>> Start loop for Theta <<<
        for (r in 1:d) {
          for (j in 1:p) {
            if(theta_mask[j]) {
              # |theta[j]| > lambda
              # REVIEW dimension drop problem
              #* diff_vec = 2 * Alpha[j, r] * Z0[, Delta[-j, r], j] %*% Alpha[-j, r][Delta[-j, r]]
              j_id_n <- length(which(Delta[-j, r]))
              if(j_id_n != 0) {
                # diff_vec = `dim<-` (Z0[, Delta[-j, r], j], c(n, j_id_n)) %*% Alpha[-j, r][Delta[-j, r]]
                if(j_id_n > 1) {
                  diff_vec =  Z0[, Delta[-j, r], j] %*%  Alpha[-j, r][Delta[-j, r]]
                } else {
                  #FIXME - diff_vec[, 1] --> diff_vec
                  diff_vec = Z0[, Delta[-j, r], j] * Alpha[-j, r][Delta[-j, r]]
                }
                diff_vec = (2 * Alpha[j, r]) * diff_vec

                bool_old = Theta_mask[j, r]
                if (bool_old) {
                  #* Case I: |Theta[j, r]| > lambda
                  a_update = sum(R^2) # t(R) %*% R
                  diff_vec = R + diff_vec
                  b_update = sum(diff_vec^2) # t(diff_vec) %*% diff_vec
                } else {
                  #* Case II: |Theta[j, r]| <= lambda
                  diff_vec = R - diff_vec
                  a_update = sum(diff_vec^2) # t(diff_vec) %*% diff_vec
                  b_update = sum(R^2) # t(R) %*% R
                }

                a_update = a_update / sig2_e0 # sig2_e0 = 2 * sig2_e
                b_update = b_update / sig2_e0

                #! sample Theta[j, r]
                Theta[j, r] = ThreeMixNorm(nm = theta[j], ns = sqrt(sig2_T), threshold, a = a_update, b = b_update)
              } else {
                # a_update = 1, b_update = 1, diff_vec = R
                #! sample Theta[j, r]
                Theta[j, r] = ThreeMixNorm(nm = theta[j], ns = sqrt(sig2_T), threshold)
              }

              #! update Theta_mask
              Theta_mask[j, r] = abs(Theta[j, r]) > lambda
              #! update Delta given theta_mask[j](old) = TRUE
              # Delta[j, r] = Theta_mask[j, r] # && theta_mask[j] <--- theta_mask[j] = TRUE
              #! update R given theta_mask[j](old) = TRUE
              #* if j_id_n == 0 ---> R = R
              #* if Theta_mask[j,r](new) != Theta_mask[j,r](old) ---> 1) case I: R = R+__; case II: R = R-__; else ---> R = R
              if (Theta_mask[j, r] != bool_old) {
                Delta[j, r] = Theta_mask[j, r]
                if (j_id_n != 0) R = diff_vec
              }
            } else {
              #* |theta[j]| <= lambda ---> a_update = b_update = t(R) %*% R / sig2_e; R = R
              #* if a_update = b_update ---> ThreeMixNorm put 1/3 weight to each component
              #! sample Theta[j, r]
              Theta[j, r] = ThreeMixNorm(nm = theta[j], ns = sqrt(sig2_T), threshold)
              # REVIEW update Theta_mask
              Theta_mask[j, r] = abs(Theta[j, r]) > lambda
              #! update Delta given theta_mask[j](old) = FALSE
              # Delta_update[j, r] = Delta_old[j, r] = FALSE <--- Delta[j, r] = FALSE = theta_mask[j] && Theta_mask[j, r] <--- theta_mask[j] = FALSE
              # Delta[j, r] = FALSE
              #! update R given theta_mask[j](old) = FALSE
              # R = R
            }
          }
        }
        # <<< End loop for Theta >>>
      # }
      #### <<< Finish updating for Theta >>> ####


      #### >>> ANCHOR Update Alpha - matrix of size p-by-d; i.e., [\alpha^1,  ..., \alpha^R] <<< ####
      # if (Alpha_include) {
        if (disp) cat("Alpha;")
        # >>> Start loop for Alpha <<<
        for (r in 1:d) {
          for (j in 1:p) {
            j_id_n <- length(which(Delta[-j, r]))
            if (Delta[j, r] && j_id_n > 0) {
              # if (Delta[j, r] && TRUE %in% Delta[-j, r]) {
              #* Case I: none-zero effect on posterior mean and variance
              double_old = Alpha[j, r]
              # REVIEW dimension drop
              # diff_vec = `dim<-` (Z0[, Delta[-j, r], j], c(n, j_id_n)) %*% Alpha[-j, r][Delta[-j, r]]
              # diff_vec = 2 * diff_vec #! symmetric property brings the factor 2
              if(j_id_n > 1) {
                diff_vec = Z0[, Delta[-j, r], j] %*% Alpha[-j, r][Delta[-j, r]]
              } else {
                #FIXME - diff_vec[, 1] --> diff_vec
                diff_vec = Z0[, Delta[-j, r], j] * Alpha[-j, r][Delta[-j, r]]
              }
              diff_vec = 2 * diff_vec #! symmetric property brings the factor 2
              diff_norm = sum(diff_vec^2) # t(diff_vec) %*% diff_vec
              # variance
              b_update = diff_norm/sig2_e + 1/sig2_A
              b_update = 1/b_update
              # mean: (R + Alpha[j,r]*diff_vec) %*% diff_vec
              a_update = t(R) %*% diff_vec + double_old * diff_norm
              a_update = b_update * (a_update/sig2_e + theta[j]/sig2_A)

              #! sample Alpha[j, r]
              Alpha[j, r] = rnorm(1, mean = a_update, sd = sqrt(b_update))
              #!update R
              R = R + (double_old - Alpha[j, r]) * diff_vec
            } else {
              #* case II: posterior distribution of Alpha[j,r] did not change; R = R
              #! sample Alpha[j, r]
              Alpha[j, r] = rnorm(1, mean = theta[j], sd = sqrt(sig2_A))
              #! update R
              # R = R
            }
          }
        }
        # <<< End loop for Alpha >>>
      # }
      #### <<< Finish updating for Alpha >>> ####
    }

    #### >>> ANCHOR Update hyperparameters <<< ####
    if (disp) cat("hyperparameters;") # DEBUG
    # >>> Update sig2_e <<<
    a_update = a_e + n0 #* n0 = n/2
    diff_norm = sum(R^2)/2
    b_update = b_e + diff_norm
    sig2_e = invgamma::rinvgamma(1, shape = a_update, rate = b_update)
    # <<< Finish sig2_e >>>


    if (beta_include) {
      # >>> Update sig2_b <<<
      a_update = a_b + p0 #* p0 = p/2
      diff2_vec = beta - gamma
      b_update = b_b + sum(diff2_vec^2)/2
      sig2_b = invgamma::rinvgamma(1, shape = a_update, rate = b_update)
      # <<< Finish sig2_b >>>
    }

    if (Alpha_include) {
      # >>> Update sig2_A <<<
      a_update = a_A + v0 #* v0 = p*d/2
      b_update = b_A + sum((Alpha - theta)^2)/2
      sig2_A = invgamma::rinvgamma(1, shape = a_update, rate = b_update)
      # <<< Finish sig2_A >>>


      # >>> Update sig2_T <<<
      a_update = a_T + v0
      b_update = b_T + sum((Theta - theta)^2)/2
      sig2_T = invgamma::rinvgamma(1, shape = a_update, rate = b_update)
      # <<< Finish sig2_T >>>
    }


    #// if (sig2V_include) {
    # >>> Update sig2_V <<<
    #NOTE exclude the case when there is no need to update sig2_V, which is trivia case for missing both X and Z
    #// a_update = a_V + p
    a_update = a_V + ifelse(joint_include, p, p0) # p0 = p/2
    # NOTE Direct updating b_update as the comment out line is faster here. But mu_update and tau2_update could be re-used for rho_update
    #// b_update = b_V + sum((V %*% Sig0) * V) / (2 * rho0) # element-wise
    # TODO feature V --> combine gamma and theta into matrix form
    #// V = cbind(gamma, theta)
    #// tau2_update = sum(V^2)
    #// mu_update = drop(V[, 1] %*% V[, 2])

    #// tau2_update = sum(gamma^2) + sum(theta^2)
    #// mu_update = sum(gamma * theta)
    #// double_old = tau2_update / (2 * rho0) - mu_update * rho / rho0
    tau2_update = ifelse(beta_include, sum(gamma^2), 0) + ifelse(Alpha_include, sum(theta^2), 0)
    if(joint_include) {
      mu_update = sum(gamma * theta)
      double_old = tau2_update / (2 * rho0) - mu_update * rho / rho0
    } else {
      double_old = tau2_update / 2
    }
    b_update = b_V + double_old
    sig2_V = invgamma::rinvgamma(1, shape = a_update, rate = b_update)
    # <<< Finish sig2_V >>>
    # }
    #### <<< Finish updating hyperparameters >>> ####


    #### >>> update rho and lambda by MH sampling <<< ####
    # ANCHOR rho and lambda

    if(rho_update) {
      if (disp) cat("rho \r\n")
      # >>> Update rho <<<
      # proposed rho
      rho_new = truncnorm::rtruncnorm(n = 1, a = a_rho, b = b_rho, mean = rho, sd = rho_ss)
      rho0_new = 1 - rho_new^2
      double_new = tau2_update / (2 * rho0_new) - mu_update * rho_new / rho0_new #* log(f(x')*sig2_V)
      # accept/reject rho
      ratio = (double_old - double_new) / sig2_V #* double_old = - double_old / sig2_V
      ratio = ratio + p0 * log(rho0/rho0_new)
      if(ratio > 0) {
        rho_accept[m + 1] = 1
        rho = rho_new
        rho0 = rho0_new
      } else {
        lu = log(runif(1, 0, 1))
        if(ratio > lu) {
          rho_accept[m + 1] = 1
          rho = rho_new
          rho0 = rho0_new
        }
      }
      # <<< Finish rho >>>
    }


    if(lambda_update) {
      if (disp) cat("lambda \r\n")
      # >>> Update lambda <<<
      # proposed lambda
      lambda_new = truncnorm::rtruncnorm(n = 1, a = a_lambda, b = b_lambda, mean = lambda, sd = lambda_ss)
      #// gamma_mask_new = abs(gamma) > lambda_new
      #// theta_mask_new = abs(theta) > lambda_new
      #// Theta_mask_new = abs(Theta) > lambda_new
      #// Delta_new <- theta_mask_new * Theta_mask_new
      if(q0 > 0) {
        diff_vec = W %*% ita
      } else {
        diff_vec = matrix(0, nrow = n)
      }
      if(beta_include) {
        gamma_mask_new = abs(gamma) > lambda_new
        diff_vec = diff_vec + X[, gamma_mask_new, drop = FALSE] %*% beta[gamma_mask_new]
      }
      if(Alpha_include) {
        theta_mask_new = abs(theta) > lambda_new
        Theta_mask_new = abs(Theta) > lambda_new
        Delta_new <- theta_mask_new * Theta_mask_new
        for (i in 1:n) {
          for (r in 1:d) {
            diff2_vec = Alpha[, r] * Delta_new[, r]
            diff_vec[i] <- diff_vec[i] + diff2_vec %*% Z[, , i] %*% diff2_vec
          }
        }
      }

      # accept/reject lambda
      diff_vec = Y - diff_vec
      ratio = diff_norm - sum(diff_vec^2)/2 #* log(f(x')/f(x_t)*sig2_e)
      ratio = ratio / sig2_e
      if(ratio > 0) {
        lambda_accept[m + 1] = 1
        lambda = lambda_new
        #// gamma_mask = gamma_mask_new
        #// theta_mask = theta_mask_new
        #// Theta_mask = Theta_mask_new
        #// Delta = Delta_new == 1 #! convert integer Delta to logic
        #// gamma_nz = length(which(gamma_mask_new))
        if(beta_include) {
          gamma_mask = gamma_mask_new
          gamma_nz = length(which(gamma_mask_new))
        }
        if(Alpha_include) {
          theta_mask = theta_mask_new
          Theta_mask = Theta_mask_new
          Delta = Delta_new == 1 #! convert integer Delta to logic
        }
        threshold = c(-lambda_new, lambda_new)
        R = diff_vec
      } else {
        lu = log(runif(1, 0, 1))
        if(ratio > lu) {
          lambda_accept[m + 1] = 1
          lambda = lambda_new
          #// gamma_mask = gamma_mask_new
          #// theta_mask = theta_mask_new
          #// Theta_mask = Theta_mask_new
          #// Delta = Delta_new == 1 #! convert integer Delta to logic
          #// gamma_nz = length(which(gamma_mask_new))
          if(beta_include) {
            gamma_mask = gamma_mask_new
            gamma_nz = length(which(gamma_mask_new))
          }
          if(Alpha_include) {
            theta_mask = theta_mask_new
            Theta_mask = Theta_mask_new
            Delta = Delta_new == 1 #! convert integer Delta to logic
          }
          threshold = c(-lambda_new, lambda_new)
          R = diff_vec
        }
      }
      # <<< Finish lambda >>>
    }
    cat("\r\n")

    # ANCHOR adjust accept rate in MH sampling
    if(m < MH_stop && m %% MH_step == 0)  { # MH_stop << burnin + M
      if(disp) cat("Update rho_ss/lambda_ss at iter:", m, "\r\n")

      if(rho_update) {
        acptrate = mean(rho_accept[(MH_flag+1):(m+1)]) # mean(rho_accept[1 + MH_flag:m])
        rho_ss = adjust_acceptance(acptrate, rho_ss, MH_acptrate)
        rhoss_save[m %/% MH_step + 1] = rho_ss
        if(disp) { # DEBUG
          cat("rho update: rho=", round(rho, 3),
              ", rho_accept=", round(acptrate, 3),
              ", rho_ss=", round(rho_ss, 3), "\r\n", sep = "")
        }
        rhoss_save[m %/% MH_step + 1] = rho_ss
      }

      if (lambda_update) {
        acptrate = mean(lambda_accept[(MH_flag+1):(m+1)]) # mean(lambda_accept[1 + MH_flag:m])
        lambda_ss = adjust_acceptance(acptrate, lambda_ss, MH_acptrate)
        if(disp) { # DEBUG
          cat("lambda update: lambda=", round(lambda, 3),
              ", lambda_accept=", round(acptrate, 3),
              ", lambda_ss=", round(lambda_ss, 3), "\r\n", sep = "")
        }
        lambdass_save[m %/% MH_step + 1] = lambda_ss
      }

      MH_flag = m + 1
    }
    #### <<< finish updating rho and lambda >>> ####



    #### >>> ANCHOR Retain values for posterior samplings <<< ####
    if (m > burnin) {
      ita_save[, m - burnin]  <- ita
      beta_save[, m - burnin] <- beta
      gamma_save[, m - burnin] <- gamma
      theta_save[, m - burnin] <- theta
      Alpha_save[, , m - burnin] <- Alpha
      Theta_save[, , m - burnin] <- Theta
      # FIXME sig2_' vectors
      sig2_e_save[m - burnin] <- sig2_e
      sig2_b_save[m - burnin] <- sig2_b
      sig2_A_save[m - burnin] <- sig2_A
      sig2_T_save[m - burnin] <- sig2_T
      sig2_V_save[m - burnin] <- sig2_V
      if(rho_update) rho_save[m - burnin] = rho
      if(lambda_update) lambda_save[m - burnin] = lambda
    }
    #### <<< Finish saving samples >>>####

  }
  #### <<< MCMC ends >>> ###

  #### >>> ANCHOR output <<< ####
  end_time <- Sys.time()
  elapse_time = end_time - start_time

  result = list(ita_save = ita_save,
                beta_save = beta_save,
                gamma_save = gamma_save,
                Alpha_save = Alpha_save,
                theta_save = theta_save,
                Theta_save = Theta_save,
                sig2_e_save = sig2_e_save,
                sig2_b_save = sig2_b_save,
                sig2_A_save = sig2_A_save,
                sig2_T_save = sig2_T_save,
                sig2_V_save = sig2_V_save,
                rho_save = rho_save,
                lambda_save = lambda_save)

  if(rho_update) result = c(result,
                            list(rhoss_save = rhoss_save,
                                 rhoacpt_save = rho_accept[(burnin+2):(burnin+M+1)]))

  if(lambda_update) result = c(result,
                               list(lambdass_save = lambdass_save,
                                    lambdaacpt_save = lambda_accept[(burnin+2):(burnin+M+1)]))

  result = c(result,
             list(elapse_time = elapse_time,
                  initial_val = Initial_list))

  return(result)

}


ThreeMixNorm <- function(nm, ns, threshold, a = 1, b = 1) {
  #'@param input
  # nm, ns mean and sd for truncated normal distribution
  # threshold is a vector of length 2, where threshold[1] is the lower boundary, and threshold[2] is the upper boundary

  #'@param output
  # y: double numeric

  ## define variables
  # w: vector of length 3
  # x: integer vector of length 3

  ## Calculate mixture proportions
  w <- rep(0, 3)
  # TODO simplified
  w[1] <- pnorm((threshold[1] - nm)/ns)  # pnorm(threshold[1], mean = nm, sd = ns)
  w[2] <- pnorm((threshold[2] - nm)/ns)
  w[3] <- 1 - w[2]
  w[2] <- w[2] - w[1]

  if(a!=b) {
    if (w[2] == 1) {
      x = 2
    } else if (w[2] == 0) {
      x = rbinom(1, 1, w[1])
      if (x == 0) x = 3
    } else {
      w[2] = exp(a-b) * w[2]
      if(is.infinite(w[2])) {
        x = 2
      } else if (w[2] == 0) {
        x = rbinom(1, 1, w[1]/(w[1]+w[3]))
        if (x == 0) x = 3
      } else {
        w = w / sum(w)
        x <- rmultinom(n = 1, size = 1, prob = w)
        x <- which(x == 1)
      }
    }
  } else {
    x <- rmultinom(n = 1, size = 1, prob = w)
    x <- which(x == 1)
  }
  # <<< added >>>

  y <- switch(x,
              truncnorm::rtruncnorm(n = 1, b = threshold[1], mean = nm, sd = ns),
              truncnorm::rtruncnorm(n = 1, a = threshold[1], b = threshold[2], mean = nm, sd = ns),
              truncnorm::rtruncnorm(n = 1, a = threshold[2], mean = nm, sd = ns)
  )
  return(y)
}


Error_Initial <- function(data_list, sd_error) {
  #! ignore ita
  q <- ncol(Data_list$W)
  q <- ifelse(is.null(q), 0, q)
  p <- nrow(data_list$Alpha)
  d <- ncol(data_list$Alpha)
  ita = data_list$ita + rnorm(q, 0, sd_error)
  beta = data_list$beta + rnorm(p, 0, sd_error)
  Alpha = data_list$Alpha + matrix(rnorm(p*d, 0, sd_error), nrow = p)
  Theta = data_list$Theta + matrix(rnorm(p*d, 0, sd_error), nrow = p)
  V = data_list$V + matrix(rnorm(p*2, 0, sd_error), nrow = p)
  n = length(data_list$sig2_vec)
  sig2_vec = data_list$sig2_vec + rnorm(n, 0, sd_error)
  for (i in 1:n) {
    while(sig2_vec[i] <= 0) data_list$sig2_vec[i] + rnorm(1, 0, sd_error)
  }
  result <- list(ita = c(data_list$b0, ita),
                 beta = beta, Alpha = Alpha, Theta = Theta, V = V,  sig2_vec = sig2_vec)
  return(result)
}

# function used to adjust the acceptance rate in the MH algorithm
adjust_acceptance <- function(accept, sgm, target) {
  y = 1 + 1000 * (accept - target)^3
  if (y < 0.9) {
    y = 0.9
  } else if (y > 1.1) {
    y = 1.1
  }
  sgm = y * sgm
  return(sgm)
}


subnet_pp <- function(X, decision = 0.5,
                      graph = FALSE, layout = on_grid()) {
  # NOTE selection summary of marginal sub-networks
  # ======================================
  #'@param input
  # X: array with 3 or 4 dimensions, where
  #    dim[1] and dim[2] is the p-by-p symmetric coefficient matrix for each sub-network predictor
  #    dim[3] indicate #(sub-networks) and dim[4] indicate length of MCMC chain if applied.
  # decision: cutoff value for MMP for edge selection in each sub-network.
  # graph: boolean to indicate whether to generate igraph data
  # layout: layout for igraph

  #'@param output
  # result: a list of MMP, selected edges and igraph data
  # ======================================

  if(dim(X)[1] != dim(X)[2])
    stop("Wrong input of X.
            Need to be the coefficient matrix with p-by-p dimension")
  # p-by-p-by-d-by-M --> p-by-p-by-d
  if(length(dim(X)) == 4) X = rowMeans(X, dims = 3)

  d = dim(X)[3]
  X_index <- as.list(seq(d))
  names(X_index) <- paste0("Subnet", seq(d))
  for(r in seq(d)) {
    X_index[[r]] <- as.data.frame(which(X[, , r] > decision, arr.ind = TRUE)) %>%
      dplyr::rename(from = row, to = col) %>%
      dplyr::filter(from > to) %>%
      dplyr::mutate(from = as.character(from),
                    to = as.character(to)) %>%
      dplyr::arrange(from, to)
  }

  result <- list()
  result$MMP <- X
  result$edges <- X_index

  if(graph) {
    X_graph <- as.list(seq(d+1))
    names(X_graph) <- c(names(X_index), "Whole")

    for(r in seq(d)) {
      X_graph[[r]] <- igraph::graph.data.frame(X_index[[r]], directed = FALSE)
    }
    X_graph[[d+1]] <- igraph::graph.data.frame(do.call(rbind, X_index) %>% arrange(from, to),
                                               directed = FALSE)
    coords <- igraph::layout_(graph = X_graph[[d+1]], layout = layout)

    for(r in seq(d)) {
      sub_coords = coords
      a <- which(!attributes(igraph::V(X_graph[[d+1]]))$names %in% attributes(igraph::V(X_graph[[r]]))$names)
      if(length(a) > 0) sub_coords = sub_coords[-a, ]
      X_graph[[r]] <- list(g = X_graph[[r]], coords = sub_coords)
      # plot(X_graph[[r]]$g, layout = X_graph[[r]]$coords,
      #      vertex.label.dist = 3, main = paste("Subnectwork", r, sep = " "))
    }
    result$graph = X_graph
  }
  return(result)
}
