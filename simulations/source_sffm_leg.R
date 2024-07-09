#' Gibbs sampling step for cumulative shrinkage processes
#' 
#' Sample the parameters of a cumulative shrinkage process (CSP)
#' for a Bayesian factor model. 
#' 
#' @details The loadings \code{Lambda} are distributed 
#' \code{Lambda[j,h]} ~ N(\code{0}, \code{theta[h]})
#' 
#' @param Lambda (p x H) matrix of loadings
#' @param z (H x 1) vector of categorical variables
#' @param alpha prior expectation for number of factors
#' @param theta_inf variance defined by the spike (near zero)
#' @param a_theta inverse-gamma shape parameter for the slab
#' @param b_theta inverse-gamma rate parameter for the slab
#' @return a list with the following elements:
#' \itemize{
#' \item \code{theta} (H x 1) vector of variance parameters
#' \item \code{z} (H x 1) vector of categorical variables
#' \item \code{pi} (H x 1) vector of probability that component h is in the spike (near zero)
#' \item \code{H_star} number of components belonging to the slab (i.e., not near zero)
#' \item \code{nu} (H x 1) vector of stick-breaking probabilities
#' \item \code{omega} (H x 1) vector of probability weights for \code{z}
#' }
#' @importFrom fGarch dstd
#' @importFrom stats dnorm rgamma sample
sampleCSP = function(Lambda,
                     z,
                     alpha = 5,
                     theta_inf = 0.05,
                     a_theta = 2,
                     b_theta = 2){
  
  # Dimensions from Lambda:
  p = nrow(Lambda); H = ncol(Lambda)
  
  # Sample the stick probabilities:
  nu = rep(0, H)
  nu[1:(H-1)] = sapply(1:(H-1), function(ell)
    rbeta(n = 1,
          shape1 = 1 + sum(z==ell),
          shape2 = alpha + sum(z > ell))
  )
  nu[H] = 1
  
  # Update the weights:
  omega = rep(0, H)
  omega[1] = nu[1]; 
  omega[2:H] = sapply(2:H, function(ell) {
    nu[ell]*prod(1 - nu[1:(ell-1)])
  })
  
  # Sample the categorical variables
  prob.ell = rep(0, H) # Pr(z[h] == ell) for ell= 1,...,H
  for(h in 1:(H-1)){
    # ell <= h
    prob.ell[1:h] = omega[1:h]*exp(sum(dnorm(Lambda[,h],
                                             mean = 0,
                                             sd = sqrt(theta_inf),
                                             log = TRUE)))
    # ell > h
    prob.ell[-(1:h)] = omega[-(1:h)]*exp(sum(dstd(Lambda[,h],
                                                  mean = 0,
                                                  sd = sqrt(b_theta/a_theta),
                                                  nu = 2*a_theta,
                                                  log = TRUE)))
    # Sample the latent variables:
    z[h] = sample(1:H, size = 1, prob = prob.ell)
  }
  # For h=H, the likelihood terms cancel:
  z[H] = sample(1:H, size = 1, prob = omega)
  
  
  # Sample theta:
  theta = sapply(1:H, function(h){
    if(z[h] <= h) {
      theta_inf
    } else {
      1/rgamma(n = 1, 
               shape = a_theta + p/2,
               rate = b_theta + 1/2*sum(Lambda[,h]^2))
    }
  })
  
  # Now return the parameters of interest:
  list(
    theta = theta, # variance parameter
    z = z,  # Categorical variable
    pi = cumsum(omega), # Probability that component h is in the spike (near zero)
    H_star = sum(z > 1:H), # Number of components belonging to the slab (i.e., not near zero)
    nu = nu, # Stick probabilities (for completeness)
    omega = omega # Weights for z (for completeness)
  )
}

# SFFM using Legramanti prior (no PX)
sffm_leg = function(Y, tau, g, K = NULL, 
                    log_prior_gamma = NULL,
                    K_0 = 5, 
                    theta_inf = 0.05, a_theta = 2, b_theta = 2,
                    gamma_init = NULL,
                    nsave = 3000, nburn = 1000, nskip = 2,
                    mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'K_0', 'Yhat', 'Ypred'),
                    Con_mat = NULL){
  # K = NULL; log_prior_gamma = NULL;
  # K_0 = 5; theta_inf = 0.05; a_theta = 2; b_theta = 2; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'Yhat', 'Ypred'); Con_mat = NULL
  #----------------------------------------------------------------------------
  # Initializations and checks for errors:
  #----------------------------------------------------------------------------
  # Convert tau to matrix, if necessary:
  tau = as.matrix(tau)
  
  # Compute the dimensions:
  n = nrow(Y);  # Number of curves
  m = ncol(Y);  # Number of observation points
  d = ncol(tau) # Dimension of the observation points
  
  # Check tau, Y dimensions:
  if(nrow(tau) != m)
    stop('nrow(tau) must be equal to ncol(Y)')
  
  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))
  
  # Rescale by observation SD (and correct parameters later):
  sdY = sd(Y, na.rm=TRUE);
  Y = Y/sdY;
  
  # Check for missingness and crudely impute:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing){
    # Indices of missing values:
    na.ind = which(is.na(Yna), arr.ind = TRUE); 
    
    # Impute using FPCA:
    Y = fdlm_init(Y, tau)$Y0
    
    # Impute using linear interpolation:
    # Y = t(apply(Y, 1, function(y) 
    #   approxfun(x = tau, y=y, rule=2)(tau)))
  }
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
  # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try) == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    g_p = function(tau, gamma) g(tau)
    
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    g_p = g
  }
  #----------------------------------------------------------------------------
  # Initialize the parametric component:
  #----------------------------------------------------------------------------
  if(is_unknown_gamma){
    # Initialize gamma, using given value or randomly:
    if(!is.null(gamma_init)){
      gamma = gamma_init
    } else gamma = runif(n = 1);
    
    # But make sure it has positive probability:
    counter = 0
    while(is.infinite(log_prior_gamma(gamma)) && counter < 1000) {
      gamma = runif(n = 1); counter = counter + 1
    } 
    if(counter == 1000) stop('Problem initializing gamma; try modifying the prior to have support on [0,1]')
    
    # And initialize the basis matrix:
    Gmat = g_p(tau, gamma)
  } else Gmat = g_p(tau)
  
  # And initialize the coefficients:
  Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  
  # Number of parametric terms:
  L = ncol(Gmat)
  #----------------------------------------------------------------------------
  # Initialize the nonparametric component:
  #----------------------------------------------------------------------------
  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Y - tcrossprod(Alpha, Gmat), tau, K); 
  Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo
  
  # Nonparmetric function matrix:
  Fmat = splineInfo$Bmat%*%Psi
  
  # Number of nonparametric terms:
  K = ncol(Fmat)
  
  # Check: expected number of active factors should not exceed K
  if(!is.null(K_0) && K_0 > K)
    stop('Expected number of factors K_0 cannot exceed the total number of factors K')
  
  # Absorb the given constraints into the basis matrix:
  if(!is.null(Con_mat)){
    if(nrow(Con_mat) != m)
      stop('The constraint matrix (Con_mat) must be m x Jc, where Jc is the number of constraints.')
    
    # Absorb the constraint:
    bup = basisConstrain(basisMat = splineInfo$Bmat,
                         penMat = splineInfo$Omega,
                         conMat = Con_mat)
    splineInfo$Bmat = bup$basisMat; splineInfo$Omega = bup$penMat
    
    # Update Psi:
    Psi = crossprod(splineInfo$Bmat, Fmat)
  }
  # Absorb the parametric constraints into the basis matrix, if known:
  if(is_unknown_gamma){
    # Unknown constraints, so add to the sampler for the FLCs
    BtCon = crossprod(splineInfo$Bmat, Gmat)
  } else {
    # Known constraints: absorb in to the basis
    bup = basisConstrain(basisMat = splineInfo$Bmat,
                         penMat = splineInfo$Omega,
                         conMat = Gmat)
    splineInfo$Bmat = bup$basisMat; splineInfo$Omega = bup$penMat
    
    # Update Psi:
    Psi = crossprod(splineInfo$Bmat, Fmat)
    
    # And no need to add to the sampler
    BtCon = NULL
  }
  # NOTE: the above code assumes that BtB is still the identity
  # We could redefine BtB (and update Psi) if necessary:
  # splineInfo$BtB = crossprod(splineInfo$Bmat)
  #----------------------------------------------------------------------------
  # Initialize the cumulative shrinkage process for nonparametric components:
  #----------------------------------------------------------------------------
  # Variance parameters:
  csp_params = sampleCSP(Lambda = Beta,
                         z = sample(1:K, K))
  theta = csp_params$theta
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Conditional mean:
  Yhat = tcrossprod(Alpha, Gmat) + tcrossprod(Beta, Fmat)
  
  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)
  
  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
  
  # SD terms for parametric components:
  sigma_alpha = apply(Alpha, 2, sd)
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, n, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('K_star', mcmc_params))) post.K_star = array(NA, c(nsave, 1))
  if(!is.na(match('K_0', mcmc_params))) post.K_0 = array(NA, c(nsave, 1))
  if(!is.na(match('Yhat', mcmc_params))) post.Yhat = array(NA, c(nsave, n, m))
  if(!is.na(match('Ypred', mcmc_params))) post.Ypred = array(NA, c(nsave, n, m))
  post_log_like_point = array(NA, c(nsave, n*m))
  #----------------------------------------------------------------------------
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    #----------------------------------------------------------------------------
    # Block 0: impute the missing data
    #----------------------------------------------------------------------------
    
    if(any.missing)
      Y[na.ind] = Yhat[na.ind] + sigma_et[na.ind[,1]]*rnorm(nrow(na.ind))
    
    #----------------------------------------------------------------------------
    # Block 1: Basis functions
    #----------------------------------------------------------------------------
    # First, sample the nonlinear parameter (if necessary)
    if(is_unknown_gamma){
      # Sample the nonlinear parameter:
      
      # Can we marginalize over alpha AND beta?
      
      # Compute the precision after marginalization:
      #prec_eps = chol2inv(chol(sigma_e^2*diag(1, m) + Fmat%*%diag(eta)%*%diag(1, K)%*%diag(eta)%*%t(Fmat)))
      #prec_eps = sigma_e^-2*(diag(1,m) - crossprod(t(Fmat)*sqrt(eta^2/(eta^2 + sigma_e^2))))
      
      # Subtract off the prior:
      #Yres = Y - tcrossprod(M*matrix(rep(eta, each = n), nrow = n), Fmat)
      Yres = Y - tcrossprod(Beta, Fmat)
      
      gamma = uni.slice(gamma, g = function(x){
        # Form the G matrix 
        G_p_x = g_p(tau, x)
        
        #  eps = tcrossprod(Alpha, G_p_x) - Yres
        # -0.5*sum(diag(tcrossprod(eps%*%prec_eps, eps)), na.rm=TRUE) + log_prior_gamma(x)
        sum(-0.5*rowSums((tcrossprod(Alpha, G_p_x) - Yres)^2, na.rm=TRUE)/sigma_et^2, na.rm=TRUE) +
          log_prior_gamma(x)
      })
      
      # Redefine the g() matrix:
      Gmat = g_p(tau, gamma)
      
      # And update the constraint matrix for the FLC sampler:
      BtCon = crossprod(splineInfo$Bmat, Gmat)
    }
    
    #----------------------------------------------------------------------------
    # Next, sample the nonparametric curves subject to orthogonality with G:
    Psi = fdlm_flc(BtY = tcrossprod(t(splineInfo$Bmat), Y - tcrossprod(Alpha, Gmat)), # BtY
                   Beta = Beta,
                   Psi = Psi,
                   BtB = diag(ncol(splineInfo$Bmat)), #splineInfo$BtB, 
                   Omega = splineInfo$Omega,
                   BtCon = BtCon,
                   lambda = tau_f_k,
                   sigmat2 = sigma_et^2)
    
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi;
    
    # Sample the smoothing parameters:
    tau_f_k = sample_lambda(tau_f_k, Psi, Omega = splineInfo$Omega, d = d, uniformPrior = TRUE, orderLambdas = FALSE)
    
    #----------------------------------------------------------------------------
    # Block 2: factors
    #----------------------------------------------------------------------------
    
    # First, sample the parametric factors Alpha:
    chQ_alpha = chol(sigma_e^-2*crossprod(Gmat) + diag(sigma_alpha^-2, L))
    lin_alpha = t(tcrossprod(Y, t(Gmat)))/sigma_e^2
    Alpha = t(backsolve(chQ_alpha, forwardsolve(t(chQ_alpha), lin_alpha) + rnorm(n*L))) 
    
    # Previous version:
    # chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    # lin_alpha = tcrossprod(Y, t(Gmat))/sigma_e^2
    # Alpha = lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L)
    
    #----------------------------------------------------------------------------
    # Next, sample the nonparametric factors Beta:
    YF = tcrossprod(Y, t(Fmat))
    
    # Previous version:
    chQ_beta = sqrt(sigma_e^-2 + rep(1/theta, each = n))
    lin_beta = YF/sigma_e^2
    Beta = lin_beta/chQ_beta^2 + 1/chQ_beta*rnorm(n*K)
    
    #----------------------------------------------------------------------------
    # Block 3: variance parameters
    #----------------------------------------------------------------------------
    
    # Update the fitted values:
    Yhat = tcrossprod(Alpha, Gmat) + tcrossprod(Beta, Fmat)
    
    # Sample the error variance:
    # Can we marginalize here?
    sigma_e = 1/sqrt(rgamma(n = 1, 
                            shape = n*m/2, 
                            rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, n)
    
    #----------------------------------------------------------------------------
    # Next, update the SD for the parametric terms:
    sigma_alpha = apply(Alpha, 2, function(x){
      1/sqrt(truncdist::rtrunc(n = 1, "gamma",
                               a = (1/100)^2, b = Inf,
                               shape = (n+1)/2,
                               rate = 1/2*sum(x^2)))
    })
    #----------------------------------------------------------------------------
    # Lastly, update the CSP parameters:
    csp_params = sampleCSP(Lambda = Beta,
                           z = csp_params$z,
                           alpha = K_0,
                           theta_inf = theta_inf,
                           a_theta = a_theta,
                           b_theta = b_theta)
    theta = csp_params$theta
    #----------------------------------------------------------------------------
    # Store the MCMC output
    #----------------------------------------------------------------------------
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples (and adjust for sdY):a
        if(!is.na(match('alpha', mcmc_params))) post.alpha[isave,,] = Alpha*sdY
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta*sdY
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma[isave,] = gamma
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e*sdY
        if(!is.na(match('K_star', mcmc_params))) post.K_star[isave,] = sum(csp_params$z > 1:K)
        if(!is.na(match('K_0', mcmc_params))) post.K_0[isave,] = K_0
        if(!is.na(match('Yhat', mcmc_params))) post.Yhat[isave,,] = Yhat*sdY
        if(!is.na(match('Ypred', mcmc_params))) post.Ypred[isave,,] = rnorm(n = n*m, mean = matrix(Yhat)*sdY, sd = rep(sigma_et,m)*sdY)
        post_log_like_point[isave,] = dnorm(matrix(Yna)*sdY, mean = matrix(Yhat)*sdY, sd = rep(sigma_et,m)*sdY, log = TRUE)
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  
  # Store the results:
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) mcmc_output$gamma = post.gamma
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('K_star', mcmc_params))) mcmc_output$K_star = post.K_star
  if(!is.na(match('K_0', mcmc_params))) mcmc_output$K_0 = post.K_0
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  
  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post_log_like_point), na.rm=TRUE)), na.rm=TRUE)
  mcmc_output$p_waic = sum(apply(post_log_like_point, 2, function(x) sd(x, na.rm=TRUE)^2), na.rm=TRUE)
  mcmc_output$WAIC = -2*(lppd - mcmc_output$p_waic)
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  
  return (mcmc_output);
}  