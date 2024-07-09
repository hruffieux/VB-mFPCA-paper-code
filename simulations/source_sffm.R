# Source Files for Semiparametric Functional Factor Models

#' MCMC algorithm for the semiparametric functional factor model
#' 
#' Run the MCMC sampling algorithm for the semiparametric functional factor model.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param K the number of (nonparametric) factors; if NULL, use SVD-based proportion of variability explained
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param K_0 hyperparameter of the cumulative shrinkage process prior: 
#' expected number of active nonparametric factors; if NULL, model as unknown 
#' with a Gamma(2,1) prior
#' @param a_1 hyperparameter of the NMIG prior: the shape parameter of the 
#' Gamma prior on the precision
#' @param a_2 hyperparameter of the NMIG prior: the rate parameter of the 
#' Gamma prior on the precision
#' @param v0 hyperparameter of the NMIG prior: the scaling for the spike component
#' of the spike-and-slab prior
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "beta" (nonparametric factors)
#' \item "fk" (nonparametric loading curves)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "K_star" (effective number of nonparametric terms)
#' \item "K_0" (prior expected number of nonparametric terms)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' }
#' @param Con_mat a \code{m x Jc} matrix of constraints for the loading curves such that
#' \code{Con_mat'fk = 0} for each loading curve \code{fk}; default is NULL for no constraints.
#' 
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{gamma}. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. For example, \code{g = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{L = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{g = function(tau, gamma) cbind(1,exp(-tau/gamma))},
#' then supply a (log) prior function via \code{log_prior_gamma} to allow for sampling of this
#' parameter. If \code{log_prior_gamma} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value, which also can be input via \code{gamma_init}. 
#'
#' @importFrom fGarch dstd
#' @export
sffm = function(Y, tau, g, K = NULL, 
                log_prior_gamma = NULL,
                K_0 = NULL, 
                a_1 = 5, a_2 = 25, v0 = 0.001, 
                gamma_init = NULL,
                nsave = 3000, nburn = 1000, nskip = 2,
                mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'K_0', 'Yhat', 'Ypred'),
                Con_mat = NULL){
  # K = NULL; log_prior_gamma = NULL;
  # K_0 = 5; a_1 = 5; a_2 = 25; v0 = 0.001; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'Yhat', 'Ypred'); Con_mat = NULL
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
    # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try)[1] == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    #g_p = function(tau, gamma) g(tau)
    g_p = function(tau, gamma) qr.Q(qr(g(tau)))
    
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    #g_p = g
    g_p = function(tau, gamma) qr.Q(qr(g(tau, gamma)))
    
  }
  
  # Check: are we sampling K_0? If so, assign hyperparameters:
  sample_K0 = is.null(K_0)
  if(sample_K0){
    # Prior: K_0 ~ Gamma(a_K, b_K)
    a_K = 2; b_K = 1;
    K_0 = ceiling(K/2) 
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
  #Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  Alpha = matrix(tcrossprod(Y, t(Gmat)), nrow = n)
  
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
  # The parameter expansion for beta:
  eta = colMeans(Beta)
  Xi = Beta/matrix(rep(eta, each = n), nrow = n) # tcrossprod(rep(1, n), eta)
  M = matrix(1, nrow = n, ncol = K); M[Xi < 0] = -1
  
  # Initialize the variance parameters:
  sigma_k = abs(eta)
  theta = rep(1, K)
  z = sample(1:K, K)
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
  sigma_alpha = apply(Alpha, 2, sd); px_sigma_alpha = rep(1, L)
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
    # chQ_alpha = chol(sigma_e^-2*crossprod(Gmat) + diag(sigma_alpha^-2, L))
    # lin_alpha = t(tcrossprod(Y, t(Gmat)))/sigma_e^2
    # Alpha = t(backsolve(chQ_alpha, forwardsolve(t(chQ_alpha), lin_alpha) + rnorm(n*L))) 
    
    chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    lin_alpha = tcrossprod(Y, t(Gmat))/sigma_e^2
    Alpha = matrix(lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L), nrow = n)
    #----------------------------------------------------------------------------
    # Next, sample the nonparametric factors (Beta) via parameter expansion:
    YF = tcrossprod(Y, t(Fmat))
    
    # (a) Sample M:
    M = matrix(1, nrow = n, ncol = K); M[runif(n = n*K) > 1/(1+ exp(-2*Xi))] = -1
    
    # (b) Sample Xi:
    chQ_xi = sqrt(rep(eta^2, each = n)/sigma_e^2 + 1)
    lin_xi = YF*rep(eta, each = n)/sigma_e^2 + M
    Xi = lin_xi/chQ_xi^2 + 1/chQ_xi*rnorm(n*K)
    
    # (c) Sample eta:
    chQ_eta = sqrt(colSums(Xi^2)/sigma_e^2 + 1/(theta*sigma_k^2))
    lin_eta = colSums(YF*Xi)/sigma_e^2
    eta = lin_eta/chQ_eta^2 + 1/chQ_eta*rnorm(K)
    
    # (d) Rescale:
    xi_scale = colMeans(abs(Xi))
    Xi = Xi/matrix(rep(xi_scale, each = n), nrow = n)
    eta = eta*xi_scale
    
    # (e) Update Beta:
    Beta = Xi*matrix(rep(eta, each = n), nrow = n) 
    
    # Previous version:
    # chQ_beta = sqrt(sigma_e^-2 + rep(sigma_beta^-2, each = n))
    # lin_beta = YF/sigma_e^2
    # Beta = lin_beta/chQ_beta^2 + 1/chQ_beta*rnorm(n*K)
    
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
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = n/2 + 1/2,
                                rate = colSums(Alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    # sigma_alpha = apply(Alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = (n+1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
    #----------------------------------------------------------------------------
    # Lastly, update the CSP parameters:
    
    # Sample the SD parameters:
    sigma_k = 1/sqrt(rgamma(n = K,
                            shape = a_1 + 0.5,
                            rate = a_2 + 0.5*eta^2/theta))
    
    # Sample the stick probabilities:
    nu = rep(0, K)
    nu[1:(K-1)] = sapply(1:(K-1), function(ell)
      rbeta(n = 1,
            shape1 = 1 + sum(z==ell),
            shape2 = K_0 + sum(z > ell))
    )
    nu[K] = 1
    
    # Sample the prior expected number of components:
    if(sample_K0){
      # Check for numerical issues: need nu[k] < 1 for k in 1:(K-1)
      if(any(nu[1:(K-1)] == 1)) nu[1:(K-1)][nu[1:(K-1)] == 1] = .999
      K_0 = rgamma(n = 1, 
                   shape = a_K + (K - 1), 
                   rate = b_K - sum(log(1 - nu[1:(K-1)])))
    }

    # Update the weights:
    omega = rep(0, K)
    omega[1] = nu[1]; 
    omega[2:K] = sapply(2:K, function(ell) {
      nu[ell]*prod(1 - nu[1:(ell-1)])
    })
    
    # Sample the categorical variables
    prob.ell = rep(0, K) # Pr(z[k] == ell) for ell= 1,...,K
    for(k in 1:(K-1)){
      # ell <= k
      prob.ell[1:k] = omega[1:k]*dstd(eta[k],
                                      mean = 0,
                                      sd = sqrt(v0*a_2/a_1),
                                      nu = 2*a_1)
      # # ell > k
      prob.ell[-(1:k)] = omega[-(1:k)]*dstd(eta[k],
                                            mean = 0,
                                            sd = sqrt(a_2/a_1),
                                            nu = 2*a_1)
      
      # Sample the latent variables:
      z[k] = sample(1:K, size = 1, prob = prob.ell)
    }
    # For h=K, the likelihood terms cancel:
    z[K] = sample(1:K, size = 1, prob = omega)
    
    # Update theta:
    theta = rep(1, K); theta[z <= 1:K] = v0
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
        if(!is.na(match('K_star', mcmc_params))) post.K_star[isave,] = sum(z > 1:K)
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
#' MCMC algorithm for the parametric functional factor model
#' 
#' Run the MCMC sampling algorithm for the parametric functional factor model.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' }
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{gamma}. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. For example, \code{g = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{L = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{g = function(tau, gamma) cbind(1,exp(-tau/gamma))},
#' then supply a (log) prior function via \code{log_prior_gamma} to allow for sampling of this
#' parameter. If \code{log_prior_gamma} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value, which also can be input via \code{gamma_init}. 
#'
#' @export
pffm = function(Y, tau, g, 
                log_prior_gamma = NULL,
                gamma_init = NULL,
                nsave = 3000, nburn = 1000, nskip = 2,
                mcmc_params = list('alpha', 'gamma', 'sigma_e', 'Yhat', 'Ypred')){
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
  # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try)[1] == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    #g_p = function(tau, gamma) g(tau)
    g_p = function(tau, gamma) qr.Q(qr(g(tau)))
    
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    #g_p = g
    g_p = function(tau, gamma) qr.Q(qr(g(tau, gamma)))
    
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
  #Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  Alpha = matrix(tcrossprod(Y, t(Gmat)), nrow = n)
  
  
  # Number of parametric terms:
  L = ncol(Gmat)
  
  # Conditional mean
  Yhat = tcrossprod(Alpha, Gmat)
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)

  # SD terms for parametric components:
  sigma_alpha = apply(Alpha, 2, sd); px_sigma_alpha = rep(1, L)
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
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

    # Impute the missing data (if necessary):
    if(any.missing)
      Y[na.ind] = Yhat[na.ind] + sigma_et[na.ind[,1]]*rnorm(nrow(na.ind))
    
    # First, sample the nonlinear parameter (if necessary)
    if(is_unknown_gamma){
      # Sample the nonlinear parameter:
      gamma = uni.slice(gamma, g = function(x){
        # Form the G matrix 
        G_p_x = g_p(tau, x)

        sum(-0.5*rowSums((tcrossprod(Alpha, G_p_x) - Y)^2, na.rm=TRUE)/sigma_et^2, na.rm=TRUE) +
          log_prior_gamma(x)  
      })
      
      # Redefine the g() matrix:
      Gmat = g_p(tau, gamma)
    }
    #----------------------------------------------------------------------------
    # Next, sample the parametric factors Alpha:
    # chQ_alpha = chol(sigma_e^-2*crossprod(Gmat) + diag(sigma_alpha^-2, L))
    # lin_alpha = t(tcrossprod(Y, t(Gmat)))/sigma_e^2
    # Alpha = t(backsolve(chQ_alpha, forwardsolve(t(chQ_alpha), lin_alpha) + rnorm(n*L))) 
    
    chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    lin_alpha = tcrossprod(Y, t(Gmat))/sigma_e^2
    Alpha = matrix(lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L), nrow = n)
    #----------------------------------------------------------------------------
    # Lastly, sample the SD parameters:
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = n/2 + 1/2,
                                rate = colSums(Alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    # sigma_alpha = apply(Alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = (n+1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
    
    
    # And update the fitted values:
    Yhat = tcrossprod(Alpha, Gmat)
    #----------------------------------------------------------------------------
    # Sample the variance:
    sigma_e = 1/sqrt(rgamma(n = 1, 
                            shape = n*m/2, 
                            rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, n)
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
        if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma[isave,] = gamma
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e*sdY
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
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) mcmc_output$gamma = post.gamma
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  
  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post_log_like_point), na.rm=TRUE)), na.rm=TRUE)
  mcmc_output$p_waic = sum(apply(post_log_like_point, 2, function(x) sd(x, na.rm=TRUE)^2), na.rm=TRUE)
  mcmc_output$WAIC = -2*(lppd - mcmc_output$p_waic)
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  
  return (mcmc_output);
}
#' MCMC algorithm for the semiparametric functional factor model w/ finite mixture priors
#' 
#' Run the MCMC sampling algorithm for the semiparametric functional factor model w/ 
#' finite mixtures instead of the stick-breaking process.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param K the number of (nonparametric) factors; if NULL, use SVD-based proportion of variability explained
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param K_0 hyperparameter of the Dirichlet prior: 
#' expected number of active nonparametric factors
#' @param a_1 hyperparameter of the NMIG prior: the shape parameter of the 
#' Gamma prior on the precision
#' @param a_2 hyperparameter of the NMIG prior: the rate parameter of the 
#' Gamma prior on the precision
#' @param v0 hyperparameter of the NMIG prior: the scaling for the spike component
#' of the spike-and-slab prior
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "beta" (nonparametric factors)
#' \item "fk" (nonparametric loading curves)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "K_star" (effective number of nonparametric terms)
#' \item "K_0" (prior expected number of nonparametric terms)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' }
#' @param Con_mat a \code{m x Jc} matrix of constraints for the loading curves such that
#' \code{Con_mat'fk = 0} for each loading curve \code{fk}; default is NULL for no constraints.
#' 
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{gamma}. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. For example, \code{g = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{L = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{g = function(tau, gamma) cbind(1,exp(-tau/gamma))},
#' then supply a (log) prior function via \code{log_prior_gamma} to allow for sampling of this
#' parameter. If \code{log_prior_gamma} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value, which also can be input via \code{gamma_init}. 
#'
#' @importFrom fGarch dstd
#' @importFrom gtools rdirichlet
#' @export
sffm_fmm = function(Y, tau, g, K = NULL, 
                log_prior_gamma = NULL,
                K_0 = 1, 
                a_1 = 5, a_2 = 25, v0 = 0.001, 
                gamma_init = NULL,
                nsave = 3000, nburn = 1000, nskip = 2,
                mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'K_0', 'Yhat', 'Ypred'),
                Con_mat = NULL){
  # K = NULL; log_prior_gamma = NULL;
  # K_0 = 5; a_1 = 5; a_2 = 25; v0 = 0.001; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'Yhat', 'Ypred'); Con_mat = NULL
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
  # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try)[1] == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    #g_p = function(tau, gamma) g(tau)
    g_p = function(tau, gamma) qr.Q(qr(g(tau)))
    
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    #g_p = g
    g_p = function(tau, gamma) qr.Q(qr(g(tau, gamma)))
    
  }
  
  # Check: are we sampling K_0? If so, assign hyperparameters:
    # NOTE: not an option here
  sample_K0 = is.null(K_0)
  if(sample_K0){
    # Prior: K_0 ~ Gamma(a_K, b_K)
    a_K = 2; b_K = 1;
    K_0 = ceiling(K/2) 
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
  #Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  Alpha = matrix(tcrossprod(Y, t(Gmat)), nrow = n)
  
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
  # The parameter expansion for beta:
  eta = colMeans(Beta)
  Xi = Beta/matrix(rep(eta, each = n), nrow = n) # tcrossprod(rep(1, n), eta)
  M = matrix(1, nrow = n, ncol = K); M[Xi < 0] = -1
  
  # Initialize the variance parameters:
  sigma_k = abs(eta)
  theta = rep(1, K)
  z = sample(1:K, K)
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
  sigma_alpha = apply(Alpha, 2, sd); px_sigma_alpha = rep(1, L)
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
    # chQ_alpha = chol(sigma_e^-2*crossprod(Gmat) + diag(sigma_alpha^-2, L))
    # lin_alpha = t(tcrossprod(Y, t(Gmat)))/sigma_e^2
    # Alpha = t(backsolve(chQ_alpha, forwardsolve(t(chQ_alpha), lin_alpha) + rnorm(n*L))) 
    
    chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    lin_alpha = tcrossprod(Y, t(Gmat))/sigma_e^2
    Alpha = matrix(lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L), nrow = n)
    #----------------------------------------------------------------------------
    # Next, sample the nonparametric factors (Beta) via parameter expansion:
    YF = tcrossprod(Y, t(Fmat))
    
    # (a) Sample M:
    M = matrix(1, nrow = n, ncol = K); M[runif(n = n*K) > 1/(1+ exp(-2*Xi))] = -1
    
    # (b) Sample Xi:
    chQ_xi = sqrt(rep(eta^2, each = n)/sigma_e^2 + 1)
    lin_xi = YF*rep(eta, each = n)/sigma_e^2 + M
    Xi = lin_xi/chQ_xi^2 + 1/chQ_xi*rnorm(n*K)
    
    # (c) Sample eta:
    chQ_eta = sqrt(colSums(Xi^2)/sigma_e^2 + 1/(theta*sigma_k^2))
    lin_eta = colSums(YF*Xi)/sigma_e^2
    eta = lin_eta/chQ_eta^2 + 1/chQ_eta*rnorm(K)
    
    # (d) Rescale:
    xi_scale = colMeans(abs(Xi))
    Xi = Xi/matrix(rep(xi_scale, each = n), nrow = n)
    eta = eta*xi_scale
    
    # (e) Update Beta:
    Beta = Xi*matrix(rep(eta, each = n), nrow = n) 
    
    # Previous version:
    # chQ_beta = sqrt(sigma_e^-2 + rep(sigma_beta^-2, each = n))
    # lin_beta = YF/sigma_e^2
    # Beta = lin_beta/chQ_beta^2 + 1/chQ_beta*rnorm(n*K)
    
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
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = n/2 + 1/2,
                                rate = colSums(Alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    # sigma_alpha = apply(Alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = (n+1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
    #----------------------------------------------------------------------------
    # Lastly, update the CSP parameters:
    
    # Sample the SD parameters:
    sigma_k = 1/sqrt(rgamma(n = K,
                            shape = a_1 + 0.5,
                            rate = a_2 + 0.5*eta^2/theta))
    
    # Sample the weights:
    n_k = sapply(1:K, function(k) sum(z==k)) # Cluster counts
    omega = as.numeric(rdirichlet(n = 1,
                       alpha = K_0/K + n_k))
    
    # Sample the prior expected number of components:
      # NOTE: not an option here
    if(sample_K0){
      stop('sampling of K0 not yet implemented for finite mixture case')
      # K_0 = rgamma(n = 1, 
      #              shape = a_K + (K - 1), 
      #              rate = b_K - sum(log(omega))/K)
    }
    
    # Sample the categorical variables
    prob.ell = rep(0, K) # Pr(z[k] == ell) for ell= 1,...,K
    for(k in 1:(K-1)){
      # ell <= k
      prob.ell[1:k] = omega[1:k]*dstd(eta[k],
                                      mean = 0,
                                      sd = sqrt(v0*a_2/a_1),
                                      nu = 2*a_1)
      # # ell > k
      prob.ell[-(1:k)] = omega[-(1:k)]*dstd(eta[k],
                                            mean = 0,
                                            sd = sqrt(a_2/a_1),
                                            nu = 2*a_1)
      
      # Sample the latent variables:
      z[k] = sample(1:K, size = 1, prob = prob.ell)
    }
    # For h=K, the likelihood terms cancel:
    z[K] = sample(1:K, size = 1, prob = omega)
    
    # Update theta:
    theta = rep(1, K); theta[z <= 1:K] = v0
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
        if(!is.na(match('K_star', mcmc_params))) post.K_star[isave,] = sum(z > 1:K)
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
#' MCMC algorithm for the spline-based semiparametric functional factor model
#' 
#' Run the MCMC sampling algorithm for the semiparametric functional factor model,
#' where the curves are modeled as splines centered at the parametric template.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param K the number of spline basis functions
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' }
#' 
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{gamma}. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. For example, \code{g = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{L = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{g = function(tau, gamma) cbind(1,exp(-tau/gamma))},
#' then supply a (log) prior function via \code{log_prior_gamma} to allow for sampling of this
#' parameter. If \code{log_prior_gamma} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value, which also can be input via \code{gamma_init}. 
#'
#' @export
sffm_spline = function(Y, tau, g, 
                   K = 15, 
                   log_prior_gamma = NULL,
                   gamma_init = NULL,
                   nsave = 3000, nburn = 1000, nskip = 2,
                   mcmc_params = list('alpha', 'gamma', 'sigma_e', 'Yhat', 'Ypred')){
  # K = 15; log_prior_gamma = NULL; gamma_init = NULL; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha', 'gamma', 'sigma_e', 'Yhat', 'Ypred')
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
  # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try)[1] == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    #g_p = function(tau, gamma) g(tau)
    g_p = function(tau, gamma) qr.Q(qr(g(tau)))
    
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    #g_p = g
    g_p = function(tau, gamma) qr.Q(qr(g(tau, gamma)))
    
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
  #Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  Alpha = matrix(tcrossprod(Y, t(Gmat)), nrow = n)
  
  # Number of parametric terms:
  L = ncol(Gmat)
  #----------------------------------------------------------------------------
  # Initialize the GP component:
  #----------------------------------------------------------------------------
  # Initialize the basis matrix, properly orthogonalized:
  bmat = getSplineInfo_d(tau01, m_eff = K-d-1, orthonormalize = FALSE)
  Fmat = qr.Q(qr(cbind(Gmat, bmat$Bmat)))[,-(1:L)] 
  K = ncol(Fmat) # update...
  
  # Coefficients:
  #Beta = matrix(tcrossprod(Y - tcrossprod(Alpha, Gmat), t(Fmat)), nrow = n)
  Beta = matrix(tcrossprod(Y, t(Fmat)), nrow = n) # faster
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Conditional mean:
  Yhat = tcrossprod(Alpha, Gmat) + tcrossprod(Beta, Fmat)
  
  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)
  
  # SD terms for parametric components:
  sigma_alpha = apply(Alpha, 2, sd); px_sigma_alpha = rep(1, L)
  
  # SD terms for spline components:
  sigma_xi = sd(Beta); px_sigma_xi = 1
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
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
      
      # And update the basis matrix:
      #BtCon = crossprod(splineInfo$Bmat, Gmat)
      Fmat = qr.Q(qr(cbind(Gmat, bmat$Bmat)))[,-(1:L)] 
    }
    #----------------------------------------------------------------------------
    # Next, sample the spline basis coefficients:
    chQ_beta = sqrt(sigma_e^-2 + sigma_xi^-2)
    lin_beta = tcrossprod(Y, t(Fmat))/sigma_e^2
    Beta = matrix(lin_beta/chQ_beta^2 + 1/chQ_beta*rnorm(n*K), nrow = n)
    #----------------------------------------------------------------------------
    # Block 2: factors
    #----------------------------------------------------------------------------
    # First, sample the parametric factors Alpha:
    # chQ_alpha = chol(sigma_e^-2*crossprod(Gmat) + diag(sigma_alpha^-2, L))
    # lin_alpha = t(tcrossprod(Y, t(Gmat)))/sigma_e^2
    # Alpha = t(backsolve(chQ_alpha, forwardsolve(t(chQ_alpha), lin_alpha) + rnorm(n*L))) 
    
    chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    lin_alpha = tcrossprod(Y, t(Gmat))/sigma_e^2
    Alpha = matrix(lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L), nrow = n)
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
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = n/2 + 1/2,
                                rate = colSums(Alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    #----------------------------------------------------------------------------
    # Update the SD for the spline basis coefficients:
    sigma_xi = 1/sqrt(rgamma(n = 1,
                             shape = n*K/2 + 1/2,
                             rate = sum(Beta^2)/2 + px_sigma_xi))
    px_sigma_xi = rgamma(n = 1, 
                         shape = 1/2 + 1/2, 
                         rate = 1/sigma_xi^2 + 1)
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
        if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma[isave,] = gamma
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e*sdY
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
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) mcmc_output$gamma = post.gamma
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  
  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post_log_like_point), na.rm=TRUE)), na.rm=TRUE)
  mcmc_output$p_waic = sum(apply(post_log_like_point, 2, function(x) sd(x, na.rm=TRUE)^2), na.rm=TRUE)
  mcmc_output$WAIC = -2*(lppd - mcmc_output$p_waic)
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  
  return (mcmc_output);
}  
#' MCMC algorithm for the semiparametric functional dynamic linear model
#' 
#' Run the MCMC sampling algorithm for the semiparametric functional dynamic linear model.
#' The parametric factors are assumed to follow independent AR(1) models, with an 
#' optional stochastic volatility model of the observation error variance.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param K the number of (nonparametric) factors; if NULL, use SVD-based proportion of variability explained
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param K_0 hyperparameter of the cumulative shrinkage process prior: 
#' expected number of active nonparametric factors; if NULL, model as unknown 
#' with a Gamma(2,1) prior
#' @param a_1 hyperparameter of the NMIG prior: the shape parameter of the 
#' Gamma prior on the precision
#' @param a_2 hyperparameter of the NMIG prior: the rate parameter of the 
#' Gamma prior on the precision
#' @param v0 hyperparameter of the NMIG prior: the scaling for the spike component
#' of the spike-and-slab prior
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "beta" (nonparametric factors)
#' \item "fk" (nonparametric loading curves)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "ar_alpha" (parametric AR coefficients)
#' \item "mu_alpha" (parametric unconditional mean)
#' \item "sigma_alpha" (parametric evolution error SD)
#' \item "K_star" (effective number of nonparametric terms)
#' \item "K_0" (prior expected number of nonparametric terms)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' \item "Yfore" (one-step forecasting distribution)
#' }
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model for the observation error variance
#' @param Con_mat a \code{m x Jc} matrix of constraints for the loading curves such that
#' \code{Con_mat'fk = 0} for each loading curve \code{fk}; default is NULL for no constraints.
#' 
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{gamma}. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. For example, \code{g = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{L = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{g = function(tau, gamma) cbind(1,exp(-tau/gamma))},
#' then supply a (log) prior function via \code{log_prior_gamma} to allow for sampling of this
#' parameter. If \code{log_prior_gamma} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value, which also can be input via \code{gamma_init}. 
#'
#' @import KFAS
#' @importFrom fGarch dstd
#' @export
sfdlm = function(Y, tau, g, K = NULL, 
                log_prior_gamma = NULL,
                K_0 = NULL, 
                a_1 = 5, a_2 = 25, v0 = 0.001,
                gamma_init = NULL,
                nsave = 3000, nburn = 1000, nskip = 2,
                mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'ar_alpha', 'mu_alpha', 'sigma_alpha','K_star', "K_0", 'Yhat', 'Ypred'),
                use_obs_SV = FALSE,
                Con_mat = NULL){
  
  # K = NULL; log_prior_gamma = NULL;
  # K_0 = 5; a_1 = 5; a_2 = 25; v0 = 0.001; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'ar_alpha', 'mu_alpha', 'sigma_alpha','K_star', 'Yhat', 'Ypred'); use_obs_SV = TRUE; Con_mat = NULL
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
  # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try)[1] == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    #g_p = function(tau, gamma) g(tau)
    g_p = function(tau, gamma) qr.Q(qr(g(tau)))
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    #g_p = g
    g_p = function(tau, gamma) qr.Q(qr(g(tau, gamma)))
  }
  
  # Check: are we sampling K_0? If so, assign hyperparameters:
  sample_K0 = is.null(K_0)
  if(sample_K0){
    # Prior: K_0 ~ Gamma(a_K, b_K)
    a_K = 2; b_K = 1;
    K_0 = ceiling(K/2) 
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
  #Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  Alpha = matrix(tcrossprod(Y, t(Gmat)), nrow = n)
  
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
  # The parameter expansion for beta:
  eta = colMeans(Beta)
  Xi = Beta/matrix(rep(eta, each = n), nrow = n) # tcrossprod(rep(1, n), eta)
  M = matrix(1, nrow = n, ncol = K); M[Xi < 0] = -1
  
  # Initialize the variance parameters:
  sigma_k = abs(eta)
  theta = rep(1, K)
  z = sample(1:K, K)
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Conditional mean:
  Yhat = tcrossprod(Alpha, Gmat) + tcrossprod(Beta, Fmat)

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Yhat)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Yhat, na.rm=TRUE)
    sigma_et = rep(sigma_e, n)
  }
  
  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
  #----------------------------------------------------------------------------
  # Now initialize the DLM parameters:
  #----------------------------------------------------------------------------
  # Unconditional means:
  mu_alpha = colMeans(Alpha); Mu_Alpha = matrix(rep(mu_alpha, each =  n), nrow = n)
  
  # AR(1) coefficients:
  ar_alpha = apply(Alpha - Mu_Alpha, 2, function(x) lm(x[-1] ~ - 1 +  x[-length(x)])$coef)
  
  # Stationarity fix:
  ar_alpha[which(abs(ar_alpha) > 0.95)] = 0.8*sign(ar_alpha[which(abs(ar_alpha) > 0.95)])
  
  # Evolution error SD::
  res_alpha = (Alpha[-1,] - Mu_Alpha[-1,]) -  
    t(ar_alpha*t(Alpha[-n,] - Mu_Alpha[-n,]))
  
  sigma_alpha = apply(res_alpha, 2, sd); px_sigma_alpha = rep(1, L)
  
  # KFAS model object:
  kfas_model = update_kfas_model(Y.dlm = Y - tcrossprod(Mu_Alpha, Gmat),
                                 Zt = array(Gmat, c(m, L, 1)),
                                 Gt = diag(ar_alpha, L),
                                 Wt = diag(sigma_alpha^2, L))
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, n, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, n))
  if(!is.na(match('ar_alpha', mcmc_params))) post.ar_alpha = array(NA, c(nsave, L))
  if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha = array(NA, c(nsave, L))
  if(!is.na(match('sigma_alpha', mcmc_params))) post.sigma_alpha = array(NA, c(nsave, L))
  if(!is.na(match('K_star', mcmc_params))) post.K_star = array(NA, c(nsave, 1))
  if(!is.na(match('K_0', mcmc_params))) post.K_0 = array(NA, c(nsave, 1))
  if(!is.na(match('Yhat', mcmc_params))) post.Yhat = array(NA, c(nsave, n, m))
  if(!is.na(match('Ypred', mcmc_params))) post.Ypred = array(NA, c(nsave, n, m))
  if(!is.na(match('Yfore', mcmc_params))) post.Yfore = array(NA, c(nsave, m))
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
      Yres = Y - tcrossprod(Beta, Fmat)
      gamma = uni.slice(gamma, g = function(x){
        # Form the G matrix 
        G_p_x = g_p(tau, x)
        
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
    # KFAS model object:
    kfas_model = update_kfas_model(Y.dlm = Y - tcrossprod(Mu_Alpha, Gmat),
                                   Zt = array(Gmat, c(m, L, 1)),
                                   sigma_et = sigma_et,
                                   Gt = diag(ar_alpha, L),
                                   Wt = diag(sigma_alpha^2, L),
                                   kfas_model = kfas_model)
    Alpha = Mu_Alpha + simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
    #----------------------------------------------------------------------------
    # Next, sample the nonparametric factors (Beta) via parameter expansion:
    YF = tcrossprod(Y, t(Fmat)); rep_sigma_et2 = rep(sigma_et^2, times = K)
    
    # (a) Sample M:
    M = matrix(1, nrow = n, ncol = K); M[runif(n = n*K) > 1/(1+ exp(-2*Xi))] = -1
    
    # (b) Sample Xi:
    chQ_xi = sqrt(rep(eta^2, each = n)/rep_sigma_et2 + 1)
    lin_xi = YF*rep(eta, each = n)/rep_sigma_et2 + M
    Xi = lin_xi/chQ_xi^2 + 1/chQ_xi*rnorm(n*K)
    
    # (c) Sample eta:
    chQ_eta = sqrt(colSums(Xi^2/rep_sigma_et2) + 1/(theta*sigma_k^2))
    lin_eta = colSums(YF*Xi/rep_sigma_et2)
    eta = lin_eta/chQ_eta^2 + 1/chQ_eta*rnorm(K)
    
    # (d) Rescale:
    xi_scale = colMeans(abs(Xi))
    Xi = Xi/matrix(rep(xi_scale, each = n), nrow = n)
    eta = eta*xi_scale
    
    # (e) Update Beta:
    Beta = Xi*matrix(rep(eta, each = n), nrow = n) 
    
    #----------------------------------------------------------------------------
    # Block 3: variance parameters
    #----------------------------------------------------------------------------
    
    # Update the fitted values:
    Yhat = tcrossprod(Alpha, Gmat) + tcrossprod(Beta, Fmat)

    # Sample the error variance:
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Yhat, svParams)
      sigma_et = svParams$sigma_et
    } else {
      sigma_e = 1/sqrt(rgamma(n = 1, 
                              shape = n*m/2, 
                              rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, n)
    }
    #----------------------------------------------------------------------------
    # Lastly, sample the AR(1) parameters:
    
    # Unconditional mean:
    mu_alpha = sampleARmu(yt = Alpha,
                          phi_j = ar_alpha,
                          sigma_tj = sigma_alpha)
    Mu_Alpha = matrix(rep(mu_alpha, each =  n), nrow = n)
    
    # AR(1) coefficients:
    ar_alpha = sampleARphi(yt = Alpha - Mu_Alpha,
                           phi_j = ar_alpha,
                           sigma_tj = sigma_alpha,
                           prior_phi = c(5,2)) #prior_phi = NULL)
    
    # Evolution error SD::
    res_alpha = (Alpha[-1,] - Mu_Alpha[-1,]) -  
      t(ar_alpha*t(Alpha[-n,] - Mu_Alpha[-n,]))
    
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = (n-1)/2 + 1/2,
                                rate = colSums(res_alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    
    # sigma_alpha = apply(res_alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = ((n-1) + 1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
    #----------------------------------------------------------------------------
    # Update the CSP parameters:
    
    # Sample the SD parameters:
    sigma_k = 1/sqrt(rgamma(n = K,
                            shape = a_1 + 0.5,
                            rate = a_2 + 0.5*eta^2/theta))
    
    # Sample the stick probabilities:
    nu = rep(0, K)
    nu[1:(K-1)] = sapply(1:(K-1), function(ell)
      rbeta(n = 1,
            shape1 = 1 + sum(z==ell),
            shape2 = K_0 + sum(z > ell))
    )
    nu[K] = 1
    
    # Sample the prior expected number of components:
    if(sample_K0){
      # Check for numerical issues: need nu[k] < 1 for k in 1:(K-1)
      if(any(nu[1:(K-1)] == 1)) nu[1:(K-1)][nu[1:(K-1)] == 1] = .999
      K_0 = rgamma(n = 1, 
                   shape = a_K + (K - 1), 
                   rate = b_K - sum(log(1 - nu[1:(K-1)])))
    }
    
    # Update the weights:
    omega = rep(0, K)
    omega[1] = nu[1]; 
    omega[2:K] = sapply(2:K, function(ell) {
      nu[ell]*prod(1 - nu[1:(ell-1)])
    })
    
    # Sample the categorical variables
    prob.ell = rep(0, K) # Pr(z[k] == ell) for ell= 1,...,K
    for(k in 1:(K-1)){
      # ell <= k
      prob.ell[1:k] = omega[1:k]*dstd(eta[k],
                                      mean = 0,
                                      sd = sqrt(v0*a_2/a_1),
                                      nu = 2*a_1)
      # # ell > k
      prob.ell[-(1:k)] = omega[-(1:k)]*dstd(eta[k],
                                            mean = 0,
                                            sd = sqrt(a_2/a_1),
                                            nu = 2*a_1)
      
      # Sample the latent variables:
      z[k] = sample(1:K, size = 1, prob = prob.ell)
    }
    # For h=K, the likelihood terms cancel:
    z[K] = sample(1:K, size = 1, prob = omega)
    
    # Update theta:
    theta = rep(1, K); theta[z <= 1:K] = v0
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
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_et*sdY
        if(!is.na(match('ar_alpha', mcmc_params))) post.ar_alpha[isave,] = ar_alpha
        if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha[isave,] = mu_alpha
        if(!is.na(match('sigma_alpha', mcmc_params))) post.sigma_alpha[isave,] = sigma_alpha
        if(!is.na(match('K_star', mcmc_params))) post.K_star[isave,] = sum(z > 1:K)
        if(!is.na(match('K_0', mcmc_params))) post.K_0[isave,] = K_0
        if(!is.na(match('Yhat', mcmc_params))) post.Yhat[isave,,] = Yhat*sdY
        if(!is.na(match('Ypred', mcmc_params))) post.Ypred[isave,,] = rnorm(n = n*m, mean = matrix(Yhat)*sdY, sd = rep(sigma_et,m)*sdY)
        post_log_like_point[isave,] = dnorm(matrix(Yna)*sdY, mean = matrix(Yhat)*sdY, sd = rep(sigma_et,m)*sdY, log = TRUE)
        
        # Forecasting distribution:
        if(!is.na(match('Yfore', mcmc_params))){
          # Parametric terms:
          alpha_fore = mu_alpha + ar_alpha*(Alpha[n,] - mu_alpha) + sigma_alpha*rnorm(n = L)
          
          # Nonparametric terms:
          m_fore = sample(c(1,-1), K, replace = TRUE, prob = c(1/2, 1/2))
          xi_fore = rnorm(n = K, mean = m_fore, sd = 1)
          eta_fore = rnorm(n = K, mean = 0, sd = sqrt(theta)*sigma_k)
          beta_fore = xi_fore*eta_fore
          
          # Observation error sd:
          sigma_fore = exp(0.5*
            (svParams$h_mu + 
               svParams$h_phi*(log(sigma_et[n]^2) - svParams$h_mu) + 
               svParams$h_sigma_eta*rnorm(n = 1))
          )
          Yfore = Gmat%*%alpha_fore + Fmat%*%beta_fore + sigma_fore*rnorm(n = m)
          post.Yfore[isave,] = Yfore*sdY 
        }
        
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
  if(!is.na(match('ar_alpha', mcmc_params))) mcmc_output$ar_alpha = post.ar_alpha
  if(!is.na(match('mu_alpha', mcmc_params))) mcmc_output$mu_alpha = post.mu_alpha
  if(!is.na(match('sigma_alpha', mcmc_params))) mcmc_output$sigma_alpha = post.sigma_alpha
  if(!is.na(match('K_star', mcmc_params))) mcmc_output$K_star = post.K_star
  if(!is.na(match('K_0', mcmc_params))) mcmc_output$K_0 = post.K_0
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  if(!is.na(match('Yfore', mcmc_params))) mcmc_output$Yfore = post.Yfore
    
  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post_log_like_point), na.rm=TRUE)), na.rm=TRUE)
  mcmc_output$p_waic = sum(apply(post_log_like_point, 2, function(x) sd(x, na.rm=TRUE)^2), na.rm=TRUE)
  mcmc_output$WAIC = -2*(lppd - mcmc_output$p_waic)
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  
  return (mcmc_output);
}  
#' MCMC algorithm for the parametric functional dynamic linear model
#' 
#' Run the MCMC sampling algorithm for the parametric functional dynamic linear model.
#' The parametric factors are assumed to follow independent AR(1) models, with an 
#' optional stochastic volatility model of the observation error variance.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "ar_alpha" (parametric AR coefficients)
#' \item "mu_alpha" (parametric unconditional mean)
#' \item "sigma_alpha" (parametric evolution error SD)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' \item "Yfore" (one-step forecasting distribution)
#' }
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model for the observation error variance
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{gamma}. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. For example, \code{g = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{L = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{g = function(tau, gamma) cbind(1,exp(-tau/gamma))},
#' then supply a (log) prior function via \code{log_prior_gamma} to allow for sampling of this
#' parameter. If \code{log_prior_gamma} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value, which also can be input via \code{gamma_init}. 
#'
#' @import KFAS
#' @export
pfdlm = function(Y, tau, g, 
                log_prior_gamma = NULL,
                gamma_init = NULL,
                nsave = 3000, nburn = 1000, nskip = 2,
                mcmc_params = list('alpha', 'gamma', 'sigma_e', 'ar_alpha', 'mu_alpha', 'sigma_alpha', 'Yhat', 'Ypred'),
                use_obs_SV = FALSE){
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
  # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try)[1] == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    #g_p = function(tau, gamma) g(tau)
    g_p = function(tau, gamma) qr.Q(qr(g(tau)))
    
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    #g_p = g
    g_p = function(tau, gamma) qr.Q(qr(g(tau, gamma)))
    
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
  #Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  Alpha = matrix(tcrossprod(Y, t(Gmat)), nrow = n)
  
  # Number of parametric terms:
  L = ncol(Gmat)
  
  # Conditional mean
  Yhat = tcrossprod(Alpha, Gmat)
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Yhat)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Yhat, na.rm=TRUE)
    sigma_et = rep(sigma_e, n)
  }
  #----------------------------------------------------------------------------
  # Now initialize the DLM parameters:
  #----------------------------------------------------------------------------
  # Unconditional means:
  mu_alpha = colMeans(Alpha); Mu_Alpha = matrix(rep(mu_alpha, each =  n), nrow = n)
  
  # AR(1) coefficients:
  ar_alpha = apply(Alpha - Mu_Alpha, 2, function(x) lm(x[-1] ~ - 1 +  x[-length(x)])$coef)
  
  # Stationarity fix:
  ar_alpha[which(abs(ar_alpha) > 0.95)] = 0.8*sign(ar_alpha[which(abs(ar_alpha) > 0.95)])

  # Evolution error SD::
  res_alpha = (Alpha[-1,] - Mu_Alpha[-1,]) -  
    t(ar_alpha*t(Alpha[-n,] - Mu_Alpha[-n,]))

  sigma_alpha = apply(res_alpha, 2, sd); px_sigma_alpha = rep(1, L)
  
  # KFAS model object:
  kfas_model = update_kfas_model(Y.dlm = Y - tcrossprod(Mu_Alpha, Gmat),
                                 Zt = array(Gmat, c(m, L, 1)),
                                 Gt = diag(ar_alpha, L),
                                 Wt = diag(sigma_alpha^2, L))
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, n))
  if(!is.na(match('ar_alpha', mcmc_params))) post.ar_alpha = array(NA, c(nsave, L))
  if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha = array(NA, c(nsave, L))
  if(!is.na(match('sigma_alpha', mcmc_params))) post.sigma_alpha = array(NA, c(nsave, L))
  if(!is.na(match('Yhat', mcmc_params))) post.Yhat = array(NA, c(nsave, n, m))
  if(!is.na(match('Ypred', mcmc_params))) post.Ypred = array(NA, c(nsave, n, m))
  if(!is.na(match('Yfore', mcmc_params))) post.Yfore = array(NA, c(nsave, m))
  post_log_like_point = array(NA, c(nsave, n*m))
  #----------------------------------------------------------------------------
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){
    
    # Impute the missing data (if necessary):
    if(any.missing)
      Y[na.ind] = Yhat[na.ind] + sigma_et[na.ind[,1]]*rnorm(nrow(na.ind))
    
    # First, sample the nonlinear parameter (if necessary)
    if(is_unknown_gamma){
      # Sample the nonlinear parameter:
      gamma = uni.slice(gamma, g = function(x){
        # Form the G matrix 
        G_p_x = g_p(tau, x)
        
        sum(-0.5*rowSums((tcrossprod(Alpha, G_p_x) - Y)^2, na.rm=TRUE)/sigma_et^2, na.rm=TRUE) +
          log_prior_gamma(x)  
      })
      
      # Redefine the g() matrix:
      Gmat = g_p(tau, gamma)
    }
    #----------------------------------------------------------------------------
    # Next, sample the parametric factors Alpha:
    # KFAS model object:
    kfas_model = update_kfas_model(Y.dlm = Y - tcrossprod(Mu_Alpha, Gmat),
                                   Zt = array(Gmat, c(m, L, 1)),
                                   sigma_et = sigma_et,
                                   Gt = diag(ar_alpha, L),
                                   Wt = diag(sigma_alpha^2, L),
                                   kfas_model = kfas_model)
    Alpha = Mu_Alpha + simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
    #----------------------------------------------------------------------------
    # Lastly, sample the AR(1) parameters:
    
    # Unconditional mean:
    mu_alpha = sampleARmu(yt = Alpha,
                      phi_j = ar_alpha,
                      sigma_tj = sigma_alpha)
    Mu_Alpha = matrix(rep(mu_alpha, each =  n), nrow = n)
    
    # AR(1) coefficients:
    ar_alpha = sampleARphi(yt = Alpha - Mu_Alpha,
                           phi_j = ar_alpha,
                           sigma_tj = sigma_alpha,
                           prior_phi = c(5,2)) #prior_phi = NULL)

    # Evolution error SD::
    res_alpha = (Alpha[-1,] - Mu_Alpha[-1,]) -  
      t(ar_alpha*t(Alpha[-n,] - Mu_Alpha[-n,]))
    
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = (n-1)/2 + 1/2,
                                rate = colSums(res_alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    # sigma_alpha = apply(res_alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = ((n-1) + 1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
    
    # And update the fitted values:
    Yhat = tcrossprod(Alpha, Gmat)
    #----------------------------------------------------------------------------
    # Sample the variance:
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Yhat, svParams)
      sigma_et = svParams$sigma_et
    } else {
      sigma_e = 1/sqrt(rgamma(n = 1, 
                              shape = n*m/2, 
                              rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, n)
    }
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
        if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma[isave,] = gamma
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_et*sdY
        if(!is.na(match('ar_alpha', mcmc_params))) post.ar_alpha[isave,] = ar_alpha
        if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha[isave,] = mu_alpha
        if(!is.na(match('sigma_alpha', mcmc_params))) post.sigma_alpha[isave,] = sigma_alpha
        if(!is.na(match('Yhat', mcmc_params))) post.Yhat[isave,,] = Yhat*sdY
        if(!is.na(match('Ypred', mcmc_params))) post.Ypred[isave,,] = rnorm(n = n*m, mean = matrix(Yhat)*sdY, sd = rep(sigma_et,m)*sdY)
        post_log_like_point[isave,] = dnorm(matrix(Yna)*sdY, mean = matrix(Yhat)*sdY, sd = rep(sigma_et,m)*sdY, log = TRUE)
        
        # Forecasting distribution:
        if(!is.na(match('Yfore', mcmc_params))){
          # Parametric terms:
          alpha_fore = mu_alpha + ar_alpha*(Alpha[n,] - mu_alpha) + sigma_alpha*rnorm(n = L)
    
          # Observation error sd:
          sigma_fore = exp(0.5*
                             (svParams$h_mu + 
                                svParams$h_phi*(log(sigma_et[n]^2) - svParams$h_mu) + 
                                svParams$h_sigma_eta*rnorm(n = 1))
          )
          Yfore = Gmat%*%alpha_fore + sigma_fore*rnorm(n = m)
          post.Yfore[isave,] = Yfore*sdY 
        }
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  
  # Store the results:
  if(!is.na(match('alpha', mcmc_params))) mcmc_output$alpha = post.alpha
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) mcmc_output$gamma = post.gamma
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('ar_alpha', mcmc_params))) mcmc_output$ar_alpha = post.ar_alpha
  if(!is.na(match('mu_alpha', mcmc_params))) mcmc_output$mu_alpha = post.mu_alpha
  if(!is.na(match('sigma_alpha', mcmc_params))) mcmc_output$sigma_alpha = post.sigma_alpha
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  if(!is.na(match('Yfore', mcmc_params))) mcmc_output$Yfore = post.Yfore
  
  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post_log_like_point), na.rm=TRUE)), na.rm=TRUE)
  mcmc_output$p_waic = sum(apply(post_log_like_point, 2, function(x) sd(x, na.rm=TRUE)^2), na.rm=TRUE)
  mcmc_output$WAIC = -2*(lppd - mcmc_output$p_waic)
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  
  return (mcmc_output);
}
#' Constrain the basis function
#' 
#' To constrain C'f = 0 for f = B*theta, 
#' we absorb the constraint into the basis matrix B. 
#' 
#' @param basisMat \code{m x J} basis matrix, 
#' where \code{m} is the number of observation points and 
#' \code{J} is the number of basis functions
#' @param penMat \code{J x J} penalty matrix
#' @param conMat \code{m x Jc} matrix of constraints
#' @return a list with the following elements:
#' \itemize{
#' \item \code{basisMat} the constrained basis matrix of dimension \code{m x (J - Jc)}
#' \item \code{penMat} the constrained penalty matrix of dimension \code{(J - Jc) x (J - Jc)}
#' \item \code{cQR} the QR decomposition of the constraint
#' }
#' @export
basisConstrain = function(basisMat, penMat, conMat){
  
  # Number of constraints:
  Jc = ncol(conMat)
  
  # Check:
  if(nrow(basisMat) != nrow(conMat))
    stop('conMat and basisMat should have the same number of rows')
  if(ncol(basisMat) != ncol(penMat) || ncol(basisMat) != nrow(penMat))
    stop('basisMat and penMat should have the same basis dimension')
  
  # Constraint:
  Con = crossprod(conMat, basisMat)
  
  # QR Decomposition:
  cQR = qr(t(Con))
  
  # Update the basis:
  basisMat = t(qr.qty(cQR, t(basisMat))[-(1:Jc),])
  
  # Update the penalty:
  penMat = qr.qty(cQR,t(qr.qty(cQR,t(penMat))))[-(1:Jc),-(1:Jc)]
  
  # And return:
  list(
    basisMat = basisMat,
    penMat = penMat,
    cQR = cQR
  )
}
#' Sample from the cumulative shrinkage prior distribution
#' 
#' Sample the (nonparametric) factors \code{Beta} (n x K) and
#' the magnitude parameters \code{eta} (K x 1) from the cumulative
#' shrinkage prior with the peNMIG parameter expansion. 
#' 
#' @param n the sample size
#' @param K the number of (nonparametric) factors
#' @param K_0 hyperparameter of the cumulative shrinkage process prior: 
#' expected number of active nonparametric factors; if NULL, model as unknown 
#' with a Gamma(2,1) prior
#' @param a_1 hyperparameter of the NMIG prior: the shape parameter of the 
#' Gamma prior on the precision
#' @param a_2 hyperparameter of the NMIG prior: the rate parameter of the 
#' Gamma prior on the precision
#' @param v0 hyperparameter of the NMIG prior: the scaling for the spike component
#' of the spike-and-slab prior
#' @return a list with the following elements:
#' \itemize{
#' \item \code{Beta} the sampled factors of dimension \code{n x K}
#' \item \code{eta} the sampled magnitudes of dimnsion \code{K x 1}
#' }
#' @export
prior_predict = function(n = 100, K = 15,  K_0 = NULL, 
                         a_1 = 5, a_2 = 25, v0 = 0.001){
  # Sample the M variables and the parameter-expanded part:
  M = sample(c(1,-1), K, replace = TRUE)
  Xi = rnorm(n = K, mean = M, sd = 1)
  
  # Sample the SD parameters:
  sigma_k = 1/sqrt(rgamma(n = K,
                          shape = a_1,
                          rate = a_2))
  
  # Prior: K_0 ~ Gamma(a_K, b_K)
  if(is.null(K_0)) 
    K_0 = rgamma(n = 1, shape = 2, rate = 1)
  
  # Sample the stick probabilities:
  nu = rep(0, K)
  nu[1:(K-1)] = sapply(1:(K-1), function(ell)
    rbeta(n = 1,
          shape1 = 1,
          shape2 = K_0)
  )
  nu[K] = 1
  
  # Update the weights:
  omega = rep(0, K)
  omega[1] = nu[1]; 
  omega[2:K] = sapply(2:K, function(ell) {
    nu[ell]*prod(1 - nu[1:(ell-1)])
  })
  
  # Sample the categorical variables
  z = sample(1:K, size = K, prob = omega)

  # Update theta:
  theta = rep(1, K); theta[z <= 1:K] = v0
  
  # And update eta:
  eta = rnorm(n = K, mean = 0, sd = sqrt(theta)*sigma_k)
  
  # Return:
  Beta = Xi*matrix(rep(eta, each = n), nrow = n) 
  
  list(
    Beta = Beta,
    eta = eta
  )
}
#' MCMC algorithm for the semiparametric functional dynamic linear model
#' with exogenous predictors.
#' 
#' Run the MCMC sampling algorithm for the semiparametric functional dynamic linear model.
#' The parametric factors are regressed on predictors with AR(1) errors. 
#' An optional stochastic volatility model is permitted for the observation error variance.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{n x p} matrix of predictors; if NULL, only include an intercept
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param K the number of (nonparametric) factors; if NULL, use SVD-based proportion of variability explained
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param K_0 hyperparameter of the cumulative shrinkage process prior: 
#' expected number of active nonparametric factors; if NULL, model as unknown 
#' with a Gamma(2,1) prior
#' @param a_1 hyperparameter of the NMIG prior: the shape parameter of the 
#' Gamma prior on the precision
#' @param a_2 hyperparameter of the NMIG prior: the rate parameter of the 
#' Gamma prior on the precision
#' @param v0 hyperparameter of the NMIG prior: the scaling for the spike component
#' of the spike-and-slab prior
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "zeta" (regression coefficients for parametric factors)
#' \item "beta" (nonparametric factors)
#' \item "fk" (nonparametric loading curves)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "ar_alpha" (parametric AR coefficients)
#' \item "mu_alpha" (parametric unconditional mean)
#' \item "sigma_alpha" (parametric evolution error SD)
#' \item "K_star" (effective number of nonparametric terms)
#' \item "K_0" (prior expected number of nonparametric terms)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' }
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model for the observation error variance
#' @param Con_mat a \code{m x Jc} matrix of constraints for the loading curves such that
#' \code{Con_mat'fk = 0} for each loading curve \code{fk}; default is NULL for no constraints.
#' 
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{gamma}. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. For example, \code{g = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{L = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{g = function(tau, gamma) cbind(1,exp(-tau/gamma))},
#' then supply a (log) prior function via \code{log_prior_gamma} to allow for sampling of this
#' parameter. If \code{log_prior_gamma} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value, which also can be input via \code{gamma_init}. 
#'
#' @import KFAS
#' @importFrom fGarch dstd
#' @export
sfosr_ar = function(Y, tau, X = NULL, g, K = NULL, 
                 log_prior_gamma = NULL,
                 K_0 = NULL, 
                 a_1 = 5, a_2 = 25, v0 = 0.001,
                 gamma_init = NULL,
                 nsave = 3000, nburn = 1000, nskip = 2,
                 mcmc_params = list('alpha', 'zeta', 'beta', 'fk', 'gamma', 'sigma_e', 'ar_alpha', 'mu_alpha', 'sigma_alpha','K_star', "K_0", 'Yhat', 'Ypred'),
                 use_obs_SV = FALSE,
                 Con_mat = NULL){
  
  # K = NULL; log_prior_gamma = NULL;
  # K_0 = NULL; a_1 = 5; a_2 = 25; v0 = 0.001; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha','zeta', 'beta', 'fk', 'gamma', 'sigma_e', 'ar_alpha', 'mu_alpha', 'sigma_alpha','K_star', 'Yhat', 'Ypred'); use_obs_SV = TRUE; Con_mat = NULL
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }

  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
  # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try)[1] == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    #g_p = function(tau, gamma) g(tau)
    g_p = function(tau, gamma) qr.Q(qr(g(tau)))
    
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    #g_p = g
    g_p = function(tau, gamma) qr.Q(qr(g(tau, gamma)))
    
  }
  
  # Check: are we sampling K_0? If so, assign hyperparameters:
  sample_K0 = is.null(K_0)
  if(sample_K0){
    # Prior: K_0 ~ Gamma(a_K, b_K)
    a_K = 2; b_K = 1;
    K_0 = ceiling(K/2) 
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
  #Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  Alpha = matrix(tcrossprod(Y, t(Gmat)), nrow = n)
  
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
  # The parameter expansion for beta:
  eta = colMeans(Beta)
  Xi = Beta/matrix(rep(eta, each = n), nrow = n) # tcrossprod(rep(1, n), eta)
  M = matrix(1, nrow = n, ncol = K); M[Xi < 0] = -1
  
  # Initialize the variance parameters:
  sigma_k = abs(eta)
  theta = rep(1, K)
  z = sample(1:K, K)
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Conditional mean:
  Yhat = tcrossprod(Alpha, Gmat) + tcrossprod(Beta, Fmat)
  
  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Yhat)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Yhat, na.rm=TRUE)
    sigma_et = rep(sigma_e, n)
  }
  
  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
  #----------------------------------------------------------------------------
  # Predictors:
  #----------------------------------------------------------------------------
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, n), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors:
  p = ncol(X)
  
  # Initialize the SSModel:
  X.arr = array(t(X), c(1, p, n))
  kfas_model = update_kfas_model(Y.dlm = as.matrix(Alpha[,1]), Zt = X.arr)
  
  # Identify all components are non-dynamic
  if(p > 1) diag(kfas_model$R[,,1])[-1] = 0
  #----------------------------------------------------------------------------
  # Now initialize the DLM parameters:
  #----------------------------------------------------------------------------
  # Unconditional means:
  mu_alpha = colMeans(Alpha); Mu_Alpha = matrix(rep(mu_alpha, each =  n), nrow = n)
  
  # AR(1) Evolution Matrix
  G_alpha = diag(p) # Replace the intercept terms as needed
  
  # AR(1) coefficients:
  ar_alpha = apply(Alpha - Mu_Alpha, 2, function(x) lm(x[-1] ~ - 1 +  x[-length(x)])$coef)
  
  # Stationarity fix:
  ar_alpha[which(abs(ar_alpha) > 0.95)] = 0.8*sign(ar_alpha[which(abs(ar_alpha) > 0.95)])
  #----------------------------------------------------------------------------
  # Initialize the regression terms:
  zeta_arr = array(0, c(n, p, L))
  for(ell in 1:L){
    # Update the evolution matrix
    G_alpha[1,1] = ar_alpha[ell]
    
    # Update the SSModel object given the new parameters
    kfas_model = update_kfas_model(Y.dlm = as.matrix(Alpha[,ell] - mu_alpha[ell]),
                                   Zt = X.arr,
                                   Gt = G_alpha,
                                   kfas_model = kfas_model)
    # Run the sampler
    zeta_arr[,,ell] = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
    
    # Conditional mean from regression equation:
    Alpha[,ell] = mu_alpha[ell] + rowSums(X*zeta_arr[,,ell])
  }

  # Evolution error SD::
  res_alpha = zeta_arr[-1,1,] - t(ar_alpha*t(zeta_arr[-n,1,]))
  
  sigma_alpha = apply(res_alpha, 2, sd); px_sigma_alpha = rep(1, L)
  
  # Storage for evolution error variance matrix:
  Wt = diag(p); W0 = diag(10^-4, p);
  #----------------------------------------------------------------------------
  # Non-intercept term:
  if(p > 1){
    # Non-dynamic setting: grab the first one (all the same) and store as (L x p-1) matrix
    zeta_reg = matrix(t(zeta_arr[1, -1, ]), nrow = L)
    
    # factor ell, predictor p:
    sigma_omega_ellp = abs(zeta_reg)
    xi_omega_ellp = matrix(1, nrow = L, ncol = p-1) # PX term
    
    # predictor p:
    lambda_omega_p = colMeans(sigma_omega_ellp)
    xi_omega_p = rep(1, (p-1)) # PX term
    
    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('zeta', mcmc_params))) post.zeta = array(NA, c(nsave, n, p, L))
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, n, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, n))
  if(!is.na(match('ar_alpha', mcmc_params))) post.ar_alpha = array(NA, c(nsave, L))
  if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha = array(NA, c(nsave, L))
  if(!is.na(match('sigma_alpha', mcmc_params))) post.sigma_alpha = array(NA, c(nsave, L))
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
      Yres = Y - tcrossprod(Beta, Fmat)
      gamma = uni.slice(gamma, g = function(x){
        # Form the G matrix 
        G_p_x = g_p(tau, x)
        
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
    YG = tcrossprod(Y, t(Gmat))
    
    # Loop over each factor ell = 1,...,L:
    for(ell in 1:L){
      
      # Update the evolution matrix
      G_alpha[1,1] = ar_alpha[ell]
      
      # Update the evolution variance:
      Wt[1,1] = W0[1,1] = sigma_alpha[ell]^2
      if(p > 1) W0[-1, -1] = diag(as.numeric(sigma_omega_ellp[ell,]^2), p - 1)
      
      # Sanity check for Wt: if variances too large, KFAS will stop running
      Wt[which(Wt > 10^6, arr.ind = TRUE)] = 10^6; W0[which(W0 > 10^6, arr.ind = TRUE)] = 10^6
      
      # Update the SSModel object given the new parameters
      kfas_model = update_kfas_model(Y.dlm = matrix(YG[,ell] - mu_alpha[ell]),
                                     Zt = X.arr,
                                     sigma_et = sigma_et,
                                     Gt = G_alpha,
                                     Wt = Wt, 
                                     W0 = W0,
                                     kfas_model = kfas_model)
      
      # Run the sampler
      zeta_arr[,,ell] = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
      
      # Conditional mean from regression equation:
      Alpha[,ell] = mu_alpha[ell] + rowSums(X*zeta_arr[,,ell])
      
    }
    #----------------------------------------------------------------------------
    # Next, sample the nonparametric factors (Beta) via parameter expansion:
    YF = tcrossprod(Y, t(Fmat)); rep_sigma_et2 = rep(sigma_et^2, times = K)
    
    # (a) Sample M:
    M = matrix(1, nrow = n, ncol = K); M[runif(n = n*K) > 1/(1+ exp(-2*Xi))] = -1
    
    # (b) Sample Xi:
    chQ_xi = sqrt(rep(eta^2, each = n)/rep_sigma_et2 + 1)
    lin_xi = YF*rep(eta, each = n)/rep_sigma_et2 + M
    Xi = lin_xi/chQ_xi^2 + 1/chQ_xi*rnorm(n*K)
    
    # (c) Sample eta:
    chQ_eta = sqrt(colSums(Xi^2/rep_sigma_et2) + 1/(theta*sigma_k^2))
    lin_eta = colSums(YF*Xi/rep_sigma_et2)
    eta = lin_eta/chQ_eta^2 + 1/chQ_eta*rnorm(K)
    
    # (d) Rescale:
    xi_scale = colMeans(abs(Xi))
    Xi = Xi/matrix(rep(xi_scale, each = n), nrow = n)
    eta = eta*xi_scale
    
    # (e) Update Beta:
    Beta = Xi*matrix(rep(eta, each = n), nrow = n) 
    
    #----------------------------------------------------------------------------
    # Block 3: variance parameters
    #----------------------------------------------------------------------------
    
    # Update the fitted values:
    Yhat = tcrossprod(Alpha, Gmat) + tcrossprod(Beta, Fmat)
    
    # Sample the error variance:
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Yhat, svParams)
      sigma_et = svParams$sigma_et
    } else {
      sigma_e = 1/sqrt(rgamma(n = 1, 
                              shape = n*m/2, 
                              rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, n)
    }
    #----------------------------------------------------------------------------
    # Sample the intercept parameters (Note: could use ASIS)
    
    # Centerend and non-centered:
    Alpha_int =  matrix(zeta_arr[,1,], nrow = n)
    
    Alpha_int_c = Alpha_int + Mu_Alpha
    
    # Sample the unconditional mean term:
    mu_alpha = sampleARmu(yt = Alpha_int_c,
                          phi_j = ar_alpha,
                          sigma_tj = sigma_alpha)
    Mu_Alpha = matrix(rep(mu_alpha, each =  n), nrow = n)

    # And update the non-centered parameter:
    Alpha_int = Alpha_int_c - Mu_Alpha
    
    # AR(1) coefficients:
    ar_alpha = sampleARphi(yt = Alpha_int,
                           phi_j = ar_alpha,
                           sigma_tj = sigma_alpha,
                           prior_phi = c(5,2)) #prior_phi = NULL)
    
    # Evolution error SD::
    res_alpha = Alpha_int[-1,] - t(ar_alpha*t(Alpha_int[-n,]))
    
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = (n-1)/2 + 1/2,
                                rate = colSums(res_alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    
    # sigma_alpha = apply(res_alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = ((n-1) + 1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
    #----------------------------------------------------------------------------
    # Sample the regression variance parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      # Non-dynamic setting: grab the first one (all the same) and store as (L x p-1) matrix
      zeta_reg = matrix(t(zeta_arr[1, -1, ]), nrow = L)
      #----------------------------------------------------------------------------
      # factor ell, predictor p:
      zeta_reg2 = zeta_reg^2; zeta_reg2 = zeta_reg2 + (zeta_reg2 < 10^-16)*10^-8
      
      sigma_omega_ellp = matrix(1/sqrt(rgamma(n = L*(p-1),
                                            shape = 1/2 + 1/2,
                                            rate = xi_omega_ellp + zeta_reg2/2)), nrow = L)
      xi_omega_ellp = matrix(rgamma(n = L*(p-1),
                                  shape = 1/2 + 1/2,
                                  rate = rep(1/lambda_omega_p^2, each = L) + 1/sigma_omega_ellp^2), nrow = L)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + L/2,
                                     rate = xi_omega_p + colSums(xi_omega_ellp)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
    #----------------------------------------------------------------------------
    # Update the CSP parameters:
    
    # Sample the SD parameters:
    sigma_k = 1/sqrt(rgamma(n = K,
                            shape = a_1 + 0.5,
                            rate = a_2 + 0.5*eta^2/theta))
    
    # Sample the stick probabilities:
    nu = rep(0, K)
    nu[1:(K-1)] = sapply(1:(K-1), function(ell)
      rbeta(n = 1,
            shape1 = 1 + sum(z==ell),
            shape2 = K_0 + sum(z > ell))
    )
    nu[K] = 1
    
    # Sample the prior expected number of components:
    if(sample_K0){
      # Check for numerical issues: need nu[k] < 1 for k in 1:(K-1)
      if(any(nu[1:(K-1)] == 1)) nu[1:(K-1)][nu[1:(K-1)] == 1] = .999
      K_0 = rgamma(n = 1, 
                   shape = a_K + (K - 1), 
                   rate = b_K - sum(log(1 - nu[1:(K-1)])))
    }
    
    # Update the weights:
    omega = rep(0, K)
    omega[1] = nu[1]; 
    omega[2:K] = sapply(2:K, function(ell) {
      nu[ell]*prod(1 - nu[1:(ell-1)])
    })
    
    # Sample the categorical variables
    prob.ell = rep(0, K) # Pr(z[k] == ell) for ell= 1,...,K
    for(k in 1:(K-1)){
      # ell <= k
      prob.ell[1:k] = omega[1:k]*dstd(eta[k],
                                      mean = 0,
                                      sd = sqrt(v0*a_2/a_1),
                                      nu = 2*a_1)
      # # ell > k
      prob.ell[-(1:k)] = omega[-(1:k)]*dstd(eta[k],
                                            mean = 0,
                                            sd = sqrt(a_2/a_1),
                                            nu = 2*a_1)
      
      # Sample the latent variables:
      z[k] = sample(1:K, size = 1, prob = prob.ell)
    }
    # For h=K, the likelihood terms cancel:
    z[K] = sample(1:K, size = 1, prob = omega)
    
    # Update theta:
    theta = rep(1, K); theta[z <= 1:K] = v0
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
        if(!is.na(match('zeta', mcmc_params))) post.zeta[isave,,,] = zeta_arr*sdY
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta*sdY
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma[isave,] = gamma
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_et*sdY
        if(!is.na(match('ar_alpha', mcmc_params))) post.ar_alpha[isave,] = ar_alpha
        if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha[isave,] = mu_alpha
        if(!is.na(match('sigma_alpha', mcmc_params))) post.sigma_alpha[isave,] = sigma_alpha
        if(!is.na(match('K_star', mcmc_params))) post.K_star[isave,] = sum(z > 1:K)
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
  if(!is.na(match('zeta', mcmc_params))) mcmc_output$zeta = post.zeta
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) mcmc_output$gamma = post.gamma
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('ar_alpha', mcmc_params))) mcmc_output$ar_alpha = post.ar_alpha
  if(!is.na(match('mu_alpha', mcmc_params))) mcmc_output$mu_alpha = post.mu_alpha
  if(!is.na(match('sigma_alpha', mcmc_params))) mcmc_output$sigma_alpha = post.sigma_alpha
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
#' MCMC algorithm for the parametric functional dynamic linear model
#' with exogenous predictors.
#' 
#' Run the MCMC sampling algorithm for the parametric functional dynamic linear model.
#' The parametric factors are regressed on predictors with AR(1) errors. 
#' An optional stochastic volatility model is permitted for the observation error variance.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param X the \code{n x p} matrix of predictors; if NULL, only include an intercept
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "zeta" (regression coefficients for parametric factors)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "ar_alpha" (parametric AR coefficients)
#' \item "mu_alpha" (parametric unconditional mean)
#' \item "sigma_alpha" (parametric evolution error SD)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' }
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model for the observation error variance
#' 
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{gamma}. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. For example, \code{g = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{L = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{g = function(tau, gamma) cbind(1,exp(-tau/gamma))},
#' then supply a (log) prior function via \code{log_prior_gamma} to allow for sampling of this
#' parameter. If \code{log_prior_gamma} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value, which also can be input via \code{gamma_init}. 
#'
#' @import KFAS
#' @export
pfosr_ar = function(Y, tau, X = NULL, g, 
                    log_prior_gamma = NULL,
                    gamma_init = NULL,
                    nsave = 3000, nburn = 1000, nskip = 2,
                    mcmc_params = list('alpha', 'zeta', 'gamma', 'sigma_e', 'ar_alpha', 'mu_alpha', 'sigma_alpha','Yhat', 'Ypred'),
                    use_obs_SV = FALSE){
  
  # log_prior_gamma = NULL; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha', 'zeta', 'gamma', 'sigma_e', 'ar_alpha', 'mu_alpha', 'sigma_alpha','Yhat', 'Ypred')
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Is the nonlinear parameter unknown?
  is_unknown_gamma = !is.null(log_prior_gamma) 
  
  # Redefine the input function to have a silent nonlinear input, if necessary:
  # Internally, we use g_p() as the function
  g_try = try(g(tau, 1), silent = TRUE)
  if(class(g_try)[1] == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    is_unknown_gamma = FALSE
    #g_p = function(tau, gamma) g(tau)
    g_p = function(tau, gamma) qr.Q(qr(g(tau)))
    
  } else {
    # g() has an unknown lambda: make sure we have a prior!
    if(is.null(log_prior_gamma))
      stop('If g() depends on unknown gamma, log_prior_gamma must be specified')
    
    #g_p = g
    g_p = function(tau, gamma) qr.Q(qr(g(tau, gamma)))
    
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
  #Alpha = tcrossprod(Y, t(Gmat))%*%chol2inv(chol(crossprod(Gmat)))
  Alpha = matrix(tcrossprod(Y, t(Gmat)), nrow = n)
  
  # Number of parametric terms:
  L = ncol(Gmat)
  
  # Conditional mean
  Yhat = tcrossprod(Alpha, Gmat)
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Yhat)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Yhat, na.rm=TRUE)
    sigma_et = rep(sigma_e, n)
  }
  #----------------------------------------------------------------------------
  # Predictors:
  #----------------------------------------------------------------------------
  if(!is.null(X)){
    # Assuming we have some predictors:
    X = as.matrix(X)
    
    # Remove any predictors which are constants/intercepts:
    const.pred = apply(X, 2, function(x) all(diff(x) == 0))
    if(any(const.pred)) X = as.matrix(X[,!const.pred])
    
    # Center and scale the (non-constant) predictors:
    # Note: may not be appropriate for intervention effects!
    #X = scale(X)
  }
  # Include an intercept:
  X = cbind(rep(1, n), X); #colnames(X)[1] = paste(intercept_model, "-Intercept", sep='')
  
  # Number of predictors:
  p = ncol(X)
  
  # Initialize the SSModel:
  X.arr = array(t(X), c(1, p, n))
  kfas_model = update_kfas_model(Y.dlm = as.matrix(Alpha[,1]), Zt = X.arr)
  
  # Identify all components are non-dynamic
  if(p > 1) diag(kfas_model$R[,,1])[-1] = 0
  #----------------------------------------------------------------------------
  # Now initialize the DLM parameters:
  #----------------------------------------------------------------------------
  # Unconditional means:
  mu_alpha = colMeans(Alpha); Mu_Alpha = matrix(rep(mu_alpha, each =  n), nrow = n)
  
  # AR(1) Evolution Matrix
  G_alpha = diag(p) # Replace the intercept terms as needed
  
  # AR(1) coefficients:
  ar_alpha = apply(Alpha - Mu_Alpha, 2, function(x) lm(x[-1] ~ - 1 +  x[-length(x)])$coef)
  
  # Stationarity fix:
  ar_alpha[which(abs(ar_alpha) > 0.95)] = 0.8*sign(ar_alpha[which(abs(ar_alpha) > 0.95)])
  #----------------------------------------------------------------------------
  # Initialize the regression terms:
  zeta_arr = array(0, c(n, p, L))
  for(ell in 1:L){
    # Update the evolution matrix
    G_alpha[1,1] = ar_alpha[ell]
    
    # Update the SSModel object given the new parameters
    kfas_model = update_kfas_model(Y.dlm = as.matrix(Alpha[,ell] - mu_alpha[ell]),
                                   Zt = X.arr,
                                   Gt = G_alpha,
                                   kfas_model = kfas_model)
    # Run the sampler
    zeta_arr[,,ell] = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
    
    # Conditional mean from regression equation:
    Alpha[,ell] = mu_alpha[ell] + rowSums(X*zeta_arr[,,ell])
  }
  
  # Evolution error SD::
  res_alpha = zeta_arr[-1,1,] - t(ar_alpha*t(zeta_arr[-n,1,]))
  
  sigma_alpha = apply(res_alpha, 2, sd); px_sigma_alpha = rep(1, L)
  
  # Storage for evolution error variance matrix:
  Wt = diag(p); W0 = diag(10^-4, p);
  #----------------------------------------------------------------------------
  # Non-intercept term:
  if(p > 1){
    # Non-dynamic setting: grab the first one (all the same) and store as (L x p-1) matrix
    zeta_reg = matrix(t(zeta_arr[1, -1, ]), nrow = L)
    
    # factor ell, predictor p:
    sigma_omega_ellp = abs(zeta_reg)
    xi_omega_ellp = matrix(1, nrow = L, ncol = p-1) # PX term
    
    # predictor p:
    lambda_omega_p = colMeans(sigma_omega_ellp)
    xi_omega_p = rep(1, (p-1)) # PX term
    
    # global:
    lambda_omega_0 = mean(lambda_omega_p)
    xi_omega_0 = 1 # PX term
  }
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('zeta', mcmc_params))) post.zeta = array(NA, c(nsave, n, p, L))
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, n))
  if(!is.na(match('ar_alpha', mcmc_params))) post.ar_alpha = array(NA, c(nsave, L))
  if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha = array(NA, c(nsave, L))
  if(!is.na(match('sigma_alpha', mcmc_params))) post.sigma_alpha = array(NA, c(nsave, L))
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
      gamma = uni.slice(gamma, g = function(x){
        # Form the G matrix 
        G_p_x = g_p(tau, x)
        
        sum(-0.5*rowSums((tcrossprod(Alpha, G_p_x) - Y)^2, na.rm=TRUE)/sigma_et^2, na.rm=TRUE) +
          log_prior_gamma(x)  
      })
      
      # Redefine the g() matrix:
      Gmat = g_p(tau, gamma)
    }
    #----------------------------------------------------------------------------
    # Next, sample the factors
    #----------------------------------------------------------------------------
    YG = tcrossprod(Y, t(Gmat))
    
    # Loop over each factor ell = 1,...,L:
    for(ell in 1:L){
      
      # Update the evolution matrix
      G_alpha[1,1] = ar_alpha[ell]
      
      # Update the evolution variance:
      Wt[1,1] = W0[1,1] = sigma_alpha[ell]^2
      if(p > 1) W0[-1, -1] = diag(as.numeric(sigma_omega_ellp[ell,]^2), p - 1)
      
      # Sanity check for Wt: if variances too large, KFAS will stop running
      Wt[which(Wt > 10^6, arr.ind = TRUE)] = 10^6; W0[which(W0 > 10^6, arr.ind = TRUE)] = 10^6
      
      # Update the SSModel object given the new parameters
      kfas_model = update_kfas_model(Y.dlm = matrix(YG[,ell] - mu_alpha[ell]),
                                     Zt = X.arr,
                                     sigma_et = sigma_et,
                                     Gt = G_alpha,
                                     Wt = Wt, W0 = W0,
                                     kfas_model = kfas_model)
      
      # Run the sampler
      zeta_arr[,,ell] = simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
      
      # Conditional mean from regression equation:
      Alpha[,ell] = mu_alpha[ell] + rowSums(X*zeta_arr[,,ell])
      
    }
    #----------------------------------------------------------------------------
    # Block 3: variance parameters
    #----------------------------------------------------------------------------
    # Update the fitted values:
    Yhat = tcrossprod(Alpha, Gmat) 
    
    # Sample the error variance:
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Yhat, svParams)
      sigma_et = svParams$sigma_et
    } else {
      sigma_e = 1/sqrt(rgamma(n = 1, 
                              shape = n*m/2, 
                              rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, n)
    }
    #----------------------------------------------------------------------------
    # Sample the intercept parameters (Note: could use ASIS)
    
    # Centerend and non-centered:
    Alpha_int =  matrix(zeta_arr[,1,], nrow = n)
    
    Alpha_int_c = Alpha_int + Mu_Alpha
    
    # Sample the unconditional mean term:
    mu_alpha = sampleARmu(yt = Alpha_int_c,
                          phi_j = ar_alpha,
                          sigma_tj = sigma_alpha)
    Mu_Alpha = matrix(rep(mu_alpha, each =  n), nrow = n)
    
    # And update the non-centered parameter:
    Alpha_int = Alpha_int_c - Mu_Alpha
    
    # AR(1) coefficients:
    ar_alpha = sampleARphi(yt = Alpha_int,
                           phi_j = ar_alpha,
                           sigma_tj = sigma_alpha,
                           prior_phi = c(5,2)) #prior_phi = NULL)
    
    # Evolution error SD::
    res_alpha = Alpha_int[-1,] - t(ar_alpha*t(Alpha_int[-n,]))
    
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = (n-1)/2 + 1/2,
                                rate = colSums(res_alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    
    # sigma_alpha = apply(res_alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = ((n-1) + 1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
    #----------------------------------------------------------------------------
    # Sample the regression variance parameters:
    #----------------------------------------------------------------------------
    # Non-intercept term:
    if(p > 1){
      # Non-dynamic setting: grab the first one (all the same) and store as (L x p-1) matrix
      zeta_reg = matrix(t(zeta_arr[1, -1, ]), nrow = L)
      #----------------------------------------------------------------------------
      # factor ell, predictor p:
      zeta_reg2 = zeta_reg^2; zeta_reg2 = zeta_reg2 + (zeta_reg2 < 10^-16)*10^-8
      sigma_omega_ellp = matrix(1/sqrt(rgamma(n = L*(p-1),
                                              shape = 1/2 + 1/2,
                                              rate = xi_omega_ellp + zeta_reg2/2)), nrow = L)
      xi_omega_ellp = matrix(rgamma(n = L*(p-1),
                                    shape = 1/2 + 1/2,
                                    rate = rep(1/lambda_omega_p^2, each = L) + 1/sigma_omega_ellp^2), nrow = L)
      #----------------------------------------------------------------------------
      # predictor p:
      lambda_omega_p = 1/sqrt(rgamma(n = p-1,
                                     shape = 1/2 + L/2,
                                     rate = xi_omega_p + colSums(xi_omega_ellp)))
      xi_omega_p = rgamma(n = p-1,
                          shape = 1/2 + 1/2,
                          rate = rep(1/lambda_omega_0^2, p-1) + 1/lambda_omega_p^2)
      #----------------------------------------------------------------------------
      # global:
      lambda_omega_0 = 1/sqrt(rgamma(n = 1,
                                     shape = 1/2 + (p-1)/2,
                                     rate = xi_omega_0 + sum(xi_omega_p)))
      xi_omega_0 = rgamma(n = 1,
                          shape = 1/2 + 1/2,
                          rate = 1 + 1/lambda_omega_0^2)
    }
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
        if(!is.na(match('zeta', mcmc_params))) post.zeta[isave,,,] = zeta_arr*sdY
        if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) post.gamma[isave,] = gamma
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_et*sdY
        if(!is.na(match('ar_alpha', mcmc_params))) post.ar_alpha[isave,] = ar_alpha
        if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha[isave,] = mu_alpha
        if(!is.na(match('sigma_alpha', mcmc_params))) post.sigma_alpha[isave,] = sigma_alpha
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
  if(!is.na(match('zeta', mcmc_params))) mcmc_output$zeta = post.zeta
  if(!is.na(match('gamma', mcmc_params)) && is_unknown_gamma) mcmc_output$gamma = post.gamma
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('ar_alpha', mcmc_params))) mcmc_output$ar_alpha = post.ar_alpha
  if(!is.na(match('mu_alpha', mcmc_params))) mcmc_output$mu_alpha = post.mu_alpha
  if(!is.na(match('sigma_alpha', mcmc_params))) mcmc_output$sigma_alpha = post.sigma_alpha
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  
  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post_log_like_point), na.rm=TRUE)), na.rm=TRUE)
  mcmc_output$p_waic = sum(apply(post_log_like_point, 2, function(x) sd(x, na.rm=TRUE)^2), na.rm=TRUE)
  mcmc_output$WAIC = -2*(lppd - mcmc_output$p_waic)
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  
  return (mcmc_output);
}  
#' MCMC algorithm for the semiparametric functional factor model with 
#' subject-specific nonlinear parameters
#' 
#' Run the MCMC sampling algorithm for the semiparametric functional factor model. 
#' Here, we allow for each curve to have its own known and unknown nonlinear parameters.
#' The unknown nonlinear parameters follow a hierarchical Gaussian model. 
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include nonlinear parameter arguments
#' @param known_params the \code{n x 1} vector of known nonlinear parameters, which appear as
#' arguments in \code{g}; default is NULL
#' @param K the number of (nonparametric) factors; if NULL, use SVD-based proportion of variability explained
#' @param K_0 hyperparameter of the cumulative shrinkage process prior: 
#' expected number of active nonparametric factors; if NULL, model as unknown 
#' with a Gamma(2,1) prior
#' @param a_1 hyperparameter of the NMIG prior: the shape parameter of the 
#' Gamma prior on the precision
#' @param a_2 hyperparameter of the NMIG prior: the rate parameter of the 
#' Gamma prior on the precision
#' @param v0 hyperparameter of the NMIG prior: the scaling for the spike component
#' of the spike-and-slab prior
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "beta" (nonparametric factors)
#' \item "fk" (nonparametric loading curves)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "K_star" (effective number of nonparametric terms)
#' \item "K_0" (prior expected number of nonparametric terms)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' }
#' @param Con_mat a \code{m x Jc} matrix of constraints for the loading curves such that
#' \code{Con_mat'fk = 0} for each loading curve \code{fk}; default is NULL for no constraints.
#' 
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include an nonlinear
#' parameter, \code{gamma}, and a known nonlinear parameter, \code{known_param}, which are
#' curve-specific. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. 
#'
#' @importFrom fGarch dstd
#' @export
sffm_sub = function(Y, tau, g, 
                    known_params = NULL,
                    K = NULL, 
                    K_0 = NULL, 
                    a_1 = 5, a_2 = 25, v0 = 0.001, 
                    gamma_init = NULL,
                    nsave = 3000, nburn = 1000, nskip = 2,
                    mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'K_0', 'Yhat', 'Ypred'),
                    Con_mat = NULL){
  
  # K = 10; K_0 = NULL; a_1 = 5; a_2 = 25; v0 = 0.001; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'Yhat', 'Ypred'); Con_mat = NULL
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Initialize:
  gamma_i = gamma_init
  mu_gamma = mean(gamma_i);
  sigma_gamma = sd(gamma_i); px_sigma_gamma = 1
    
  # Define the ORTHOGONALIZED functions:
  g_p = function(tau, gamma, known_params){
    qr.Q(qr(
      g(tau, gamma, known_params)
    ))
  } 
  
  # Check: are we sampling K_0? If so, assign hyperparameters:
  sample_K0 = is.null(K_0)
  if(sample_K0){
    # Prior: K_0 ~ Gamma(a_K, b_K)
    a_K = 2; b_K = 1;
    K_0 = ceiling(K/2) 
  }
  #----------------------------------------------------------------------------
  # Initialize the parametric component:
  #----------------------------------------------------------------------------
  # List of basis matrices:
  Glist = lapply(1:n, function(i) g_p(tau, gamma_i[i], known_params[i]))
  
  # Number of parametric terms:
  L = ncol(Glist[[1]])
  
  # Sample mean of basis matrices:
  Gavg = matrix(rowMeans(matrix(unlist(Glist), m*L, n)), m, L)
  #Gavg = 0; for(i in 1:n) Gavg = Gavg + 1/n*Glist[[i]]

  # And initialize the coefficients:
  Alpha = matrix(t(sapply(1:n, function(i) 
    crossprod(Glist[[i]], Y[i,]))), nrow = n)
  
  # Fitted parametric part:
  GAlpha = t(sapply(1:n, function(i)
    tcrossprod(Alpha[i,], Glist[[i]])))
  #----------------------------------------------------------------------------
  # Initialize the nonparametric component:
  #----------------------------------------------------------------------------
  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Y - GAlpha, tau, K); 
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
  BtCon = crossprod(splineInfo$Bmat, Gavg)
  # NOTE: the above code assumes that BtB is still the identity
  # We could redefine BtB (and update Psi) if necessary:
  # splineInfo$BtB = crossprod(splineInfo$Bmat)
  #----------------------------------------------------------------------------
  # Initialize the cumulative shrinkage process for nonparametric components:
  #----------------------------------------------------------------------------
  # The parameter expansion for beta:
  eta = colMeans(Beta)
  Xi = Beta/matrix(rep(eta, each = n), nrow = n) # tcrossprod(rep(1, n), eta)
  M = matrix(1, nrow = n, ncol = K); M[Xi < 0] = -1
  
  # Initialize the variance parameters:
  sigma_k = abs(eta)
  theta = rep(1, K)
  z = sample(1:K, K)
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Conditional mean:
  Yhat = GAlpha + tcrossprod(Beta, Fmat)
  
  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)
  
  # Initialize the FLC smoothing parameters (conditional MLE):
  tau_f_k = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - (d+1))/crossprod(x, splineInfo$Omega)%*%x)
  
  # SD terms for parametric components:
  sigma_alpha = apply(Alpha, 2, sd); px_sigma_alpha = rep(1, L)
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, n, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('gamma', mcmc_params))) post.gamma = array(NA, c(nsave, n))
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
    # First, sample the nonlinear parameters
    
    # Subtract off the NP part:
    Yres = Y - tcrossprod(Beta, Fmat)
    
    # Sample the nonlinear parameters for each subject:
    gamma_i = sapply(1:n, function(i){
      uni.slice(gamma_i[i], g = function(x){
        # Form the G matrix 
        G_p_x = g_p(tau, x, known_params[i])
        
        # log-likelihood + log-prior
        -0.5*sum(((tcrossprod(Alpha[i,], G_p_x) - Yres[i,])/sigma_et[i])^2, na.rm = TRUE) +
          dnorm(x, mean = mu_gamma, sd = sigma_gamma, log = TRUE)
      })
    })
    
    # Sample the mean parameter:
    chQ_mu = sqrt(n/sigma_gamma^2 + 1/10)
    lin_mu = sum(gamma_i)/sigma_gamma^2
    mu_gamma = lin_mu/chQ_mu^2 + 1/chQ_mu*rnorm(n=1)
    
    # Sample the sd parameter:
    sigma_gamma = 1/sqrt(rgamma(n = 1,
                                shape = n/2 + 1/2,
                                rate = sum((gamma_i - mu_gamma)^2)/2 + px_sigma_gamma))
    px_sigma_gamma = rgamma(n = 1, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_gamma^2 + 1)
    
    
    # Redefine the g() list:
    #Glist = lapply(gamma_i, function(gamma) g_p(tau, gamma))
    Glist = lapply(1:n, function(i) g_p(tau, gamma_i[i], known_params[i]))
    
    # Sample mean of basis matrices:
    Gavg = matrix(rowMeans(matrix(unlist(Glist), m*L, n)), m, L)
    
    # And update the constraint matrix for the FLC sampler:
    BtCon = crossprod(splineInfo$Bmat, Gavg)
    
    # And update the additional terms:
    GAlpha = t(sapply(1:n, function(i)
      tcrossprod(Alpha[i,], Glist[[i]])))
    #----------------------------------------------------------------------------
    # Next, sample the nonparametric curves subject to orthogonality with G:
    Psi = fdlm_flc(BtY = tcrossprod(t(splineInfo$Bmat), Y - GAlpha), # BtY
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
    chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    lin_alpha = t(sapply(1:n, function(i) 
      crossprod(Glist[[i]], Y[i,] - Fmat%*%Beta[i,])/sigma_e^2))
    Alpha = matrix(lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L), nrow = n)
    
    # And update:
    GAlpha = t(sapply(1:n, function(i)
      tcrossprod(Alpha[i,], Glist[[i]])))
    #----------------------------------------------------------------------------
    # Next, sample the nonparametric factors (Beta) via parameter expansion:
    #YF = tcrossprod(Y, t(Fmat))
    YF = tcrossprod(Y - GAlpha, t(Fmat))
    
    # (a) Sample M:
    M = matrix(1, nrow = n, ncol = K); M[runif(n = n*K) > 1/(1+ exp(-2*Xi))] = -1
    
    # (b) Sample Xi:
    chQ_xi = sqrt(rep(eta^2, each = n)/sigma_e^2 + 1)
    lin_xi = YF*rep(eta, each = n)/sigma_e^2 + M
    Xi = lin_xi/chQ_xi^2 + 1/chQ_xi*rnorm(n*K)
    
    # (c) Sample eta:
    chQ_eta = sqrt(colSums(Xi^2)/sigma_e^2 + 1/(theta*sigma_k^2))
    lin_eta = colSums(YF*Xi)/sigma_e^2
    eta = lin_eta/chQ_eta^2 + 1/chQ_eta*rnorm(K)
    
    # (d) Rescale:
    xi_scale = colMeans(abs(Xi))
    Xi = Xi/matrix(rep(xi_scale, each = n), nrow = n)
    eta = eta*xi_scale
    
    # (e) Update Beta:
    Beta = Xi*matrix(rep(eta, each = n), nrow = n) 
    
    #----------------------------------------------------------------------------
    # Block 3: variance parameters
    #----------------------------------------------------------------------------
    
    # Update the fitted values:
    Yhat = GAlpha + tcrossprod(Beta, Fmat)
    
    # Sample the error variance:
    # Can we marginalize here?
    sigma_e = 1/sqrt(rgamma(n = 1, 
                            shape = n*m/2, 
                            rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, n)
    
    #----------------------------------------------------------------------------
    # Next, update the SD for the parametric terms:
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = n/2 + 1/2,
                                rate = colSums(Alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    # sigma_alpha = apply(Alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = (n+1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
    #----------------------------------------------------------------------------
    # Lastly, update the CSP parameters:
    
    # Sample the SD parameters:
    sigma_k = 1/sqrt(rgamma(n = K,
                            shape = a_1 + 0.5,
                            rate = a_2 + 0.5*eta^2/theta))
    
    # Sample the stick probabilities:
    nu = rep(0, K)
    nu[1:(K-1)] = sapply(1:(K-1), function(ell)
      rbeta(n = 1,
            shape1 = 1 + sum(z==ell),
            shape2 = K_0 + sum(z > ell))
    )
    nu[K] = 1
    
    # Sample the prior expected number of components:
    if(sample_K0){
      # Check for numerical issues: need nu[k] < 1 for k in 1:(K-1)
      if(any(nu[1:(K-1)] == 1)) nu[1:(K-1)][nu[1:(K-1)] == 1] = .999
      K_0 = rgamma(n = 1, 
                   shape = a_K + (K - 1), 
                   rate = b_K - sum(log(1 - nu[1:(K-1)])))
    }
    
    # Update the weights:
    omega = rep(0, K)
    omega[1] = nu[1]; 
    omega[2:K] = sapply(2:K, function(ell) {
      nu[ell]*prod(1 - nu[1:(ell-1)])
    })
    
    # Sample the categorical variables
    prob.ell = rep(0, K) # Pr(z[k] == ell) for ell= 1,...,K
    for(k in 1:(K-1)){
      # ell <= k
      prob.ell[1:k] = omega[1:k]*dstd(eta[k],
                                      mean = 0,
                                      sd = sqrt(v0*a_2/a_1),
                                      nu = 2*a_1)
      # # ell > k
      prob.ell[-(1:k)] = omega[-(1:k)]*dstd(eta[k],
                                            mean = 0,
                                            sd = sqrt(a_2/a_1),
                                            nu = 2*a_1)
      
      # Sample the latent variables:
      z[k] = sample(1:K, size = 1, prob = prob.ell)
    }
    # For h=K, the likelihood terms cancel:
    z[K] = sample(1:K, size = 1, prob = omega)
    
    # Update theta:
    theta = rep(1, K); theta[z <= 1:K] = v0
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
        if(!is.na(match('gamma', mcmc_params))) post.gamma[isave,] = gamma_i
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e*sdY
        if(!is.na(match('K_star', mcmc_params))) post.K_star[isave,] = sum(z > 1:K)
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
  if(!is.na(match('gamma', mcmc_params))) mcmc_output$gamma = post.gamma
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
#' MCMC algorithm for the parametric functional factor model with 
#' subject-specific nonlinear parameters
#' 
#' Run the MCMC sampling algorithm for the parametric functional factor model.
#' Here, we allow for each curve to have its own known and unknown nonlinear parameters.
#' The unknown nonlinear parameters follow a hierarchical Gaussian model. 
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param known_params the \code{n x 1} vector of known nonlinear parameters, which appear as
#' arguments in \code{g}; default is NULL
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "alpha" (parametric factors)
#' \item "gamma" (parametric function parameter)
#' \item "sigma_e" (observation error SD)
#' \item "Yhat" (fitted values)
#' \item "Ypred" (posterior predictive values)
#' }
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' 
#' @details  The parametric function \code{g} should input an \code{m x d} matrix
#' of observation points, \code{tau}, and may include an nonlinear
#' parameter, \code{gamma}, and a known nonlinear parameter, \code{known_param}, which are
#' curve-specific. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. 
#'
#' @export
pffm_sub = function(Y, tau, g, 
                    known_params = NULL,
                    gamma_init = NULL,
                    nsave = 3000, nburn = 1000, nskip = 2,
                    mcmc_params = list('alpha', 'gamma', 'sigma_e', 'Yhat', 'Ypred')){
  
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
    #Y = fdlm_init(Y, tau)$Y0
    
    # Impute using column, row, and overall mean:
    Y = apply(Y, 2, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y})
    if(any(is.na(Y))) Y = t(apply(Y, 1, function(y){y[is.na(y)] = mean(y, na.rm=TRUE); y}))
    if(any(is.na(Y))) Y[is.na(Y)] = mean(Y, na.rm=TRUE)
  }
  #----------------------------------------------------------------------------
  # Define the nonlinear components:
  #----------------------------------------------------------------------------
  # Initialize:
  gamma_i = gamma_init
  mu_gamma = mean(gamma_i);
  sigma_gamma = sd(gamma_i); px_sigma_gamma = 1
  
  # Define the ORTHOGONALIZED functions:
  g_p = function(tau, gamma, known_params){
    qr.Q(qr(
      g(tau, gamma, known_params)
    ))
  } 
  #----------------------------------------------------------------------------
  # Initialize the parametric component:
  #----------------------------------------------------------------------------
  # List of basis matrices:
  #Glist = lapply(gamma_i, function(gamma) g_p(tau, gamma))
  Glist = lapply(1:n, function(i) g_p(tau, gamma_i[i], known_params[i]))
  
  # Number of parametric terms:
  L = ncol(Glist[[1]])

  # And initialize the coefficients:
  Alpha = matrix(t(sapply(1:n, function(i) 
    crossprod(Glist[[i]], Y[i,]))), nrow = n)
  
  # Fitted parametric part:
  Yhat = t(sapply(1:n, function(i)
    tcrossprod(Alpha[i,], Glist[[i]])))
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)
  
  # SD terms for parametric components:
  sigma_alpha = apply(Alpha, 2, sd); px_sigma_alpha = rep(1, L)
  #----------------------------------------------------------------------------
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('alpha', mcmc_params))) post.alpha = array(NA, c(nsave, n, L))
  if(!is.na(match('gamma', mcmc_params))) post.gamma = array(NA, c(nsave, n))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsave, 1))
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
    # First, sample the nonlinear parameters
    
    # Sample the nonlinear parameters for each subject:
    gamma_i = sapply(1:n, function(i){
      uni.slice(gamma_i[i], g = function(x){
        # Form the G matrix 
        G_p_x = g_p(tau, x, known_params[i])
        
        # log-likelihood + log-prior
        -0.5*sum(((tcrossprod(Alpha[i,], G_p_x) - Y[i,])/sigma_et[i])^2, na.rm = TRUE) +
          dnorm(x, mean = mu_gamma, sd = sigma_gamma, log = TRUE)
      })
    })
    
    # Sample the mean parameter:
    chQ_mu = sqrt(n/sigma_gamma^2 + 1/10)
    lin_mu = sum(gamma_i)/sigma_gamma^2
    mu_gamma = lin_mu/chQ_mu^2 + 1/chQ_mu*rnorm(n=1)
    
    # Sample the sd parameter:
    sigma_gamma = 1/sqrt(rgamma(n = 1,
                                shape = n/2 + 1/2,
                                rate = sum((gamma_i - mu_gamma)^2)/2 + px_sigma_gamma))
    px_sigma_gamma = rgamma(n = 1, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_gamma^2 + 1)
    
    # Redefine the g() list:
    #Glist = lapply(gamma_i, function(gamma) g_p(tau, gamma))
    Glist = lapply(1:n, function(i) g_p(tau, gamma_i[i], known_params[i]))
    #----------------------------------------------------------------------------
    # Block 2: factors
    #----------------------------------------------------------------------------
    # First, sample the parametric factors Alpha:
    chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    lin_alpha = t(sapply(1:n, function(i) 
      crossprod(Glist[[i]], Y[i,])/sigma_e^2))
    Alpha = matrix(lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L), nrow = n)
    #----------------------------------------------------------------------------
    # Block 3: variance parameters
    #----------------------------------------------------------------------------
    
    # Update the fitted values:
    Yhat = t(sapply(1:n, function(i)
      tcrossprod(Alpha[i,], Glist[[i]])))
    
    # Sample the error variance:
    # Can we marginalize here?
    sigma_e = 1/sqrt(rgamma(n = 1, 
                            shape = n*m/2, 
                            rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, n)
    
    #----------------------------------------------------------------------------
    # Next, update the SD for the parametric terms:
    sigma_alpha = 1/sqrt(rgamma(n = L,
                                shape = n/2 + 1/2,
                                rate = colSums(Alpha^2)/2 + px_sigma_alpha))
    px_sigma_alpha = rgamma(n = L, 
                            shape = 1/2 + 1/2, 
                            rate = 1/sigma_alpha^2 + 1)
    # sigma_alpha = apply(Alpha, 2, function(x){
    #   1/sqrt(truncdist::rtrunc(n = 1, "gamma",
    #                            a = (1/100)^2, b = Inf,
    #                            shape = (n+1)/2,
    #                            rate = 1/2*sum(x^2)))
    # })
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
        if(!is.na(match('gamma', mcmc_params))) post.gamma[isave,] = gamma_i
        if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[isave,] = sigma_e*sdY
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
  if(!is.na(match('gamma', mcmc_params))) mcmc_output$gamma = post.gamma
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.na(match('Ypred', mcmc_params))) mcmc_output$Ypred = post.Ypred
  
  # Compute WAIC:
  lppd = sum(log(colMeans(exp(post_log_like_point), na.rm=TRUE)), na.rm=TRUE)
  mcmc_output$p_waic = sum(apply(post_log_like_point, 2, function(x) sd(x, na.rm=TRUE)^2), na.rm=TRUE)
  mcmc_output$WAIC = -2*(lppd - mcmc_output$p_waic)
  
  print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
  
  return (mcmc_output);
}  
