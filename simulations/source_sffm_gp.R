#' MCMC algorithm for the Gaussian process semiparametric functional factor model
#' 
#' Run the MCMC sampling algorithm for the semiparametric functional factor model,
#' where the curves are modeled as Gaussian processes centered at the parametric template.
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param log_prior_gamma a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param rho the correlation parameter for the Gaussian process
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
#' @import FastGP
sffm_gp = function(Y, tau, g, 
                   log_prior_gamma = NULL,
                   gamma_init = NULL,
                   rho = 0.5,
                   nsave = 3000, nburn = 1000, nskip = 2,
                   mcmc_params = list('alpha', 'gamma', 'sigma_e', 'Yhat', 'Ypred')){
  # log_prior_gamma = NULL; gamma_init = NULL; rho = 0.5; nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha', 'gamma', 'sigma_e', 'Yhat', 'Ypred')
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
  
  # matrix of absolute distances between tau's (for GP)
  tauDiffs = getTDiffs(tau01)						        
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
  Rinv = chol2inv(chol(corrFun(tauDiffs, rho)))
  sigma_h = sigma_e = 1
  chQ_h = chol(1/sigma_h^2*Rinv + diag(1/sigma_e^2, m))
  lin_h = 1/sigma_e^2*(Y - tcrossprod(Alpha, Gmat))
  h_all = t(backsolve(chQ_h, forwardsolve(t(chQ_h), t(lin_h)) + rnorm(n*m))) 
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Conditional mean:
  Yhat = tcrossprod(Alpha, Gmat) + h_all
  
  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)
  
  # SD terms for parametric components:
  sigma_alpha = apply(Alpha, 2, sd); px_sigma_alpha = rep(1, L)
  
  # SD terms for GP components:
  sigma_h = sqrt(sum(diag(crossprod(tcrossprod(h_all, Rinv), h_all)))/(n*m))
  px_sigma_h = 1
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
      Yres = Y - h_all
      
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
    }
    
    # And the GPs:
    #prec_lik = chol2inv(chol(Gmat%*%diag(sigma_alpha^2)%*%t(Gmat) + diag(sigma_e^2, m)))
    prec_lik = 1/sigma_e^2*(
      diag(1,m) - crossprod(t(Gmat)*sqrt(1/(sigma_e^2/sigma_alpha^2 + 1)))
      #diag(1,m) - Gmat%*%diag(1/(sigma_e^2/sigma_alpha^2 + 1))%*%t(Gmat)
    )
    chQ_h = chol(1/sigma_h^2*Rinv + prec_lik)
    lin_h = tcrossprod(Y, prec_lik)
    h_all = t(backsolve(chQ_h, forwardsolve(t(chQ_h), t(lin_h)) + rnorm(n*m)))   
    #----------------------------------------------------------------------------
    # Block 2: factors
    #----------------------------------------------------------------------------
    # First, sample the parametric factors Alpha:
    # chQ_alpha = chol(sigma_e^-2*crossprod(Gmat) + diag(sigma_alpha^-2, L))
    # lin_alpha = t(tcrossprod(Y, t(Gmat)))/sigma_e^2
    # Alpha = t(backsolve(chQ_alpha, forwardsolve(t(chQ_alpha), lin_alpha) + rnorm(n*L))) 
    
    chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    lin_alpha = tcrossprod(Y - h_all, t(Gmat))/sigma_e^2
    Alpha = matrix(lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L), nrow = n)
    #----------------------------------------------------------------------------
    # Block 3: variance parameters
    #----------------------------------------------------------------------------
    
    # Update the fitted values:
    Yhat = tcrossprod(Alpha, Gmat) + h_all
    
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
    # Update the SD for the GPs:
    sigma_h = 1/sqrt(rgamma(n = 1,
                             shape = n*m/2 + 1/2,
                             rate = sum(diag(crossprod(tcrossprod(h_all, Rinv), h_all)))/2 + px_sigma_h))
    px_sigma_h = rgamma(n = 1, 
                         shape = 1/2 + 1/2, 
                         rate = 1/sigma_h^2 + 1)
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
#' MCMC algorithm for the Gaussian process semiparametric functional factor model with 
#' subject-specific nonlinear parameters
#' 
#' Run the MCMC sampling algorithm for the Gaussian process semiparametric functional factor model. 
#' Here, we allow for each curve to have its own known and unknown nonlinear parameters.
#' The unknown nonlinear parameters follow a hierarchical Gaussian model. 
#' 
#' @param Y the \code{n x m} data observation matrix, where \code{n} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau the \code{m x d} matrix of coordinates of observation points
#' @param g a function to compute the parametric component, which must return a \code{m x L} matrix
#' for \code{L} the number of parametric curves; may include nonlinear parameter arguments
#' @param known_params the \code{n x 1} vector of known nonlinear parameters, which appear as
#' arguments in \code{g}; default is NULL
#' @param gamma_init initial value for gamma; if NULL, initialize randomly
#' @param rho the correlation parameter for the Gaussian process
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
#' of observation points, \code{tau}, and may include an nonlinear
#' parameter, \code{gamma}, and a known nonlinear parameter, \code{known_param}, which are
#' curve-specific. The function should return a \code{m x L} matrix, where \code{L} is the
#' number of parametric functions. 
#'
#' @importFrom fGarch dstd
#' @export
sffm_gp_sub = function(Y, tau, g, 
                    known_params = NULL,
                    gamma_init = NULL,
                    rho = 0.5,
                    nsave = 3000, nburn = 1000, nskip = 2,
                    mcmc_params = list('alpha', 'beta', 'fk', 'gamma', 'sigma_e', 'K_star', 'K_0', 'Yhat', 'Ypred'),
                    Con_mat = NULL){
  
  # nsave = 1000; nburn = 1000; nskip = 0; mcmc_params = list('alpha',  'gamma', 'sigma_e', 'K_star', 'Yhat', 'Ypred'); Con_mat = NULL
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
  
  # matrix of absolute distances between tau's (for GP)
  tauDiffs = getTDiffs(tau01)		
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
  # Initialize the GP component:
  #----------------------------------------------------------------------------
  Rinv = chol2inv(chol(corrFun(tauDiffs, rho)))
  sigma_h = sigma_e = 1
  chQ_h = chol(1/sigma_h^2*Rinv + diag(1/sigma_e^2, m))
  lin_h = 1/sigma_e^2*(Y - GAlpha)
  h_all = t(backsolve(chQ_h, forwardsolve(t(chQ_h), t(lin_h)) + rnorm(n*m))) 
  #----------------------------------------------------------------------------
  # Initialize the remaining terms:
  #----------------------------------------------------------------------------
  # Conditional mean:
  Yhat = GAlpha + h_all
  
  # Initialize the observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE); sigma_et = rep(sigma_e, n)
  
  # SD terms for parametric components:
  sigma_alpha = apply(Alpha, 2, sd); px_sigma_alpha = rep(1, L)
  
  # SD terms for GP components:
  sigma_h = sqrt(sum(diag(crossprod(tcrossprod(h_all, Rinv), h_all)))/(n*m))
  px_sigma_h = 1
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
    
    # Subtract off the NP part:
    Yres = Y - h_all
    
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

    
    # And update the additional terms:
    GAlpha = t(sapply(1:n, function(i)
      tcrossprod(Alpha[i,], Glist[[i]])))
    #----------------------------------------------------------------------------
    # And the GPs:
    for(i in 1:n){
      prec_lik = 1/sigma_e^2*(
        diag(1,m) - crossprod(t(Glist[[i]])*sqrt(1/(sigma_e^2/sigma_alpha^2 + 1)))
      )
      chQ_h = chol(1/sigma_h^2*Rinv + prec_lik)
      lin_h = crossprod(prec_lik, Y[i,])
      h_all[i,] = t(backsolve(chQ_h, forwardsolve(t(chQ_h), lin_h) + rnorm(m)))   
    }
    #----------------------------------------------------------------------------
    # Block 2: factors
    #----------------------------------------------------------------------------
    # First, sample the parametric factors Alpha:
    chQ_alpha = sqrt(sigma_e^-2 + rep(sigma_alpha^-2, each = n))
    lin_alpha = t(sapply(1:n, function(i) 
      crossprod(Glist[[i]], Y[i,] - h_all[i,])/sigma_e^2))
    Alpha = matrix(lin_alpha/chQ_alpha^2 + 1/chQ_alpha*rnorm(n*L), nrow = n)
    
    # And update:
    GAlpha = t(sapply(1:n, function(i)
      tcrossprod(Alpha[i,], Glist[[i]])))
    #----------------------------------------------------------------------------
    # Block 3: variance parameters
    #----------------------------------------------------------------------------
    
    # Update the fitted values:
    Yhat = GAlpha + h_all
    
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
    # Update the SD for the GPs:
    sigma_h = 1/sqrt(rgamma(n = 1,
                            shape = n*m/2 + 1/2,
                            rate = sum(diag(crossprod(tcrossprod(h_all, Rinv), h_all)))/2 + px_sigma_h))
    px_sigma_h = rgamma(n = 1, 
                        shape = 1/2 + 1/2, 
                        rate = 1/sigma_h^2 + 1)
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
############################################################
# Correlation function:
############################################################
corrFun = function(d,rho) (1 + sqrt(5)*d/rho + 5*d^2/(3*rho^2))*exp(-sqrt(5)*d/rho)
############################################################
# Compute absolute difference of "times" points in a matrix
############################################################
getTDiffs = function(times){
  tMat = matrix(rep(times,length(times)), nrow=length(times), byrow=FALSE)
  abs(tMat - t(tMat))
}
