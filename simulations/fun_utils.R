require(ellipse)

create_named_list <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}

get_list_corr_Zeta_univ_outlier <- function(N, L, p, vec_sd_zeta = NULL,
                                            list_vec_rho = rep(list(1/(2:(L+1))^0.2), p)) {
  
  if (is.null(vec_sd_zeta)) {
    vec_sd_zeta <- 1/(1:L)
  } else {
    stopifnot(length(vec_sd_zeta) == L)
  }
  
  list_Zeta <- vector("list", p)
  for (l in 1:L) {
    R <- matrix(sapply(list_vec_rho, function(ll) ll[l]), nrow = p,  ncol = p)
    diag(R) <- 1
    L_mat <- t(chol(R))
    tZ <- matrix(sapply(1:p, function(j) rnorm(N, 0, vec_sd_zeta[l])),
                 ncol = N, byrow = TRUE)
    Zeta_l <- as.matrix(t(L_mat %*% tZ))
    
    list_Zeta <- lapply(1:p, function(j) cbind(list_Zeta[[j]], Zeta_l[,j]))
  }
  
  list_Zeta
}


gauss_mfpca_data_outlier <- function(N, p, n, L, n_g, vec_sd_eps, mu_func, Psi_func,
                             time_obs = NULL, generate_from_univ = FALSE,
                             vec_sd_zeta = NULL, vec_rho_Zeta = NULL) {
  
  if (is.null(time_obs)) {
    time_obs <- vector("list", length = N)
    for(i in 1:N) {
      
      time_obs[[i]] <- vector("list", length = p)
      for(j in 1:p) {
        
        time_obs[[i]][[j]] <- sort(runif(n[i, j]))
      }
    }
  }
  
  time_g <- seq(0, 1, length.out = n_g)
  
  mu_g <- vector("list", length = p)
  Psi_g <- vector("list", length = p)
  for(j in 1:p) {
    mu_g[[j]] <- mu_func(time_g, j = j)
    Psi_g[[j]] <- Psi_func(time_g, j = j, p = p)
  }
  
  if (generate_from_univ) {
    
    if (is.null(vec_rho_Zeta) | all(unlist(vec_rho_Zeta) == 0)) {
      corr_Zeta <- FALSE
      list_Zeta <- NULL
    } else {
      corr_Zeta <- TRUE
      list_Zeta <- get_list_corr_Zeta_univ_outlier(N, L, p, vec_sd_zeta = vec_sd_zeta,
                                                   list_vec_rho = vec_rho_Zeta)
    }
    
    for (j in 1:p) {
      
      if (corr_Zeta) {
        Zeta_univ <- list_Zeta[[j]]
      } else {
        Zeta_univ <- get_Zeta(N, L, vec_sd_zeta = vec_sd_zeta)
        list_Zeta <- append(list_Zeta, list(Zeta_univ))
      }
      
      Y_univ <- get_Y_outlier(N, n, p, time_obs, Zeta_univ, vec_sd_eps, mu_func, Psi_func)
      if (j == 1) {
        Y <- lapply(Y_univ, function(ll) ll[1])
      } else {
        Y <- lapply(seq_along(Y), function(ii) append(Y[[ii]], Y_univ[[ii]][j]))
      }
    }
    Zeta <- list_Zeta
    
  } else {
    
    Zeta <- get_Zeta(N, L, vec_sd_zeta = vec_sd_zeta)
    
    Y <- get_Y_outlier(N, n, p, time_obs, Zeta, vec_sd_eps, mu_func, Psi_func)
    
  }
  
  create_named_list(time_obs, Zeta, time_g, mu_g, Psi_g, Y)
  
}

get_Y_outlier <- function(N, n, p, time_obs, Zeta, vec_sd_eps, mu_func, Psi_func) {
  
  Y <- vector("list", length = N)
  for(i in 1:N) {
    
    Y[[i]] <- vector("list", length = p)
    for(j in 1:p) {
      
      resid_vec <- rnorm(n[i, j], 0, vec_sd_eps[j])
      mean_vec <- mu_func(time_obs[[i]][[j]], j = j) +
        Psi_func(time_obs[[i]][[j]], j = j, p = p) %*% Zeta[i, ]
      Y[[i]][[j]] <- as.vector(mean_vec + resid_vec)
    }
  }
  
  Y
}

convert_data_for_happ <- function(time_obs, Y) {
  
  require(MFPCA)
  
  list_funData <- NULL
  for (j in 1:p) {
    
    time_obs_j <- sapply(time_obs, function(time_obs_i) time_obs_i[[j]])
    Y_j <- sapply(Y, function(Y_i) Y_i[[j]])
    ifd_j <- irregFunData(time_obs_j, Y_j)
    fd_j <- as.funData(ifd_j)
    list_funData <- append(list_funData, list(fd_j))
    
  }
  names(list_funData) <- paste0("Variable_", 1:p)
  mfd <- multiFunData(list_funData)
  mfd
  
}

cprod <- function(x, y) {
  
  if(missing(y)) {
    
    if(!is.vector(x)) {
      
      stop("Use the crossprod function for matrix inner products")
    }
    
    y <- x
  }
  
  if(!is.vector(y) & !is.vector(x)) {
    
    stop("Use the crossprod function for matrix inner products")
  }
  
  ans <- as.vector(crossprod(x, y))
  return(ans)
}

display_fit_list_happ <- function(p_sample, N_sample, 
                                  time_obs, time_g, time_g_happ, 
                                  Y,
                                  Y_hat, Y_low, Y_upp,
                                  Y_hat_happ, # Y_low_happ = NULL, Y_upp_happ = NULL, 
                                  offset = 0.1,
                                  col_data = "grey55", col = "black", col_happ = "blue",
                                  lwd = 1.2, lwd_happ = 1.2) {
  
  p_sample <- sort(p_sample)
  N_sample <- sort(N_sample)
  
  list_ylim <- list()
  for (j in p_sample) {
    
    vec_lim <- c(min(unlist(sapply(Y[N_sample], function(Y_i)Y_i[[j]])),
                     unlist(sapply(Y_low[N_sample], function(Y_i)Y_i[[j]]))
                     # ,
                     # unlist(sapply(Y_low_happ[N_sample], function(Y_i)Y_i[[j]]))
    ),
    max(unlist(sapply(Y[N_sample], function(Y_i)Y_i[[j]])),
        unlist(sapply(Y_upp[N_sample], function(Y_i)Y_i[[j]]))
        #                  ,
        #                  unlist(sapply(Y_low_happ[N_sample], function(Y_i)Y_i[[j]]))
    ))
    
    list_ylim <- append(list_ylim, list(c(vec_lim[1]-offset, vec_lim[2]+offset)))
  }
  
  
  par(mfrow = c(length(N_sample), length(p_sample)))
  for (i in N_sample) {
    if (i == N_sample[length(N_sample)]) {
      par(mar = c(4,4.5,1.5,1))
    } else if (i == N_sample[1]){
      par(mar = c(1.5,4.5,4,1))
    } else {
      par(mar = c(2.6,4.5,2.6,1))
    }
    jj <- 1
    for (j in p_sample) {
      plot(time_obs[[i]][[j]], Y[[i]][[j]],
           main = ifelse(i == N_sample[1], paste0("Variable ", j, "\n"), ""),
           xlab = ifelse(i == N_sample[length(N_sample)], "time", ""),
           ylab = ifelse(j == p_sample[1], parse(text=paste0("Y[", i, "]")), ""),
           pch = 20, col = col_data,
           xlim = c(0, 1),
           ylim = list_ylim[[jj]])
      
      lines(time_g, Y_hat[[i]][[j]], col=col, lwd = lwd)
      lines(time_g, Y_low[[i]][[j]], col=col,lwd = lwd,lty = 2)
      lines(time_g, Y_upp[[i]][[j]], col=col,lwd = lwd,lty = 2)
      
      lines(list_time_g_happ[[j]], Y_hat_happ[[i]][[j]], col=col_happ, lwd = lwd_happ)
      # lines(time_g, Y_low_happ[[i]][[j]], col=col_happ,lwd = lwd_happ,lty = 2)
      # lines(time_g, Y_upp_happ[[i]][[j]], col=col_happ,lwd = lwd_happ,lty = 2)
      
      jj <- jj +1
    }
    
  }
  
}


display_eigenfunctions_happ <- function(L, time_g, list_time_g_happ, 
                                        mu_g, Psi_g, 
                                        mu_hat, list_Psi_hat, 
                                        mu_hat_happ, list_Psi_hat_happ,
                                        lwd = 2, data_col = "red",
                                        vec_col_add = c("black", "blue"), 
                                        add_fct_lwd = rep(1.2, 2),
                                        vec_flip_happ = rep(1, L)) {
  
  par(mfcol = c(1+L, p))
  
  for (j in 1:p) {
    ylim_m <- c(min(c(mu_g[[j]], mu_hat[, j], mu_hat_happ[[j]])),
                max(c(mu_g[[j]], mu_hat[, j], mu_hat_happ[[j]])))
    
    plot(time_g, mu_g[[j]], type = "l", xlab = "time", 
         ylab = paste0("j = ", j), main = "Mean function",
         col = data_col[1], lwd = 2,
         ylim = ylim_m)
    lines(time_g, mu_hat[, j], col = vec_col_add[1], lwd = add_fct_lwd[1])
    lines(list_time_g_happ[[j]], mu_hat_happ[[j]], col = vec_col_add[2], lwd = add_fct_lwd[2])
    
    for (l in 1:L) {
      Psi_g_l <- sapply(Psi_g, function(Psi_g_j) Psi_g_j[,l])
      Psi_hat_l <- list_Psi_hat[[l]]
      Psi_hat_happ_l <- list_Psi_hat_happ[[l]]
      ylim_l <- c(min(c(Psi_g_l[,j], Psi_hat_l[,j], Psi_hat_happ_l[[j]])),
                  max(c(Psi_g_l[,j], Psi_hat_l[,j], Psi_hat_happ_l[[j]])))
      
      plot(time_g, Psi_g_l[,j], type = "l", xlab = "time", 
           ylab = paste0("j = ", j), main = paste0("Eigenfunction l = ", l),
           col = data_col[1], lwd = 2,
           ylim = ylim_l)
      lines(time_g, Psi_hat_l[, j], col = vec_col_add[1], lwd = add_fct_lwd[1])
      lines(list_time_g_happ[[j]], vec_flip_happ[l] * Psi_hat_happ_l[[j]], 
            col = vec_col_add[2], lwd = add_fct_lwd[2])
      
    }
    
  }
  
}



get_C <- function(time_obs, N, p, K, n_g) {
  
  # Set up fixed parameters
  
  time_vec <- unlist(time_obs)
  t_min <- 1.01*min(time_vec) - 0.01*max(time_vec)
  t_max <- 1.01*max(time_vec) - 0.01*min(time_vec)
  int_knots <- quantile(unique(time_vec), seq(0, 1, length = K)[-c(1, K)])
  
  C <- vector("list", length = N)
  for(i in 1:N) {
    
    C[[i]] <- vector("list", length = p)
    for(j in 1:p) {
      
      X <- X_design(time_obs[[i]][[j]])
      Z <- ZOSull(time_obs[[i]][[j]], range.x = c(0, 1), intKnots = int_knots)
      C[[i]][[j]] <- cbind(X, Z)
    }
  }
  
  # Set up plotting grid
  
  time_g <- seq(0, 1, length.out = n_g)
  
  X_g <- X_design(time_g)
  Z_g <- ZOSull(time_g, range.x = c(t_min, t_max), intKnots = int_knots)
  C_g <- cbind(X_g, Z_g)
  
  create_named_list(C, time_g, C_g)
}



# Post-inference orthonormalisation procedure for MCMC inference [deprecated].
#
# This function is used to orthonormalise the eigenfunctions and scores
# inferred using R Stan - deprecated, will be removed.
#
# @param stan_obj Object obtained from R Stan code for FPCA inference.
# @param C_g Design matrix C(t) constructed from the set of K spline functions
#            based on the dense time grid.
# @param Psi_g Reference eigenfunctions (if available, e.g., in simulations)
#              used to flip the sign of the resulting scores and eigenfunctions.
#
# @return An object containing the orthnormalised eigenfunctions and
#         uncorrelated scores.
#  [Note: not exported anymore, copied in VMP_FPCA/simulations/fun_utils.R ]
#
summarise_mcmc_multivariate <- function(stan_obj, C_g, Psi_g, L_sim = NULL, pred_interval = TRUE) {
  
  mcmc_samples <- rstan::extract(stan_obj, permuted=FALSE)
  n_mcmc <- dim(mcmc_samples)[1]
  time_g <- C_g[,2]
  n_g <- dim(C_g)[1]
  p <- length(Psi_g)
  
  fpca_params <- dimnames(mcmc_samples)$parameters
  L <- length(fpca_params[grep("beta_psi", fpca_params, fixed=TRUE)])/2/p
  if (is.null(L_sim)) {
    L_sim <- L
  }
  N <- length(fpca_params[grep("zeta", fpca_params, fixed=TRUE)])/L
  
  mu_g_mcmc <- Psi_g_mcmc <- vector("list", length=p)
  if (pred_interval) {
    sigma_eps_mcmc <- vector("list", length=p)
  }
  for (j in 1:p) {
    
    beta_mu_cols <- fpca_params[grep(paste0("beta_mu[", j), fpca_params, fixed=TRUE)]
    beta_mu_mcmc <- mcmc_samples[, 1, beta_mu_cols]
    
    u_mu_cols <- fpca_params[grep(paste0("u_mu[", j), fpca_params, fixed=TRUE)]
    u_mu_mcmc <- mcmc_samples[, 1, u_mu_cols]
    
    if (pred_interval) {
      sigma_eps_cols <- fpca_params[grep(paste0("sigma_eps[", j), fpca_params, fixed=TRUE)]
      sigma_eps_mcmc[[j]] <- mcmc_samples[, 1,  sigma_eps_cols]
    }
    
    nu_mu_mcmc <- t(cbind(beta_mu_mcmc, u_mu_mcmc))
    mu_g_mcmc[[j]] <- C_g%*%nu_mu_mcmc
    
    Psi_g_mcmc[[j]] <- vector("list", length=L)
    for(l in 1:L) {
      
      beta_psi_l_cols <- fpca_params[grep(paste0("beta_psi[", j, ",", l, ","), 
                                          fpca_params, fixed=TRUE)]
      beta_psi_mcmc <- mcmc_samples[, 1, beta_psi_l_cols]
      
      u_psi_l_cols <- fpca_params[grep(paste0("u_psi[", j, ",", l, ","), 
                                       fpca_params, fixed=TRUE)]
      u_psi_mcmc <- mcmc_samples[, 1, u_psi_l_cols]
      
      nu_psi_mcmc <- t(cbind(beta_psi_mcmc, u_psi_mcmc))
      Psi_g_mcmc[[j]][[l]] <- C_g%*%nu_psi_mcmc
    }
    
  }
  
  zeta_mcmc <- vector("list", length=N)
  for(i in 1:N) {
    
    zeta_i <- paste("zeta[", i, ",", sep="")
    zeta_i_cols <- fpca_params[grep(zeta_i, fpca_params, fixed=TRUE)]
    zeta_mcmc[[i]] <- mcmc_samples[, 1, zeta_i_cols]
  }
  
  one_N <- rep(1, N)
  mu_hat <- vector("list", length=n_mcmc)
  Psi_hat <- vector("list", length=n_mcmc)
  Zeta_hat <- vector("list", length=n_mcmc)
  if (pred_interval) {
    sigma_eps_hat <- vector("list", length=n_mcmc)
  }
  for(k in 1:n_mcmc) {
    
    Psi <- vector("list", length=p)
    for (j in 1:p) {
      
      Psi[[j]] <- matrix(NA, n_g, L)
      for(l in 1:L) {
        
        Psi[[j]][,l] <- Psi_g_mcmc[[j]][[l]][,k]
      }
      
    }
    
    Zeta <- matrix(NA, N, L)
    for(i in 1:N) {
      
      Zeta[i,] <- zeta_mcmc[[i]][k,]
    }
    
    # Orthogonalisation:
    #
    mu_g_mcmc_k <- Reduce(c, lapply(mu_g_mcmc, function(ll) ll[,k]))
    Psi_k <- Reduce(rbind, Psi)
    if (pred_interval) {
      sigma_eps_hat[[k]] <- Reduce(c, lapply(sigma_eps_mcmc, function(ll) ll[k]))
    }
    
    svd_Psi <- svd(Psi_k)
    U_psi <- svd_Psi$u
    D_psi <- diag(svd_Psi$d)
    V_psi <- svd_Psi$v
    
    zeta_rotn <- t(Zeta %*% V_psi %*% D_psi)
    C_zeta <- cov(t(zeta_rotn))
    eigen_C <- eigen(C_zeta)
    Q <- eigen_C$vectors
    Lambda <- diag(eigen_C$values)
    
    Psi_tilde <- U_psi %*% Q %*% sqrt(Lambda)
    Zeta_tilde <- crossprod(zeta_rotn, Q %*% solve(sqrt(Lambda)))
    
    mu_hat[[k]] <- split(mu_g_mcmc_k, rep(1:p, each = n_g))
    
    Psi_hat[[k]] <- matrix(NA, p*n_g, L)
    Zeta_hat[[k]] <- matrix(NA, N, L)
    norm_vec <- rep(NA, L)
    time_int_vec <- seq(0, p, length = n_g*p)
    for(l in 1:L) {
      
      norm_vec[l] <- sqrt(trapint(time_int_vec, Psi_tilde[, l]^2))
      Psi_hat[[k]][, l] <- Psi_tilde[, l]/norm_vec[l]
      Zeta_hat[[k]][, l] <- norm_vec[l]*Zeta_tilde[, l]
    }
    
    for(l in 1:L_sim) {
      Psi_g_comb <- vector("list", length = p)
      for(j in 1:p) {
        
        Psi_g_comb[[j]] <- Psi_g[[j]][, l]
      }
      Psi_g_comb <- Reduce(c, Psi_g_comb)
      
      inner_prod_sign <- sign(cprod(Psi_g_comb, Psi_hat[[k]][, l]))
      if(inner_prod_sign == -1) {
        
        Psi_hat[[k]][, l] <- -Psi_hat[[k]][, l]
        Zeta_hat[[k]][, l] <- -Zeta_hat[[k]][, l]
      }
    }
    Psi_hat[[k]] <- lapply(split(Psi_hat[[k]], rep(1:p, each = n_g)), 
                           matrix, nrow = n_g, ncol = L)
  }
  
  
  # Summarise the MCMC outputs:
  #
  Y_g_mcmc_summary <- vector("list", length=N)
  for(i in 1:N) {
    
    Y_g_mcmc_summary[[i]] <- vector("list", length=p)
    
    for (j in 1:p) {
      
      Y_g_mcmc <- mean_mcmc <- matrix(NA, n_g, n_mcmc)
      for(k in 1:n_mcmc) {
        
        mean_mcmc[,k] <- mu_hat[[k]][[j]] + Psi_hat[[k]][[j]]%*%Zeta_hat[[k]][i,]
        
        if (pred_interval) {
          Y_g_mcmc[,k] <- mean_mcmc[,k] + rnorm(1, mean = 0, sd = sigma_eps_hat[[k]][[j]]) 
        } else {
          Y_g_mcmc[,k] <- mean_mcmc[,k]
        }
        
      }
      
      Y_g_mcmc_summary[[i]][[j]] <- matrix(NA, nrow=n_g, ncol=3)
      Y_g_mcmc_summary[[i]][[j]][,1] <- apply(Y_g_mcmc, 1, quantile, 0.025)
      Y_g_mcmc_summary[[i]][[j]][,2] <- apply(Y_g_mcmc, 1, mean)
      Y_g_mcmc_summary[[i]][[j]][,3] <- apply(Y_g_mcmc, 1, quantile, 0.975)
      
    }
    
    
  }
  
  zeta_mcmc_summary <- vector("list", length=N)
  for(i in 1:N) {
    
    zeta_mcmc_i <- matrix(NA, n_mcmc, 2)
    for(k in 1:n_mcmc) {
      
      zeta_mcmc_i[k,] <- Zeta_hat[[k]][i,1:2]
    }
    
    zeta_mcmc_mean <- apply(zeta_mcmc_i, 2, mean)
    zeta_mcmc_cov <- cov(zeta_mcmc_i)
    
    zeta_mcmc_ellipse <- ellipse(
      zeta_mcmc_cov,
      centre=zeta_mcmc_mean,
      level=0.95
    )
    
    zeta_mcmc_summary[[i]] <- list(zeta_mcmc_mean, zeta_mcmc_ellipse)
    names(zeta_mcmc_summary[[i]]) <- c("mean", "credible boundary")
  }
  
  gbl_mcmc_summary <- vector("list", length = L + 1)
  gbl_mcmc_summary[[1]] <- Reduce(cbind, lapply(1:p, function(j) 
    apply(Reduce(cbind, lapply(mu_hat, function(mu_hat_k) mu_hat_k[[j]])), 1, mean)))
  
  for(l in 1:L) {
    gbl_mcmc_summary[[l+1]] <- matrix(NA, n_g, p)
    for(j in 1:p) {
      gbl_mcmc_summary[[l+1]][, j] <- Reduce("+", lapply(Psi_hat, function(Psi_hat_k) Psi_hat_k[[j]][,l]))/n_mcmc
    }
  }
  
  
  gbl_ci_mcmc_summary <- vector("list", length = L + 1)
  gbl_ci_mcmc_summary[[1]] <- vector("list", length = p)
  for(j in 1:p) {
    gbl_ci_mcmc_summary[[1]][[j]] <- t(apply(simplify2array(lapply(mu_hat, function(mu_hat_k) mu_hat_k[[j]])), 1,
                                             quantile, prob = c(0.025, 0.975)))
  }
  
  for(l in 1:L) {
    gbl_ci_mcmc_summary[[l+1]] <- vector("list", length = p) # matrix(NA, n_g, p)
    for(j in 1:p) {
      gbl_ci_mcmc_summary[[l+1]][[j]] <- t(apply(simplify2array(lapply(Psi_hat, function(Psi_hat_k) Psi_hat_k[[j]][,l])), 1,
                                                 quantile, prob = c(0.025, 0.975)))
    }
  }
  
  # Summary outputs:
  
  outputs <- list(Y_g_mcmc_summary, gbl_mcmc_summary, gbl_ci_mcmc_summary, zeta_mcmc_summary)
  names(outputs) <- c("Y_g_mcmc_summary", "gbl_mcmc_summary", "gbl_ci_mcmc_summary", "zeta_mcmc_summary")
  
  return(outputs)
}


display_scores_custom_axes <- function(l_ind, N_sample, Zeta, Zeta_hat, zeta_ellipse,
                                       Zeta_hat_add = NULL, zeta_ellipse_add = NULL,
                                       vec_col = c("black", "blue"), data_col = "red",
                                       mfrow = NULL) {
  
  n_sample <- length(N_sample)
  
  par(mfrow = c(1, 1))
  zeta_labels <- vector("list", length = n_sample)
  zeta_id <- rep(NA, n_sample)
  
  if(is.null(mfrow)) {
    mfrow <- c(floor(sqrt(n_sample)), ceiling(sqrt(n_sample)))
  }
  
  for(i in 1:n_sample) {
    
    N_i <- N_sample[i]
    
    zeta_id[i] <- parse(text=paste("zeta[", N_i, "]", sep=""))
    zeta_val <- eval(bquote(expression(zeta[.(N_i)])))
    zeta_labels[[i]] <- rep(zeta_val, nrow(zeta_ellipse[[N_i]]))
  }
  zeta_labels <- do.call(c, zeta_labels)
  zeta_labels <- factor(zeta_labels, levels=zeta_id)
  
  strip.math <- function(
    which.given, which.panel, var.name, factor.levels, ...
  ) {
    
    fl <- zeta_id
    
    strip.default(which.given,which.panel,var.name,fl,...)
  }
  
  zeta_ellipse_mat <- Reduce(rbind, zeta_ellipse[N_sample])
  zeta_ellipse_x <- zeta_ellipse_mat[,1]
  zeta_ellipse_y <- zeta_ellipse_mat[,2]
  
  
  if (!is.null(Zeta_hat_add)) {
    zeta_ellipse_add_mat <- Reduce(rbind, zeta_ellipse_add[N_sample])
    zeta_ellipse_add_x <- zeta_ellipse_add_mat[,1]
    zeta_ellipse_add_y <- zeta_ellipse_add_mat[,2]
  }
  
  
  if (is.list(Zeta)) {
    xtmp_z <- as.vector(sapply(N_sample, function(ns) c(zeta_ellipse[[ns]][,1],
                                                        zeta_ellipse[[ns]][,3],
                                                        unlist(lapply(Zeta, function(zz) zz[ns, 1])),
                                                        unlist(lapply(Zeta, function(zz) zz[ns, 3])),
                                                        Zeta_hat[ns,1],
                                                        Zeta_hat[ns,3])))
    ytmp_z <- as.vector(sapply(N_sample, function(ns) c(zeta_ellipse[[ns]][,2],
                                                        zeta_ellipse[[ns]][,4],
                                                        unlist(lapply(Zeta, function(zz) zz[ns, 2])),
                                                        unlist(lapply(Zeta, function(zz) zz[ns, 4])),
                                                        Zeta_hat[ns,2],
                                                        Zeta_hat[ns,4])))
    
  } else {
    xtmp_z <- as.vector(sapply(N_sample, function(ns) c(zeta_ellipse[[ns]][,1],
                                                        Zeta[ns, 1],
                                                        Zeta_hat[ns,1])))
    ytmp_z <- as.vector(sapply(N_sample, function(ns) c(zeta_ellipse[[ns]][,2],
                                                        Zeta[ns, 2],
                                                        Zeta_hat[ns,2])))
  }
  xlim_z <- c(min(xtmp_z), max(xtmp_z))
  ylim_z <- c(min(ytmp_z), max(ytmp_z))
  
  Zeta <- Zeta[,l_ind]
  Zeta_hat <- Zeta_hat[,l_ind]
  Zeta_hat_add <- Zeta_hat_add[,l_ind]
  
  
  score_plots <- xyplot(
    zeta_ellipse_y ~ zeta_ellipse_x | zeta_labels, groups = zeta_labels,
    data=data.frame(
      zeta_ellipse_x = zeta_ellipse_x, zeta_ellipse_y = zeta_ellipse_y,
      zeta_labels = zeta_labels
    ),
    layout=mfrow, main="",
    strip=strip.math,
    xlim = 1.1*xlim_z,
    ylim = 1.1*ylim_z,
    par.strip.text=list(cex=0.8),
    par.settings = list(layout.heights = list(strip = 1),
                        strip.background=list(col="white")),
    # key=list(space="top",
    #          lines=list(col=c("grey55","blue"), lty=1, lwd=1),
    #          text=list(c("Simulated scores", "Estimated scores"))),
    xlab=paste0("Scores FPC ", l_ind[1]),
    ylab=paste0("Scores FPC ", l_ind[2]),
    as.table=TRUE,
    panel=function(x, y, subscripts, groups) {
      
      iPan <- panel.number()
      i <- rep(1:n_sample, each=1)[iPan]
      # panel.grid()
      panel.xyplot(
        zeta_ellipse[[N_sample[i]]][,1], zeta_ellipse[[N_sample[i]]][,2],
        col=vec_col[1], type="l", lwd=1.5
      )
      panel.xyplot(
        Zeta_hat[N_sample[i],1], Zeta_hat[N_sample[i],2],
        col= vec_col[1], type="p", pch=16, cex=0.7
      )
      
      if (is.list(Zeta)) {
        for (j in 1:length(Zeta)) {
          panel.xyplot(
            Zeta[[j]][N_sample[i], 1], Zeta[[j]][N_sample[i], 2],
            # col=grDevices::adjustcolor(data_col, alpha.f = 1/j^0.9),
            col = data_col, type="p", pch=16, cex=0.7
          )
        }
        
      } else {
        panel.xyplot(
          Zeta[N_sample[i], 1], Zeta[N_sample[i], 2],
          col=data_col, type="p", pch=16, cex=0.7
        )
      }
      
      if (!is.null(Zeta_hat_add)) {
        panel.xyplot(
          zeta_ellipse_add[[N_sample[i]]][,1], zeta_ellipse_add[[N_sample[i]]][,2],
          col=vec_col[2], type="l", lwd=1.5
        )
        panel.xyplot(
          Zeta_hat_add[N_sample[i],1], Zeta_hat_add[N_sample[i],2],
          col= vec_col[2], type="p", pch=16, cex=0.7
        )
      }
    }
  )
  
  print(score_plots)
}

