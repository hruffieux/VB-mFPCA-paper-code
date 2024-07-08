rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

bool_cluster <- T
if (bool_cluster) {
  
  #!/usr/bin/env Rscript
  args = commandArgs(trailingOnly=TRUE)
  job_id <- as.integer(args[1]) # 1 to 6
  
  CORE_DIR_MRC_BSU <- Sys.getenv("CORE_DIR_MRC_BSU")
  out_dir <- file.path(CORE_DIR_MRC_BSU, "mFPCA_output/")
} else {
  
  job_id <- 2 # 1 to 8
  
  CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")
  out_dir <- file.path(CORE_DIR_ICLOUD, "mFPCA_output/")
}

main_dir <- file.path(CORE_DIR, "bayesian-mFPCA-paper-code/simulations/")
setwd(main_dir)

require(bayesFPCA)

source("fun_utils.R")

seed <- 1
set.seed(seed)

# Establish simulation variables:

p <- 3                                        # number of responses
N <- 100                                      # number of curves
N_t_min <- 15                                 # minimum number of time observations for each curve
N_t_max <- 25                                 # maximum number of time observations for each curve
n <- matrix(sample(N_t_min:N_t_max, N*p,
                   replace = TRUE), N, p)     # number of time observations

L_max <- 10
lambda_L <- 1

choose_both <- T #
if (choose_both) {                            # model choice for both K and L
  choose_K_pve <- T 
  K <-  NULL
  K_min <- 5
  K_max <- 20
  unif_dist_K <- T 
  if (!unif_dist_K) {
    lambda_K <- 10                            # truncated Poisson distribution
  } else {
    lambda_K <- NULL                          # discrete uniform distribution on Unif[K_min:K_max]
  }
} else {
  choose_K_pve <- F
  bool_K_Ruppert <- T
  if (bool_K_Ruppert) {
    K <- NULL                                 # default based on the rule in Ruppert 2002
  } else {
    n_int_knots <- 5                          # number of interior knots
    K <- n_int_knots + 2                      # number of spline basis functions
  }
  K_min <- K_max <- NULL
  lambda_K <- NULL
}


L_sim_max <- 8
L_sim <- c(1:L_sim_max)[job_id] 
stopifnot(L_sim <= L_max)

tol  <- 1e-5                                  
maxit <- 1000                             

n_g <- 1000                               

sd_eps <- 1
vec_sd_eps <- rep(sd_eps, p)

ind_exp <- 3
exponent_sd_zeta <- c(-1/2, -1/4, -1/8, -1/16)[ind_exp]
vec_sd_zeta <- (1:L_sim)^exponent_sd_zeta 

bool_orth_splines <- T

generate_from_univ <- T

if (generate_from_univ) {
  id_strength <- 6 # up to 6
  rho_Zeta <- seq(0, 1, by = 0.2)[id_strength]
  vec_rho_Zeta <- rep(rho_Zeta, L_sim)
  cat("Correlation of simulated scores across variables for eigenfunctions 1 to L_sim: ",
      vec_rho_Zeta, "\n")
} else {
  vec_rho_Zeta <- NULL
}


if (rho_Zeta == 1) { # perfect correlation, so we can directly simulate from the multivariate model
  generate_from_univ <- FALSE 
  vec_rho_Zeta <- NULL
}


n_repl <- 100
n_cpus <- 50

bool_save <- T
if (bool_save) {
  res_dir <- paste0(out_dir, "/model_choice_n_repl_", n_repl, 
                    ifelse(unif_dist_K, "_unif_dist_K", ""), "_", 
                    ifelse(choose_both, "L_and_K", "L"), 
                    ifelse(choose_K_pve, "_choose_K_pve", ""), "_orth_Bsplines_",  
                    ifelse(bool_orth_splines, "T", "F"), "_p_", p, "_N_", N, 
                    "_Nt_min_", N_t_min, "_max_", N_t_max, "_K_min_", K_min, 
                    "_K_max_", K_max, "_L_max_", L_max, 
                    ifelse(choose_both, paste0("_lambda_K_", lambda_K), 
                           paste0("_K_", K)), "_lambda_L_", lambda_L, 
                    "_exponent_sd_zeta_", 
                    paste0(format(exponent_sd_zeta , digits = 2), collapse = "-"),
                    "_sd_eps_", unique(vec_sd_eps), "_tol_", tol, 
                    "_maxit_", maxit, "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}

f_mu <- function(t, j) (-1)^j*2*sin((2*pi+j)*t)

f_psi_1 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(pi*t)
f_psi_2 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(pi*t)
f_psi_3 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(2*pi*t)
f_psi_4 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(2*pi*t)
f_psi_5 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(3*pi*t)
f_psi_6 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(3*pi*t)
f_psi_7 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(4*pi*t)
f_psi_8 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(4*pi*t)
f_psi_9 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(5*pi*t)
f_psi_10 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(5*pi*t)


f_Psi_periodic <- function(time_obs, j, p) {
  ans <- cbind(f_psi_1(time_obs, j, p), 
               f_psi_2(time_obs, j, p),
               f_psi_3(time_obs, j, p),
               f_psi_4(time_obs, j, p),
               f_psi_5(time_obs, j, p),
               f_psi_6(time_obs, j, p),
               f_psi_7(time_obs, j, p),
               f_psi_8(time_obs, j, p),
               f_psi_9(time_obs, j, p),
               f_psi_10(time_obs, j, p))[, (L_max-L_sim+1):L_max, drop = FALSE] 
  return(ans)
}


f_Psi_orthBsplines <- function(time_obs, j, p) {
  
  require(Splinets) 
  k <- 2 # order 
  n_knots <- max(L_sim_max + 1, 3) 
  knots <- seq(0, 1, length.out = n_knots + 2)
  so <- splinet(knots, smorder = k, norm = TRUE) 
  
  splines <- evspline(so$os, x = time_obs)[,-1, drop = FALSE] 
  stopifnot(ncol(splines)>=L_sim) 
  
  ans <- (-1)^j * sqrt(1/p) * splines[, (L_sim_max-L_sim+1):L_sim_max, drop = FALSE] # renormalise 
  
  return(ans)
  
}

f_Psi <- ifelse(bool_orth_splines, f_Psi_orthBsplines, f_Psi_periodic)

vec_L <- 1:L_max
if (choose_both) {
  vec_K <- K_min:K_max # 2:K_max # 16 cores if between 5 and 20
} else {
  vec_K <- NULL
}

out <- parallel::mclapply(1:n_repl, function(repl) {
  
  set.seed(repl)
  
  # Set up the data:
  #
  data <- generate_fpca_data(N, p, n, L_sim, n_g, vec_sd_eps, f_mu, f_Psi,
                             vec_sd_zeta = vec_sd_zeta,
                             generate_from_univ = generate_from_univ,
                             vec_rho_Zeta = vec_rho_Zeta)
  
  time_obs <- data$time_obs
  Zeta <- data$Zeta
  Y <- data$Y
  
  time_g <- data$time_g
  mu_g <- data$mu_g
  Psi_g <- data$Psi_g
  
  if (choose_both) {
    all_cum_pve <- all_rmse <- all_ise <- replicate(length(vec_L), vector("list", length(vec_K)), simplify = FALSE)
    all_elbos <- all_log_unnorm_p_model_given_y <- matrix(NA, nrow = length(vec_L), ncol = length(vec_K))
    names(all_cum_pve) <- names(all_rmse) <- names(all_ise) <- paste0("L_", vec_L)
    all_cum_pve <- all_rmse <- all_ise <- lapply(all_cum_pve, function(ll) {names(ll) <-  paste0("K_", vec_K); ll})
    rownames(all_elbos) <- rownames(all_log_unnorm_p_model_given_y) <- paste0("L_", vec_L)
    colnames(all_elbos) <- colnames(all_log_unnorm_p_model_given_y) <- paste0("K_", vec_K)
  } else {
    all_cum_pve <- all_rmse <- all_ise <- vector("list", length(vec_L))
    all_elbos <- all_log_unnorm_p_model_given_y <- rep(NA, length(vec_L))
    names(all_cum_pve) <- names(all_rmse) <- names(all_ise) <- names(all_elbos) <- names(all_log_unnorm_p_model_given_y) <- paste0("L_", vec_L)
  }
  
  for (ind_ll in seq_along(vec_L)) {
    
    ll <- vec_L[ind_ll]
    
    if (choose_both) {
      
      for(ind_kk in seq_along(vec_K)) {
        
        kk <- vec_K[ind_kk]
        
        mfvb_res <- run_mfvb_fpca(time_obs, Y, L = ll, K = kk, n_g = NULL, 
                                  time_g = time_g, 
                                  tol = tol, maxit = maxit,
                                  plot_elbo = F, Psi_g = NULL)
      
        eigenvalues <- apply(mfvb_res$Zeta_hat, 2, function(vv) var(vv)) 
        all_cum_pve[[ind_ll]][[ind_kk]] <- cumsum(eigenvalues / sum(eigenvalues) * 100)
        
        all_elbos[ind_ll, ind_kk] <- mfvb_res$elbo
        
        if (unif_dist_K) {
          all_log_unnorm_p_model_given_y[ind_ll, ind_kk] <- mfvb_res$elbo + ll * log(lambda_L) - lfactorial(ll) - lambda_L - log(length(vec_K)) 
        } else {
          all_log_unnorm_p_model_given_y[ind_ll, ind_kk] <- mfvb_res$elbo + ll * log(lambda_L) - lfactorial(ll) - lambda_L + kk * log(lambda_K) - lfactorial(kk) - lambda_K
        }
        
        rmse_mfpca <- apply(Zeta[,1:min(ll, L_sim), drop = F] - mfvb_res$Zeta[,1:min(ll, L_sim), drop = F], 2, function(x) sqrt(mean(x^2))) # one error per FPC (vector of length L)
        if (ll < L_sim) {
          rmse_mfpca <- c(rmse_mfpca, rep(NA, L_sim - ll))
        }
        names(rmse_mfpca) <- paste0("FPC_", 1:L_sim)
        
        all_rmse[[ind_ll]][[ind_kk]] <- rmse_mfpca
        
        
        ise_mfpca <- matrix(NA, L_sim+1, p)
        for (l in 1:(min(ll, L_sim)+1)) {
          for(j in 1:p) {
            if (l == 1) {
              ise_mfpca[l, j] <- trapint(time_g, (mfvb_res$mu_hat[,j] - mu_g[[j]])^2)
            } else {
              ise_mfpca[l, j] <- trapint(time_g, (mfvb_res$list_Psi_hat[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
            }
          }
        }
        rownames(ise_mfpca) <- c("mu", paste0("Psi_", 1:L_sim))
        colnames(ise_mfpca) <- paste0("Variable_", 1:p)
        
        
        all_ise[[ind_ll]][[ind_kk]] <- ise_mfpca

      }

    } else {
      
      mfvb_res <- run_mfvb_fpca(time_obs, Y, L = ll, K = K, n_g = NULL, 
                                time_g = time_g, # here we use the same time grid as that used to simulate the true mean and eigen- functions
                                tol = tol, maxit = maxit,
                                plot_elbo = F, Psi_g = NULL)
      
      eigenvalues <- apply(mfvb_res$Zeta_hat, 2, function(vv) var(vv)) 
      all_cum_pve[[ind_ll]] <- cumsum(eigenvalues / sum(eigenvalues) * 100)
      
      all_elbos[ind_ll] <- mfvb_res$elbo
      all_log_unnorm_p_model_given_y[ind_ll] <- mfvb_res$elbo + ll * log(lambda_L) - lfactorial(ll) - lambda_L 
      
      rmse_mfpca <- apply(Zeta[,1:min(ll, L_sim), drop = F] - mfvb_res$Zeta[,1:min(ll, L_sim), drop = F], 2, function(x) sqrt(mean(x^2))) # one error per FPC (vector of length L)
      if (ll < L_sim) {
        rmse_mfpca <- c(rmse_mfpca, rep(NA, L_sim - ll))
      }
      names(rmse_mfpca) <- paste0("FPC_", 1:L_sim)
      
      all_rmse[[ind_ll]] <- rmse_mfpca
      
      ise_mfpca <- matrix(NA, L_sim+1, p)
      for (l in 1:(min(ll, L_sim)+1)) {
        for(j in 1:p) {
          if (l == 1) {
            ise_mfpca[l, j] <- trapint(time_g, (mfvb_res$mu_hat[,j] - mu_g[[j]])^2)
          } else {
            ise_mfpca[l, j] <- trapint(time_g, (mfvb_res$list_Psi_hat[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
          }
        }
      }
      rownames(ise_mfpca) <- c("mu", paste0("Psi_", 1:L_sim))
      colnames(ise_mfpca) <- paste0("Variable_", 1:p)
      
      all_ise[[ind_ll]] <- ise_mfpca
      
    }
  
  }
  
  
  all_unnorm_p_model_given_y <- exp(all_log_unnorm_p_model_given_y - max(all_log_unnorm_p_model_given_y)) 
  
  all_p_model_given_y <- all_unnorm_p_model_given_y / sum(all_unnorm_p_model_given_y)

  
  if (choose_both && choose_K_pve) {
    all_unnorm_p_model_given_y_L_10 <- exp(all_log_unnorm_p_model_given_y[paste0("L_", L_max), ] - max(all_log_unnorm_p_model_given_y[paste0("L_", L_max),])) 
    
    all_p_model_given_y_L_10 <- all_unnorm_p_model_given_y_L_10 / sum(all_unnorm_p_model_given_y_L_10)
    
    K_model_choice_L_10 <- vec_K[which.max(all_p_model_given_y_L_10)]
    all_cum_pve_model_choice_K <- all_cum_pve[[paste0("L_", L_max)]][[paste0("K_", K_model_choice_L_10)]]
    
    
    L_model_choice <- vec_L[which(all_p_model_given_y== max(all_p_model_given_y),  # Extract row name of min
                                  arr.ind = TRUE)[ , 1]]
    K_model_choice <- vec_K[which(all_p_model_given_y== max(all_p_model_given_y),  # Extract row name of min
                                  arr.ind = TRUE)[ , 2]]
    
  } else {
    all_cum_pve_model_choice_K <- K_model_choice_L_10 <- L_model_choice <- K_model_choice <- NULL
  }
  
  create_named_list(all_p_model_given_y, all_elbos, 
                    all_cum_pve, all_cum_pve_model_choice_K, K_model_choice_L_10,
                    all_rmse, all_ise, L_model_choice, K_model_choice)

  
}, mc.cores = n_cpus)

names(out) <- paste0("repl_", 1:n_repl)

list_all_p_model_given_y <- lapply(out, "[[", "all_p_model_given_y")
list_all_cum_pve <- lapply(out, "[[", "all_cum_pve")
list_all_cum_pve_model_choice_K <- lapply(out, "[[", "all_cum_pve_model_choice_K")
list_all_elbos <- lapply(out, "[[", "all_elbos")
list_all_rmse <- lapply(out, "[[", "all_rmse")
list_all_ise <- lapply(out, "[[", "all_ise")
vec_K_model_choice <- sapply(out, "[[", "K_model_choice")
vec_K_model_choice_L_10 <- sapply(out, "[[", "K_model_choice_L_10")
vec_L_model_choice <- sapply(out, "[[", "L_model_choice")

if (choose_both) {
  
  mean_all_p_model_given_y <- apply(simplify2array(list_all_p_model_given_y), 1:2, mean)
  sd_all_p_model_given_y <- apply(simplify2array(list_all_p_model_given_y), 1:2, sd)
  
  if (!choose_K_pve) {
    K_pve <- median(sapply(1:p, function(j) max(round(min(median(n[,j]/4), 40)), 7)))
    
    mat_all_cum_pve <- do.call(rbind, lapply(lapply(list_all_cum_pve, "[[", paste0("L_", L_max)), "[[", paste0("K_", K_pve)))
    
    add_label <- paste0(" with K = ", K_pve)
    
  } else {
    mat_all_cum_pve <- do.call(rbind, list_all_cum_pve_model_choice_K)
    
    add_label <- " with model choice for K"
    
    if (bool_save) {
      pdf(paste0(res_dir, "/boxplot_chosen_K_for_L_10_across_replicates_for_L_sim_", L_sim, ".pdf"),
          width = 5, height = 5, paper='special')
    }
    
    boxplot(vec_K_model_choice_L_10, main ="Chosen K for L_10 (pve) across replicates")
    
    if (bool_save) {
      dev.off()
    }
  }
  
  
  mean_list_rmse <- sd_list_rmse <- vector("list", L_sim)
  for (l_sim in 1:L_sim) {
    mean_list_rmse[[l_sim]] <- sd_list_rmse[[l_sim]] <- matrix(NA, nrow = length(vec_L), ncol = length(vec_K))
    rownames(mean_list_rmse[[l_sim]]) <- rownames(sd_list_rmse[[l_sim]]) <- paste0("L_", vec_L)
    colnames(mean_list_rmse[[l_sim]]) <- colnames(sd_list_rmse[[l_sim]]) <- paste0("K_", vec_K)
    for (ind_ll in seq_along(vec_L)) {
      for (ind_kk in seq_along(vec_K)) {
        mean_list_rmse[[l_sim]][ind_ll, ind_kk] <- mean(sapply(list_all_rmse, function(all_rmse) all_rmse[[ind_ll]][[ind_kk]][l_sim]))
        sd_list_rmse[[l_sim]][ind_ll, ind_kk] <- sd(sapply(list_all_rmse, function(all_rmse) all_rmse[[ind_ll]][[ind_kk]][l_sim]))
      } 
    }
  }
  names(mean_list_rmse) <- names(sd_list_rmse) <- paste0("FPC_", 1:L_sim)
  
  
  mean_list_ise <- sd_list_ise <- vector("list", L_sim+1)
  for (ind in 1:(L_sim+1)) {
    mean_list_ise[[ind]] <- sd_list_ise[[ind]] <- matrix(NA, nrow = length(vec_L), ncol = length(vec_K))
    rownames(mean_list_ise[[ind]]) <- rownames(sd_list_ise[[ind]]) <- paste0("L_", vec_L)
    colnames(mean_list_ise[[ind]]) <- colnames(sd_list_ise[[ind]]) <- paste0("K_", vec_K)
    for (ind_ll in seq_along(vec_L)) {
      for (ind_kk in seq_along(vec_K)) {
        mean_list_ise[[ind]][ind_ll, ind_kk] <- mean(sapply(list_all_ise, function(all_ise) mean(all_ise[[ind_ll]][[ind_kk]][ind,])))
        sd_list_ise[[ind]][ind_ll, ind_kk] <- sd(sapply(list_all_ise, function(all_ise) mean(all_ise[[ind_ll]][[ind_kk]][ind,])))
      } 
    }
  }
  names(mean_list_ise) <- names(sd_list_ise) <- c("mu", paste0("Psi_", 1:L_sim))

  
} else {
  
  mat_all_cum_pve <- do.call(rbind, lapply(list_all_cum_pve, "[[", paste0("L_", L_max)))
  
  mean_all_cum_pve <- apply(simplify2array(list_cum_pve), 1, mean)
  
  mat_all_p_model_given_y <- do.call(rbind, list_all_p_model_given_y)
  mean_all_p_model_given_y <- apply(simplify2array(list_all_p_model_given_y), 1, mean)
  sd_all_p_model_given_y <- apply(simplify2array(list_all_p_model_given_y), 1, sd)
  
  print(paste0("Simulated L_sim: ", L_sim))
  print("Inferred: ")
  names(mean_all_p_model_given_y)[which(mean_all_p_model_given_y == max(mean_all_p_model_given_y))]
  
  add_label <- NULL
  
  
  mean_list_rmse <- sd_list_rmse <- vector("list", L_sim)
  for (l_sim in 1:L_sim) {
    mean_list_rmse[[l_sim]] <- sd_list_rmse[[l_sim]] <- rep(NA, length(vec_L))
    names(mean_list_rmse[[l_sim]]) <- names(sd_list_rmse[[l_sim]]) <- paste0("L_", vec_L)
    for (ind_ll in seq_along(vec_L)) {
        mean_list_rmse[[l_sim]][ind_ll] <- mean(sapply(list_all_rmse, function(all_rmse) all_rmse[[ind_ll]][l_sim]))
        sd_list_rmse[[l_sim]][ind_ll] <- sd(sapply(list_all_rmse, function(all_rmse) all_rmse[[ind_ll]][l_sim]))
    }
  }
  names(mean_list_rmse) <- names(sd_list_rmse) <- paste0("FPC_", 1:L_sim)
  
  
  mean_list_ise <- sd_list_ise <- vector("list", L_sim+1)
  for (ind in 1:(L_sim+1)) {
    mean_list_ise[[ind]] <- sd_list_ise[[ind]] <- rep(NA, length(vec_L))
    names(mean_list_ise[[ind]]) <- names(sd_list_ise[[ind]]) <- paste0("L_", vec_L)
    for (ind_ll in seq_along(vec_L)) {
      for (ind_kk in seq_along(vec_K)) {
        mean_list_ise[[ind]][ind_ll] <- mean(sapply(list_all_ise, function(all_ise) mean(all_ise[[ind_ll]][ind,])))
        sd_list_ise[[ind]][ind_ll] <- sd(sapply(list_all_ise, function(all_ise) mean(all_ise[[ind_ll]][ind,])))
      } 
    }
  }
  names(mean_list_ise) <- names(sd_list_ise) <- c("mu", paste0("Psi_", 1:L_sim))
}



if (choose_both) {
  
require(gplots)

  # Define a custom color palette transitioning from white to black
  my_palette <- colorRampPalette(c("white", "black"))(n = 500)
  #initiate cols with all black
  
  cols <- rep('black', nrow(mean_all_p_model_given_y))
  #turn red the specified rows in tf
  cols[rownames(mean_all_p_model_given_y) == paste0("L_", L_sim)] <- 'red'
  
  if (bool_save) {
    pdf(paste0(res_dir, "/heatmap_prob_L_sim_", L_sim, ".pdf"),
        width = 8, height = 7, paper='special')
  }
  # Adjust the size of the heatmap
  par(mar = c(5, 4, 4, 8)) # Adjust the margins as needed
  
  # Create the heatmap with the modified color palette
  heatmap.2(mean_all_p_model_given_y,
            dendrogram='none',
            Rowv=FALSE,
            Colv=FALSE,
            trace='none',
            main = "Model choice probabilities",
            sepwidth=c(0.01, 0.004),  # width of the borders
            sepcolor='black',    
            colRow = cols,
            colsep = c(0:ncol(mean_all_p_model_given_y)),
            rowsep = c(0:nrow(mean_all_p_model_given_y)),
            col=my_palette,  # Set the custom color palette
            key=TRUE,       # Show the color key
            keysize = 1.0,  # Adjust the size of the color key
            cexRow = 1.0,   # Adjust the size of row labels
            cexCol = 1.0,   # Adjust the size of column labels
            scale="none")   # Turn off scaling
  
   # Reset the margins to default after plotting
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  if (bool_save) {
    dev.off()
  }
  
  
  for (ll in 1:L_sim) {
    if (bool_save) {
      pdf(paste0(res_dir, "/heatmap_L_sim_", L_sim, "_RMSE_FPC_", ll, ".pdf"),
          width = 8, height = 7, paper='special')
    }
    # Adjust the size of the heatmap
    par(mar = c(5, 4, 4, 8)) # Adjust the margins as needed
    
    # Create the heatmap with the modified color palette
    heatmap.2(mean_list_rmse[[ll]],
              dendrogram='none',
              Rowv=FALSE,
              Colv=FALSE,
              trace='none',
              main = paste0("RMSE FPC ", ll),
              sepwidth=c(0.01, 0.004),  # width of the borders
              sepcolor='black',    
              colRow = cols,
              na.color = grDevices::adjustcolor("red", alpha.f = 0.2),
              colsep = c(0:ncol(mean_list_rmse[[ll]])),
              rowsep = c(0:nrow(mean_list_rmse[[ll]])),
              col=my_palette, # Set the custom color palette
              key=TRUE,       # Show the color key
              keysize = 1.0,  # Adjust the size of the color key
              cexRow = 1.0,   # Adjust the size of row labels
              cexCol = 1.0,   # Adjust the size of column labels
              scale="none")   # Turn off scaling
    
    # Reset the margins to default after plotting
    par(mar = c(5, 4, 4, 2) + 0.1)
    
    if (bool_save) {
      dev.off()
    }
  }
  
  if (bool_save) {
    pdf(paste0(res_dir, "/heatmap_L_sim_", L_sim, "_ISE_mu.pdf"),
        width = 8, height = 7, paper='special')
  }
  # Adjust the size of the heatmap
  par(mar = c(5, 4, 4, 8)) # Adjust the margins as needed
  
  # Create the heatmap with the modified color palette
  heatmap.2(mean_list_ise$mu,
            dendrogram='none',
            Rowv=FALSE,
            Colv=FALSE,
            trace='none',
            main = paste0("ISE mu(t)"),
            sepwidth=c(0.01, 0.004),  # width of the borders
            sepcolor='black',    
            colRow = cols,
            na.color = grDevices::adjustcolor("red", alpha.f = 0.2),
            colsep = c(0:ncol(mean_list_ise$mu)),
            rowsep = c(0:nrow(mean_list_ise$mu)),
            col=my_palette,  # Set the custom color palette
            key=TRUE,       # Show the color key
            keysize = 1.0,  # Adjust the size of the color key
            cexRow = 1.0,   # Adjust the size of row labels
            cexCol = 1.0,   # Adjust the size of column labels
            scale="none")   # Turn off scaling
  
  # Reset the margins to default after plotting
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  if (bool_save) {
    dev.off()
  }
  
  for (ll in 2:(L_sim+1)) {
    if (bool_save) {
      pdf(paste0(res_dir, "/heatmap_L_sim_", L_sim, "_ISE_Psi_", ll-1, ".pdf"),
          width = 8, height = 7, paper='special')
    }
    # Adjust the size of the heatmap
    par(mar = c(5, 4, 4, 8)) # Adjust the margins as needed
    
    # Create the heatmap with the modified color palette
    heatmap.2(mean_list_ise[[ll]],
              dendrogram='none',
              Rowv=FALSE,
              Colv=FALSE,
              trace='none',
              main = paste0("ISE Psi_", ll -1, "(t)"),
              sepwidth=c(0.01, 0.004),  # width of the borders
              sepcolor='black',    
              colRow = cols,
              na.color = grDevices::adjustcolor("red", alpha.f = 0.2),
              colsep = c(0:ncol(mean_list_ise[[ll]])),
              rowsep = c(0:nrow(mean_list_ise[[ll]])),
              col=my_palette,  # Set the custom color palette
              key=TRUE,       # Show the color key
              keysize = 1.0,  # Adjust the size of the color key
              cexRow = 1.0,   # Adjust the size of row labels
              cexCol = 1.0,   # Adjust the size of column labels
              scale="none")   # Turn off scaling
    
    # Reset the margins to default after plotting
    par(mar = c(5, 4, 4, 2) + 0.1)
    
    if (bool_save) {
      dev.off()
    }
  }
  
  
  # recover K corresponding to Ruppert's rule of thumb
  K_pve <- median(sapply(1:p, function(j) max(round(min(median(n[,j]/4), 40)), 7)))
  list_rmse <- vector("list", L_sim)
  list_ise <- vector("list", (L_sim)+1)
  names(list_rmse) <- paste0("FPC_", 1:L_sim)
  names(list_ise) <- c("mu", paste0("Psi_", 1:L_sim))
  
  for (ind in 1:L_sim) { # K_model_choice and L_model_choice depends on the replicate.
    list_rmse[[ind]] <- matrix(NA, nrow = n_repl, ncol = 5)
    list_rmse[[ind]][,1] <- sapply(seq_along(list_all_rmse), function(repl) list_all_rmse[[repl]][[paste0("L_", vec_L_model_choice[repl])]][[paste0("K_", vec_K_model_choice[repl])]][ind]) #  L_model_choice, K_model_choice
    list_rmse[[ind]][,2] <- sapply(seq_along(list_all_rmse), function(repl) list_all_rmse[[repl]][[paste0("L_", L_max)]][[paste0("K_", vec_K_model_choice_L_10[repl])]][ind]) #  L_10, K_model_choice
    list_rmse[[ind]][,3] <- sapply(seq_along(list_all_rmse), function(repl) list_all_rmse[[repl]][[paste0("L_", L_max)]][[paste0("K_", K_pve)]][ind]) #  L_10, K_Ruppert # but this is taking the median across all variables, so not exactly Ruppert rule
    list_rmse[[ind]][,4] <- sapply(seq_along(list_all_rmse), function(repl) list_all_rmse[[repl]][[paste0("L_", L_max)]][[paste0("K_", sample(vec_K, 1))]][ind]) # pick one K at random, for L = 10
    list_rmse[[ind]][,5] <- sapply(seq_along(list_all_rmse), function(repl) list_all_rmse[[repl]][[paste0("L_", sample(setdiff(vec_L, 1:(ind-1)), 1))]][[paste0("K_", sample(vec_K, 1))]][ind]) # pick one K and one L at random, setdiff(vec_L, 1:(ind-1)) is to exclude the case with ind (simulated L) > L and so we can't compute the error for L = ind
    
    rownames(list_rmse[[ind]]) <- paste0("repl_", 1:n_repl)
    colnames(list_rmse[[ind]]) <- c("MC for K & L", "MC for K, L = 10", "ART for K, L = 10", "any K, L = 10", "any K & any L") # ART = approx rule of thumb Ruppoert as taking the median across variables to avoid re-running the algorithms
  
    if (bool_save) {
      pdf(paste0(res_dir, "/boxplots_L_sim_", L_sim, "_RMSE_FPC_", ind, ".pdf"),
          width = 10, height = 5, paper='special')
    }
    
    print(boxplot(list_rmse[[ind]], main =paste0("RMSE FPC ", ind), las = 1, cex.axis=.8))

    if (bool_save) {
      dev.off()
    }
  }

  
  for (ind in 1:(L_sim+1)) { # K_model_choice and L_model_choice depends on the replicate.
    list_ise[[ind]] <- matrix(NA, nrow = n_repl, ncol = 5)
    list_ise[[ind]][,1] <- sapply(seq_along(list_all_ise), function(repl) mean(list_all_ise[[repl]][[paste0("L_", vec_L_model_choice[repl])]][[paste0("K_", vec_K_model_choice[repl])]][ind,])) #  L_model_choice, K_model_choice
    list_ise[[ind]][,2] <- sapply(seq_along(list_all_ise), function(repl) mean(list_all_ise[[repl]][[paste0("L_", L_max)]][[paste0("K_", vec_K_model_choice_L_10[repl])]][ind,])) #  L_10, K_model_choice
    list_ise[[ind]][,3] <- sapply(seq_along(list_all_ise), function(repl) mean(list_all_ise[[repl]][[paste0("L_", L_max)]][[paste0("K_", K_pve)]][ind,])) #  L_10, K_Ruppert # but this is taking the median across all variables, so not exactly Ruppert rule
    list_ise[[ind]][,4] <- sapply(seq_along(list_all_ise), function(repl) mean(list_all_ise[[repl]][[paste0("L_", L_max)]][[paste0("K_", sample(vec_K, 1))]][ind,])) # pick one K at random, for L = 10
    list_ise[[ind]][,5] <- sapply(seq_along(list_all_ise), function(repl) mean(list_all_ise[[repl]][[paste0("L_", sample(setdiff(vec_L, 1:(ind-1)), 1))]][[paste0("K_", sample(vec_K, 1))]][ind,])) # pick one K and one L at random, setdiff(vec_L, 1:(ind-1)) is to exclude the case with ind (simulated L) > L and so we can't compute the error for L = ind
    
    rownames(list_ise[[ind]]) <- paste0("repl_", 1:n_repl)
    colnames(list_ise[[ind]]) <- c("MC for K & L", "MC for K, L = 10", "ART for K, L = 10", "any K, L = 10", "any K & any L")  # ART = approx rule of thumb Ruppoert as taking the median across variables to avoid re-running the algorithms
 
    if (bool_save) {
      pdf(paste0(res_dir, "/boxplots_L_sim_", L_sim, "_ISE_",  names(list_ise)[ind], ".pdf"),
          width = 10, height = 5, paper='special')
    }
    
    print(boxplot(list_ise[[ind]], main =paste0("ISE ", names(list_ise)[ind]), las = 1, cex.axis=.8))
    
    if (bool_save) {
      dev.off()
    }  
  }

} else {
  
  if (bool_save) {
    pdf(paste0(res_dir, "/boxplots_prob_L_sim_", L_sim, ".pdf"),
             width = 10, height = 5, paper='special')
  }
  
  boxplot(mat_all_p_model_given_y, main ="Model choice probabilities")
  abline(v = which(colnames(mat_all_p_model_given_y) == paste0("L_", L_sim)), col = "red")
  
  if (bool_save) {
    dev.off()
  }
  
}

if (bool_save) {
  pdf(paste0(res_dir, "/boxplots_cum_pve_using_L_10_for_L_sim_", L_sim, ".pdf"),
      width = 7, height = 5, paper='special')
}
colnames(mat_all_cum_pve) <- gsub("FPC_", "L_", colnames(mat_all_cum_pve))
boxplot(mat_all_cum_pve, main = paste0("Cumulated proportion of variance explained", add_label), ylim = c(0,100))
abline(h = 95, col = "gray60", lty = 2)
abline(v = which(colnames(mat_all_cum_pve) == paste0("L_", L_sim)), col = "red")

if (bool_save) {
  dev.off()
}

if (bool_save) {
  save(list_all_rmse, list_all_ise, 
       list_rmse, list_ise,
       mean_list_rmse, mean_list_ise,
       sd_list_rmse, sd_list_ise,
       L_sim, mean_all_p_model_given_y, sd_all_p_model_given_y, 
       list_all_p_model_given_y, list_all_cum_pve, list_all_elbos, 
       vec_K_model_choice_L_10, vec_K_model_choice, vec_L_model_choice,
       list_all_cum_pve_model_choice_K, 
       n_repl, choose_both, unif_dist_K, 
       choose_K_pve, bool_orth_splines, 
       p, N, N_t_min, N_t_max, L_max, lambda_L, 
       K, K_min, K_max, lambda_K, 
       exponent_sd_zeta, vec_sd_eps, tol, maxit, seed, n,
       mat_all_cum_pve, mean_all_p_model_given_y, n_g,                          
       file = paste0(res_dir, "/output_L_sim_", L_sim, ".RData"))
}

