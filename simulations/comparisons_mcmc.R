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
  
  job_id <- 3 # 1 to 6
  
  CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")
  out_dir <- file.path(CORE_DIR_ICLOUD, "mFPCA_output/")
}

main_dir <- file.path(CORE_DIR, "bayesian-mFPCA-paper-code/simulations/")
setwd(main_dir)

source("fun_utils.R") 

require(scales)
require(bayesFPCA)
require(rstan)
rstan_options(auto_write = TRUE)
require(ggplot2)


seed <- 123
set.seed(seed)

p <- 3                                        # number of variables
N <- c(50, 100, 200, 300, 400, 500)[job_id]   # number of curves
N_t_min <- 10                                 # minimum number of time observations for each curve
N_t_max <- 20                                 # maximum number of time observations for each curve
T_vec <- sample(N_t_min:N_t_max,              # number of time observations for each curve
                N, replace = T)               
n_time_obs <- sum(T_vec)                      # total number of observations (i.e., for all multivariate curves)
n <- matrix(rep(T_vec, p), ncol = p)          # matrix whose entries are the number of observations per subject and variable

K_Ruppert <- T
if (K_Ruppert) {
  K <- NULL
} else {
  n_int_knots <- 5                            # number of interior knots
  K <- n_int_knots + 2                        # number of spline basis functions
}

L_sim <- 2                                    # guess of the number of FPCA basis functions
stopifnot(L_sim == 2)                         # if different L_sim wanted, then need to modify the function 
                                              # summarise_mcmc as outputs the first two sets of scores
L <- 10
L_mcmc <- 10
data_col <- "grey50"                          # colour of the data in the plots


mfvb_col <- "darkgreen"                       # colour of the VB lines in the plots
vmp_col <- "black"                            


tol <- 1e-5
maxit <- 1000

n_burnin <- 1000                              # length of burn-in.
n_mcmc <- 1000                                # number of mcmc iterations
n_thin <- 1                                   # thinning factor.
mcmc_col <- "blue"                            # colour of the MCMC lines in the plots

sigma_eps <- rep(1, p)                        # sd of the residuals

n_g <- 1000                                   # length of plotting grid

bool_multiple_repl <- T                       # whether or not to perform multiple VB simulations to estimate the mean and eigen-functions
n_repl <- 100                                 # number of replicates when bool_multiple_repl = T
n_cpus <- 10                                  # number of cpus used for parallel execution when bool_multiple_repl = T

bool_mcmc_repl <- T

bool_save <- T
if (bool_save) {
  res_dir <- paste0(out_dir, "/comparison_mcmc_corr_metrics", 
                    ifelse(K_Ruppert, "_K_Ruppert", ""),
                    ifelse(bool_mcmc_repl, "", "_no_repl_mcmc"), "_p_", p, 
                    "_N_", N, "_Nt_min_", N_t_min, "_max_", N_t_max, "_K_", K, 
                    "_L_sim_", L_sim, "_L_", L, "_L_mcmc_", L_mcmc, "_sigeps_", 
                    unique(sigma_eps), "_nitermcmc_", n_mcmc, "_tol_", tol, 
                    "_maxit_", maxit, "_n_repl_", n_repl, "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}


mu_beta <- 0
sigsq_beta <- 1e10
sigma_beta <- sqrt(sigsq_beta)
A <- 1e5

sigma_zeta <- 1/(1:L_sim)       # for simulating the data only.
sigma_zeta_init <- 1/(1:L_mcmc) # for the initialisation of the mcmc algorithm


f_mu <- function(t, j) (-1)^j*2*sin((2*pi+j)*t)

f_psi_1 <- function(t, j, p) (-1)^j*sqrt(2/p)*cos(2*pi*t)
f_psi_2 <- function(t, j, p) (-1)^j*sqrt(2/p)*sin(2*pi*t)


f_Psi <- function(time_obs, j, p) {
  ans <- cbind(f_psi_1(time_obs, j, p),
               f_psi_2(time_obs, j, p))
  return(ans)
}

data <- generate_fpca_data(N, p, n, L = L_sim, n_g, sigma_eps,
                           f_mu, f_Psi,
                           vec_sd_zeta = sigma_zeta,
                           generate_from_univ = F)

time_obs <- data$time_obs
zeta <- data$Zeta
Y <- data$Y

time_g <- data$time_g
mu_g <- data$mu_g
Psi_g <- data$Psi_g

if (K_Ruppert) {
  K <- unique(sapply(1:p, function(j) max(round(min(median(sapply(time_obs, function(time_obs_i) length(time_obs_i[[j]]))/4), 40)), 7)))
}
list_grid <- get_grid_objects(time_obs, K, time_g = time_g)

C <- list_grid$C

X <- lapply(C, function(C_i) lapply(C_i, function(C_i_j) C_i_j[, 1:2]))
Z <- lapply(C, function(C_i) lapply(C_i, function(C_i_j) C_i_j[, -c(1:2)]))

C_g <- list_grid$C_g
C_g_same_grid <- unique(C_g)[[1]]

list_X <- lapply(1:p, function(j) Reduce(rbind, sapply(1:N, function(i) X[[i]][[j]])))
X_array <- array(NA_real_, dim = c(p, n_time_obs, 2))
for(j in 1:p) X_array[j,,] <- list_X[[j]]

list_Z <- lapply(1:p, function(j) Reduce(rbind, sapply(1:N, function(i) Z[[i]][[j]])))
Z_array <- array(NA_real_, dim = c(p, n_time_obs, K))
for(j in 1:p) Z_array[j,,] <- list_Z[[j]]

Y_mat <- t(sapply(1:p, function(j) Reduce(c, sapply(1:N, function(i) Y[[i]][[j]]))))

mfpca_model <- "

	data {

		int<lower=1> p;                // number of variables
		int<lower=1> N;                // number of curves
		int<lower=N> n_time_obs;       // total number of time observations # HR: for now, same number of obs for all p variables
		int<lower=1> K;                // number of splines
		int<lower=1> L_mcmc;           // number of basis functions
		real<lower=0> sigma_beta;      // fixed effects prior standard deviation
		real<lower=0> A;               // cauchy hyperparameter
    vector<lower=0>[L_mcmc] sigma_zeta; // prior standard deviation of the scores
		matrix[n_time_obs, 2] X[p];    // array [p, n_time_obs, 2], X[j] is an rbind of all design matrices for variable j
		matrix[n_time_obs, K] Z[p];    // array [p, n_time_obs, K], Z[j] is an rbind of all spline design matrices for variable j
		int<lower=1> T_vec[N];         // vector of time observations for each curve
		vector[n_time_obs] Y[p];       // multivariate case: matrix of dim p x n_time_obs, where the jth row are the obs for all subjects (stacked) for variable j
	}

	parameters {

		matrix[N,L_mcmc] zeta;

		vector<lower=0>[p] sigma_eps;

		vector[2] beta_mu[p];
		vector[K] u_mu[p];
		vector<lower=0>[p] sigma_mu;

		matrix[L_mcmc, 2] beta_psi[p];
		matrix[L_mcmc, K] u_psi[p];
		vector<lower=0>[L_mcmc] sigma_psi[p];
	}

	transformed parameters {

    vector[n_time_obs] mu[p];
	  matrix[L_mcmc, n_time_obs] psi[p];

		for (j in 1:p) {

			mu[j] = X[j]*beta_mu[j] + Z[j]*u_mu[j];

  		for (l in 1:L_mcmc) {
  			psi[j, l] = beta_psi[j, l]*X[j]' + u_psi[j, l]*Z[j]';
  		}

  	}

	}

	model {

		int pos;

		for (j in 1:p) {

			pos = 1;

  		for(i in 1:N) {

  			// Temporary vectors
  			vector[T_vec[i]] mu_i;
  			matrix[L_mcmc, T_vec[i]] psi_i;
  			vector[T_vec[i]] Y_i_hat;

  			mu_i = segment(mu[j], pos, T_vec[i]);
  			psi_i = block(psi[j], 1, pos, L_mcmc, T_vec[i]);
  			Y_i_hat = mu_i + to_vector(zeta[i]*psi_i);

  			segment(Y[j], pos, T_vec[i]) ~ normal(Y_i_hat, sigma_eps[j]);

  			pos = pos + T_vec[i];

  			for (l in 1:L_mcmc) {
  			  zeta[i,l] ~ normal(0, sigma_zeta[l]);
  			}
  		}


  		sigma_eps[j] ~ cauchy(0, A);

  		beta_mu[j] ~ normal(0, sigma_beta);
  		u_mu[j] ~ normal(0, sigma_mu[j]);
  		sigma_mu[j] ~ cauchy(0, A);

  		for(l in 1:L_mcmc) {

  			beta_psi[j, l] ~ normal(0, sigma_beta);
  			u_psi[j, l] ~ normal(0, sigma_psi[j, l]);
  			sigma_psi[j, l] ~ cauchy(0, A);

  		}
		}
	}
"

all_data <- list(p=p,
                 N=N, n_time_obs=n_time_obs, K=K, L_mcmc=L_mcmc,
                 sigma_beta=sigma_beta, A=A,
                 sigma_zeta=sigma_zeta_init,
                 X=X_array,
                 Z=Z_array,
                 T_vec=T_vec,
                 Y=Y_mat
)

time_mcmc <- system.time({
  
  compile_obj <- stan(
    model_code=mfpca_model, data=all_data,
    iter=1, chains=1
  )
  
  stan_obj <- stan(
    model_code=mfpca_model, data=all_data, warmup=n_burnin,
    iter=(n_burnin+n_mcmc), chains=1, thin=n_thin,
    refresh=100, fit=compile_obj
  )
  
  mcmc_summary <- summarise_mcmc_multivariate(stan_obj, C_g_same_grid, Psi_g, 
                                              L_sim = L_sim, 
                                              pred_interval = T)
})

print(time_mcmc)

Y_mcmc_summary <- mcmc_summary$Y_g_mcmc_summary
gbl_mcmc_summary <- mcmc_summary$gbl_mcmc_summary
gbl_ci_mcmc_summary <- mcmc_summary$gbl_ci_mcmc_summary
zeta_mcmc_summary <- mcmc_summary$zeta_mcmc_summary

Y_mcmc_low <- lapply(Y_mcmc_summary, function(Y_i) lapply(Y_i, function(Y_i_j) Y_i_j[,1]))
Y_mcmc_hat <- lapply(Y_mcmc_summary, function(Y_i) lapply(Y_i, function(Y_i_j) Y_i_j[,2]))
Y_mcmc_upp <- lapply(Y_mcmc_summary, function(Y_i) lapply(Y_i, function(Y_i_j) Y_i_j[,3]))



set.seed(seed)

list_hyper <- set_hyper(sigma_zeta = 1, 
                        sigma_beta = sigma_beta, A = A)

vmp_res <- run_vmp_fpca(time_obs, Y, L = L, K = K, n_g = NULL, time_g = time_g, 
                        tol = tol, maxit = maxit,
                        plot_elbo = F, Psi_g = Psi_g,
                        list_hyper = list_hyper)

Y_hat <- vmp_res$Y_hat
Y_low <- vmp_res$Y_low
Y_upp <- vmp_res$Y_upp

mu_hat <- vmp_res$mu_hat
list_Psi_hat <- vmp_res$list_Psi_hat
Zeta_hat <- vmp_res$Zeta_hat
list_zeta_ellipse <- vmp_res$list_zeta_ellipse


bool_mfvb <- T
if (bool_mfvb) {
  
  set.seed(seed)
  
  mfvb_res <- run_mfvb_fpca(time_obs, Y, L = L, K = K, tol = tol, maxit = maxit, 
                            n_g = NULL, time_g = time_g, Psi_g = Psi_g, 
                            list_hyper = list_hyper)
  
  
  Y_hat_mfvb <- mfvb_res$Y_hat
  Y_low_mfvb <- mfvb_res$Y_low
  Y_upp_mfvb <- mfvb_res$Y_upp
  
  mu_hat_mfvb <- mfvb_res$mu_hat
  list_Psi_hat_mfvb <- mfvb_res$list_Psi_hat
  Zeta_hat_mfvb <- mfvb_res$Zeta_hat
  list_zeta_ellipse_mfvb <- mfvb_res$list_zeta_ellipse
  
} else {
  
  mfvb_res <- NULL
  
}


n_plots_fit <- 10   # number of plots produced (each with a random set of subjects)
n_sample_fit <- 3   # number of subjects shown
p_sample <- 1:p     # variables to display

for (ff in 1:n_plots_fit) {
  
  set.seed(ff)
  N_sample_ff <- sort(sample(1:N, n_sample_fit))
  if (bool_save) {
    pdf(paste0(res_dir, "/data_samples_", paste0(N_sample_ff, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit_list(p_sample, N_sample_ff, time_obs, time_g, Y, Y_hat = NULL,
                   Y_low = NULL, Y_upp = NULL)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/fits_samples_", paste0(N_sample_ff, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit_list(p_sample, N_sample_ff, time_obs, time_g,
                   Y, Y_hat, Y_low, Y_upp, # vmp
                   Y_hat_add = Y_mcmc_hat, # mfvb
                   Y_low_add = Y_mcmc_low,
                   Y_upp_add = Y_mcmc_upp)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/fits_samples_mfvb_", paste0(N_sample_ff, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit_list(p_sample, N_sample_ff, time_obs, time_g,
                   Y, Y_hat_mfvb, Y_low_mfvb, Y_upp_mfvb, # mfvb
                   Y_hat_add = Y_mcmc_hat, 
                   Y_low_add = Y_mcmc_low,
                   Y_upp_add = Y_mcmc_upp)
  if (bool_save) {
    dev.off()
  }
  
}


if (bool_save) {
  pdf(paste0(res_dir, "/eigenfunctions.pdf"),
      width = 6.3, height = 5, paper='special')
}

display_eigenfunctions(L_sim, time_g, mu_g, Psi_g,
                       mu_hat, list_Psi_hat, # vmp = grey
                       mu_hat_add = gbl_mcmc_summary[[1]], # mcmc = blue
                       list_Psi_hat_add =gbl_mcmc_summary[2:3],
                       vec_col_add = c(vmp_col, mcmc_col),
                       vec_lwd = c(2,1))
if (bool_save) {
  dev.off()
}

n_plots_scores <- 4
n_sample_ps <- 16

for (ps in 1:n_plots_scores) {
  
  set.seed(ps)
  N_sample_ps <- sort(sample(1:N, n_sample_ps))
  
  if (bool_save) { 
    pdf(paste0(res_dir, "/scores_samples_", paste0(N_sample_ps, collapse = "-"), ".pdf"),
        width = 6, height = 6, paper='special')
  }
  display_scores(N_sample_ps, zeta, Zeta_hat, list_zeta_ellipse,
                 Zeta_hat_add =  t(sapply(zeta_mcmc_summary, "[[", "mean")),
                 zeta_ellipse_add =  lapply(zeta_mcmc_summary, "[[", "credible boundary"),
                 mfrow = c(ceiling(sqrt(n_sample_ps)), tail(sqrt(n_sample_ps))))
  if (bool_save) {
    dev.off()
  }
  
  
}



rmse_mcmc <- apply(zeta - t(sapply(zeta_mcmc_summary, "[[", "mean")), 2, function(x) sqrt(mean(x^2)))
rmse_vmp <- apply(zeta - Zeta_hat[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
rmse_mfvb <- apply(zeta - Zeta_hat_mfvb[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
names(rmse_mcmc) <- names(rmse_vmp) <- names(rmse_mfvb) <- paste0("FPC_", 1:L_sim)

ise_mcmc <- ise_vmp <- ise_mfvb <- matrix(NA, L_sim+1, p)

for (l in 1:(L_sim+1)) {
  
  for(j in 1:p) {
    
    if (l == 1) {
      ise_mcmc[l, j] <- trapint(time_g, (gbl_mcmc_summary[[1]][,j] - mu_g[[j]])^2)
      ise_vmp[l, j] <- trapint(time_g, (mu_hat[,j] - mu_g[[j]])^2)
      ise_mfvb[l, j] <- trapint(time_g, (mu_hat_mfvb[,j] - mu_g[[j]])^2)
    } else {
      
      ise_mcmc[l, j] <- trapint(time_g, (gbl_mcmc_summary[[l]][,j]- Psi_g[[j]][,l-1])^2)
      ise_vmp[l, j] <- trapint(time_g, (list_Psi_hat[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
      ise_mfvb[l, j] <- trapint(time_g, (list_Psi_hat_mfvb[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
    }
    
  }
  
}
rownames(ise_mcmc) <- rownames(ise_vmp) <- rownames(ise_mfvb) <- c("mu", paste0("Psi_", 1:L_sim))
colnames(ise_mcmc) <- colnames(ise_vmp) <- colnames(ise_mfvb) <- paste0("Variable_", 1:p)




if (bool_multiple_repl) {
  
  list_res <- parallel::mclapply(1:n_repl, function(repl) {
    
    set.seed(repl)
    
    Y_repl <- vector("list", length = N)
    
    zeta_repl <- matrix(NA, N, L_sim)
    for(i in 1:N) {
      
      zeta_repl[i,] <- MASS::mvrnorm(1, rep(0, L_sim), diag(sigma_zeta^2))
      
      Y_repl[[i]] <- vector("list", length = p)
      for(j in 1:p) {
        
        resid_vec <- rnorm(n[i, j], 0, sigma_eps[j])
        mean_vec <- f_mu(time_obs[[i]][[j]], j) + f_Psi(time_obs[[i]][[j]], j, p = p) %*% zeta_repl[i, ]
        Y_repl[[i]][[j]] <- as.vector(mean_vec + resid_vec)
      }
    }
    
    
    time_vmp_repl <- system.time({vmp_res_repl <- run_vmp_fpca(time_obs, Y_repl, L = L, K = K, 
                                                               n_g = NULL, time_g = time_g, 
                                                               tol = tol, maxit = maxit,
                                                               plot_elbo = FALSE, Psi_g = Psi_g,
                                                               list_hyper = list_hyper)})["elapsed"]
    
    Y_hat <- vmp_res_repl$Y_hat
    Y_low <- vmp_res_repl$Y_low
    Y_upp <- vmp_res_repl$Y_upp
    
    mu_hat <- vmp_res_repl$mu_hat
    list_Psi_hat <- vmp_res_repl$list_Psi_hat
    Zeta_hat <- vmp_res_repl$Zeta_hat
    list_zeta_ellipse <- vmp_res_repl$list_zeta_ellipse
    
    
    # do not compare time because a fixed number of iterations.
    time_mfvb_repl <- system.time({mfvb_res_repl <- run_mfvb_fpca(time_obs, Y_repl, L = L, K = K, 
                                                                  n_g = NULL, time_g = time_g, 
                                                                  tol = tol, maxit = maxit, 
                                                                  Psi_g = Psi_g, list_hyper = list_hyper)})["elapsed"]
    
    Y_hat_mfvb <- mfvb_res_repl$Y_hat
    Y_low_mfvb <- mfvb_res_repl$Y_low
    Y_upp_mfvb <- mfvb_res_repl$Y_upp
    
    mu_hat_mfvb <- mfvb_res_repl$mu_hat
    list_Psi_hat_mfvb <- mfvb_res_repl$list_Psi_hat
    Zeta_hat_mfvb <- mfvb_res_repl$Zeta_hat
    list_zeta_ellipse_mfvb <- mfvb_res_repl$list_zeta_ellipse
    
    
    if (bool_mcmc_repl) {

      Y_repl_mat <- t(sapply(1:p, function(j) Reduce(c, sapply(1:N, function(i) Y_repl[[i]][[j]]))))
      
      
      all_data_repl <- list(p=p,
                            N=N, n_time_obs=n_time_obs, K=K, L_mcmc=L_mcmc,
                            sigma_beta=sigma_beta, A=A,
                            sigma_zeta=sigma_zeta_init,
                            X=X_array,
                            Z=Z_array,
                            T_vec=T_vec,
                            Y=Y_repl_mat
      )
      
      time_mcmc_repl <- system.time({
        compile_obj_repl <- stan(
          model_code=mfpca_model, data=all_data_repl,
          iter=1, chains=1
        )
        
        stan_obj_repl <- stan(
          model_code=mfpca_model, data=all_data_repl, warmup=n_burnin,
          iter=(n_burnin+n_mcmc), chains=1, thin=n_thin,
          refresh=100, fit=compile_obj_repl
        )
        
        mcmc_summary_repl <- summarise_mcmc_multivariate(stan_obj_repl, C_g_same_grid, Psi_g, L_sim = L_sim, pred_interval = T)
      })["elapsed"]
      
      
      gbl_mcmc_summary <- mcmc_summary_repl$gbl_mcmc_summary
      mu_hat_mcmc <- gbl_mcmc_summary[[1]]
      list_Psi_hat_mcmc <- gbl_mcmc_summary[-1]
      zeta_mcmc_summary <- mcmc_summary_repl$zeta_mcmc_summary 
      
      rmse_mcmc <- apply(zeta_repl - t(sapply(zeta_mcmc_summary, "[[", "mean")), 2, function(x) sqrt(mean(x^2)))
      names(rmse_mcmc) <- paste0("FPC_", 1:L_sim)
      
      ise_mcmc <- matrix(NA, L_sim+1, p)
      
      for (l in 1:(L_sim+1)) {
        
        for(j in 1:p) {
          
          if (l == 1) {
            ise_mcmc[l, j] <- trapint(time_g, (mu_hat_mcmc[,j] - mu_g[[j]])^2)
          } else {
            ise_mcmc[l, j] <- trapint(time_g, (list_Psi_hat_mcmc[[l-1]][,j]- Psi_g[[j]][,l-1])^2)
          }
          
        }
        
      }
      rownames(ise_mcmc) <- c("mu", paste0("Psi_", 1:L_sim))
      colnames(ise_mcmc) <- paste0("Variable_", 1:p)
      
    } else{
      mu_hat_mcmc <- list_Psi_hat_mcmc <- rmse_mcmc <- ise_mcmc <- time_mcmc_repl <- NULL
    }
    
    
    rmse_vmp <- apply(zeta_repl - Zeta_hat[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
    rmse_mfvb <- apply(zeta_repl - Zeta_hat_mfvb[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
    names(rmse_vmp) <- names(rmse_mfvb) <- paste0("FPC_", 1:L_sim)
    
    ise_vmp <- ise_mfvb <- matrix(NA, L_sim+1, p)
    
    for (l in 1:(L_sim+1)) {
      
      for(j in 1:p) {
        
        if (l == 1) {
          ise_vmp[l, j] <- trapint(time_g, (mu_hat[,j] - mu_g[[j]])^2)
          ise_mfvb[l, j] <- trapint(time_g, (mu_hat_mfvb[,j] - mu_g[[j]])^2)
        } else {
          ise_vmp[l, j] <- trapint(time_g, (list_Psi_hat[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
          ise_mfvb[l, j] <- trapint(time_g, (list_Psi_hat_mfvb[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
        }
        
      }
      
    }
    rownames(ise_vmp) <- rownames(ise_mfvb) <- c("mu", paste0("Psi_", 1:L_sim))
    colnames(ise_vmp) <- colnames(ise_mfvb) <- paste0("Variable_", 1:p)
    
    
    create_named_list(mu_hat, list_Psi_hat,
                      mu_hat_mfvb, list_Psi_hat_mfvb,
                      mu_hat_mcmc, list_Psi_hat_mcmc,
                      rmse_vmp, rmse_mfvb, rmse_mcmc,
                      ise_vmp, ise_mfvb, ise_mcmc,
                      time_vmp_repl, time_mfvb_repl, time_mcmc_repl)
  }, mc.cores = n_cpus)
  
  list_mu_hat <- lapply(list_res, "[[", "mu_hat")
  list_list_Psi_hat  <- lapply(list_res, "[[", "list_Psi_hat")
  
  list_mu_hat_mfvb <- lapply(list_res, "[[", "mu_hat_mfvb")
  list_list_Psi_hat_mfvb  <- lapply(list_res, "[[", "list_Psi_hat_mfvb")
  
  if (bool_mcmc_repl) {
    list_mu_hat_mcmc <- lapply(list_res, "[[", "mu_hat_mcmc")
    list_list_Psi_hat_mcmc  <- lapply(list_res, "[[", "list_Psi_hat_mcmc")
  }
  
  mat_rmse_vmp <- sapply(list_res, "[[", "rmse_vmp")
  list_ise_vmp  <- lapply(list_res, "[[", "ise_vmp")
  
  mat_rmse_mfvb <- sapply(list_res, "[[", "rmse_mfvb")
  list_ise_mfvb  <- lapply(list_res, "[[", "ise_mfvb")
  
  colnames(mat_rmse_vmp) <- colnames(mat_rmse_mfvb) <- names(list_ise_vmp) <- names(list_ise_mfvb) <- paste0("repl_", 1:n_repl)
  
  if (bool_mcmc_repl) {
    mat_rmse_mcmc <- sapply(list_res, "[[", "rmse_mcmc")
    list_ise_mcmc  <- lapply(list_res, "[[", "ise_mcmc")
    colnames(mat_rmse_mcmc) <- names(list_ise_mcmc) <- paste0("repl_", 1:n_repl)
  }
  
  avg_ise_vmp <- lapply(list_ise_vmp, function(mm) rowSums(mm))  
  avg_ise_mfvb <- lapply(list_ise_mfvb, function(mm) rowSums(mm))  
  
  meths <- c("VMP", "MFVB")
  ise_mu_repl <- cbind(sapply(avg_ise_vmp, "[[", "mu"), sapply(avg_ise_mfvb, "[[", "mu"))
  ise_Psi_1_repl <- cbind(sapply(avg_ise_vmp, "[[", "Psi_1"), sapply(avg_ise_mfvb, "[[", "Psi_1"))
  ise_Psi_2_repl <- cbind(sapply(avg_ise_vmp, "[[", "Psi_2"), sapply(avg_ise_mfvb, "[[", "Psi_2"))
  colnames(ise_mu_repl) <- colnames(ise_Psi_1_repl) <- colnames(ise_Psi_2_repl) <- meths
  
  rmse_FPC_1_repl <-  cbind(mat_rmse_vmp["FPC_1",], mat_rmse_mfvb["FPC_1",])
  rmse_FPC_2_repl <-  cbind(mat_rmse_vmp["FPC_2",], mat_rmse_mfvb["FPC_2",])
  colnames(rmse_FPC_1_repl) <- colnames(rmse_FPC_2_repl) <- meths
  
  if (bool_mcmc_repl) {
    avg_ise_mcmc <- lapply(list_ise_mcmc, function(mm) rowSums(mm))  
    ise_mu_repl <- cbind(ise_mu_repl, sapply(avg_ise_mcmc, "[[", "mu"))
    ise_Psi_1_repl <- cbind(ise_Psi_1_repl, sapply(avg_ise_mcmc, "[[", "Psi_1"))
    ise_Psi_2_repl <- cbind(ise_Psi_2_repl, sapply(avg_ise_mcmc, "[[", "Psi_2"))
    colnames(ise_mu_repl)[3] <- colnames(ise_Psi_1_repl)[3] <- colnames(ise_Psi_2_repl)[3] <- "MCMC"
    
    rmse_FPC_1_repl <-  cbind(rmse_FPC_1_repl, mat_rmse_mcmc["FPC_1",])
    rmse_FPC_2_repl <-  cbind(rmse_FPC_2_repl, mat_rmse_mcmc["FPC_2",])
    colnames(rmse_FPC_1_repl)[3] <- colnames(rmse_FPC_2_repl)[3] <- "MCMC"
    meths <- c(meths, "MCMC")
  }
  
  vec_time_vmp <- sapply(list_res, "[[", "time_vmp_repl")
  vec_mfvb_vmp <- sapply(list_res, "[[", "time_mfvb_repl")
  vec_mcmc_vmp <- sapply(list_res, "[[", "time_mcmc_repl")
  mat_time <- cbind(vec_time_vmp, vec_mfvb_vmp, vec_mcmc_vmp)
  colnames(mat_time) <- meths
  
  if (bool_save) {
    pdf(paste0(res_dir, "/runtime_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  df_runtime <- reshape2::melt(mat_time)
  df_runtime$Var2 <- factor(as.character(df_runtime$Var2), levels = meths)
  pl <- ggplot(df_runtime, aes(x = Var2,
                               y = value)) +
    geom_boxplot() + labs(x="",
                          y="Runtime") +
    ggtitle("Runtime")
  
  pl <- pl + theme_classic() 
  
  print(pl)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/log_runtime_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  df_runtime <- reshape2::melt(mat_time)
  df_runtime$Var2 <- factor(as.character(df_runtime$Var2), levels = meths)
  pl <- ggplot(df_runtime, aes(x = Var2,
                               y = value)) +
    scale_y_continuous(trans=log10_trans(), 
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() + labs(x="",
                          y="Runtime (log-scale)") +
    ggtitle("Runtime")
  
  pl <- pl + theme_classic() 
  
  print(pl)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/runtime_no_mfvb_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  df_runtime <- reshape2::melt(mat_time[,-2])
  df_runtime$Var2 <- factor(as.character(df_runtime$Var2), levels = meths)
  pl <- ggplot(df_runtime, aes(x = Var2,
                               y = value)) +
    geom_boxplot() + labs(x="",
                          y="Runtime") +
    ggtitle("Runtime")
  
  pl <- pl + theme_classic() 
  
  print(pl)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/log_runtime_no_mfvb_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  df_runtime <- reshape2::melt(mat_time[,-2])
  df_runtime$Var2 <- factor(as.character(df_runtime$Var2), levels = meths)
  pl <- ggplot(df_runtime, aes(x = Var2,
                               y = value)) +
    scale_y_continuous(trans=log10_trans(), 
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() + labs(x="",
                          y="Runtime (log-scale)") +
    ggtitle("Runtime")
  
  pl <- pl + theme_classic() 
  
  print(pl)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/ise_mu_avg_across_variables_vmp_mfvb_mcmc_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  vec_col <- c(vmp_col, mfvb_col)
  df_ise_mu <- reshape2::melt(ise_mu_repl)
  df_ise_mu$Var2 <- factor(as.character(df_ise_mu$Var2), levels = meths)
  pl <- ggplot(df_ise_mu, aes(x = Var2,
                              y = value)) +
    scale_y_continuous(trans=log10_trans(), 
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() + labs(x="",
                          y="ISE of estimated mu(t) (log-scale)", fill = "Estimation") +
    ggtitle("Integrated squared error of mean function")
  
  pl <- pl + theme_classic() 
  
  if (!bool_mcmc_repl) {
    indiv_ise_initial_run_mu <- data.frame(name=meths, value=c(mean(ise_vmp["mu",]), mean(ise_mfvb["mu",])))
    indiv_ise_initial_run_mu$name <- factor(indiv_ise_initial_run_mu$name, levels = meths)
    
    pl <- pl + geom_point(data= indiv_ise_initial_run_mu, aes(x=name, y=value, color=name), 
                          position=position_jitter(0), color=vec_col, size=5, pch=20)
    pl <- pl + geom_hline(yintercept=mean(ise_mcmc["mu",]), linetype="dashed", color = mcmc_col, size = 0.6) +
      theme(plot.title = element_text(size = 7, face = "bold"))
  }
  
  print(pl)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/ise_Psi_1_avg_across_variables_vmp_mfvb_mcmc_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  df_ise_Psi_1 <- reshape2::melt(ise_Psi_1_repl)
  df_ise_Psi_1$Var2 <- factor(as.character(df_ise_Psi_1$Var2), levels = meths)
  pl <- ggplot(df_ise_Psi_1, aes(x = Var2,
                                 y = value)) +
    scale_y_continuous(trans=log10_trans(), 
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() + labs(x="",
                          y="ISE of estimated Psi_1(t) (log-scale)", fill = "Estimation") +
    ggtitle("Integrated squared error of 1st eigenfunction")
  
  pl <- pl + theme_classic() 
  
  if (!bool_mcmc_repl) {
    indiv_ise_initial_run_Psi_1 <- data.frame(name=meths, value=c(mean(ise_vmp["Psi_1",]), mean(ise_mfvb["Psi_1",])))
    indiv_ise_initial_run_Psi_1$name <- factor(indiv_ise_initial_run_Psi_1$name, levels = meths)
    
    pl <- pl + geom_point(data= indiv_ise_initial_run_Psi_1, aes(x=name, y=value, color=name), 
                          position=position_jitter(0), color=vec_col, size=5, pch=20)
    pl <- pl + geom_hline(yintercept=mean(ise_mcmc["Psi_1",]), linetype="dashed", color = mcmc_col, size = 0.6) +
      theme(plot.title = element_text(size = 7, face = "bold"))
  }
  print(pl)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/ise_Psi_2_avg_across_variables_vmp_mfvb_mcmc_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  df_ise_Psi_2 <- reshape2::melt(ise_Psi_2_repl)
  df_ise_Psi_2$Var2 <- factor(as.character(df_ise_Psi_2$Var2), levels = meths)
  pl <- ggplot(df_ise_Psi_2, aes(x = Var2,
                                 y = value)) +
    scale_y_continuous(trans=log10_trans(), 
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    geom_boxplot() + labs(x="",
                          y="ISE of estimated Psi_2(t) (log-scale)", fill = "Estimation") +
    ggtitle("Integrated squared error of 2nd eigenfunction")
  
  pl <- pl + theme_classic() 
  
  if (!bool_mcmc_repl) {
    indiv_ise_initial_run_Psi_2 <- data.frame(name=meths, value=c(mean(ise_vmp["Psi_2",]), mean(ise_mfvb["Psi_2",])))
    indiv_ise_initial_run_Psi_2$name <- factor(indiv_ise_initial_run_Psi_2$name, levels = meths)
    
    pl <- pl + geom_point(data= indiv_ise_initial_run_Psi_2, aes(x=name, y=value, color=name), 
                          position=position_jitter(0), color=vec_col, size=5, pch=20)
    pl <- pl + geom_hline(yintercept=mean(ise_mcmc["Psi_2",]), linetype="dashed", color = mcmc_col, size = 0.6) +
      theme(plot.title = element_text(size = 7, face = "bold"))
  }
  print(pl)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/rmse_FPC_1_vmp_mfvb_mcmc_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  df_ise_FPC_1 <- reshape2::melt(rmse_FPC_1_repl)
  df_ise_FPC_1$Var2 <- factor(as.character(df_ise_FPC_1$Var2), levels = meths)
  pl <- ggplot(df_ise_FPC_1, aes(x = Var2,
                                 y = value)) +
    geom_boxplot() + labs(x="",
                          y="RMSE of estimated FPC 1 zeta", fill = "Estimation") +
    ggtitle("RMSE of estimated FPC 1 scores")
  
  pl <- pl + theme_classic() 
  
  if (!bool_mcmc_repl) {
    indiv_rmse_initial_run_FPC1 <- data.frame(name=meths, value=c(rmse_vmp["FPC_1"], rmse_mfvb["FPC_1"]))
    indiv_rmse_initial_run_FPC1 $name <- factor(indiv_rmse_initial_run_FPC1 $name, levels = meths)
    
    pl <- pl + geom_point(data= indiv_rmse_initial_run_FPC1 , aes(x=name, y=value, color=name), 
                          position=position_jitter(0), color=vec_col, size=5, pch=20)
    pl <- pl + geom_hline(yintercept=rmse_mcmc["FPC_1"], linetype="dashed", color = mcmc_col, size = 0.6) +
      theme(plot.title = element_text(size = 7, face = "bold"))
  }
  print(pl)
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/rmse_FPC_2_vmp_mfvb_mcmc_n_repl_", n_repl, ".pdf"),
        width = 3, height = 5, paper='special')
  }
  df_ise_FPC_2 <- reshape2::melt(rmse_FPC_2_repl)
  df_ise_FPC_2$Var2 <- factor(as.character(df_ise_FPC_2$Var2), levels = meths)
  pl <- ggplot(df_ise_FPC_2, aes(x = Var2,
                                 y = value)) +
    geom_boxplot() + labs(x="",
                          y="RMSE of estimated FPC 2 zeta", fill = "Estimation") +
    ggtitle("RMSE of estimated FPC 2 scores")
  
  pl <- pl + theme_classic() 
  
  if (!bool_mcmc_repl) {
    indiv_rmse_initial_run_FPC1 <- data.frame(name=meths, value=c(rmse_vmp["FPC_2"], rmse_mfvb["FPC_2"]))
    indiv_rmse_initial_run_FPC1 $name <- factor(indiv_rmse_initial_run_FPC1 $name, levels = meths)
    
    pl <- pl + geom_point(data= indiv_rmse_initial_run_FPC1 , aes(x=name, y=value, color=name), 
                          position=position_jitter(0), color=vec_col, size=5, pch=20)
    pl <- pl + geom_hline(yintercept=rmse_mcmc["FPC_2"], linetype="dashed", color = mcmc_col, size = 0.6) +
      theme(plot.title = element_text(size = 7, face = "bold"))
  }
  print(pl)
  if (bool_save) {
    dev.off()
  }
  
  
  
  if (bool_save) {
    pdf(paste0(res_dir, "/eigenfunctions_vmp_mcmc_n_repl_", n_repl, ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  
  display_eigenfunctions(L_sim, time_g, mu_g, Psi_g,
                         mu_hat = gbl_mcmc_summary[[1]],
                         list_Psi_hat = gbl_mcmc_summary[2:3],
                         mu_hat_add = list_mu_hat,
                         list_Psi_hat_add = list_list_Psi_hat,
                         mu_hat_ci = gbl_ci_mcmc_summary[[1]],
                         list_Psi_hat_ci = gbl_ci_mcmc_summary[2:3],
                         vec_col_add = c(mcmc_col, vmp_col),
                         vec_lwd = c(1.2, 0.5))
  
  if (bool_save) {
    dev.off()
    pdf(paste0(res_dir, "/eigenfunctions_mfvb_mcmc_n_repl_", n_repl, ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  
  display_eigenfunctions(L_sim, time_g, mu_g, Psi_g,
                         mu_hat = gbl_mcmc_summary[[1]],
                         list_Psi_hat = gbl_mcmc_summary[2:3],
                         mu_hat_add = list_mu_hat_mfvb,
                         list_Psi_hat_add = list_list_Psi_hat_mfvb,
                         mu_hat_ci = gbl_ci_mcmc_summary[[1]],
                         list_Psi_hat_ci = gbl_ci_mcmc_summary[2:3],
                         vec_col_add = c(mcmc_col, vmp_col),
                         vec_lwd = c(1.2, 0.5))
  
  if (bool_save) {
    dev.off()
  }
  
}

rm(compile_obj)
rm(mcmc_summary)
rm(list_res)
rm(mfvb_res)
rm(vmp_res)
rm(stan_obj)
if (bool_save) {
  save.image(file = file.path(res_dir, "output.RData"))
}
