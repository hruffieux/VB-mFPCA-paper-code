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

  job_id <- 1 # 1 to 6

  CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")
  out_dir <- file.path(CORE_DIR_ICLOUD, "mFPCA_output/")
}

main_dir <- file.path(CORE_DIR, "VB-mFPCA-paper-code/simulations/")
setwd(main_dir)

source("fun_utils.R")

require(scales)
require(bayesFPCA)
require(ggplot2)

seed <- 123
set.seed(seed)

p <- c(1, 5, 10, 25, 50, 75)[job_id]        # number of variables

N <- 200                                    # number of subjects

N_t_ind <- 1 
N_t_min <- c(2, 5, 10, 15)[N_t_ind]         # minimum number of time observations for each curve
N_t_max <- c(8, 15, 20, 25)[N_t_ind]        # maximum number of time observations for each curve

n <- matrix(sample(N_t_min:N_t_max, N*p,
                   replace = TRUE), N, p)   # number of time observations

K_Ruppert <- T
if (K_Ruppert) { # adaptive choice of K (rule of thumb adapted from Ruppert 2002)
  K <- NULL
} else {
  n_int_knots <- 5                          # number of interior knots
  K <- n_int_knots + 2                      # number of spline basis functions
}

L_sim <- 2                                  # number of simulated components
L <- 10                                     # number of estimated components

sigma_zeta <- 1/(1:L_sim) 
sigma_eps <- rep(1, p)  

tol <- 1e-5
maxit <- 5000
n_g <- 1000                                 # length of plotting grid

n_repl <- 50
n_cpus <- 50

bool_save <- T
if (bool_save) {
  res_dir <- paste0(out_dir, "/runtime_vs_p_", 
                    ifelse(K_Ruppert, "_K_Ruppert", ""), "_p_", p, "_N_", N, 
                    "_Nt_min_", N_t_min, "_max_", N_t_max, "_K_", K, "_L_sim_", 
                    L_sim, "_L_", L, "_sigeps_", unique(sigma_eps), "_tol_", tol, 
                    "_maxit_", maxit, "_n_repl_", n_repl, "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}

# Choose the simulated mean function and the FPCA basis functions:
#
f_mu <- function(t, j) (-1)^j*2/sqrt(p)*sin((2*pi)*t)
f_psi_1 <- function(t, j, p) (-1)^j*sqrt(2/p)*cos(2*pi*t)
f_psi_2 <- function(t, j, p) (-1)^j*sqrt(2/p)*sin(2*pi*t)


f_Psi <- function(time_obs, j, p) {
  ans <- cbind(f_psi_1(time_obs, j, p),
               f_psi_2(time_obs, j, p))
  return(ans)
}

list_res <- parallel::mclapply(1:n_repl, function(repl) {

  data <- generate_fpca_data(N, p, n, L_sim, n_g, sigma_eps,
                             f_mu, f_Psi,
                             vec_sd_zeta = sigma_zeta,
                             generate_from_univ = F)

  time_obs <- data$time_obs
  zeta <- data$Zeta
  Y <- data$Y

  time_g <- data$time_g
  mu_g <- data$mu_g
  Psi_g <- data$Psi_g

  set.seed(seed)

  time_vmp <- system.time({vmp_res <- run_vmp_fpca(time_obs, Y, L = L, K = NULL, 
                                                   n_g = NULL, time_g = time_g, 
                                                   tol = tol, maxit = maxit,
                                                   plot_elbo = F, Psi_g = Psi_g)})["elapsed"]

  Y_hat <- vmp_res$Y_hat
  Y_low <- vmp_res$Y_low
  Y_upp <- vmp_res$Y_upp

  mu_hat <- vmp_res$mu_hat
  list_Psi_hat <- vmp_res$list_Psi_hat
  Zeta_hat <- vmp_res$Zeta_hat
  list_zeta_ellipse <- vmp_res$list_zeta_ellipse
  
  n_iter_vmp <- vmp_res$n_iter


  set.seed(seed)

  time_mfvb <- system.time({mfvb_res <- run_mfvb_fpca(time_obs, Y, L = L, 
                                                      K = NULL, tol = tol, 
                                                      maxit = maxit, n_g = NULL,
                                                      time_g = time_g, 
                                                      Psi_g = Psi_g, 
                                                      check_elbo = FALSE)})["elapsed"]

  Y_hat_mfvb <- mfvb_res$Y_hat
  Y_low_mfvb <- mfvb_res$Y_low
  Y_upp_mfvb <- mfvb_res$Y_upp

  mu_hat_mfvb <- mfvb_res$mu_hat
  list_Psi_hat_mfvb <- mfvb_res$list_Psi_hat
  Zeta_hat_mfvb <- mfvb_res$Zeta_hat
  list_zeta_ellipse_mfvb <- mfvb_res$list_zeta_ellipse

  n_iter_mfvb <- mfvb_res$n_iter
  
  if (repl == 1) {

    n_plots_fit <- 10                        # number of plots produced (each with a random set of subjects)
    n_sample_fit <- 3                        # number of subjects shown
    p_sample <- sort(sample(1:p, min(p, 5))) # variables to display

    for (ff in 1:n_plots_fit) {

      set.seed(ff)
      N_sample_ff <- sort(sample(1:N, n_sample_fit))
      if (bool_save) {
        pdf(paste0(res_dir, "/data_samples_", paste0(N_sample_ff, collapse = "-"), ".pdf"),
            width = 6.3, height = 5, paper='special')
      }
      print(display_fit_list(p_sample, N_sample_ff, time_obs, time_g, Y, Y_hat = NULL,
                       Y_low = NULL, Y_upp = NULL))
      if (bool_save) {
        dev.off()
        pdf(paste0(res_dir, "/fits_samples_vmp_", paste0(N_sample_ff, collapse = "-"), ".pdf"),
            width = 6.3, height = 5, paper='special')
      }
      print(display_fit_list(p_sample, N_sample_ff, time_obs, time_g,
                       Y, Y_hat, Y_low, Y_upp))
      if (bool_save) {
        dev.off()
        pdf(paste0(res_dir, "/fits_samples_mfvb_", paste0(N_sample_ff, collapse = "-"), ".pdf"),
            width = 6.3, height = 5, paper='special')
      }
      print(display_fit_list(p_sample, N_sample_ff, time_obs, time_g,
                       Y, Y_hat_mfvb, Y_low_mfvb, Y_upp_mfvb))
      if (bool_save) {
        dev.off()
      }
    }


    if (bool_save) {
      pdf(paste0(res_dir, "/eigenfunctions_vmp.pdf"),
          width = 6.3, height = 5, paper='special')
    }

    print(display_eigenfunctions(L = L_sim, time_g,
                                 mu_g, Psi_g, mu_hat, list_Psi_hat, p_sample = p_sample))
    if (bool_save) {
      dev.off()
      pdf(paste0(res_dir, "/eigenfunctions_mfvb.pdf"),
          width = 6.3, height = 5, paper='special')
    }

    print(display_eigenfunctions(L = L_sim, time_g,
                                 mu_g, Psi_g, mu_hat_mfvb, list_Psi_hat_mfvb, p_sample = p_sample))
    if (bool_save) {
      dev.off()
    }

    n_plots_scores <- 4
    n_sample_ps <- 16

    for (ps in 1:n_plots_scores) {

      set.seed(ps)
      N_sample_ps <- sort(sample(1:N, n_sample_ps))

      if (bool_save) {
        pdf(paste0(res_dir, "/scores_vmp_samples_", paste0(N_sample_ps, collapse = "-"), ".pdf"),
            width = 6, height = 6, paper='special')
      }
      print(display_scores(N_sample_ps, zeta, Zeta_hat, list_zeta_ellipse,
                     mfrow = c(ceiling(sqrt(n_sample_ps)), tail(sqrt(n_sample_ps)))))
      if (bool_save) {
        dev.off()
        pdf(paste0(res_dir, "/scores_samples_mfvb_", paste0(N_sample_ps, collapse = "-"), ".pdf"),
            width = 6, height = 6, paper='special')
      }
      print(display_scores(N_sample_ps, zeta, Zeta_hat_mfvb, list_zeta_ellipse_mfvb,
                     mfrow = c(ceiling(sqrt(n_sample_ps)), tail(sqrt(n_sample_ps)))))
      if (bool_save) {
        dev.off()
      }
    }
  }

  rmse_vmp <- apply(zeta - Zeta_hat[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
  rmse_mfvb <- apply(zeta - Zeta_hat_mfvb[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
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

  create_named_list(ise_vmp, ise_mfvb, rmse_vmp, rmse_mfvb, time_vmp, time_mfvb, n_iter_vmp, n_iter_mfvb)
}, mc.cores = n_cpus)


list_ise_vmp <- lapply(list_res, "[[", "ise_vmp")
list_rmse_vmp <- lapply(list_res, "[[", "rmse_vmp")
list_ise_mfvb <- lapply(list_res, "[[", "ise_mfvb")
list_rmse_mfvb <- lapply(list_res, "[[", "rmse_mfvb")
vec_runtime_vmp <- sapply(list_res, "[[", "time_vmp")
vec_runtime_mfvb <- sapply(list_res, "[[", "time_mfvb")
vec_n_iter_vmp <- sapply(list_res, "[[", "n_iter_vmp")
vec_n_iter_mfvb <- sapply(list_res, "[[", "n_iter_mfvb")

names(list_ise_vmp) <- names(list_rmse_vmp) <- names(list_ise_mfvb) <- 
  names(list_rmse_mfvb) <- names(vec_runtime_vmp) <- names(vec_runtime_mfvb) <- 
  names(vec_n_iter_vmp) <- names(vec_n_iter_mfvb) <- paste0("repl_", 1:n_repl)


if (bool_save) {
  save(list_ise_vmp, list_ise_mfvb,
       list_rmse_vmp, list_rmse_mfvb,
       vec_runtime_vmp, vec_runtime_mfvb,
       vec_n_iter_vmp, vec_n_iter_mfvb,
       file = file.path(res_dir, "output.RData"))
}

print("AVG RUNTIME VMP")
print(mean(vec_runtime_vmp))

print("AVG RUNTIME MFVB")
print(mean(vec_runtime_mfvb))

if (bool_save) {
  pdf(paste0(res_dir, "/runtime.pdf"), width = 6, height = 3.3, paper='special')
}
mat_runtime <- cbind(vec_runtime_mfvb / 60, vec_runtime_vmp / 60)
colnames(mat_runtime) <- c("MFVB", "VMP")
boxplot(mat_runtime, main = "Runtime in minutes")
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/number_of_iterations_tol_", tol, ".pdf"), width = 6, height = 3.3, paper='special')
}
mat_n_iter <- cbind(vec_n_iter_mfvb, vec_n_iter_vmp)
colnames(mat_n_iter) <- c("MFVB", "VMP")
boxplot(mat_n_iter, main = "Total number of iterations before convergence")
if (bool_save) {
  dev.off()
}
