rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

out_dir <- file.path(CORE_DIR, "VB-mFPCA-paper-code/output/")

main_dir <- file.path(CORE_DIR, "VB-mFPCA-paper-code/simulations/")
setwd(main_dir)

require(bayesFPCA)

source("fun_utils.R")

seed <- 1
set.seed(seed)

p <- 3                                        # number of responses
N <- 200                                      # number of curves
n_obs <- list(30:70, 40:70, 30:50)
n <- matrix(Reduce(                           # number of time observations
  cbind, lapply(
    n_obs[1:p], sample,
    size = N, replace = TRUE
  )
), nrow = N)

K <- NULL                                     # number of spline basis functions (set internally according to Ruppert 2002)
L_true <- 4                                   # number of simulated FPCA basis functions
L <- 10                                       # number of inferred FPCA basis functions (conservative upper bound)

sigma_zeta_vec <- 1/(1:L_true)                # sd of simulated scores
sigma_eps <- rep(1, p)                        # sd of simulated residuals

col_fixed <- "darkseagreen3"
col_est <- "#FFB90F"
col_sim <- "black"
n_g <- 1000                                   # length of the plotting grid

bool_save <- T

mu_func <- function(time_obs, j) {
  ans <- ((-1)^j)*sin(pi*time_obs)
  return(ans)
}

# Setting where functions explaining less variation are more wiggly
psi_func <- function(time_obs, j, l) {
  ans <- ((-1)^j)*sqrt(2/p)*sin((2*l-1)*pi*time_obs)
  return(ans)
}

Psi_func <- function(time_obs, j, p) {
  ans <- sapply(1:L, function(l) psi_func(time_obs, j, l))
  return(ans)
}


mfpca_data <- generate_fpca_data(N, p, n, L_true, n_g, sigma_eps, mu_func,
                                 Psi_func, vec_sd_zeta = sigma_zeta_vec)

time_obs <- mfpca_data$time_obs
Zeta <- mfpca_data$Zeta
mu_g <- mfpca_data$mu_g
Psi_g <- mfpca_data$Psi_g
Y <- mfpca_data$Y


rt_mfpca <- system.time(mfpca_res <- run_mfvb_fpca(time_obs, Y, L = L,
                                                   n_g = n_g, Psi_g = Psi_g,
                                                   seed = seed))


rt_mfpca_est <- system.time(mfpca_res_est <- run_mfvb_fpca(time_obs, Y, L = L,
                                                           n_g = n_g, Psi_g = Psi_g,
                                                           fixed_score_variance = F))

time_g <- mfpca_res$time_g
Y_hat <- mfpca_res$Y_hat
Y_low <- mfpca_res$Y_low
Y_upp <- mfpca_res$Y_upp
mu_hat <- mfpca_res$mu_hat
list_Psi_hat <- mfpca_res$list_Psi_hat
Zeta_hat <- mfpca_res$Zeta_hat
list_zeta_ellipse <- mfpca_res$list_zeta_ellipse

Y_hat_est <- mfpca_res_est$Y_hat
Y_low_est <- mfpca_res_est$Y_low
Y_upp_est <- mfpca_res_est$Y_upp
mu_hat_est <- mfpca_res_est$mu_hat
list_Psi_hat_est <- mfpca_res_est$list_Psi_hat
Zeta_hat_est <- mfpca_res_est$Zeta_hat
list_zeta_ellipse_est <- mfpca_res_est$list_zeta_ellipse

set.seed(seed)
n_sample <- 6                                 # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))       # specific curves for the plots

display_fit_list(1:p, N_sample, time_obs, time_g, Y, Y_hat, Y_low, Y_upp,
                 Y_hat_add = Y_hat_est, Y_low_add = Y_low_est, Y_upp_add = Y_upp_est,
                 col = col_fixed,
                 col_add = col_est)

if (bool_save) {
  pdf(paste0(out_dir, "/latent_functions.pdf"), width = 7.5, height = 7.5, paper='special')
}
display_eigenfunctions(L_true, time_g, mu_g, Psi_g, mu_hat, list_Psi_hat,
                       mu_hat_add = mu_hat_est, list_Psi_hat_add = list_Psi_hat_est,
                       data_col = col_sim, p_sample = NULL,
                       vec_col_add = c(col_fixed, col_est))

if (bool_save) {
  dev.off()
}

if (L > 1) { # scores for the first two components

  set.seed(seed)
  n_sample <- 12                                # number of curves for the plots
  N_sample <- sort(sample(1:N, n_sample))       # specific curves for the plots

  # scores FPC 1 vs 2
  display_scores(N_sample, Zeta, Zeta_hat, list_zeta_ellipse,
                 Zeta_hat_add = Zeta_hat_est, zeta_ellipse_add = list_zeta_ellipse_est,
                 data_col = col_sim, vec_col = c(col_fixed, col_est))


  list_zeta_ellipse_l_3_4 <- list_zeta_ellipse_l_3_4_est <- vector("list", length=N)
  for(i in 1:N) {

    zeta_mean <- Zeta_hat[i,][3:4]
    zeta_mean_est <- Zeta_hat_est[i,][3:4]

    zeta_ellipse <- ellipse(
      mfpca_res$Cov_zeta_hat[[i]][3:4, 3:4],
      centre = zeta_mean,
      level = 0.95
    )

    zeta_ellipse_est <- ellipse(
      mfpca_res_est$Cov_zeta_hat[[i]][3:4, 3:4],
      centre = zeta_mean_est,
      level = 0.95
    )

    list_zeta_ellipse_l_3_4[[i]] <- zeta_ellipse
    list_zeta_ellipse_l_3_4_est[[i]] <- zeta_ellipse_est

  }

  if (bool_save) {
    pdf(paste0(out_dir, "/scores_FPC_3_vs_4.pdf"), width = 7.5, height = 7.5, paper='special')
  }
  # custom function to display scores FPC 3 vs 4
  display_scores_custom_axes(l_ind = c(3, 4), N_sample, Zeta, Zeta_hat,
                             list_zeta_ellipse_l_3_4,
                             Zeta_hat_add = Zeta_hat_est,
                             zeta_ellipse_add = list_zeta_ellipse_l_3_4_est,
                             data_col = col_sim,
                             vec_col = c(col_fixed, col_est))
  if (bool_save) {
    dev.off()
  }

}


cat(paste0("Number of iterations -- bayesFPCA: ", mfpca_res$n_iter,
           ", bayesFPCA est: ", mfpca_res_est$n_iter, "\n",
           "Total runtime (sec) -- bayesFPCA: ", format(rt_mfpca["elapsed"]),
           ", bayesFPCA est: ", format(rt_mfpca_est["elapsed"]), "\n"))
