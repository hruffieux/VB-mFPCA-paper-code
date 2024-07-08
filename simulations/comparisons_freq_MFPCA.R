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
  
  job_id <- 4 # 1 to 4
  
  CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")
  out_dir <- file.path(CORE_DIR_ICLOUD, "mFPCA_output/")
}

main_dir <- file.path(CORE_DIR, "bayesian-mFPCA-paper-code/simulations/")
setwd(main_dir)

require(bayesFPCA)
require(MFPCA)

source("fun_utils.R")

seed <- 1
set.seed(seed)

p <- 3                                         # number of variables
N <- 100                                       # number of curves
N_t_min <- seq(5, 65, by = 20)[job_id]
N_t_max <- seq(35, 95, by = 20)[job_id]
n <- matrix(sample(N_t_min:N_t_max, N*p,
                   replace = TRUE), N, p)      # number of time observations

model_choice_K <- T
if (model_choice_K) {
  K <- NULL
  n_cpus <- 4 
} else {
  bool_K_Ruppert <- T
  if (bool_K_Ruppert) {
    K <- NULL                                   # default based on the rule in Ruppert 2002
  } else {
    n_int_knots <- 5                            # number of interior knots
    K <- n_int_knots + 2                        # number of spline basis functions
  }
}

L_sim <- 2
L <- 10.                                        # number of FPCA basis functions

tol  <- 1e-5                                    # convergence tolerance
maxit <- 1000                                   # maximum number of mfvb iterations

n_g <- 1000                                     # length of the plotting grid

vec_sd_zeta <- 1/(1:L_sim)
vec_sd_eps <- rep(1, p)                         # sd of the residuals (here same for all variables)


generate_from_univ <- T

if (generate_from_univ) {
  id_strength <- 6 
  rho_Zeta <- seq(0, 1, by = 0.2)[id_strength]
  vec_rho_Zeta <- rep(rho_Zeta, L)
  cat("Correlation of simulated scores across variables for eigenfunctions 1 to L: ",
      vec_rho_Zeta)
} else {
  vec_rho_Zeta <- NULL
}

bool_save <- T
if (bool_save) {
  res_dir <- paste0(out_dir, "/sim_Happ_", ifelse(generate_from_univ,
                           paste0("gen_from_univ_corr_",
                                  paste0(format(rho_Zeta, digits = 2),
                                         collapse = "-"), "_"), ""),
                    "p_", p, "_N_", N, "_Nt_min_", N_t_min, "_max_", N_t_max, 
                    "_L_sim_", L_sim, "_L_", L, ifelse(model_choice_K, 
                                                       "_model_choice_K", 
                                                       paste0("K_", K)),
                    "_sd_zeta_", paste0(format(vec_sd_zeta , digits = 2), collapse = "-"),
                    "_sd_eps_", unique(vec_sd_eps), "_tol_", tol, "_maxit_", 
                    maxit, "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}


####################################################
#
#  SIMULATE  THE  DATA
#
####################################################

  
if (rho_Zeta == 1) { # perfect correlation, so we can directly simulate from the multivariate model
  generate_from_univ <- FALSE 
  vec_rho_Zeta <- NULL
}

f_mu <- function(t, j) (-1)^j*2*sin((2*pi+j)*t)
f_psi_1 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(2*pi*t)
f_psi_2 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(2*pi*t)


f_Psi <- function(time_obs, j, p) {
  ans <- cbind(f_psi_1(time_obs, j, p), 
               f_psi_2(time_obs, j, p))
  return(ans)
}

data <- generate_fpca_data(N, p, n, L = L_sim, n_g, vec_sd_eps, f_mu, f_Psi,
                           vec_sd_zeta = vec_sd_zeta,
                           generate_from_univ = generate_from_univ,
                           vec_rho_Zeta = vec_rho_Zeta)

time_obs <- data$time_obs
Zeta <- data$Zeta
Y <- data$Y

time_g <- data$time_g
mu_g <- data$mu_g
Psi_g <- data$Psi_g

mfd <- convert_data_for_happ(time_obs, Y)


list_type_j <- list(type = "uFPCA")
list_type <- rep(list(list_type_j), p)
bool_ci <- F # pointwise bootstrap confidence bands are calculated for eigenvalues and for Psi_l
             # does not give the confidence bands for the fit itself??
if(bool_ci) {
  nb_boot <- 100
} else {
  nb_boot <- NULL
}

set.seed(seed)

run_time_happ <- system.time({ 
  happ_res <- MFPCA(mfd, 
                    M = L_sim, 
                    uniExpansions = list_type, 
                    fit = TRUE, 
                    bootstrap = bool_ci,
                    nBootstrap = nb_boot
                    )
  }) 
summary(happ_res)

print(run_time_happ)

list_time_g_happ <- lapply(happ_res$meanFunction, function(fd_j) as.data.frame(fd_j)$argvals1)
sapply(list_time_g_happ, length) # evaluation grid (differs for each variable...)

mu_hat_happ <- lapply(happ_res$meanFunction, function(fd_j) as.data.frame(fd_j)$X)
sapply(mu_hat_happ, length)

list_Psi_hat_happ <- lapply(1:L_sim, function(l) 
  lapply(happ_res$functions, function(fd_j) { 
    df <- as.data.frame(fd_j); df$X[df$obs == l] 
    }))

Zeta_hat_happ <- happ_res$scores

Y_hat_happ <- lapply(1:N, function(i) lapply(1:p, function(j) {
  df_j <- as.data.frame(happ_res$fit[[j]])
  df_j$X[df_j$obs == i]
}))

jj <- 1
Psi_g_happ_j <- f_Psi(list_time_g_happ[[jj]], j = jj, p = p)  

for (l in 1:L_sim) {
  cprod_sign <- sign(cprod(list_Psi_hat_happ[[l]][[jj]], Psi_g_happ_j[,l]))
  
  for (j in 1:p) {
    list_Psi_hat_happ[[l]][[j]] <- cprod_sign*list_Psi_hat_happ[[l]][[j]]
  }

  Zeta_hat_happ[,l] <- cprod_sign*Zeta_hat_happ[,l]
}
# }

set.seed(seed)

if (model_choice_K) {
  run_time_mfvb <- system.time({
    mfvb_res <- run_mfvb_fpca_model_choice(time_obs, Y, L = L, n_g = NULL, time_g = time_g,
                              tol = tol, maxit = maxit, Psi_g = Psi_g, n_cpus = n_cpus)
  })
} else {
  run_time_mfvb <- system.time({
    mfvb_res <- run_mfvb_fpca(time_obs, Y, L = L, K = K, n_g = NULL, time_g = time_g, 
                              tol = tol, maxit = maxit, Psi_g = Psi_g)
  })
}

print(run_time_mfvb)

Y_hat <- mfvb_res$Y_hat
Y_low <- mfvb_res$Y_low
Y_upp <- mfvb_res$Y_upp
mu_hat <- mfvb_res$mu_hat
list_Psi_hat <- mfvb_res$list_Psi_hat
Zeta_hat <- mfvb_res$Zeta_hat
list_zeta_ellipse <- mfvb_res$list_zeta_ellipse

vec_col_add <- c("dodgerblue1", "darkgrey")

n_plots_fit <- 10   # number of plots produced (each with a random set of subjects)
n_sample_fit <- 3   # number of subjects shown
p_sample <- 1:p     # variables to display

for (ff in 1:n_plots_fit) {
  
  set.seed(ff)
  N_sample_ff <- sort(sample(1:N, n_sample_fit))

  if (bool_save) {
    pdf(paste0(res_dir, "/fits_samples_", paste0(N_sample_ff, collapse = "-"), ".pdf"),
        width = 6.3, height = 5, paper='special')
  }
  display_fit_list_happ(p_sample, N_sample_ff, time_obs, time_g, time_g_happ, 
                        Y, Y_hat, Y_low, Y_upp, Y_hat_happ, 
                        col = vec_col_add[1], col_happ = vec_col_add[2])
  if (bool_save) {
    dev.off()
  }
  
}
             

if (bool_save) {
  pdf(paste0(res_dir, "eigenfunctions.pdf"), width = 7, height = 7, paper='special')
}
vec_flip_happ <- c(1, 1) 
display_eigenfunctions_happ(L = L_sim, time_g, list_time_g_happ, 
                            mu_g, Psi_g, 
                            mu_hat, list_Psi_hat, 
                            mu_hat_happ, list_Psi_hat_happ, 
                            vec_col_add = vec_col_add,
                            vec_flip_happ = vec_flip_happ, 
                            data_col = "black")

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
  display_scores(N_sample_ps, Zeta, Zeta_hat[,1:L_sim], list_zeta_ellipse,
                 Zeta_hat_add = Zeta_hat_happ,
                 zeta_ellipse_add = NULL,
                 vec_col = vec_col_add, data_col = "black",
                 mfrow = c(ceiling(sqrt(n_sample_ps)), tail(sqrt(n_sample_ps))))
  if (bool_save) {
    dev.off()
  }
  
}


if (!generate_from_univ) {
  rmse_mfpca <- apply(Zeta - Zeta_hat[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
  rmse_happ <- apply(Zeta - Zeta_hat_happ, 2, function(x) sqrt(mean(x^2)))
  
  cat("mFPCA rmse is:", rmse_mfpca, "\n")
  cat("Happ rmse is:", rmse_happ, "\n")
}

ise_mfpca <- ise_happ <- matrix(NA, L_sim+1, p)

for (l in 1:(L_sim+1)) {
  
  for(j in 1:p) {
    
      mu_g_happ_j <- f_mu(list_time_g_happ[[j]], j = j)
      Psi_g_happ_j <- f_Psi(list_time_g_happ[[j]], j = j, p = p)
    
    if (l == 1) {
      ise_mfpca[l, j] <- trapint(time_g, (mu_hat[,j] - mu_g[[j]])^2)
      ise_happ[l, j] <- trapint(list_time_g_happ[[j]], (mu_hat_happ[[j]] - mu_g_happ_j)^2) 
    } else {
      ise_mfpca[l, j] <- trapint(time_g, (list_Psi_hat[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
      ise_happ[l, j] <- trapint(list_time_g_happ[[j]], (list_Psi_hat_happ[[l-1]][[j]] - Psi_g_happ_j[,l-1])^2) 
    }
    
  }
  
}
rownames(ise_mfpca) <- rownames(ise_happ) <- c("mu", paste0("Psi_", 1:L_sim))
colnames(ise_mfpca) <- colnames(ise_happ) <- paste0("Variable_", 1:p)


add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

if (bool_save) {
  pdf(paste0(res_dir, "/ISE_one_replicate.pdf"), width = 6.5, height = 5, 
      paper='special')
}
par(mfrow = c(1, 1))
plot(1:p, log(ise_mfpca["mu",]), col = "black", type = "b", lwd = 2, xaxt='n',
     xlab = "", ylab = "Log ISE", 
     main = "Integrated squared error mean and eigenfunctions",
     ylim = c(min(c(as.vector(log(ise_mfpca)), as.vector(log(ise_happ)))), 
              max(c(as.vector(log(ise_mfpca)), as.vector(log(ise_happ))))))
axis(1, at = 1:p, labels = paste0("Variable ", 1:p))
lines(1:p, log(ise_happ["mu",]), col = "black", type = "b", lwd = 2, lty = 2)
lines(1:p, log(ise_mfpca["Psi_1",]), col = "grey50", type = "b", lwd = 2)
lines(1:p, log(ise_happ["Psi_1",]), col = "grey50", type = "b", lwd = 2, lty = 2)
lines(1:p, log(ise_mfpca["Psi_2",]), col = "grey80", type = "b", lwd = 2)
lines(1:p, log(ise_happ["Psi_2",]), col = "grey80", type = "b", lwd = 2, lty = 2)
add_legend("bottomright", c("mu", "Psi_1", "Psi_2", "mFPCA    ", "Happ"), horiz = T,
           lwd = 2, col = c("black", "grey50", "grey80", "black", "black"), 
           lty = c(rep(NA, 3), 1, 2), pch = c(rep(20, 3), rep(NA, 2)), bty = "n")
if (bool_save) {
  dev.off()
  save.image(file = file.path(res_dir, "output.RData"))
}
