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
  
  job_id <- 5 # 1 to 6
  
  CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")
  out_dir <- file.path(CORE_DIR_ICLOUD, "mFPCA_output/")
}

main_dir <- file.path(CORE_DIR, "bayesian-mFPCA-paper-code/simulations/")
setwd(main_dir)

require(bayesFPCA)

source("fun_utils.R")

seed <- 1
set.seed(seed)

p <- 6                                          # number of responses
N <- c(50, 75, 100, 150, 200)[job_id]           # number of curves
N_t_min <- 50
N_t_max <- 75
n <- matrix(sample(N_t_min:N_t_max, N*(p-1),
                   replace = TRUE), N, p-1)     # number of time observations

N_t_min_small <- 5
N_t_max_small <- 10
n <- cbind(sample(N_t_min_small:N_t_max_small, N,replace = TRUE), n)

model_choice_K <- TRUE
if (model_choice_K) {
  K <- NULL
} else {
  bool_K_Ruppert <- F
  if (bool_K_Ruppert) {
    K <- NULL                                   # default based on the rule in Ruppert 2002
  } else {
    n_int_knots <- 5                            # number of interior knots
    K <- n_int_knots + 2                        # number of spline basis functions
  }
}

L_sim <- 2
L <- 10                                        # number of FPCA basis functions

tol  <- 1e-5                                   # convergence tolerance 
maxit <- 1000                                  # maximum number of vmp iterations

n_g <- 1000                                    # length of the plotting grid

n_cpus <- 16

bool_multiple_repl_eigen <- FALSE              

if (bool_multiple_repl_eigen) {
  n_repl_eigen <- 4
  n_cpus_eigen <- 4
}

ind_exp <- 1
exponent_sd_zeta <- c(-1, -1/2, -1/4, -1/8, -1/16)[ind_exp]
vec_sd_zeta <- (1:L_sim)^exponent_sd_zeta
sd_eps <- 2.5
vec_sd_eps <- rep(sd_eps, p)


id_strength <- 6 # up to 6
rho_Zeta <- seq(0, 1, by = 0.2)[id_strength]
vec_rho_Zeta <- rep(rho_Zeta, L_sim)
cat("Correlation of simulated scores across variables for eigenfunctions 1 to L_sim: ",
    vec_rho_Zeta)


bool_save <- T
if (bool_save) {
  res_dir <- paste0(out_dir, "/comparison_sparse_",
                    ifelse(model_choice_K, "model_choice_K_", ""),
                    paste0(format(rho_Zeta, digits = 2), collapse = "-"),
                    "p_", p, "_N_", N, "_Nt_min_small_", N_t_min_small, 
                    "_max_small_", N_t_max_small, "_Nt_min_", N_t_min, "_max_", 
                    N_t_max, "_K_", K, "_L_sim_", L_sim, "_L_", L, "_sd_zeta_", 
                    paste0(format(vec_sd_zeta , digits = 2), collapse = "-"),
                    "_sd_eps_", unique(vec_sd_eps), "_tol_", tol, "_maxit_", 
                    maxit, "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}

if (rho_Zeta == 1) { # perfect correlation, so we can directly simulate from the multivariate model
  vec_rho_Zeta <- NULL
}

f_mu <- function(t, j) (-1)^j*2*sin((2*pi+j)*t)

f_psi_1 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(2*pi*t)
f_psi_2 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(2*pi*t)
f_psi_3 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(3*pi*t)
f_psi_4 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(3*pi*t)
f_psi_5 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(4*pi*t)
f_psi_6 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(4*pi*t)
f_psi_7 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(5*pi*t)
f_psi_8 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(5*pi*t)


f_Psi <- function(time_obs, j, p) {
  ans <- cbind(f_psi_1(time_obs, j, p),
               f_psi_2(time_obs, j, p),
               f_psi_3(time_obs, j, p),
               f_psi_4(time_obs, j, p),
               f_psi_5(time_obs, j, p), 
               f_psi_6(time_obs, j, p),
               f_psi_7(time_obs, j, p),
               f_psi_8(time_obs, j, p))[,1:L_sim, drop = F]
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


set.seed(seed)

if (model_choice_K) {
  system.time(res <- run_mfvb_fpca_model_choice(time_obs, Y, L = L, n_g = NULL,
                                   time_g = time_g, 
                                   tol = tol, maxit = maxit, Psi_g = Psi_g, n_cpus = n_cpus))
} else {
  system.time(res <- run_mfvb_fpca(time_obs, Y, L = L, K = K, n_g = NULL,
                                   time_g = time_g, 
                                   tol = tol, maxit = maxit,
                                   plot_elbo = FALSE, Psi_g = Psi_g))
}

Y_hat <- res$Y_hat
Y_low <- res$Y_low
Y_upp <- res$Y_upp
mu_hat <- res$mu_hat
list_Psi_hat <- res$list_Psi_hat
Zeta_hat <- res$Zeta_hat
list_zeta_ellipse <- res$list_zeta_ellipse

cumulated_pve <- res$cumulated_pve
L_thres_95 <- sum(cumulated_pve < 95) + 1
L_thres_99 <- sum(cumulated_pve < 99) + 1

K <- unique(res$K)


bool_univ <- T

if (bool_univ) { 

  list_time_obs_univ <- lapply(1:p, function(j) lapply(time_obs, function(time_obs_i) time_obs_i[j]))
  list_Y_univ <- lapply(1:p, function(j) lapply(Y, function(Y_i) Y_i[j]))
  
  if (model_choice_K) {
  list_res_univ <- lapply(seq_along(list_Y_univ), function(j)
    run_mfvb_fpca_model_choice(list_time_obs_univ[[j]], list_Y_univ[[j]], L = L, n_g = NULL,
                 time_g = time_g, tol = tol, maxit = maxit,
                 Psi_g = Psi_g[j], n_cpus = n_cpus))
  } else {
    list_res_univ <- parallel::mclapply(seq_along(list_Y_univ), function(j)
      run_mfvb_fpca(list_time_obs_univ[[j]], list_Y_univ[[j]], L = L, K = K, n_g = NULL,
                    time_g = time_g, tol = tol, maxit = maxit,
                    plot_elbo = FALSE, Psi_g = Psi_g[j]),
      mc.cores = min(n_cpus, p)) 
  }
  names(list_res_univ) <- paste0("run_univ_Y", 1:p)

  Y_hat_univ <- lapply(1:N, function(i) Reduce(cbind, sapply(list_res_univ, "[[", "Y_hat")[i,]))
  Y_low_univ <- lapply(1:N, function(i) Reduce(cbind, sapply(list_res_univ, "[[", "Y_low")[i,]))
  Y_upp_univ <- lapply(1:N, function(i) Reduce(cbind, sapply(list_res_univ, "[[", "Y_upp")[i,]))

  mu_hat_univ <- Reduce(cbind, lapply(1:p, function(j) list_res_univ[[j]]$mu_hat))
  list_Psi_hat_univ <- lapply(1:L, function(l) Reduce(cbind, lapply(1:p, function(j) list_res_univ[[j]]$list_Psi_hat[[l]])))

  cumulated_pve_univ <- sapply(list_res_univ, "[[", "cumulated_pve")
  L_thres_95_univ <- apply(cumulated_pve_univ, 2, function(cum_pve_j) sum(cum_pve_j < 95) + 1)
  L_thres_99_univ <- apply(cumulated_pve_univ, 2, function(cum_pve_j) sum(cum_pve_j < 99) + 1)
  
  K_univ <- sapply(list_res_univ, "[[", "K")
  
  ise_mfpca <- ise_fpca <- matrix(NA, L_sim+1, p)

  for (l in 1:(L_sim+1)) {

    for(j in 1:p) {

      if (l == 1) {
        ise_mfpca[l, j] <- trapint(time_g, (mu_hat[,j] - mu_g[[j]])^2)
        ise_fpca[l, j] <- trapint(time_g, (mu_hat_univ[,j] - mu_g[[j]])^2)
      } else {
        ise_mfpca[l, j] <- trapint(time_g, (list_Psi_hat[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
        ise_fpca[l, j] <- trapint(time_g, (list_Psi_hat_univ[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
      }

    }

  }
  rownames(ise_mfpca) <- rownames(ise_fpca) <- c("mu", paste0("Psi_", 1:L_sim))
  colnames(ise_mfpca) <- colnames(ise_fpca) <- paste0("Variable_", 1:p)

  if (bool_save) {
    save(ise_mfpca, ise_fpca, res, 
         file = file.path(out_dir, "output.RData"))
  }

} else if (bool_save) {
  save(res, file = file.path(res_dir, "output.RData"))
}



Zeta_hat_mFPCA <- res$Zeta_hat
sd_Zeta_hat_mFPCA <- cbind(sqrt(sapply(res$Cov_zeta_hat, function(cov_ii) cov_ii[1, 1])),
                           sqrt(sapply(res$Cov_zeta_hat, function(cov_ii) cov_ii[2, 2])))
colnames(sd_Zeta_hat_mFPCA) <- c("FPC_1", "FPC_2")
rownames(sd_Zeta_hat_mFPCA) <- rownames(Zeta_hat_mFPCA)

list_Zeta_hat <- lapply(list_res_univ, function(ll) sqrt(p)*ll$Zeta_hat) # adjust as not orthonormal in L2 but multivariate space
list_sd_Zeta_hat <- lapply(list_res_univ, function(ll) { sd_Zeta_hat <- cbind(sqrt(sapply(ll$Cov_zeta_hat, function(cov_ii) p*cov_ii[1, 1])),
                                                          sqrt(sapply(ll$Cov_zeta_hat, function(cov_ii) p*cov_ii[2, 2])))
                                                          rownames(sd_Zeta_hat) <- rownames(Zeta_hat_mFPCA)
                                                          colnames(sd_Zeta_hat) <- c("FPC_1", "FPC_2")
                                                          sd_Zeta_hat
                                                          }) 


list_Zeta_hat <- append(list_Zeta_hat, list(Zeta_hat_mFPCA))
list_sd_Zeta_hat <- append(list_sd_Zeta_hat, list(sd_Zeta_hat_mFPCA))


rownames(Zeta) <- rownames(Zeta_hat_mFPCA)
colnames(Zeta) <- colnames(Zeta_hat_mFPCA)[1:L_sim]
list_Zeta_hat <- append(list_Zeta_hat, list(Zeta))

sd_Zeta <- sd_Zeta_hat_mFPCA
sd_Zeta[] <- 0
list_sd_Zeta_hat <- append(list_sd_Zeta_hat, list(sd_Zeta))

names(list_Zeta_hat) <- names(list_sd_Zeta_hat)<- c(paste0("FPCA Y", 1:p), "mFPCA", "True")


require(tidyverse)
require(ggforestplot)
require(dplyr)
require(reshape2)

list_tb <- lapply(seq_along(list_Zeta_hat), function(id) { df <- data.frame(rep(names(list_Zeta_hat)[id], nrow(list_Zeta_hat[[id]])),
                                                                            melt(list_Zeta_hat[[id]]),
                                                                            melt(list_sd_Zeta_hat[[id]])[, "value"])
names(df) <- c("Type", "Patient", "Component", "estimate", "sd")
df
})
names(list_tb) <- names(list_Zeta_hat)

tb <- do.call(rbind, list_tb) %>% arrange(Patient)  %>% arrange(Component)
rownames(tb) <- NULL
tb$Type <- factor(tb$Type, levels = names(list_tb)[c(p+2, p+1, p:1)])

tb$Patient <- gsub("_", " ", tb$Patient)

n_plots_scores <- 50
n_sample_ps <- 8

for (ps in 1:n_plots_scores) {
  
  set.seed(ps)
  patients <- paste0("subj ", sample(1:N, n_sample_ps))

  tb_fpc1 <- as_tibble(tb[ tb$Patient %in% patients &
      tb$Component == "FPC_1", ])
  p1 <- forestplot(
    df = tb_fpc1,
    name = Patient,
    estimate = estimate,
    se = sd,
    logodds = FALSE,
    colour = Type,
    shape = Type,
    title = "Scores FPC 1",
    xlab = "zeta 1 with 95% credible interval"
  ) + ggplot2::scale_shape_manual(
    values = c(24L, 23L, rep(21L, p)),
    labels = names(list_tb)[c(p+2,p+1, p:1)]
  )  #+ theme(legend.position = "none")
  
  tb_fpc2 <- as_tibble(tb[ tb$Patient %in% patients &
                             tb$Component == "FPC_2", ])
  
  p2 <- forestplot(
    df = tb_fpc2,
    name = Patient,
    estimate = estimate,
    se = sd,
    logodds = FALSE,
    colour = Type,
    shape = Type,
    title = "Scores FPC 2",
    xlab = "zeta 2 with 95% credible interval"
  ) + ggplot2::scale_shape_manual(
    values = c(24L, 23L, rep(21L, p)),
    labels = names(list_tb)[c(p+2,p+1, p:1)]
  )
  
  if (bool_save) { 
    pdf(paste0(res_dir, "/scores_univ_multiv_forest_samples_", paste0(patients, collapse = "-"), ".pdf"),
        width = 11, height = 7, paper='special')
  }
  gridExtra::grid.arrange(p1, p2, nrow = 1)
  if (bool_save) {
    dev.off()
  }

}
