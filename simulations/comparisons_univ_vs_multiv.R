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

main_dir <- file.path(CORE_DIR, "bayesian-mFPCA-paper-code/simulations/")
setwd(main_dir)

require(bayesFPCA)

source("fun_utils.R")

seed <- 123
set.seed(seed)

p <- 6                                        # number of responses
N <- 100                                      # number of curves
N_t_min <- 5 
N_t_max <- 10
n <- matrix(sample(N_t_min:N_t_max, N*p,
                   replace = TRUE), N, p)     # number of time observations

model_choice_K <- T
if (model_choice_K) {
  K <- NULL
} else {
  bool_K_Ruppert <- T
  if (bool_K_Ruppert) {
    K <- NULL                                 # default based on the rule in Ruppert 2002
  } else {
    n_int_knots <- 5                          # number of interior knots
    K <- n_int_knots + 2                      # number of spline basis functions
  }
}

L_sim <- 2                                    # number of FPCA basis functions
L <- 10

tol  <- 1e-5                                  # convergence tolerance 
maxit_mfvb <- 1000                            # maximum number of vmp iterations


n_g <- 1000                                   # length of the plotting grid

n_repl <- 100
n_cpus <- 16

ind_exp <- 4
exponent_sd_zeta <- c(-1, -1/2, -1/4, -1/8, -1/16)[ind_exp]
vec_sd_zeta <- (1:L_sim)^exponent_sd_zeta
vec_sd_eps <- rep(1, p)                       

generate_from_univ <- T
rho_Zeta <- seq(0, 1, by = 0.2)[job_id]
vec_rho_Zeta <- rep(rho_Zeta, L_sim)
cat("Correlation of simulated scores across variables for eigenfunctions 1 to L_sim: ",
    vec_rho_Zeta)


bool_save <- T
if (bool_save) {
  res_dir <- paste0(out_dir, "/comparison_univ_vs_multiv_n_repl_", n_repl, "_",
                    ifelse(model_choice_K, "model_choice_K_", ""),
                                  paste0(format(rho_Zeta, digits = 2), ""),
                    "p_", p, "_N_", N, "_Nt_min_", N_t_min, "_max_", N_t_max, 
                    "_K_", K, "_L_sim_", L_sim, "_L_", L, "_sd_zeta_", 
                    paste0(format(vec_sd_zeta , digits = 2), collapse = "-"),
                    "_sd_eps_", unique(vec_sd_eps), "_tol_", tol, "_maxit_", 
                    maxit_mfvb, "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}

if (rho_Zeta == 1) { # perfect correlation, so we can directly simulate from the multivariate model
  generate_from_univ <- FALSE
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


list_res <- parallel::mclapply(1:n_repl, function(repl) {

  data <- generate_fpca_data(N, p, n, L = L_sim, n_g, vec_sd_eps,
                             f_mu, f_Psi,
                             vec_sd_zeta = vec_sd_zeta,
                             generate_from_univ = generate_from_univ,
                             vec_rho_Zeta = vec_rho_Zeta)

  time_obs <- data$time_obs
  Zeta <- data$Zeta
  Y <- data$Y

  time_g <- data$time_g
  mu_g <- data$mu_g
  Psi_g <- data$Psi_g

  if (model_choice_K) {
    mfvb_res <- run_mfvb_fpca_model_choice(time_obs, Y, L = L, n_g = NULL,
                                      time_g = time_g, 
                                      tol = tol, maxit = maxit_mfvb, Psi_g = Psi_g, n_cpus = 1)
  } else {
    mfvb_res <- run_mfvb_fpca(time_obs, Y, L = L, K = K, n_g = NULL,
                         time_g = time_g, 
                         tol = tol, maxit = maxit_mfvb,
                         plot_elbo = FALSE, Psi_g = Psi_g)
  }

  Y_hat <- mfvb_res$Y_hat
  Y_low <- mfvb_res$Y_low
  Y_upp <- mfvb_res$Y_upp

  mu_hat <- mfvb_res$mu_hat
  list_Psi_hat <- mfvb_res$list_Psi_hat
  Zeta_hat <- mfvb_res$Zeta_hat
  list_zeta_ellipse <- mfvb_res$list_zeta_ellipse

  cumulated_pve <- mfvb_res$cumulated_pve
  L_thres_95 <- sum(cumulated_pve < 95) + 1
  L_thres_99 <- sum(cumulated_pve < 99) + 1
  
  K <- unique(mfvb_res$K) 
  
  if (generate_from_univ) {

    rmse_mfpca <- matrix(NA, L_sim, p)
    for (j in 1:p) { 
      rmse_mfpca[,j] <- apply(Zeta[[j]] - mfvb_res$Zeta[,1:L_sim, drop = F], 2, function(x) sqrt(mean(x^2))) # one error per FPC
    }
    rownames(rmse_mfpca) <- paste0("FPC_", 1:L_sim)
    colnames(rmse_mfpca) <- paste0("Variable_", 1:p)
    
  } else {
    
    rmse_mfpca <- apply(Zeta - mfvb_res$Zeta[,1:L_sim, drop = F], 2, function(x) sqrt(mean(x^2))) # one error per FPC (vector of length L)
    names(rmse_mfpca) <- paste0("FPC_", 1:L_sim)
    
  }

  list_time_obs_univ <- lapply(1:p, function(j) lapply(time_obs, function(time_obs_i) time_obs_i[j]))
  list_Y_univ <- lapply(1:p, function(j) lapply(Y, function(Y_i) Y_i[j]))

  list_mfvb_res_univ <- lapply(seq_along(list_Y_univ), function(j) {

    if (model_choice_K) {
      run_mfvb_fpca_model_choice(list_time_obs_univ[[j]], list_Y_univ[[j]], L = L,
                                 n_g = NULL, time_g = time_g, tol = tol, 
                                 maxit = maxit_mfvb, Psi_g = Psi_g[j], n_cpus = 1)
      
    } else {
      run_mfvb_fpca(list_time_obs_univ[[j]], list_Y_univ[[j]], L = L, K = K, 
                    n_g = NULL, time_g = time_g, tol = tol, maxit = maxit_mfvb,
                    plot_elbo = FALSE, Psi_g = Psi_g[j])
    }
    
  }
  )
  names(list_mfvb_res_univ) <- paste0("run_univ_Y", 1:p)

  cumulated_pve_univ <- sapply(list_mfvb_res_univ, "[[", "cumulated_pve")
  L_thres_95_univ <- apply(cumulated_pve_univ, 2, function(cum_pve_j) sum(cum_pve_j < 95) + 1)
  L_thres_99_univ <- apply(cumulated_pve_univ, 2, function(cum_pve_j) sum(cum_pve_j < 99) + 1)
  
  K_univ <- sapply(list_mfvb_res_univ, "[[", "K")
  
  list_Zeta_univ <- lapply(1:p, function(j)  sqrt(p)*list_mfvb_res_univ[[j]]$Zeta_hat) 

  if (generate_from_univ) {

    rmse_fpca <- matrix(NA, L_sim, p)
    for (j in 1:p) { 
      rmse_fpca[,j] <- apply(Zeta[[j]] - list_Zeta_univ[[j]][, 1:L_sim, drop = F], 2, function(x) sqrt(mean(x^2))) # one error per FPC
    }
   
  } else {

    rmse_fpca <- matrix(NA, L_sim, p)
    for (j in 1:p) {
      rmse_fpca[,j] <- apply(Zeta - list_Zeta_univ[[j]][, 1:L_sim, drop = F], 2, function(x) sqrt(mean(x^2))) # one error per FPC
      
    }

  }
  rownames(rmse_fpca) <- paste0("FPC_", 1:L_sim)
  colnames(rmse_fpca) <- paste0("Variable_", 1:p)

  mu_hat_univ <- Reduce(cbind, lapply(1:p, function(j) list_mfvb_res_univ[[j]]$mu_hat))
  list_Psi_hat_univ <- lapply(1:L, function(l) Reduce(cbind, lapply(1:p, function(j) list_mfvb_res_univ[[j]]$list_Psi_hat[[l]]/sqrt(p))))  # adjust as not orthonormal in the multivariate space

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

  create_named_list(mfvb_res, K,
                    cumulated_pve, L_thres_95, L_thres_99,
                    list_mfvb_res_univ, K_univ,
                    cumulated_pve_univ, L_thres_95_univ, L_thres_99_univ,
                    rmse_mfpca, rmse_fpca,
                    ise_mfpca, ise_fpca)
  

}, mc.cores = n_cpus)
names(list_res) <- paste0("repl_", 1:n_repl)

list_rmse_mfpca <- lapply(list_res, "[[", "rmse_mfpca")
list_rmse_fpca <- lapply(list_res, "[[", "rmse_fpca")

list_ise_mfpca <- lapply(list_res, "[[", "ise_mfpca")
list_ise_fpca <- lapply(list_res, "[[", "ise_fpca")

list_cumulated_pve <- lapply(list_res, "[[", "cumulated_pve")
list_cumulated_pve_univ <- lapply(list_res, "[[", "cumulated_pve_univ")

list_L_thres_95 <- lapply(list_res, "[[", "L_thres_95")
list_L_thres_95_univ <- lapply(list_res, "[[", "L_thres_95_univ")

list_L_thres_99 <- lapply(list_res, "[[", "L_thres_99")
list_L_thres_99_univ <- lapply(list_res, "[[", "L_thres_99_univ")

list_K <- lapply(list_res, "[[", "K")
list_K_univ <- lapply(list_res, "[[", "K_univ")

estimated_L_thres_95 <- sapply(1:p, function(k) sapply(list_L_thres_95_univ, "[[", paste0("run_univ_Y", k)))
estimated_L_thres_95 <- cbind(estimated_L_thres_95, unlist(list_L_thres_95))
colnames(estimated_L_thres_95) <- c(paste0("FPCA Y", 1:p),  "mFPCA")

estimated_L_thres_99 <- sapply(1:p, function(k) sapply(list_L_thres_99_univ, "[[", paste0("run_univ_Y", k)))
estimated_L_thres_99 <- cbind(estimated_L_thres_99, unlist(list_L_thres_99))
colnames(estimated_L_thres_99) <- c(paste0("FPCA Y", 1:p),  "mFPCA")

all_K <- sapply(1:p, function(k) sapply(list_K_univ, "[[", paste0("run_univ_Y", k)))
all_K <- cbind(all_K, unlist(list_K))
colnames(all_K) <- c(paste0("FPCA Y", 1:p),  "mFPCA")


if (bool_save) {
  save(list_rmse_fpca, list_rmse_mfpca, list_ise_fpca, list_ise_mfpca,
       list_cumulated_pve, list_cumulated_pve_univ,
       list_L_thres_95, list_L_thres_95_univ,
       list_L_thres_99, list_L_thres_99_univ,
       list_K, list_K_univ,
       file = file.path(res_dir, "output.RData"))
}

if (generate_from_univ) {

  list_mat_rmse_FPC <- NULL
  for (l in 1:L_sim) {
    mat_rmse_FPC <- cbind(t(sapply(list_rmse_mfpca, function(ll) ll[paste0("FPC_", l),])),
                               t(sapply(list_rmse_fpca, function(ll) ll[paste0("FPC_", l),])))
    
    colnames(mat_rmse_FPC) <- c(paste0("mFPCA \n variable ", 1:p),
                                paste0("FPCA \n variable ", 1:p))
    
    list_mat_rmse_FPC <- append(list_mat_rmse_FPC, list(mat_rmse_FPC))
  }
  names(list_mat_rmse_FPC) <- paste0("FPC_", 1:L_sim)

} else {
  
  list_mat_rmse_FPC <- NULL
  for (l in 1:L_sim) {
    mat_rmse_FPC <-cbind(sapply(list_rmse_mfpca, function(ll) ll[paste0("FPC_", l)]),
                           t(sapply(list_rmse_fpca, function(ll) ll[paste0("FPC_", l),])))
    
    colnames(mat_rmse_FPC) <- c("mFPCA\n",
                                paste0("FPCA \n variable ", 1:p))
    rownames(mat_rmse_FPC) <- paste0("repl_", 1:n_repl)
    
    list_mat_rmse_FPC <- append(list_mat_rmse_FPC, list(mat_rmse_FPC))
  }
  names(list_mat_rmse_FPC) <- paste0("FPC_", 1:L_sim)
  
}


mat_ise_mu <- cbind(t(sapply(list_ise_mfpca, function(ll) ll["mu",])),
                    t(sapply(list_ise_fpca, function(ll) ll["mu",])))

colnames(mat_ise_mu) <- c(paste0("mFPCA \n variable ", 1:p),
                          paste0("FPCA \n variable ", 1:p))


list_mat_ise_Psi <- NULL
for (l in 1:L_sim) {
  mat_ise_Psi <- cbind(t(sapply(list_ise_mfpca, function(ll) ll[paste0("Psi_", l),])),
                       t(sapply(list_ise_fpca, function(ll) ll[paste0("Psi_", l),])))
  
  colnames(mat_ise_Psi) <- c(paste0("mFPCA \n variable ", 1:p),
                             paste0("FPCA \n variable ", 1:p))
  
  list_mat_ise_Psi <- append(list_mat_ise_Psi, list(mat_ise_Psi))
}
names(list_mat_rmse_FPC) <- paste0("FPC_", 1:L_sim)


for (l in 1:L_sim) {
  if (bool_save) {
    pdf(paste0(res_dir, "/rmse_scores_FPC_", l, ".pdf"), width = 7.5, height = 5.5, paper='special')
  }
  boxplot(list_mat_rmse_FPC[[l]], main = paste0("Root-mean-square error for estimated scores FPC ", l),
          sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
          xaxt = "n")
  axis(side = 1, at = 1:ncol(list_mat_rmse_FPC[[l]]), labels = colnames(list_mat_rmse_FPC[[l]]), tick = FALSE)
  if (bool_save) {
    dev.off()
  }
}

if (bool_save) {
  pdf(paste0(res_dir, "/ise_mu.pdf"), width = 7.5, height = 5.5, paper='special')
}
boxplot(log(mat_ise_mu),
        main = "Integrated square error for estimated mean functions",
        sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
        ylab = "log ISE",
        xaxt = "n")
axis(side = 1, at = 1:ncol(mat_ise_mu), labels = colnames(mat_ise_mu), tick = FALSE)
if (bool_save) {
  dev.off()
}

for (l in 1:L_sim) {
  if (bool_save) {
    pdf(paste0(res_dir, "/ise_Psi_", l, ".pdf"), width = 7.5, height = 5.5,paper='special')
  }
  boxplot(log(list_mat_ise_Psi[[l]]),
          main = paste0("Integrated square error for estimated eigenfunctions FPC ", l),
          sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
          ylab = "log ISE",
          xaxt = "n")
  axis(side = 1, at = 1:ncol(list_mat_ise_Psi[[l]]), labels = colnames(list_mat_ise_Psi[[l]]), tick = FALSE)
  if (bool_save) {
    dev.off()
  }
}


if (generate_from_univ) {
  
  list_mat_rmse_avg_FPC <- NULL
  for (l in 1:L_sim) {
    mat_rmse_avg_FPC <- cbind(colMeans(sapply(list_rmse_mfpca, function(ll) ll[paste0("FPC_", l),])),
                              colMeans(sapply(list_rmse_fpca, function(ll) ll[paste0("FPC_", l),])))
    colnames(mat_rmse_avg_FPC) <- c("mFPCA", "FPCA")
    
    list_mat_rmse_avg_FPC <- append(list_mat_rmse_avg_FPC, list(mat_rmse_avg_FPC))
  }
  names(list_mat_rmse_avg_FPC) <- paste0("FPC_", 1:L_sim)
  
} else {
  
  list_mat_rmse_avg_FPC <- NULL
  for (l in 1:L_sim) {
    mat_rmse_avg_FPC <- cbind(sapply(list_rmse_mfpca, function(ll) ll[paste0("FPC_", l)]),
                              colMeans(sapply(list_rmse_fpca, function(ll) ll[paste0("FPC_", l),])))
    colnames(mat_rmse_avg_FPC) <- c("mFPCA", "FPCA")
    
    list_mat_rmse_avg_FPC <- append(list_mat_rmse_avg_FPC, list(mat_rmse_avg_FPC))
  }
  names(list_mat_rmse_avg_FPC) <- paste0("FPC_", 1:L_sim)
  
}

mat_ise_mu_avg <- cbind(colMeans(sapply(list_ise_mfpca, function(ll) ll["mu",])),
                        colMeans(sapply(list_ise_fpca, function(ll) ll["mu",])))
colnames(mat_ise_mu_avg ) <- c("mFPCA", "FPCA")

list_mat_ise_Psi_avg <- NULL
for (l in 1:L_sim) {
  mat_ise_Psi_avg <- cbind(colMeans(sapply(list_ise_mfpca, function(ll) ll[paste0("Psi_", l),])),
                             colMeans(sapply(list_ise_fpca, function(ll) ll[paste0("Psi_", l),])))
  colnames(mat_ise_Psi_avg) <- c("mFPCA", "FPCA")
  
  list_mat_ise_Psi_avg <- append(list_mat_ise_Psi_avg, list(mat_ise_Psi_avg))
}
names(list_mat_rmse_avg_FPC) <- paste0("FPC_", 1:L_sim)


for (l in 1:L_sim) {
  if (bool_save) {
    pdf(paste0(res_dir, "/rmse_scores_avg_FPC_", l, ".pdf"), width = 7.5, height = 5.5, paper='special')
  }
  boxplot(list_mat_rmse_avg_FPC[[l]], main = paste0("Root-mean-square error for estimated scores FPC ", l),
          sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
          xaxt = "n")
  axis(side = 1, at = 1:ncol(list_mat_rmse_avg_FPC[[l]]), 
       labels = colnames(list_mat_rmse_avg_FPC[[l]]), tick = FALSE)
  if (bool_save) {
    dev.off()
  }
}
if (bool_save) {
  pdf(paste0(res_dir, "/ise_mu_avg.pdf"), width = 7.5, height = 5.5, paper='special')
}
boxplot(log(mat_ise_mu_avg),
        main = "Integrated square error for estimated mean functions",
        sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
        ylab = "log ISE",
        xaxt = "n")
axis(side = 1, at = 1:ncol(mat_ise_mu_avg), labels = colnames(mat_ise_mu_avg), tick = FALSE)
if (bool_save) {
  dev.off()
}

for (l in 1:L_sim) {
  if (bool_save) {
    pdf(paste0(res_dir, "/ise_Psi_", l, "_avg.pdf"), width = 7.5, height = 5.5,paper='special')
  }
  boxplot(log(list_mat_ise_Psi_avg[[l]]),
          main = paste0("Integrated square error for estimated eigenfunctions FPC ", l),
          sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
          ylab = "log ISE",
          xaxt = "n")
  axis(side = 1, at = 1:ncol(list_mat_ise_Psi_avg[[l]]), 
       labels = colnames(list_mat_ise_Psi_avg[[l]]), tick = FALSE)
  if (bool_save) {
    dev.off()
  }
}

if (bool_save) {
  pdf(paste0(res_dir, "/estimated_L_thres_95.pdf"),
      width = 7.8, height = 5.8, paper='special')
}
boxplot(estimated_L_thres_95, main = "Learnt number of eigenfunctions", 
        ylab = "Learnt L (PVE thres 95%)") 
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/estimated_L_thres_99.pdf"),
      width = 7.8, height = 5.8, paper='special')
}
boxplot(estimated_L_thres_99, main = "Learnt number of eigenfunctions", 
        ylab = "Learnt L (PVE thres 99%)") 
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/all_K.pdf"),
      width = 7.8, height = 5.8, paper='special')
}
boxplot(all_K, main = ifelse(model_choice_K, "Learnt number of splines", 
                             "Number of splines (rule of thumb)"), 
        ylab = "K") 
if (bool_save) {
  dev.off()
}


if (bool_save) {
  save(list_rmse_fpca, list_rmse_mfpca, 
       list_mat_rmse_FPC, list_mat_rmse_avg_FPC, 
       list_ise_fpca, list_ise_mfpca,
       mat_ise_mu, list_mat_ise_Psi, 
       mat_ise_mu_avg, list_mat_ise_Psi_avg, 
       list_cumulated_pve, list_cumulated_pve_univ,
       list_L_thres_95, list_L_thres_95_univ,
       list_L_thres_99, list_L_thres_99_univ,
       list_K, list_K_univ,
       file = file.path(res_dir, "output.RData"))
}

