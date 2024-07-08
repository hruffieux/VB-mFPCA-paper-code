rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

bool_cluster <- TRUE
if (bool_cluster) {

  #!/usr/bin/env Rscript
  args = commandArgs(trailingOnly=TRUE)
  job_id <- as.integer(args[1]) # 1 to 9

  CORE_DIR_MRC_BSU <- Sys.getenv("CORE_DIR_MRC_BSU")
  out_dir <- file.path(CORE_DIR_MRC_BSU, "mFPCA_output/")
} else {

  job_id <- 1 

  CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")
  out_dir <- file.path(CORE_DIR_ICLOUD, "mFPCA_output/")
}

main_dir <- file.path(CORE_DIR, "bayesian-mFPCA-paper-code/simulations/")
setwd(main_dir)

require(bayesFPCA)
require(MFPCA)

source("fun_utils.R")

seed <- 123
set.seed(seed)

p <- 3
N <- 100

bool_mfpca_only <- FALSE
if (bool_mfpca_only) {
  N_t_min <- seq(85, 245, by = 20)[job_id]
  N_t_max <- seq(115, 275, by = 20)[job_id]
} else {
  N_t_min <- seq(5, 225, by = 20)[job_id]
  N_t_max <- seq(35, 255, by = 20)[job_id]
}
n <- matrix(sample(N_t_min:N_t_max, N*p,
                   replace = TRUE), N, p)      # number of time observations

model_choice_K <- T
if (model_choice_K) {
  K <- NULL
} else {
  bool_K_Ruppert <- T
  if (bool_K_Ruppert) {
    K <- NULL                                  # default based on the rule in Ruppert 2002
  } else {
    n_int_knots <- 5                           # number of interior knots
    K <- n_int_knots + 2                       # number of spline basis functions
  }
}

L_sim <- 2
L <- 10                                        # number of FPCA basis functions

tol  <- 1e-5                                   # convergence tolerance 
maxit <- 1000                                  # maximum number of mfvb iterations

n_g <- 1000                                    # length of the plotting grid

vec_sd_zeta <- 1/(1:L_sim)
vec_sd_eps <- rep(1, p)                        # sd of the residuals (here same for all variables)

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

n_repl <- 200
n_cpus <- 25

bool_save <- T
if (bool_save) {
  res_dir <- paste0(out_dir, "/comp_Happ_", 
                    ifelse(bool_mfpca_only, "mfpca_only_", ""), 
                    "n_repl_", n_repl, "_", 
                           ifelse(generate_from_univ,
                                  paste0("gen_from_univ_corr_",
                                         paste0(format(rho_Zeta, digits = 2),
                                                collapse = "-"), "_"), ""),
                    "p_", p, "_N_", N, "_Nt_min_",
                    N_t_min, "_max_", N_t_max, "_L_sim_", L_sim, "_L_", L,
                    ifelse(model_choice_K, "_model_choice_K", paste0("K_", K)),
                    "_sd_zeta_", paste0(format(vec_sd_zeta , digits = 2), collapse = "-"),
                    "_sd_eps_", unique(vec_sd_eps), "_tol_", tol, 
                    "_maxit_", maxit, "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}

if (rho_Zeta == 1) { 
  generate_from_univ <- FALSE
  vec_rho_Zeta <- NULL
}

list_res <- parallel::mclapply(1:n_repl, function(repl) {

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


  # Set up the data:
  #
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


  if (!bool_mfpca_only) {
    
    list_type_j <- list(type = "uFPCA")
    list_type <- rep(list(list_type_j), p)
    bool_ci <- F 
  
    if(bool_ci) {
      nb_boot <- 100
    } else {
      nb_boot <- NULL
    }
    
   # set.seed(seed)
    
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
    sapply(list_time_g_happ, length) 
    
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
    
    # flip sign
    # for (j in 1:p) { # we can base the sign on the first variable j = 1 and apply the corresponding swap to all others as same sign for a given l
    jj <- 1
    Psi_g_happ_j <- f_Psi(list_time_g_happ[[jj]], j = jj, p = p)  # not the same grid for Happ!
    
    for (l in 1:L_sim) {
      cprod_sign <- sign(cprod(list_Psi_hat_happ[[l]][[jj]], Psi_g_happ_j[,l]))
      
      for (j in 1:p) {
        list_Psi_hat_happ[[l]][[j]] <- cprod_sign*list_Psi_hat_happ[[l]][[j]]
      }
      
      Zeta_hat_happ[,l] <- cprod_sign*Zeta_hat_happ[,l]
    }
    # }
    
    
  } else {
    happ_res <- run_time_happ <- NULL
  }
  
  # set.seed(seed)

  if (model_choice_K) {
    run_time_mfvb <- system.time({
      mfvb_res <- run_mfvb_fpca_model_choice(time_obs, Y, L = L, n_g = NULL, time_g = time_g, 
                                             tol = tol, maxit = maxit, Psi_g = Psi_g, n_cpus = 1)
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
  cumulated_pve <- mfvb_res$cumulated_pve
  L_thres_95 <- sum(cumulated_pve < 95) + 1
  L_thres_99 <- sum(cumulated_pve < 99) + 1


  if (generate_from_univ) {
    
    rmse_mfpca <- rmse_happ <- matrix(NA, L_sim, p)
    for (j in 1:p) { 
      rmse_mfpca[,j] <- apply(Zeta[[j]] - Zeta_hat[,1:L_sim], 2, function(x) sqrt(mean(x^2))) 
      if (!bool_mfpca_only) rmse_happ[,j] <- apply(Zeta[[j]] - Zeta_hat_happ, 2, function(x) sqrt(mean(x^2)))
    }
    rownames(rmse_mfpca) <- rownames(rmse_happ) <- paste0("FPC_", 1:L_sim)
    colnames(rmse_mfpca) <- colnames(rmse_happ) <-paste0("Variable_", 1:p)
  } else {
    rmse_mfpca <- apply(Zeta - Zeta_hat[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
    rmse_mfpca_reswap <- apply(Zeta + Zeta_hat[,1:L_sim], 2, function(x) sqrt(mean(x^2)))
    for (l in 1:L_sim) {
      if (rmse_mfpca_reswap[l] < rmse_mfpca[l]) {
        warning("Wrong score / eigenfunction sign for mFPCA still likely, re-swap.")
        Zeta_hat[, l] <- -Zeta_hat[,l]
        list_Psi_hat[[l]] <- -list_Psi_hat[[l]]
      }
    }
    if (!bool_mfpca_only) {
      rmse_happ <- apply(Zeta - Zeta_hat_happ, 2, function(x) sqrt(mean(x^2)))
      rmse_happ_reswap <- apply(Zeta + Zeta_hat_happ, 2, function(x) sqrt(mean(x^2)))
      for (l in 1:L_sim) {
        if (rmse_happ_reswap[l] < rmse_happ[l]) {
          warning("Wrong score / eigenfunction sign for Happ still likely, re-swap.")
          Zeta_hat_happ[, l] <- -Zeta_hat_happ[,l]
          for (j in 1:p) {
            list_Psi_hat_happ[[l]][[j]] <- -list_Psi_hat_happ[[l]][[j]]
          }
        }
      }
    } else {
      rmse_happ <- rep(NA, L_sim)
    }
    names(rmse_happ) <- paste0("FPC_", 1:L_sim)
  }
  
  ise_mfpca <- ise_happ <- matrix(NA, L_sim+1, p)

  for (l in 1:(L_sim+1)) {

    for(j in 1:p) {

      if (!bool_mfpca_only) mu_g_happ_j <- f_mu(list_time_g_happ[[j]], j = j)
      if (!bool_mfpca_only) Psi_g_happ_j <- f_Psi(list_time_g_happ[[j]], j = j, p = p)

      if (l == 1) {
        ise_mfpca[l, j] <- trapint(time_g, (mu_hat[,j] - mu_g[[j]])^2)
        if (!bool_mfpca_only) ise_happ[l, j] <- trapint(list_time_g_happ[[j]], (mu_hat_happ[[j]] - mu_g_happ_j)^2) 
      } else {
        ise_mfpca[l, j] <- trapint(time_g, (list_Psi_hat[[l-1]][,j] - Psi_g[[j]][,l-1])^2)
        if (!bool_mfpca_only) ise_happ[l, j] <- trapint(list_time_g_happ[[j]], (list_Psi_hat_happ[[l-1]][[j]] - Psi_g_happ_j[,l-1])^2) 
      }

    }

  }
  rownames(ise_mfpca) <- rownames(ise_happ) <- c("mu", paste0("Psi_", 1:L_sim))
  colnames(ise_mfpca) <- colnames(ise_happ) <- paste0("Variable_", 1:p)

  create_named_list(mfvb_res, run_time_mfvb,
                    happ_res, run_time_happ,
                    rmse_mfpca, rmse_happ,
                    ise_mfpca, ise_happ,
                    L_thres_95, L_thres_99)

}, mc.cores = n_cpus)
names(list_res) <- paste0("repl_", 1:n_repl)


list_runtime_mfpca <- lapply(list_res, "[[", "run_time_mfvb")
list_runtime_happ <- lapply(list_res, "[[", "run_time_happ")

list_rmse_mfpca <- lapply(list_res, "[[", "rmse_mfpca")
list_rmse_happ <- lapply(list_res, "[[", "rmse_happ")

list_ise_mfpca <- lapply(list_res, "[[", "ise_mfpca")
list_ise_happ <- lapply(list_res, "[[", "ise_happ")

vec_L_thres_95 <- sapply(list_res, "[[", "L_thres_95")
vec_L_thres_99 <- sapply(list_res, "[[", "L_thres_99")


if (bool_save) {
  save(list_rmse_happ, list_rmse_mfpca, list_ise_happ, list_ise_mfpca,
       vec_L_thres_95, vec_L_thres_99,
       file = file.path(res_dir, "output.RData"))
}


if (generate_from_univ) {

  list_mat_rmse_FPC <- vector("list", L_sim)
  for (l in 1:L_sim) {
    list_mat_rmse_FPC[[l]] <- cbind(t(sapply(list_rmse_mfpca, function(ll) ll[paste0("FPC_", l),])),
                                    t(sapply(list_rmse_happ, function(ll) ll[paste0("FPC_", l),])))
    
    colnames(list_mat_rmse_FPC[[l]]) <- c(paste0("mFPCA \n variable ", 1:p),
                                          paste0("Happ \n variable ", 1:p))
    rownames(list_mat_rmse_FPC[[l]]) <- paste0("repl_", 1:n_repl)
  }

} else {

  list_mat_rmse_FPC <- vector("list", L_sim)
  for (l in 1:L_sim) {
    list_mat_rmse_FPC[[l]] <- cbind(sapply(list_rmse_mfpca, function(ll) ll[paste0("FPC_", l)]),
                                    sapply(list_rmse_happ, function(ll) ll[paste0("FPC_", l)]))
    
    colnames(list_mat_rmse_FPC[[l]]) <- c("mFPCA\n",
                                          paste0("Happ \n"))
    rownames(list_mat_rmse_FPC[[l]]) <- paste0("repl_", 1:n_repl)
  }

}

mat_runtime <-cbind(sapply(list_runtime_mfpca, function(ll) ll["elapsed"]),
                       sapply(list_runtime_happ, function(ll) ll["elapsed"]))


if (bool_mfpca_only) {
  mat_runtime_tmp <- matrix(NA, nrow = n_repl, ncol = 2)
  mat_runtime_tmp[,1] <- unlist(mat_runtime[,1])
  mat_runtime_tmp[,2] <- NA
  mat_runtime <- mat_runtime_tmp
  rm(mat_runtime_tmp)
}

colnames(mat_runtime) <- c("mFPCA\n", paste0("Happ \n"))
rownames(mat_runtime) <-  paste0("repl_", 1:n_repl)

mat_ise_mu <- cbind(t(sapply(list_ise_mfpca, function(ll) ll["mu",])),
                    t(sapply(list_ise_happ, function(ll) ll["mu",])))

list_mat_ise_Psi <- vector("list", L_sim)
for (l in 1:L_sim) {
  list_mat_ise_Psi[[l]] <- cbind(t(sapply(list_ise_mfpca, function(ll) ll[paste0("Psi_", l),])),
                                 t(sapply(list_ise_happ, function(ll) ll[paste0("Psi_", l),])))
  
  colnames(list_mat_ise_Psi[[l]]) <- c(paste0("mFPCA \n variable ", 1:p),
                                       paste0("Happ \n variable ", 1:p))
}

colnames(mat_ise_mu) <- c(paste0("mFPCA \n variable ", 1:p),
                          paste0("Happ \n variable ", 1:p))

if (bool_save) {
  pdf(paste0(res_dir, "/runtime.pdf"), width = 7.5, height = 5.5, paper='special')
}
boxplot(mat_runtime, main = "Runtime", xaxt = "n")
axis(side = 1, at = 1:ncol(mat_runtime), labels = colnames(mat_runtime), tick = FALSE)
if (bool_save) {
  dev.off()
}

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
          main = paste0("Integrated square error for estimated eigenfunctions ", l),
          sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
          ylab = "log ISE",
          xaxt = "n")
  axis(side = 1, at = 1:ncol(list_mat_ise_Psi[[l]]), labels = colnames(list_mat_ise_Psi[[l]]), tick = FALSE)
  if (bool_save) {
    dev.off()
  }
}
if (bool_save) {
  pdf(paste0(res_dir, "/mat_L_thres_95.pdf"), width = 7.5, height = 5.5, paper='special')
}
boxplot(vec_L_thres_95,
        main = "Learnt number of eigenfunction (PVE > 95%)",
        ylab = "Learnt L",
        xaxt = "n")
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/mat_L_thres_99.pdf"), width = 7.5, height = 5.5, paper='special')
}
boxplot(vec_L_thres_99,
        main = "Learnt number of eigenfunction (PVE > 99%)",
        ylab = "Learnt L",
        xaxt = "n")
if (bool_save) {
  dev.off()
}

if (generate_from_univ) {

  list_mat_rmse_avg_FPC <- vector("list", L_sim)
  for (l in 1:L_sim) {
    list_mat_rmse_avg_FPC[[l]] <- cbind(colMeans(sapply(list_rmse_mfpca, function(ll) ll[paste0("FPC_", l),])),
                                        colMeans(sapply(list_rmse_happ, function(ll) ll[paste0("FPC_", l),])))

    if (bool_save) {
      pdf(paste0(res_dir, "/rmse_scores_avg_FPC_", l, ".pdf"), width = 7.5, height = 5.5, paper='special')
    }
    boxplot(list_mat_rmse_avg_FPC[[l]], main = paste0("Root-mean-square error for estimated scores FPC ", l),
            sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
            xaxt = "n")
    axis(side = 1, at = 1:ncol(list_mat_rmse_avg_FPC[[l]]), labels = colnames(list_mat_rmse_avg_FPC[[l]]), tick = FALSE)
    if (bool_save) {
      dev.off()
    }
  }
} else {
  list_mat_rmse_avg_FPC <- list_mat_rmse_FPC
}

mat_ise_mu_avg <- cbind(colMeans(sapply(list_ise_mfpca, function(ll) ll["mu",])),
                        colMeans(sapply(list_ise_happ, function(ll) ll["mu",])))

colnames(mat_ise_mu_avg) <- c("mFPCA", "Happ")
  
list_mat_ise_Psi_avg <- vector("list", L_sim)
for (l in 1:L_sim) {
  list_mat_ise_Psi_avg[[l]] <- cbind(colMeans(sapply(list_ise_mfpca, function(ll) ll[paste0("Psi_", l),])),
                                     colMeans(sapply(list_ise_happ, function(ll) ll[paste0("Psi_", l),])))

  colnames(list_mat_ise_Psi_avg[[l]]) <- c("mFPCA", "Happ")
  colnames(list_mat_rmse_avg_FPC[[l]]) <- c("mFPCA", "Happ")
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
          main = paste0("Integrated square error for estimated eigenfunctions ", l),
          sub = paste0("Score correlation: ", unique(vec_rho_Zeta)),
          ylab = "log ISE",
          xaxt = "n")
  axis(side = 1, at = 1:ncol(list_mat_ise_Psi_avg[[l]]), labels = colnames(list_mat_ise_Psi_avg[[l]]), tick = FALSE)
  if (bool_save) {
    dev.off()
  }
}

if (bool_save) {
  save(list_rmse_happ, list_rmse_mfpca,
       list_mat_rmse_FPC, 
       list_mat_rmse_avg_FPC, 
       list_ise_happ, list_ise_mfpca,
       mat_ise_mu, list_mat_ise_Psi, 
       mat_ise_mu_avg, list_mat_ise_Psi_avg, 
       list_runtime_happ, list_runtime_mfpca, mat_runtime,
       vec_L_thres_95, vec_L_thres_99,
       file = file.path(res_dir, "output.RData"))
}

