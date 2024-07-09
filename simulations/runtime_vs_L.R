rm(list = ls())

bool_cluster <- T
if (bool_cluster) {
  CORE_DIR_MRC_BSU <- Sys.getenv("CORE_DIR_MRC_BSU")
  out_dir <- file.path(CORE_DIR_MRC_BSU, "mFPCA_output/")
} else {
  CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")
  out_dir <- file.path(CORE_DIR_ICLOUD, "mFPCA_output/")
}

require(bayesFPCA)
require(ggplot2)
require(parallel)

seed <- 1
set.seed(seed)

p <- 3                                        # number of responses
N <- 100                                      # number of subjects
N_t_min <- 5                                  # minimum number of time observations for each curve
N_t_max <- 10                                 # maximum number of time observations for each curve
n <- matrix(sample(N_t_min:N_t_max, N*p,
                   replace = TRUE), N, p)     # number of time observations

model_choice_K <- F
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
tol <- 1e-5                                   # convergence criterion
maxit <- 1000                                 # maximum number of iterations

sigma_zeta_vec <- 1/(1:L_sim)                 # sd for first and second scores
sigma_eps <- rep(1, p)                        # sd of the residuals

n_g <- 1000                                   # length of the plotting grid

mu_func <- function(time_obs, j) {
  ans <- ((-1)^j)*3*sin(pi*j*time_obs) - 3/2
  return(ans)
}

psi_1 <- function(time_obs, j) {
  ans <- sqrt(2/p)*sin(2*pi*j*time_obs)
  return(ans)
}

psi_2 <- function(time_obs, j) {
  ans <- sqrt(2/p)*cos(2*pi*j*time_obs)
  return(ans)
}

Psi_func <- function(time_obs, j, p) {
  ans <- cbind(psi_1(time_obs, j), psi_2(time_obs, j))
  return(ans)
}

if (L_sim > 2) {
  stop("Data generation not implemented for L_sim > 2 with the above Psi functions.")
}

n_repl <- 100
vec_L <- 2:10

n_cpus <- 50

bool_save <- TRUE
if (bool_save) {
  
  res_dir <- paste0(out_dir, "/runtime_vs_L_p_", p, ifelse(model_choice_K,
                                                           "_model_choice_K",
                                                           paste0("K_", K)), 
                    "_N_", N, "_N_t_min_", N_t_min, "_N_t_max_", N_t_max, 
                    "_L_sim_", L_sim, "_tol_", tol, "_maxit_", maxit, 
                    "_sigma_zeta_vec_", paste(sigma_zeta_vec, collapse = "_"), 
                    "_sigma_eps_", paste(sigma_eps, collapse = "_"), 
                    "_n_repl_", n_repl, "_vec_L_", paste(vec_L, collapse = "_"), 
                    "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}


run_fpca_replicate <- function(rep) {
  set.seed(seed + rep) 
  mfpca_data <- generate_fpca_data(N, p, n, L_sim, n_g, sigma_eps, mu_func, 
                                   Psi_func, vec_sd_zeta = sigma_zeta_vec)
  time_obs <- mfpca_data$time_obs
  Y <- mfpca_data$Y
  Psi_g <- mfpca_data$Psi_g
  
  replicate_results <- list()
  
  for (L in vec_L) {
    mfvb_runtime <- system.time({
      mfpca_res <- run_mfvb_fpca(time_obs, Y, L, K = K, n_g = n_g, tol = tol, 
                                 maxit = maxit, Psi_g = Psi_g)
    })["elapsed"]
    
    vmp_runtime <- system.time({
      vmp_res <- run_vmp_fpca(time_obs, Y, L, K = K, n_g = n_g, tol = tol, 
                              maxit = maxit, Psi_g = Psi_g)
    })["elapsed"]
    
    replicate_results[[L]] <- c(mfvb_runtime, vmp_runtime)
  }
  
  return(replicate_results)
}

runtime_results <- mclapply(1:n_repl, run_fpca_replicate, mc.cores = n_cpus)

runtime_df <- do.call(rbind, lapply(1:n_repl, function(rep) {
  do.call(rbind, lapply(vec_L, function(L) {
    data.frame(Replicate = rep, L = L, Method = c("MFVB", "VMP"), 
               Runtime = runtime_results[[rep]][[L]])
  }))
}))

if (bool_save) {
  pdf(paste0(res_dir, "/runtime_vs_L.pdf"), width = 6.7, height = 5)
}
ggplot(runtime_df, aes(x = factor(L), y = Runtime, fill = Method)) +
  geom_boxplot() +
  labs(title = "Runtime profiling wrt number of inferred components L",
       x = "L",
       y = "Runtime (seconds)") +
  theme_classic() +
  scale_fill_manual(values = c("MFVB" = "darkseagreen3", "VMP" = "dodgerblue1"))

if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/log_runtime.pdf"), width =5.8, height = 3.6, paper='special') 
}
pl <- ggplot(runtime_df,  aes(x = factor(L), y = Runtime, fill = Method)) +
  geom_boxplot() +
  labs(title = "Runtime profiling wrt number of inferred components L",
       x = "L",
       y = "Runtime (seconds, log-scale)") +
  scale_y_continuous(trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  ggtitle(paste0("Runtime as a function of L (p = ", p, ", n = ", N, ", n_i = ", 
                 (N_t_max - N_t_min)/2 +N_t_min, ")")) + 
  geom_boxplot(position = position_dodge2(preserve = "single", padding = 0.2))

pl+scale_fill_manual(values = c("MFVB" = "darkseagreen3", "VMP" = "dodgerblue1")) + 
  theme_classic() + geom_boxplot(lwd=0.3, fatten = 2, outlier.size = 0.5)+
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line.x = element_line(color = "black", size = 0.1),
    axis.line.y = element_line(color = "black", size = 0.1)
  )

if (bool_save) {
  dev.off()
  save.image(file = file.path(res_dir, "output.RData"))
}
