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
  
  job_id <- 5 # 1 to 5
  
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
  bool_K_Ruppert <- T
  if (bool_K_Ruppert) {
    K <- NULL                                   # default based on the rule in Ruppert 2002
  } else {
    n_int_knots <- 5                            # number of interior knots
    K <- n_int_knots + 2                        # number of spline basis functions
  }
}

L_sim <- 2
L <- 10                                         # number of FPCA basis functions

tol  <- 1e-5                                    # convergence tolerance 
maxit <- 1000                                   # maximum number of vmp iterations

n_g <- 1000                                     # length of the plotting grid

ind_exp <- 1
exponent_sd_zeta <- c(-1, -1/2, -1/4, -1/8, -1/16)[ind_exp]
vec_sd_zeta <- (1:L_sim)^exponent_sd_zeta
sd_eps <- 2.5
vec_sd_eps <- rep(sd_eps, p)                    # sd of the residuals (here same for all variables)

generate_from_univ <- T

if (generate_from_univ) {
  id_strength <- 6 
  rho_Zeta <- seq(0, 1, by = 0.2)[id_strength]
  vec_rho_Zeta <- rep(rho_Zeta, L_sim)
  cat("Correlation of simulated scores across variables for eigenfunctions 1 to L_sim: ",
      vec_rho_Zeta)
} else {
  vec_rho_Zeta <- NULL
}

n_cpus <- 50
n_repl <- 500

bool_save <- TRUE
if (bool_save) {
  res_dir <- paste0(out_dir, "/coverage_n_repl_", n_repl, "_sim_pulling_",
                    ifelse(model_choice_K, "model_choice_K_", ""),
                    ifelse(generate_from_univ, paste0("gen_from_univ_corr_",
                                  paste0(format(rho_Zeta, digits = 2),
                                         collapse = "-"), "_"), ""),
                    "p_", p, "_N_", N, "_Nt_min_small_",
                    N_t_min_small, "_max_small_", N_t_max_small, "_Nt_min_",
                    N_t_min, "_max_", N_t_max, "_K_", K, "_L_sim_", L_sim, 
                    "_L_", L, "_sd_zeta_", 
                    paste0(format(vec_sd_zeta , digits = 2), collapse = "-"),
                    "_sd_eps_", unique(vec_sd_eps), "_tol_", tol, "_maxit_", 
                    maxit, "_seed_", seed, "/")
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

f_mu <- function(t, j) (-1)^j*2*sin((2*pi+j)*t)
f_psi_1 <- function(t, j, p) (-1)^j * sqrt(2/p)*cos(2*pi*t)
f_psi_2 <- function(t, j, p) (-1)^j * sqrt(2/p)*sin(2*pi*t)


f_Psi <- function(time_obs, j, p) {
  ans <- cbind(f_psi_1(time_obs, j, p),
               f_psi_2(time_obs, j, p))
  return(ans)
}


list_res <- parallel::mclapply(1:n_repl, function(repl) {
  
  set.seed(repl)
  
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
  

  if (model_choice_K) {
    res <- run_mfvb_fpca_model_choice(time_obs, Y, L = L, n_g = NULL,
                         time_g = time_g, 
                         tol = tol, maxit = maxit, Psi_g = Psi_g, n_cpus = 1)
  } else {
    res <- run_mfvb_fpca(time_obs, Y, L = L, K = K, n_g = NULL,
                         time_g = time_g, 
                         tol = tol, maxit = maxit,
                         plot_elbo = FALSE, Psi_g = Psi_g)
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
  
  list_time_obs_univ <- lapply(1:p, function(j) lapply(time_obs, function(time_obs_i) time_obs_i[j]))
  list_Y_univ <- lapply(1:p, function(j) lapply(Y, function(Y_i) Y_i[j]))
  
  if (model_choice_K) {
    list_res_univ <- lapply(seq_along(list_Y_univ), function(j)
      run_mfvb_fpca_model_choice(list_time_obs_univ[[j]], list_Y_univ[[j]], L = L, n_g = NULL,
                    time_g = time_g, tol = tol, maxit = maxit,
                    Psi_g = Psi_g[j], n_cpus = 1)) 
  } else {
    list_res_univ <- lapply(seq_along(list_Y_univ), function(j)
      run_mfvb_fpca(list_time_obs_univ[[j]], list_Y_univ[[j]], L = L, K = K, n_g = NULL,
                    time_g = time_g, tol = tol, maxit = maxit,
                    plot_elbo = FALSE, Psi_g = Psi_g[j]))
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
  
  Zeta_hat_mFPCA <- res$Zeta_hat
  sd_Zeta_hat_mFPCA <- cbind(sqrt(sapply(res$Cov_zeta_hat, function(cov_ii) cov_ii[1, 1])),
                             sqrt(sapply(res$Cov_zeta_hat, function(cov_ii) cov_ii[2, 2])))
  colnames(sd_Zeta_hat_mFPCA) <- c("FPC_1", "FPC_2")
  rownames(sd_Zeta_hat_mFPCA) <- rownames(Zeta_hat_mFPCA)
  
  list_Zeta_hat <- lapply(list_res_univ, function(ll) sqrt(p)*ll$Zeta_hat) 
  list_sd_Zeta_hat <- lapply(list_res_univ, function(ll) { sd_Zeta_hat <- cbind(sqrt(sapply(ll$Cov_zeta_hat, function(cov_ii) p*cov_ii[1, 1])),
                                                                                sqrt(sapply(ll$Cov_zeta_hat, function(cov_ii) p*cov_ii[2, 2])))
  rownames(sd_Zeta_hat) <- rownames(Zeta_hat_mFPCA)
  colnames(sd_Zeta_hat) <- c("FPC_1", "FPC_2")
  sd_Zeta_hat
  })
  
  list_Zeta_hat <- append(list_Zeta_hat, list(Zeta_hat_mFPCA))
  list_sd_Zeta_hat <- append(list_sd_Zeta_hat, list(sd_Zeta_hat_mFPCA))
  
  names(list_Zeta_hat) <- names(list_sd_Zeta_hat)<- c(paste0("FPCA Y", 1:p), "mFPCA")
  
  rownames(Zeta) <- rownames(Zeta_hat_mFPCA)
  colnames(Zeta) <- colnames(Zeta_hat_mFPCA)[1:L_sim]
  
  list_credible_interval_length <- lapply(list_sd_Zeta_hat, function(ll) 2*1.96*ll)
  list_coverage <- lapply(seq_along(list_Zeta_hat), function(k) Zeta >= list_Zeta_hat[[k]][, 1:L_sim] - 1.96*list_sd_Zeta_hat[[k]] & Zeta <= list_Zeta_hat[[k]][, 1:L_sim] + 1.96*list_sd_Zeta_hat[[k]] )
  names(list_coverage) <- names(list_Zeta_hat)
  
  list_mean_credible_interval_length <- lapply(list_credible_interval_length, function(ll) colMeans(ll))
  list_mean_coverage <- lapply(list_coverage, function(ll) colMeans(ll))
  
  create_named_list(list_credible_interval_length, list_coverage, 
                    list_mean_credible_interval_length, list_mean_coverage, 
                    cumulated_pve, cumulated_pve_univ,
                    L_thres_95, L_thres_95_univ, 
                    L_thres_99, L_thres_99_univ,
                    K, K_univ)
  
  # does the variance need to change as well?
  
  
}, mc.cores = n_cpus)
names(list_res) <- paste0("repl_", 1:n_repl)

list_credible_interval_length <- lapply(list_res, "[[", "list_credible_interval_length")
list_coverage <- lapply(list_res, "[[", "list_coverage")

list_mean_credible_interval_length <- lapply(list_res, "[[", "list_mean_credible_interval_length")
list_mean_coverage <- lapply(list_res, "[[", "list_mean_coverage")

list_cumulated_pve <- lapply(list_res, "[[", "cumulated_pve")
list_cumulated_pve_univ <- lapply(list_res, "[[", "cumulated_pve_univ")

list_L_thres_95 <- lapply(list_res, "[[", "L_thres_95")
list_L_thres_95_univ <- lapply(list_res, "[[", "L_thres_95_univ")

list_L_thres_99 <- lapply(list_res, "[[", "L_thres_99")
list_L_thres_99_univ <- lapply(list_res, "[[", "L_thres_99_univ")

list_K <- lapply(list_res, "[[", "K")
list_K_univ <- lapply(list_res, "[[", "K_univ")


mean_coverage_FPC_1 <- sapply(1:p, function(k) sapply(list_mean_coverage, "[[", paste0("FPCA Y", k))["FPC_1",])
mean_coverage_FPC_1 <- cbind(mean_coverage_FPC_1, sapply(list_mean_coverage, "[[", "mFPCA")["FPC_1",])
colnames(mean_coverage_FPC_1) <- c(paste0("FPCA Y", 1:p),  "mFPCA")

mean_coverage_FPC_2 <- sapply(1:p, function(k) sapply(list_mean_coverage, "[[", paste0("FPCA Y", k))["FPC_2",])
mean_coverage_FPC_2 <- cbind(mean_coverage_FPC_2, sapply(list_mean_coverage, "[[", "mFPCA")["FPC_2",])
colnames(mean_coverage_FPC_2) <- c(paste0("FPCA Y", 1:p),  "mFPCA")

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
  pdf(paste0(res_dir, "/coverage_scores_FPC_1.pdf"),
      width = 7.8, height = 5.8, paper='special')
}
boxplot(mean_coverage_FPC_1, main = "Mean coverage scores FPC 1", 
        ylim = c(0, 1), ylab = "Mean coverage of 95% CI for Zeta 1") 
abline(h = 0.95, lty = 3)
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/coverage_scores_FPC_2.pdf"),
      width = 7.8, height = 5.8, paper='special')
}
boxplot(mean_coverage_FPC_2, main = "Mean coverage scores FPC 2", 
        ylim = c(0, 1), ylab = "Mean coverage of 95% CI for Zeta 2")
abline(h = 0.95, lty = 3)
if (bool_save) {
  dev.off()
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


require(ggplot2)
require(tidyr)
require(dplyr)

colnames(mean_coverage_FPC_1) <- colnames(mean_coverage_FPC_2) <- c("FPCA X1", "FPCA X2", "FPCA X3", "FPCA X4", "FPCA X5", "FPCA X6", "mFPCA")

# Convert matrices to data frames
df_FPC_1 <- as.data.frame(mean_coverage_FPC_1)
df_FPC_2 <- as.data.frame(mean_coverage_FPC_2)

# Add an identifier for the group
df_FPC_1$Scores <- "FPC_1"
df_FPC_2$Scores <- "FPC_2"

# Reshape the data frames to long format, excluding the Scores column
df_FPC_1_long <- df_FPC_1 %>%
  pivot_longer(cols = -Scores, names_to = "Component", values_to = "Value")

df_FPC_2_long <- df_FPC_2 %>%
  pivot_longer(cols = -Scores, names_to = "Component", values_to = "Value")

# Combine the long format data frames
combined_data <- bind_rows(df_FPC_1_long, df_FPC_2_long)

# Set custom colors for the components
component_colors <- c("FPC_1" = "grey", "FPC_2" = "white")

pl <- ggplot(combined_data, aes(x = Component, y = Value, fill = Scores)) +
  geom_boxplot(position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_fill_manual(values = component_colors) +
  labs(
    title = "Average coverage for FPC score credible intervals",
    y = "Coverage of 95% credible intervals\n(average across i = 1, ..., n)",
    x = ""
  ) +
  scale_x_discrete(labels = c("FPCA X1" = bquote("FPCA" ~ x^.(1)), 
                              "FPCA X2" = bquote("FPCA" ~ x^.(2)), 
                              "FPCA X3" = bquote("FPCA" ~ x^.(3)), 
                              "FPCA X4" = bquote("FPCA" ~ x^.(4)), 
                              "FPCA X5" = bquote("FPCA" ~ x^.(5)), 
                              "FPCA X6" = bquote("FPCA" ~ x^.(6)), 
                              "mFPCA" = bquote("mFPCA" ~ "[" * x^.(1) * ",..., "* x^.(6) * "]"))) +
  theme_classic() +  theme(title = element_text(face="bold"),
                           axis.title = element_text(face="plain"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line.x = element_line(color = "black", size = 0.1),
    axis.line.y = element_line(color = "black", size = 0.1))

if (bool_save) {
  pdf(paste0(res_dir, "/coverage_v2.pdf"),
      width = 8, height = 3.85, paper='special')
}
# Display the plot
print(pl)

if (bool_save) {
  dev.off()
}

# Create the ggplot boxplot with custom fill colors and LaTeX-style labels
combined_data$Scores <- factor(combined_data$Scores, levels = c("FPC_2", "FPC_1"))
combined_data$Component <- factor(combined_data$Component, levels = c("mFPCA", paste0("FPCA X", p:1)))

pl <- ggplot(combined_data, aes(x = Component, y = Value, fill = Scores)) +
  geom_boxplot(position = position_dodge(0.9), fatten = 2, outlier.size = 0.9) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_fill_manual(values = component_colors) +
  labs(
    title = "Average coverage\nfor FPC score credible intervals",
    y = "Coverage of 95% credible intervals\n (average across i = 1, ..., n)",
    x = ""
  ) +
  scale_x_discrete(labels = c("FPCA X1" = bquote("FPCA" ~ x^.(1)), 
                              "FPCA X2" = bquote("FPCA" ~ x^.(2)), 
                              "FPCA X3" = bquote("FPCA" ~ x^.(3)), 
                              "FPCA X4" = bquote("FPCA" ~ x^.(4)), 
                              "FPCA X5" = bquote("FPCA" ~ x^.(5)), 
                              "FPCA X6" = bquote("FPCA" ~ x^.(6)), 
                              "mFPCA" = bquote("mFPCA" ~ "[" * x^.(1) * ",...," * x^.(6) * "]"))) +
  theme_classic() +  theme(title = element_text(face="bold"),
                           axis.title = element_text(face="plain"),
                           panel.border = element_rect(color = "black", fill = NA, size = 1),
                           axis.line.x = element_line(color = "black", size = 0.1),
                           axis.line.y = element_line(color = "black", size = 0.1))+ coord_flip()
if (bool_save) {
  pdf(paste0(res_dir, "/coverage_poster.pdf"),
      width = 6, height = 6, paper='special')
}
# Display the plot
print(pl)

if (bool_save) {
  dev.off()
}


mean_credible_interval_length_FPC_1 <- sapply(1:p, function(k) sapply(list_mean_credible_interval_length, "[[", paste0("FPCA Y", k))["FPC_1",])
mean_credible_interval_length_FPC_1 <- cbind(mean_credible_interval_length_FPC_1, sapply(list_mean_credible_interval_length, "[[", "mFPCA")["FPC_1",])
colnames(mean_credible_interval_length_FPC_1) <- c(paste0("FPCA Y", 1:p),  "mFPCA")

mean_credible_interval_length_FPC_2 <- sapply(1:p, function(k) sapply(list_mean_credible_interval_length, "[[", paste0("FPCA Y", k))["FPC_2",])
mean_credible_interval_length_FPC_2 <- cbind(mean_credible_interval_length_FPC_2, sapply(list_mean_credible_interval_length, "[[", "mFPCA")["FPC_2",])
colnames(mean_credible_interval_length_FPC_2) <- c(paste0("FPCA Y", 1:p),  "mFPCA")

if (bool_save) {
  pdf(paste0(res_dir, "/CI_length_scores_FPC_1.pdf"),
      width = 7.8, height = 5.8, paper='special')
}
boxplot(mean_credible_interval_length_FPC_1, main = "Mean credible interval length scores FPC 1", ylab = "Mean length of 95% CI for Zeta 1")
if (bool_save) {
  dev.off()
  pdf(paste0(res_dir, "/CI_length_scores_FPC_2.pdf"),
      width = 7.8, height = 5.8, paper='special')
}
boxplot(mean_credible_interval_length_FPC_2, main = "Mean credible interval length scores FPC 2",  ylab = "Mean length of 95% CI for Zeta 2")
if (bool_save) {
  dev.off()
}

# Assuming these are the column names
colnames(mean_credible_interval_length_FPC_1) <- colnames(mean_credible_interval_length_FPC_2) <-
  c("FPCA X1", "FPCA X2", "FPCA X3", "FPCA X4", "FPCA X5", "FPCA X6", "mFPCA")

# Convert matrices to data frames
df_FPC_1 <- as.data.frame(mean_credible_interval_length_FPC_1)
df_FPC_2 <- as.data.frame(mean_credible_interval_length_FPC_2)

# Add an identifier for the group
df_FPC_1$Scores <- "FPC_1"
df_FPC_2$Scores <- "FPC_2"

# Reshape the data frames to long format, excluding the Scores column
df_FPC_1_long <- df_FPC_1 %>%
  pivot_longer(cols = -Scores, names_to = "Component", values_to = "Value")

df_FPC_2_long <- df_FPC_2 %>%
  pivot_longer(cols = -Scores, names_to = "Component", values_to = "Value")

# Combine the long format data frames
combined_data <- bind_rows(df_FPC_1_long, df_FPC_2_long)

# Set custom colors for the components
component_colors <- c("FPC_1" = "grey", "FPC_2" = "white")

# Create the ggplot boxplot with custom fill colors and LaTeX-style labels
pl <- ggplot(combined_data, aes(x = Component, y = Value, fill = Scores)) +
  geom_boxplot(position = position_dodge(0.9)) +
  scale_fill_manual(values = component_colors) +
  labs(
    title = "Average length of FPC score credible intervals",
    y = "Length of 95% credible intervals\n(average across i = 1, ..., n)",
    x = ""
  ) +
  scale_x_discrete(labels = c("FPCA X1" = bquote("FPCA" ~ x^.(1)),
                              "FPCA X2" = bquote("FPCA" ~ x^.(2)),
                              "FPCA X3" = bquote("FPCA" ~ x^.(3)),
                              "FPCA X4" = bquote("FPCA" ~ x^.(4)),
                              "FPCA X5" = bquote("FPCA" ~ x^.(5)),
                              "FPCA X6" = bquote("FPCA" ~ x^.(6)),
                              "mFPCA" = bquote("mFPCA" ~ "[" * x^.(1) * ",...," * x^.(6) * "]"))) +
  theme_classic() +  theme(title = element_text(face="bold"),
                           axis.title = element_text(face="plain"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line.x = element_line(color = "black", size = 0.1),
    axis.line.y = element_line(color = "black", size = 0.1))


if (bool_save) {
  pdf(paste0(res_dir, "/interval_length_v2.pdf"),
      width =  8, height = 3.85, paper='special')
}
# Display the plot
print(pl)
if (bool_save) {
  dev.off()
}


if (bool_save) {
  save(list_res, file = file.path(res_dir, "output.RData"))
}

