rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

bool_cluster <- TRUE
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

seed <- 1
set.seed(seed)

# Number of simulations:
Nsims = 500

# Number of functional observations:
n = 200

# Number of observation points
m = 15

# Parametric function:
ind_type <- 1
par_fun = c("linear", "nelson-siegel")[ind_type] 

# TRUE number of nonparametric factors:
L_nonparam_true = c(0, 1, 3)[job_id]  

# Root-signal-to-noise ratio
RSNR = 3 

# for bayesFCPA inference only
L <- 10 
tol <- 1e-3 

if (par_fun == "linear") {
  L_sim <- L_nonparam_true + 2
} else if (par_fun == "nelson-siegel"){
  L_sim <- L_nonparam_true + 3
}

model_choice_K <- T
if (model_choice_K) {
  K <- NULL
} else {
  bool_K_Ruppert <- T
  if (bool_K_Ruppert) {
    K <- NULL                     # default based on the rule in Ruppert 2002
  } else {
    n_int_knots <- 5              # number of interior knots
    K <- n_int_knots + 2          # number of spline basis functions
  }
}

n_cpus <- 100

bool_save <- TRUE
if (bool_save) {
  res_dir <- paste0(out_dir, "/comparison_SFFM_PFFM_p_1_", 
                    ifelse(model_choice_K, "_model_choice_K", paste0("K_", K)),
                    "_N_", n, "_Nt_", m, "_L_sim_", L_sim, "_param_fct_", par_fun, 
                    "_L_nonparam_true_", L_nonparam_true, "_RSNR_", RSNR, 
                    "_bayesFPCA_tol_", tol, "_seed_", seed, "_n_repl_", Nsims, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}

# Known or unknown parameters?
unknown_gamma = FALSE # Only relevant for Nelson-Siegel basis 

# MCMC iterations:
nsave = 5000; nburn = 1000; nskip = 0

# Upper bound on number of factors:
L_nonparam = ifelse(L_nonparam_true < 8, 
                    10, # for most cases 
                    15) # when L_nonparam_true is larger
#----------------------------------------------------------------------------
# Load in the necessary files:
library(ggplot2)
library(dfosr); library(fGarch); 
library(coda); library(fields); library(gtools)
require(bayesFPCA)
source('SFFM_scripts/helper_functions.R')
source('SFFM_scripts/source_sffm.R')
source('SFFM_scripts/source_sffm_gp.R')
#----------------------------------------------------------------------------
# Define the parametric functions and the priors:
#----------------------------------------------------------------------------
# Observation points:
tau = seq(0, 1, length.out = m)
#----------------------------------------------------------------------------
if(par_fun == "linear"){
  unknown_gamma = FALSE
  g = function(tau) qr.Q(qr(cbind(1, tau)))
  log_prior_gamma = NULL
}
#----------------------------------------------------------------------------
if(par_fun == "nelson-siegel"){
  gamma_true = 0.0609; 
  if(unknown_gamma){
    g = function(tau, gamma){tau = tau*100 + 1; 
    qr.Q(qr(cbind(1, (1 - exp(-tau*gamma))/(tau*gamma), (1 - exp(-tau*gamma))/(tau*gamma) - exp(-tau*gamma))))}
    
    # Center the prior at the Diebold-Li estimate with weakly informative variance:
    mu_gamma = 0.0609; var_gamma = 1/2
    log_prior_gamma = function(x) dgamma(x, 
                                         shape = mu_gamma^2/var_gamma, 
                                         rate = mu_gamma/var_gamma, log = TRUE)
    
  } else {
    gamma = gamma_true; g = function(tau){tau = tau*100 + 1; 
    qr.Q(qr(cbind(1, (1 - exp(-tau*gamma))/(tau*gamma), (1 - exp(-tau*gamma))/(tau*gamma) - exp(-tau*gamma))))}
    log_prior_gamma = NULL
  }
}

# Parametric function:
if(unknown_gamma){G_true = g(tau, gamma_true)} else G_true = g(tau); 
L_G_true = ncol(G_true)
#----------------------------------------------------------------------------
# Storage:
#----------------------------------------------------------------------------
e_names = c('PFFM', "PFFM+gp", "SFFM(fmm)", "SFFM", "bayesFPCA")
rmse_y_true = array(NA, c(Nsims, length(e_names)), dimnames = list(1:Nsims, e_names))
pi_cover = pi_width = array(NA, c(Nsims, length(e_names)), dimnames = list(1:Nsims, e_names))
rmse_alpha = array(NA, c(Nsims, length(e_names[e_names != "bayesFPCA"])), 
                   dimnames = list(1:Nsims, e_names[e_names != "bayesFPCA"]))
ci_alpha_cover = ci_alpha_width = array(NA, c(Nsims, length(e_names[e_names != "bayesFPCA"])), 
                                        dimnames = list(1:Nsims, e_names[e_names != "bayesFPCA"]))
prob_L_nonparam = prob_L_nonparam_fmm = array(NA, c(Nsims, L_nonparam+1))

mat_cum_pve_bayesFPCA <- array(NA, c(Nsims, L), dimnames = list(1:Nsims,  1:L))

for(ni in 1:Nsims) {
  set.seed(ni)
  #----------------------------------------------------------------------------
  # Simulate the parametric part:
  Alpha_true = matrix(rnorm(n = n*L_G_true), ncol = L_G_true)
  Y_true_p = tcrossprod(Alpha_true, G_true) 
  
  # Simulate the nonparametric part:
  if(L_nonparam_true > 0){
    # Orthogonal polynomial, then constrain to be orthogonal to G_true:
    F_true = poly(tau, L_nonparam_true + 1)[,-1]
    F_true = qr.Q(qr(cbind(G_true, F_true)))[,-(1:L_G_true)]
    
    Beta_true = matrix(sapply(1:L_nonparam_true, 
                              function(h) rnorm(n = n, sd = 1/(h + 1))), ncol=L_nonparam_true)
    Y_true_np = tcrossprod(Beta_true, F_true) 
  } else Y_true_np = 0  
  
  # Simulate the data:
  Y_true = Y_true_p + Y_true_np
  Y = Y_true + sd(Y_true)/RSNR*rnorm(n*m)
  
  # Test data:
  Y_test = Y_true + sd(Y_true)/RSNR*rnorm(n*m)
  
  # Plot a sampled curve:
  if(ni == 1){
    i = sample(1:n, 1); plot(tau, Y[i,]); lines(tau, Y_true[i,])
  }
  #----------------------------------------------------------------------------
  # Parametric functional factor model
  #----------------------------------------------------------------------------
  fit_p = pffm(Y = Y, tau = tau, g = g, log_prior_gamma = log_prior_gamma,
               nsave = nsave, nburn = nburn, nskip = nskip)
  #----------------------------------------------------------------------------
  # Parametric functional factor model centered at GP
  #----------------------------------------------------------------------------
  fit_gp = sffm_gp(Y = Y, tau = tau, g = g, log_prior_gamma = log_prior_gamma,
                   rho = 0.5,
                   nsave = nsave, nburn = nburn, nskip = nskip)
  #----------------------------------------------------------------------------
  # Semiparametric functional factor model w/ finite mixture prior
  #----------------------------------------------------------------------------
  fit_fmm = sffm_fmm(Y = Y, tau = tau, g = g, log_prior_gamma = log_prior_gamma,
                     K = L_nonparam, K_0 = 1, a_1 = 5, a_2 = 25, v0 = 0.001,
                     nsave = nsave, nburn = nburn, nskip = nskip)
  #----------------------------------------------------------------------------
  # Semiparametric functional factor model (proposed)
  #----------------------------------------------------------------------------
  fit_sp = sffm(Y = Y, tau = tau, g = g, log_prior_gamma = log_prior_gamma,
                K = L_nonparam, K_0 = NULL, a_1 = 5, a_2 = 25, v0 = 0.001,
                nsave = nsave, nburn = nburn, nskip = nskip)
  #----------------------------------------------------------------------------
  # bayesFPCA model (univariate VMP)
  #----------------------------------------------------------------------------
  
  time_obs <- lapply(1:n, function(i) tau)
  Y_list <- lapply(1:n, function(i) Y[i,])
  
  
  if (model_choice_K) {
    rt <- system.time(fit_bayesfpca <- run_mfvb_fpca_model_choice(time_obs = time_obs, 
                                                                  Y = Y_list, L = L, tol = tol, n_cpus = 1)) 
  } else {
    rt <- system.time(fit_bayesfpca <- run_mfvb_fpca(time_obs = time_obs, Y = Y_list, L = L, K = K, tol = tol)) 
  }
  print(paste0("bayesFPCA runtime: ", rt["elapsed"]))
  
  mat_cum_pve_bayesFPCA[ni,] <- fit_bayesfpca$cumulated_pve
  
  #----------------------------------------------------------------------------
  # Results:
  #----------------------------------------------------------------------------
  # RMSE for Y_true:
  rmse_y_true[ni,1] = sqrt(mean((Y_true - colMeans(fit_p$Yhat))^2))
  rmse_y_true[ni,2] = sqrt(mean((Y_true - colMeans(fit_gp$Yhat))^2))
  rmse_y_true[ni,3] = sqrt(mean((Y_true - colMeans(fit_spline$Yhat))^2))
  rmse_y_true[ni,4] = sqrt(mean((Y_true - colMeans(fit_fmm$Yhat))^2))
  rmse_y_true[ni,5] = sqrt(mean((Y_true - colMeans(fit_sp$Yhat))^2))
  
  ind_time_obs_rec <- sapply(tau, function(time_obs_i) which.min(abs(time_obs_i - fit_bayesfpca$time_g)))
  Y_hat_obs_rec <- t(sapply(fit_bayesfpca$Y_hat, function(ll) ll[ind_time_obs_rec]))
  rmse_y_true[ni,6] = sqrt(mean((Y_true - Y_hat_obs_rec)^2))
  #----------------------------------------------------------------------------
  # Prediction intervals:
  pi_p = pi_gp = pi_spline = pi_fmm = pi_sp = pi_bayesfpca = array(NA, c(n,m,2))
  for(i in 1:n){
    pi_p[i,,] = t(apply(fit_p$Ypred[,i,], 2, quantile, c(.025, .975)))  
    pi_gp[i,,] = t(apply(fit_gp$Ypred[,i,], 2, quantile, c(.025, .975)))
    pi_spline[i,,] = t(apply(fit_spline$Ypred[,i,], 2, quantile, c(.025, .975)))
    pi_fmm[i,,] = t(apply(fit_fmm$Ypred[,i,], 2, quantile, c(.025, .975)))
    pi_sp[i,,] = t(apply(fit_sp$Ypred[,i,], 2, quantile, c(.025, .975)))
    
    pi_bayesfpca[i,,] <- cbind(fit_bayesfpca$Y_low[[i]][ind_time_obs_rec], fit_bayesfpca$Y_upp[[i]][ind_time_obs_rec]) 
  }
  # Coverage:
  pi_cover[ni,1] = mean(pi_p[,,1] <= Y_test & pi_p[,,2] >= Y_test)
  pi_cover[ni,2] = mean(pi_gp[,,1] <= Y_test & pi_gp[,,2] >= Y_test)
  pi_cover[ni,3] = mean(pi_spline[,,1] <= Y_test & pi_spline[,,2] >= Y_test)
  pi_cover[ni,4] = mean(pi_fmm[,,1] <= Y_test & pi_fmm[,,2] >= Y_test)
  pi_cover[ni,5] = mean(pi_sp[,,1] <= Y_test & pi_sp[,,2] >= Y_test)
  pi_cover[ni,6] = mean(pi_bayesfpca[,,1] <= Y_test & pi_bayesfpca[,,2] >= Y_test)
  
  # Width
  pi_width[ni,1] = mean(pi_p[,,2] - pi_p[,,1])
  pi_width[ni,2] = mean(pi_gp[,,2] - pi_gp[,,1])
  pi_width[ni,3] = mean(pi_spline[,,2] - pi_spline[,,1])
  pi_width[ni,4] = mean(pi_fmm[,,2] - pi_fmm[,,1])
  pi_width[ni,5] = mean(pi_sp[,,2] - pi_sp[,,1])
  pi_width[ni,6] = mean(pi_bayesfpca[,,2] - pi_bayesfpca[,,1])
  #----------------------------------------------------------------------------
  # K* inference from the SFFM and fmm version:
  prob_L_nonparam[ni,] = sapply(0:L_nonparam, function(js) mean(js == fit_sp$K_star, na.rm=TRUE)); 
  prob_L_nonparam_fmm[ni,] = sapply(0:L_nonparam, function(js) mean(js == fit_fmm$K_star, na.rm=TRUE)); 
  #----------------------------------------------------------------------------
  rm(fit_p, fit_gp, fit_spline, fit_fmm, fit_sp)
  #----------------------------------------------------------------------------
  # Progress:
  #----------------------------------------------------------------------------
  print('----------------------------------------------------------------------------')
  print(ni)  
  print('----------------------------------------------------------------------------')
}

if (bool_save) {
  save.image(file = file.path(res_dir, "output.RData"))
}


if (bool_save) {
  pdf(paste0(res_dir, "/proportion_variance_explained_L_sim_", L_sim, ".pdf"),
      width = 6.7, height = 5)
}
par(mfrow = c(1,1), mar = c(5.1, 5.1, 4.1, 2.1), las = 1)
boxplot(mat_cum_pve_bayesFPCA, type = "o", pch = 20, ylim = c(0, 100),
        col = "white",
        main = "Cumulated PVE across components",
        xlab = "Number of components",
        ylab = "Cumulated PVE",
        cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, notch = FALSE, outline = FALSE)
abline(h = 95, col = "gray20", lty = 2, lwd = 2)
abline(v = L_sim, col = "red", lwd = 2.5, lty = 3)
legend("bottomright", 
       leg = c("Simulated nb of latent functions", 
               "Cumulated PVE threshold"), col = c("red", "gray20"),
       lty = c(3, 2), lwd = c(2.5, 2), bty = "n")
if (bool_save) {
  dev.off()
}


fname = paste("temp", # modify w/ your directory
              par_fun,"_L_nonparamtrue=", L_nonparam_true, sep="")
use_inds = seq_along(e_names)# match(c("PFFM", "PFFM+gp", "SFFM"), e_names)
#----------------------------------------------------------------------------
if(bool_save) pdf(file = paste0(res_dir, "/", fname, "_rmse-y.pdf"), width = 10, height = 6)
par(mai = c(1, 2.5, 1, 0.5), las = 1);
bp = boxplot(rmse_y_true[,use_inds], 
             main = paste('Prediction errors'),
             horizontal = TRUE,
             xlab = 'Root mean squared prediction error',
             ylim = c(0.0, 0.2),
             col = 'white',
             cex.lab = 2, cex.main = 2, cex.axis = 2, notch = FALSE, outline = FALSE);
if(bool_save) dev.off()
#----------------------------------------------------------------------------
# Mean predictive interval width:
if(bool_save) pdf(file = paste0(res_dir, "/", fname, "_pi.pdf"), width = 10, height = 6)
par(mai = c(1, 2.5, 1, 0.5), las = 1);
bp0 = boxplot(pi_width[,use_inds], plot = FALSE)
bp <- boxplot(
  pi_width[, use_inds],
  main = "Mean prediction interval widths", 
  horizontal = TRUE,
  col = 'white',
  xlab = "Interval width",
  ylim = c(0, 1.3),
  cex.lab = 2, cex.main = 2, cex.axis = 2,
  notch = FALSE, outline = FALSE
)
text(bp$stats[5,], 1:length(use_inds), 
     labels = paste(' ', round(100*apply(pi_cover[,use_inds], 2, mean, na.rm=TRUE)), '%', sep=''), 
     cex = 2, adj = c(0,NA), col='red')
if(bool_save) dev.off()
#----------------------------------------------------------------------------
if(bool_save) pdf(file = paste0(res_dir, "/", fname, "_L_nonparamprob", ".pdf"), width = 3, height = 6)
prob_L_nonparam_all = cbind(prob_L_nonparam_fmm[,L_nonparam_true + 1],
                            prob_L_nonparam[,L_nonparam_true + 1]
); colnames(prob_L_nonparam_all) = c("SFFM-fmm", "SFFM")
par(las = 1, las = 2, mar = c(7.9, 4.6, 4.1, 2.1));
bp = boxplot(prob_L_nonparam_all, main = paste("Posterior probability \n for L (non-param)"), 
             horizontal = F,
             col = 'white',
             ylab = "Posterior probability",
             ylim = c(0, 1),
             cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, notch = F, outline = FALSE);
if(bool_save) dev.off()
#----------------------------------------------------------------------------
# CRPS:
library(scoringRules)
if(bool_save) pdf(file = paste0(res_dir, "/", fname, "_L_nonparamcrps.pdf"), width = 10, height = 3.7)
crps_L_nonparam_all = cbind(
  crps_sample(rep(1+L_nonparam_true, Nsims), 
              prob_L_nonparam_fmm),
  crps_sample(rep(1+L_nonparam_true, Nsims), 
              prob_L_nonparam)
); colnames(crps_L_nonparam_all) = c("SFFM-fmm", "SFFM")
par(mai = c(.5, 3.25, 1, 0.5), las = 1);
bp = boxplot(crps_L_nonparam_all, main = paste("Ranked probability score for L (non-param)"), 
             horizontal = TRUE,
             col = 'white',
             cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5, notch = F, outline = FALSE);
if(bool_save) dev.off()
# }
#----------------------------------------------------------------------------
# Probability of overestimating L_nonparam_true:
round(mean(rowSums(prob_L_nonparam[,-(1:(L_nonparam_true+1))]), na.rm=TRUE), 3)
round(mean(rowSums(prob_L_nonparam_fmm[,-(1:(L_nonparam_true+1))]), na.rm=TRUE), 3)

# Probability exceed zero
round(summary(1 - prob_L_nonparam[,1]), 3)
round(summary(1 - prob_L_nonparam_fmm[,1]), 3)

# Probability mass within one of the true value:
round(summary(rowSums(prob_L_nonparam[, max(L_nonparam_true - 1, 0):(L_nonparam_true + 1) + 1])), 3)
round(summary(rowSums(prob_L_nonparam_fmm[, max(L_nonparam_true - 1, 0):(L_nonparam_true + 1) + 1])), 3)
#----------------------------------------------------------------------------
