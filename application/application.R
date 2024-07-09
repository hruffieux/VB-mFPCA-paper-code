rm(list = ls())

CORE_DIR <- Sys.getenv("CORE_DIR")

bool_cluster <- F
if (bool_cluster) {
  
  #!/usr/bin/env Rscript
  args = commandArgs(trailingOnly=TRUE)
  job_id <- as.integer(args[1]) # 1 to 20
  
  CORE_DIR_MRC_BSU <- Sys.getenv("CORE_DIR_MRC_BSU")
  out_dir <- file.path(CORE_DIR_MRC_BSU, "mFPCA_output/")
} else {
  
  job_id <- 1 # 1 to 6
  
  CORE_DIR_ICLOUD <- Sys.getenv("CORE_DIR_ICLOUD")
  out_dir <- file.path(CORE_DIR_ICLOUD, "mFPCA_output/")
}

data_dir <- file.path(CORE_DIR, "covid-19-metabo/data/")
main_dir <- file.path(CORE_DIR, "bayesian-mFPCA-paper-code/application/")
setwd(main_dir)


require(bayesFPCA)
require(dplyr)

source("fun_application.R")

if (bool_cluster) {
  max_job <- 20
  seed <- seq(1, max_job, by = 1)[job_id] 
} else {
  seed <- 12345 
}
set.seed(seed)

version_ms <- c("20201009", "11DEC2020", "20210311", "20210614", "20210630")[5]
version_nmr <- "201120"
version_glyc <- "010721"
version_lipid <- "20210614"
version_ques <- "20210505"
version_info <- c("20201123", "20210218", "20210427", "20210427_fu")[4]
version_flow <- c("201221", "HC2_210429")[2]
version_cytof <- "201005"

keep_replicates <- FALSE

load(paste0(data_dir,
            "preprocessed_data/metabolites_ms_v", version_ms,
            "_nmr_v", version_nmr,
            "_glyc_v", version_glyc,
            "_lipids_v", version_lipid,
            "_clin_v", version_info,
            "_ques_v", version_ques,
            ifelse(keep_replicates, "_with_replicates", ""),
            ".RData"))

load(paste0(data_dir, "preprocessed_data/cell_types_flow_v", version_flow,
            "_cytof_v", version_cytof,
            "_clin_v", version_info, ".RData"))

data_id <- 3

data_name <- c("NMR_SM", # 1
               "NMR_LP_abundances", # 2
               "MS", # 3
               "lipids", # 4
               "cell_types", # 5
               "cytokines", # 6
               "complements", # 7
               "spike_ig", # 8
               "glyc", # 9
               "ratios", # 10
               "log_ratios", # 11
               "clinical", # 12
               "bleed", # 13
               "log_bleed", # 14
               "cell_types_main", # 15
               "inflammation" # 16
)[data_id]


df_comb <- list(df_nmr_sm_comb,
                df_nmr_lp_comb,
                df_ms_comb,
                df_lipid_comb,
                df_ct_comb,
                df_cytokine_comb,
                df_cplt_comb,
                df_auc_comb,
                df_glyc_comb,
                df_all_ratios_comb,
                df_all_ratios_comb,
                df_info,
                df_bleed_avg_comb,
                df_bleed_avg_log_comb,
                df_flow_main_comb,
                df_info
                # df_cytokine_comb
)[[data_id]]

all_var <- sort(names(list(df_nmr_sm,
                           df_nmr_lp,
                           df_ms,
                           df_lipid,
                           df_ct,
                           df_cytokine,
                           df_cplt,
                           df_auc,
                           df_glyc,
                           df_all_ratios,
                           df_all_ratios,
                           df_info,
                           df_bleed_avg,
                           df_bleed_avg_log,
                           df_flow_main,
                           df_info)[[data_id]]))

if (data_id == 6) {
  all_var <- setdiff(all_var, c("IFNg", "TNFa"))
} else if (data_id == 3) {
  
  all_var <- c("log_CRP",
               "IL10",
               "GlycB",
               "Quinolinic acid",
               "Tryptophan") 
  vec_flip <- c(1, 1)
} else if (data_id == 9) {
  all_var <- c("GlycA", "GlycB")
} else if (data_id == 16) {
  all_var <- c("log_CRP") 
  vec_flip <- c(1, -1)
} else if (data_id == 5) {
  all_var <- c("NK",  "IL10", "log_CRP")
} else if (data_id == 7) {
  all_var <- names(df_cplt)
} else if (data_id %in% c(13, 14)) {

  all_var <- list(
    c("CRP", "LY#"),
    c("CRP", "NLR", "Platelet count", "Lymphocytes",
      "Fibrinogen", "D-Dimers", "LACTATE"),
    c("CRP", "NLR", "Platelet count", "Lymphocytes",
      "Fibrinogen", "D-DIMER", "LACTATE"),
    c("CRP", "NLR", "Platelet count", "LY#",
      "Fibrinogen", "D-Dimers", "LACTATE"),
    c("CRP", "NLR", "Platelet count", "LY#",
      "Fibrinogen", "D-DIMER", "LACTATE")
  )[[1]]
}

mess_version <- c(paste0("_nmr_v", version_nmr),
                  paste0("_nmr_v", version_nmr),
                  paste0("_ms_v", version_ms),
                  paste0("_lipid_v", version_lipid),
                  paste0("_flow_v", version_flow, "_cytof_v", version_cytof),
                  rep("", 3),
                  paste0("_glyc_v", version_glyc),
                  paste0("_ms_v", version_ms),
                  paste0("_ms_v", version_ms),
                  rep("", 3),
                  paste0("_flow_v", version_flow),
                  "")[data_id]

bool_save <- F

days_thres <- 50
tnew <- 0:(days_thres-1) # timepoints to evaluate the prediction at

min_timepoints <- 2 

selected_severity_groups <- c("B", "C", "D", "E")

df_comb_no_HC <- choose_severity_groups(df_comb, selected_severity_groups)
df_comb_no_HC <- df_comb_no_HC[df_comb_no_HC$days_from_sx_or_swab_imputed < days_thres, ]

subset_subj <- F
all_var_mFPCA <- c("log_CRP",
                   "IL10",
                   "GlycB",
                   "Quinolinic acid",
                   "Tryptophan") # OK

bool_both <- T 
if (subset_subj & !isTRUE(all.equal(all_var, all_var_mFPCA))) {

  load(paste0(out_dir, "COVID_application/surv_irreg_grid_long_covid_",gsub(" ", "_",
                                                                               paste0(all_var_mFPCA, collapse = "_")),
              ifelse(subset_subj, "_subset_subj", ""),
              ifelse(bool_both, "_joint_surv", ""), "_seed_", seed, ".RData"))

  df_comb_no_HC <- df_comb_no_HC[df_comb_no_HC$subject_id %in% subject_ids, ]

  rm(pval_surv_EF1, pval_surv_EF2, n, subject_ids)

}


bool_rescale_time <- TRUE
data <- prepare_data_irregular(df_comb_no_HC, all_var, min_nb_timepoints = min_timepoints,
                               days_thres = days_thres, bool_rescale_time = bool_rescale_time)

Y <- data$Y
p <- data$p
N <- data$N
time_obs <- data$time_obs
subject_ids <- names(Y)

print(paste0("Number of subjects: ", N))

model_choice_K <- T

if (model_choice_K) {
  K <- NULL
} else {
  bool_K_Ruppert <- F
  if (bool_K_Ruppert) {
    K <- NULL                       # default based on the rule in Ruppert 2002
  } else {
    n_int_knots <- 5                # number of interior knots
    K <- n_int_knots + 2            # number of spline basis functions
  }
}
 
L <- 10                             # number of FPCA basis functions

tol  <- 1e-4                        # convergence tolerance 
maxit_mfvb <- 1000                  # maximum number of mfvb iterations

n_g <- 1000                         # length of the plotting grid

n_cpus <- 4

if (bool_save) {
  res_dir <- paste0(out_dir, "/FPCA_application_", mess_version,
                    "_clin_v", version_info, "_",
                    paste0(selected_severity_groups, collapse = "-"),
                    "_days_thres_", days_thres, "_",
                    ifelse(model_choice_K, "model_choice_K_", ""),
                    ifelse(length(all_var)>10,
                           paste0(substring(all_var, 1, 7), collapse = "-"),
                           paste0(all_var, collapse = "-")),
                    "_min_timepoints_", min_timepoints,
                    ifelse(subset_subj, "_subset_subj", ""), "_L_", L,
                    "_tol_", tol,
                    "_seed_", seed, "/")
  dir.create(res_dir)
  sink(paste(res_dir, "out.txt", sep = ""), append = F, split = T, type = "output")
  sink(file(paste(res_dir, "err.txt", sep = ""), open = "wt"), type = "message")
} else {
  res_dir <- NULL
}


####################################################
#
#  MFPCA RUN
#
####################################################


if (model_choice_K) {
  run_time <- system.time({ 
    res <- run_mfvb_fpca_model_choice(time_obs, Y, L = L, n_g = n_g, 
                        tol = tol, maxit = maxit_mfvb,
                        n_cpus = n_cpus, seed = seed)
  })
} else {
  run_time <- system.time({ 
    res <- run_mfvb_fpca(time_obs, Y, L = L, K = K, n_g = n_g,
                        tol = tol, maxit = maxit_mfvb,
                        plot_elbo = F, seed = seed)
  })
}

cumulated_pve <- res$cumulated_pve

print(cumulated_pve)
L_thres_95 <- sum(cumulated_pve < 95) + 1
L_thres_99 <- sum(cumulated_pve < 99) + 1

K <- unique(res$K)

print(format(cumulated_pve, digits = 4, scientific = F)) 
var_expl <-  res$cumulated_pve[1:2]
var_expl[2] <- var_expl[2] - var_expl[1]
time_g <- res$time_g


if (bool_save) {
  pdf(paste0(res_dir, "/prop_var_explained.pdf"),
      width = 4, height = 2.9, paper='special')
}

library(ggplot2)

# Assuming cumulated_pve is a numeric vector
component_index <- seq_along(cumulated_pve) # create component index

# Create a data frame for ggplot
df <- data.frame(
  ComponentIndex = component_index,
  PVE = cumulated_pve
)

# Plot using ggplot2
ggplot(df, aes(x = ComponentIndex, y = PVE)) +
  geom_point(shape = 20, size = 2) +
  geom_line() +
  labs(
    x = "Component index, l",
    y = "PVE",
    title = "Cumulated PVE across components"
  ) +
  theme_classic() +
  scale_x_continuous(breaks = 1:10) +
  ylim(0, 100) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.line.x.bottom = element_line(color = "black", size = 0.1),
    axis.line.y.left = element_line(color = "black", size = 0.1),
    axis.line.x.top = element_line(color = "black", size = 0.1),
    axis.line.y.right = element_line(color = "black", size = 0.1),
    plot.title = element_text(face = "bold")
  )

if (bool_save) {
  dev.off()
}

mu_hat <- res$mu_hat
list_Psi_hat <- res$list_Psi_hat

Y_hat <- res$Y_hat
Y_low <- res$Y_low
Y_upp <- res$Y_upp
Zeta_hat <- res$Zeta_hat
zeta_ellipse <- res$list_zeta_ellipse
Cov_zeta_hat <- res$Cov_zeta_hat

out_flip <- flip_sign(vec_flip, list_Psi_hat, Zeta_hat, zeta_ellipse)

list_Psi_hat <- out_flip$list_Psi_hat
Zeta_hat <- out_flip$Zeta_hat
zeta_ellipse <- out_flip$zeta_ellipse


scores <- get_df_scores(subject_ids, Zeta_hat, df_info, vec_col)

if (bool_save) {
  save.image(file = file.path(res_dir, "objects.RData"))
}

# Select subjects to display
#
vec_subj <- c(scores$subject_id[which.min(scores$EF1)],
              scores$subject_id[which.max(scores$EF1)],
              scores$subject_id[which.min(scores$EF2)],
              scores$subject_id[which.max(scores$EF2)]) 

for (vv in all_var) {
  if (bool_save) {
    pdf(paste0(res_dir, "/predicted_trajectories_", vv, "_extreme_subjects_trunc.pdf"),
        width = 6.3, height = 6.3, paper='special')
  }

  plot_trajectories_irregular_grid(vec_subj, vv, Y, Y_hat, Y_upp, Y_low,
                                   days_thres, time_obs, time_g, df_comb,
                                   scores, bool_rescale_time,
                                   ylim_offset = 2.6)

  if (bool_save) {
    dev.off()
  }

}

if (bool_save) {
  pdf(paste0(res_dir, "/predicted_trajectories_all_var_extreme_subjects_trunc_large_offset.pdf"),
      width = 12.5, height = 9.5, paper='special')
}
par(mfrow = c(length(vec_subj), length(all_var)), mar = c(4, 4.1, 2, 2.1))
for (subj in vec_subj) {
  for (vv in all_var) {
    if (subj == vec_subj[1]) {
      main <- gsub("log_", "", vv)
    } else {
      main <- ""
    }

    if (subj == vec_subj[length(vec_subj)]) {
      xlab <- "Days from symptom onset"
    } else {
      xlab <- ""
    }

    if (vv == all_var[1]) {
      ylim_offset <- 2.6
    } else {
      ylim_offset <- 1.75
    }

    if (any(df_comb$severity == "HC")) {
      lo <- quantile(df_comb[df_comb$severity == "HC", vv], probs = 0.25, na.rm = TRUE)
      up <- quantile(df_comb[df_comb$severity == "HC", vv], probs = 0.75, na.rm = TRUE)

      mmax <- max(sapply(Y[vec_subj], function(y_subj) max(y_subj[[vv]], na.rm = T)), up, na.rm = T)
      mmin <- min(sapply(Y[vec_subj], function(y_subj) min(y_subj[[vv]], na.rm = T)), up, na.rm = T)
    } else {
      mmax <- max(sapply(Y[vec_subj], function(y_subj) max(y_subj[[vv]], na.rm = T)))
      mmin <- min(sapply(Y[vec_subj], function(y_subj) min(y_subj[[vv]], na.rm = T)))
    }

    plot_trajectories_irregular_grid(subj, vv, Y, Y_hat, Y_upp, Y_low,
                                     days_thres, time_obs, time_g, df_comb,
                                     scores, bool_rescale_time,
                                     ylim_offset = ylim_offset, main = main,
                                     xlab = xlab, mfrow = F, bool_sub = F,
                                     mmin = mmin, mmax = mmax)
  }

}
if (bool_save) {
  dev.off()
}

vec_col_grey <- vec_col
vec_col_grey["A"] <- "grey90"
vec_col_grey["B"] <- "grey80"
vec_col_grey["C"] <- "grey65"
vec_col_grey["D"] <- "grey45"
vec_col_grey["E"] <- "grey5"


if (bool_save) {
  pdf(paste0(res_dir, "/scatterplot_scores_col_by_severity_grey.pdf"),
      width = 5, height = 5, paper='special')
}
scatterplot_scores(scores, var_expl = var_expl,  vec_col = vec_col_grey)
if (bool_save) {
  dev.off()
}

subj_dead <- scores$subject_id[scores$hospital_outcome %in% "dead"]
n_dead <- length(subj_dead)
set.seed(seed)
n_repl <- 1e5
p_emp_EF1 <- (sum(unlist(parallel::mclapply(1:n_repl, function(repl) {mean(sample(scores$EF1, n_dead)) >= mean(scores$EF1[scores$subject_id %in% subj_dead])},
                                            mc.cores = n_cpus))) + 1) / (n_repl + 1)
print(p_emp_EF1)

p_emp_EF2 <- (sum(unlist(parallel::mclapply(1:n_repl, function(repl) {mean(sample(scores$EF2, n_dead)) <= mean(scores$EF2[scores$subject_id %in% subj_dead])},
                                            mc.cores = n_cpus))) + 1) / (n_repl + 1)
print(p_emp_EF2)


for (vv in all_var) {
  if (bool_save) {
    pdf(paste0(res_dir, "/eigenfunctions_", vv, "_bw.pdf"),
        width =  4.65,
        height = 5.1, paper='special')
  }

  plot_eigenfunctions_irregular_grid(time_g, list_Psi_hat, vv, L = 2, var_expl = NULL, days_thres, bool_rescale_time = TRUE)
  if (bool_save) {
    dev.off()
  }
}


t_col <- function(color, percent = 50, name = NULL) {
 
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)

  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)

  ## Save the color
  invisible(t.col)
}

if (bool_save) {
  pdf(paste0(res_dir, "/eigenfunctions_all_var.pdf"),
      width =  13.9,
      height = 3.15, paper='special')
}
par(mfrow = c(1, length(all_var)))

for (vv in all_var) {
  if (vv == all_var[1]) {
    ylab <- "Eigenfunctions"
    bool_leg <- T
  } else {
    ylab <- ""
    bool_leg <- F
  }
  plot_eigenfunctions_irregular_grid(time_g, list_Psi_hat, vv, L = 2, var_expl = var_expl, bool_leg = bool_leg,
                                     days_thres, bool_rescale_time = TRUE, ylab = ylab,
                                     vec_col = c("dodgerblue1", t_col("dodgerblue1", percent = 70)), vec_lty = c(1, 1))
}
if (bool_save) {
  dev.off()
}



method <- c("t.test", "anova", "wilcox.test")[1]
bool_log10_display <- F
res_tests <- plot_boxplots(scores,
                           vec_var = "EF1",
                           group_var = "severity",
                           vec_var_to_log = NULL,
                           method = method,
                           main = "Scores FPC 1 by severity class",
                           xlab = "Severity class",
                           ylab = "Scores FPC 1",
                           vec_col = vec_col_grey,
                           save_path = res_dir,
                           bool_log10_disp = bool_log10_display,
                           bool_overall = F,
                           width = 4.6,
                           height = 3.7)


bool_ques_analysis <- T
if (bool_ques_analysis) {
  require(dplyr)
  # Obtain long-covid composite scores
  #
  df_ques_sub <-  unique(subset(df_ques_wide,
                                select = c("subject_id", "recorded_time")))
  df_info_sub <- unique(subset(df_info,
                               select = c("subject_id", "date_sx_or_swab")))
  df_merge <- left_join(df_ques_sub, df_info_sub, by = "subject_id")
  df_merge <- inner_join(df_merge, scores)

  df_merge$date_diff <- as.numeric(as.Date(format(as.POSIXct(as.character(df_merge$recorded_time), format='%d/%m/%Y %H:%M'), format='%Y-%m-%d')) - #as.character(df_merge$recorded_time), format='%m/%d/%Y %H:%M')-
                                     as.Date(as.character(df_merge$date_sx_or_swab), format="%Y-%m-%d")) / 365 * 12

  rownames(df_merge) <- paste0(df_merge$subject_id, "_",
                               format(as.POSIXct(as.character(df_merge$recorded_time), format='%d/%m/%Y %H:%M'),
                                      format='%d/%m/%Y'))

  df_ques_num_all_mth <- df_ques_num_all_mth[rownames(df_ques_num_all_mth) %in% rownames(df_merge),, drop = FALSE]
  mat_ques <- t(df_ques_num_all_mth)

  table(table(df_merge$subject_id)) # up to 2 questionnaires per subject
  hist(as.numeric(df_merge$date_diff))

  min(as.numeric(df_merge$date_diff))
  max(as.numeric(df_merge$date_diff))
  mean(as.numeric(df_merge$date_diff))

  ques_subject_ids <- sub("_[^_]+$", "", colnames(mat_ques))

  stopifnot(all.equal(rownames(mat_ques), df_ques_names$QuPos[match(rownames(mat_ques),
                                                                    df_ques_names$QuPos)]))

  rownames(mat_ques) <- df_ques_names$disp_name[match(rownames(mat_ques), df_ques_names$QuPos)]
  mat_ques <- mat_ques[order(rownames(mat_ques)),]
  mat_ques <- mat_ques[!(rownames(mat_ques) %in% "Other (specify in comment)"),]

  id_perc_ques <- which(rownames(mat_ques) == "If you have not made a full recovery at what % of your normal function are you currently at?")
  id_full_rec_ques <- which(rownames(mat_ques) == "Have you made a full physical and mental recovery form Covid-19?")

  mat_ques[id_perc_ques, mat_ques[id_full_rec_ques,] %in% 5 & is.na(mat_ques[id_perc_ques,])] <- 5
  rownames(mat_ques)[rownames(mat_ques) ==  "If you have not made a full recovery at what % of your normal function are you currently at?"] <- "Q0 physical and mental recovery"
  bool_sub_ques <- TRUE
  if (bool_sub_ques) {
    mat_ques <- mat_ques[rownames(mat_ques) %in% "Q0 physical and mental recovery"
                         | startsWith(rownames(mat_ques), "Q"), ]
    mat_ques <- mat_ques[order(as.numeric(gsub("Q", "", gsub("a)", "", gsub("b)", "", gsub("c)", "", gsub("d)", "", gsub( " .*$", "", rownames(mat_ques) )))))))), ]
  }

  mat_ques <- mat_ques[,order(ordered(sub("_[^_]+$", "", colnames(mat_ques)),
                                      levels = unique(scores$subject_id)))]

  ques_subject_ids <- sub("_[^_]+$", "", colnames(mat_ques))


  question_short <- rownames(mat_ques)
  names(question_short) <- rownames(mat_ques)

  id_questions <- c(1:4, 6:17) # excluding persisiting fever as no subject had it
  question_short[id_questions] <- c("Physical and mental recovery",
                                    "Dyspnoea",
                                    "Cough",
                                    "Chest Pain on exertion, palpitations or swollen ankles",
                                    # "Persisting fever (2 months or more)",
                                    "New leg swelling in one leg or shortness of breath with chest pain",
                                    "New skin rashes or sores",
                                    "Voice alteration",
                                    "Difficulties eating, drinking or swallowing", # (cough, choking, food avoidance)",
                                    "Constant noisy breathing or throat whistling",
                                    "Anosmia or dysgeusia",
                                    "Difficulty to gain or maintain weight, loss of appetite",
                                    "New neurology in one or more limbs",
                                    "New pain in one or more parts of the body",
                                    "General muscle weakness, balance or range of movement of joints",
                                    "Fatigue",
                                    "Cognition: memory, concentration and thinking skills")
  questions <- rownames(mat_ques)[id_questions]
  question_short <- question_short[id_questions]


  df_all <- as.data.frame(t(mat_ques))
  df_all <- abs(round(df_all - 5, digits=3)) # now large score is bad
  df_all$subject_id <- sub("_[^_]+$", "", rownames(df_all))


  df_all <- aggregate(df_all[,!(colnames(df_all) %in% "subject_id")],
                      list(df_all$subject_id),
                      mean)
  colnames(df_all)[colnames(df_all)== "Group.1"] <- "subject_id"

  df_all <- full_join(df_all, scores[!duplicated(scores$subject_id),])
  df_all$ques_score <- rowSums(df_all[, questions]) / length(questions)

  df_all_comp <- df_all
  df_all_comp$ques_score <- rowSums(df_all_comp[, questions]) / length(questions)

  vec_pval_EF1 <- vec_pval_EF2 <- vec_cor_EF1  <- vec_cor_EF2 <- NULL
  for (qq in names(question_short)) {

    print(qq)
    print("===================")
    print("EF1")
    getOption("na.action") # omits na.
    ct_qq_EF1 <- cor.test(df_all_comp$EF1,
                          df_all_comp[,qq],
                          method = "kendall")
    if (!is.na(ct_qq_EF1$p.value) && ct_qq_EF1$p.value < 0.05) print(ct_qq_EF1$p.value)
    vec_pval_EF1 <- c(vec_pval_EF1, ct_qq_EF1$p.value)
    vec_cor_EF1 <- c(vec_cor_EF1, ct_qq_EF1$estimate)


    print("EF2")
    ct_qq_EF2 <- cor.test(df_all_comp$EF2,
                          df_all_comp[,qq],
                          method = "kendall")
    if (!is.na(ct_qq_EF2$p.value) && ct_qq_EF2$p.value < 0.05) print(ct_qq_EF2$p.value)
    vec_pval_EF2 <- c(vec_pval_EF2, ct_qq_EF2$p.value)
    vec_cor_EF2 <- c(vec_cor_EF2, ct_qq_EF2$estimate)
  }
  names(vec_pval_EF1) <- names(vec_pval_EF2) <- names(vec_cor_EF1) <- names(vec_cor_EF2) <- question_short
  vec_cor_EF1 <- vec_cor_EF1[-1]
  vec_cor_EF2 <- vec_cor_EF2[-1]

  vec_adj_pval_EF1 <- p.adjust(vec_pval_EF1[-1], method = "BH") 
  vec_adj_pval_EF2 <- p.adjust(vec_pval_EF2[-1], method = "BH")

  vec_pval_EF1[1]
  vec_pval_EF2[1]

  tb_signif_sympt <- apply(cbind(vec_adj_pval_EF1, vec_adj_pval_EF2), 2, function(cc) sapply(cc, add_signif_label))
  tb_cor_sympt <- cbind(round(vec_cor_EF1,2), tb_signif_sympt[,1],
                        round(vec_cor_EF2, 2), tb_signif_sympt[,2])
  
  require(xtable)
  print(xtable(tb_cor_sympt))

  which(vec_adj_pval_EF1 < 0.05) 
  which(vec_adj_pval_EF2 < 0.05) 


  plot(df_all_comp$EF2, df_all_comp$`Q10 Any new pain in one or more parts of the body`,
       pch = 20, col = df_all_comp$col_severity,
       main = format(cor(df_all_comp$EF2, df_all_comp$`Q10 Any new pain in one or more parts of the body`,
                         use =  "complete.obs", method = "kendall"), digits = 3),
       xlab = "Score FPC 2", #  large - is bad
       ylab = "New pain in one or more parts of the body") 

  plot(df_all_comp$EF2, df_all_comp$`Q12 Fatigue`,
       pch = 20, col = df_all_comp$col_severity,
       main = format(cor(df_all_comp$EF2, df_all_comp$`Q12 Fatigue`,
                         use =  "complete.obs", method = "kendall"), digits = 3),
       xlab = "Score FPC 2", #  large - is bad
       ylab = "Fatigue") 


  ct_EF1 <- cor.test(df_all_comp$EF1,
                     df_all_comp$`Q0 physical and mental recovery`,
                     method = "kendall")#,
  ct_EF2 <- cor.test(df_all_comp$EF2,
                     df_all_comp$`Q0 physical and mental recovery`,
                     method = "kendall")

  n_complete_obs <- sum(!is.na(df_all_comp$EF2) & !is.na(df_all_comp$`Q0 physical and mental recovery`))
  cor_fpc_lc <- c(ct_EF1$estimate, ct_EF2$estimate)
  pval_fpc_lc <- c(ct_EF1$p.value, ct_EF2$p.value)
  names(cor_fpc_lc) <- names(pval_fpc_lc) <- c("Severity scores", "Recovery scores")


  if (bool_save) {
    save(cor_fpc_lc, pval_fpc_lc, n_complete_obs, file = paste0(out_dir,
                                                                "COVID_application/corr_irreg_grid_long_covid_",
                                                                gsub(" ", "_",
                                                                     paste0(all_var, collapse = "_")),
                                                                ifelse(subset_subj, "_subset_subj", ""), "_seed_", seed, ".RData"))
  }

  bool_table <- T
  if (bool_table) {

    vec_var <- c("log_CRP", "IL10", "GlycB", "Quinolinic acid", "Tryptophan")
    tb_cor <- tb_pval <- vec_n_complete_obs <- NULL
    for (vv in vec_var) {
      load(paste0(out_dir, "COVID_application/corr_irreg_grid_long_covid_",  gsub(" ", "_", vv),
                  ifelse(subset_subj, "_subset_subj", ""), "_seed_", seed, ".RData"))
      tb_cor <- rbind(tb_cor, cor_fpc_lc)
      tb_pval <- rbind(tb_pval, pval_fpc_lc)
      vec_n_complete_obs <- c(vec_n_complete_obs, n_complete_obs)
    }
    colnames(tb_cor) <- colnames(tb_pval) <- names(cor_fpc_lc)
    rownames(tb_cor) <- rownames(tb_pval) <- paste0(gsub("log_", "", vec_var), " FPCA, n = ", vec_n_complete_obs)

    load(paste0(out_dir, "COVID_application/corr_irreg_grid_long_covid_",gsub(" ", "_",
                                                                                 paste0(all_var, collapse = "_")),
                ifelse(subset_subj, "_subset_subj", ""), "_seed_", seed, ".RData"))
    tb_cor <- rbind(tb_cor, cor_fpc_lc)
    tb_pval <- rbind(tb_pval, pval_fpc_lc)
    vec_n_complete_obs <- c(vec_n_complete_obs, n_complete_obs)
    rownames(tb_cor)[length(vec_var)+1] <- rownames(tb_pval)[length(vec_var)+1] <- paste0("mFPCA, n = ", n_complete_obs)

    tb_signif <- apply(tb_pval, 2, function(cc) sapply(cc, add_signif_label))

    require(xtable)
    tb_cor_overall <- round(tb_cor,2)
    mat_overall <- cbind(tb_cor_overall[,1], tb_signif[,1],
                         tb_cor_overall[,2], tb_signif[,2])
    colnames(mat_overall) <- c("Severity scores", "", "Recovery scores", "")
    print(xtable(mat_overall))
  
  }
}

bool_surv <- T
bool_bayesian_surv <- T
if (bool_surv) { # try bayesian: indeptCoxph from library(spBayesSurv)

  df_surv <- right_join(df_comb, scores)

  df_surv <- get_early_or_late_samples(df_surv,
                                       vec_var = NULL,
                                       bool_early = FALSE,
                                       single_sample_per_subject = TRUE)$df_comb
  df_surv$time <- df_surv$hospital_outcome_days_from_start_sx_or_first_pos_swab
  df_surv$time[!(df_surv$hospital_outcome %in% "dead")] <- df_surv$last_observation_days_from_start_sx_or_first_pos_swab[!(df_surv$hospital_outcome %in% "dead")]

  df_surv$hospital_outcome[df_surv$severity %in% c("HC", "A", "B")] <- "alive"

  df_surv$status <- ifelse(df_surv$hospital_outcome == "dead", 1, 0) # 0 = censored, 1 = dead


  require(survival)
  require(survminer)

  aa <- table(df_surv$status, df_surv$time)
  time_death <- as.numeric(colnames(aa)[aa[2,]==1])

  if (bool_both) {
    
    if (bool_bayesian_surv) {
      require(brms)
      set.seed(seed) #123
      res_cox_both <- brm(time | cens(1 - status) ~ EF1 + EF2, # status should be coded 0 1 here not 2 1
                      data = df_surv, family = cox())
      res_sum_both <- summary(res_cox_both)
      print(res_sum_both)
      
    } else {
      res_cox_both <- coxph(Surv(time, status) ~ EF1 + EF2, data = df_surv) 
      res_sum_both <- summary(res_cox_both)
      print(res_sum_both)
      pval_surv_EF1 <- res_sum_both$coefficients[1, 5]
      pval_surv_EF2 <- res_sum_both$coefficients[2, 5]
      print(pval_surv_EF1)
      print(pval_surv_EF2)
    }


  } else {
    

    if (bool_bayesian_surv) {
      require(brms)
      set.seed(seed)
      
      res_cox1 <- brm(time | cens(1 - status) ~ EF1, # 1 + 
                     data = df_surv, family = cox()) # 2 -status and not 1 - status, because status coded 1 and 2, not 0 and 1
      res_sum1 <- summary(res_cox1)
      print(res_sum1)
      
      res_cox2 <- brm(time | cens(1 - status) ~ EF2, # 1 + 
                      data = df_surv, family = cox())#brmsfamily("cox"))
      res_sum2 <- summary(res_cox2)
      print(res_sum2)
   
    } else {
      res_cox1 <- coxph(Surv(time, status) ~ EF1, data = df_surv) 
      print(res_cox1)
      
      res_sum1 <- summary(res_cox1)
      pval_surv_EF1 <- res_sum1$coefficients[5]
      
      res_cox2 <- coxph(Surv(time, status) ~ EF2, data = df_surv)
      print(res_cox2)
      
      res_sum2 <- summary(res_cox2)
      pval_surv_EF2 <- res_sum2$coefficients[5]
      
      print(pval_surv_EF2)      
    }

  }

  n <- nrow(df_surv)

  if (bool_save & !bool_bayesian_surv) {
    save(pval_surv_EF1, pval_surv_EF2, n, subject_ids, file = paste0(out_dir,
                                                                     "COVID_application/surv_irreg_grid_long_covid_",
                                                                     gsub(" ", "_",
                                                                          paste0(all_var, collapse = "_")),
                                                                     ifelse(subset_subj, "_subset_subj", ""),
                                                                     ifelse(bool_both, "_joint_surv", ""),
                                                                     "_seed_", seed, ".RData"))
  }

  bool_table <- F
  if (bool_table) {

    vec_var <- c("log_CRP", "IL10", "GlycB", "Quinolinic acid", "Tryptophan")
    tb_surv_minlog10_pval <- tb_surv_pval <- vec_n <- NULL
    for (vv in vec_var) {
      load(paste0(out_dir, "COVID_application/surv_irreg_grid_long_covid_",  gsub(" ", "_", vv),
                  ifelse(subset_subj, "_subset_subj", ""),
                  ifelse(bool_both, "_joint_surv", ""), "_seed_", seed, ".RData"))

      tb_surv_pval <- rbind(tb_surv_pval, c(pval_surv_EF1, pval_surv_EF2))
      tb_surv_minlog10_pval <- rbind(tb_surv_minlog10_pval, c(-log10(pval_surv_EF1), -log10(pval_surv_EF2)))
      vec_n <- c(vec_n, n)
    }
    colnames(tb_surv_minlog10_pval) <- colnames(tb_surv_pval) <- c("Severity scores", "Recovery scores")
    rownames(tb_surv_minlog10_pval) <- rownames(tb_surv_pval) <- paste0(gsub("log_", "", vec_var), " FPCA, n = ", vec_n)

    load(paste0(out_dir, "COVID_application/surv_irreg_grid_long_covid_",gsub(" ", "_",
                                                                                 paste0(all_var, collapse = "_")),
                ifelse(subset_subj, "_subset_subj", ""), ifelse(bool_both, "_joint_surv", ""), "_seed_", seed, ".RData"))
    tb_surv_pval <- rbind(tb_surv_pval, c(pval_surv_EF1, pval_surv_EF2))
    tb_surv_minlog10_pval <- rbind(tb_surv_minlog10_pval, c(-log10(pval_surv_EF1), -log10(pval_surv_EF2)))
    vec_n <- c(vec_n, n)
    rownames(tb_surv_pval)[length(vec_var)+1] <- rownames(tb_surv_minlog10_pval)[length(vec_var)+1] <- paste0("mFPCA, n = ", n)

    tb_surv_signif <- apply(tb_surv_pval, 2, function(cc) sapply(cc, add_signif_label))

    require(gplots)
    require(RColorBrewer)
    pal <- colorRampPalette(c("blue", "white", "red"))(1000)
    br <- seq(-1, 1, length.out=1001)

    if (bool_save) {
      pdf(paste0(res_dir, "heatmap_correlation_with_survival", ifelse(subset_subj, "_subset_subj", ""), ".pdf"),
          width = 8, height = 6, paper='special')
    }
    par(cex.main=1)
    main_hp <- heatmap.2(as.matrix(tb_surv_minlog10_pval), #as.matrix(-log10(tb_surv_pval)), #as.matrix(tb_surv_cor),
                         dendrogram='none',
                         Rowv=FALSE, Colv=FALSE,
                         trace='none', scale="none",
                         xlab = "",
                         ylab = "",
                         cexCol = 1.2,
                         cexRow = 1,
                         key=T, keysize=0.4,
                         margins = c(10, 20),
                         key.title="",
                         cellnote = tb_surv_signif, #format(tb_surv_pval, digits = 1),# tb_surv_signif,
                         notecex=1.25,
                         notecol="black",
                         main = "Kendall's tau correlation between FPCA scores and long COVID composite score",
                         key.xlab="", key.ylab="", col=pal, density.info="none",
                         breaks=br,na.color="grey90")

    if (bool_save) {
      dev.off()
    }

    require(xtable)
    xtable(format(tb_surv_pval, digits = 1, scientific = F))

  }
}

