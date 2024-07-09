create_named_list <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}

get_early_or_late_samples <- function(df_comb, vec_var, bool_early, 
                                      trunc_days = NULL, 
                                      drop_HC = FALSE, nzf_thres = NULL,
                                      single_sample_per_subject = TRUE,
                                      enforce_subj_w_samples_avail_after_trunc_days = FALSE) {
  if (!bool_early) {
    stopifnot(!enforce_subj_w_samples_avail_after_trunc_days) # not used for late samples.
  }
  
  # get first sample from each subjects, if it is < trunc_days
  # if two samples from a same subject < trunc_days, only the first is kept.
  if (!is.null(vec_var)) {
    stopifnot(all(vec_var %in% colnames(df_comb)))
  } else {
    stopifnot(is.null(nzf_thres))
  }
  
  stopifnot(all(!is.na(df_comb$days_from_sx_or_swab_imputed[df_comb$severity != "HC"])))
  
  
  if (bool_early & enforce_subj_w_samples_avail_after_trunc_days) { 
    
    df_comb <- filter_out_subj_with_no_late_samples(df_comb, trunc_days)
    
  }
  
  if (bool_early) {
    df_comb <- df_comb[order(df_comb$days_from_sx_or_swab_imputed),]
  } else {
    df_comb <- df_comb[order(df_comb$days_from_sx_or_swab_imputed, decreasing = TRUE),]
  }
  
  if (single_sample_per_subject) {
    df_comb <- df_comb[!duplicated(df_comb$subject_id),, drop = FALSE]
  }
  
  if (!bool_early) { # reorder increasingly
    df_comb <- df_comb[order(df_comb$days_from_sx_or_swab_imputed),]
  }
  
  if (!is.null(trunc_days)) {
    
    if (bool_early) {
      bool_keep <- df_comb$days_from_sx_or_swab_imputed <= trunc_days
    } else {
      bool_keep <- df_comb$days_from_sx_or_swab_imputed > trunc_days
    }
    
  } else {
    bool_keep <- rep(TRUE, length(df_comb$days_from_sx_or_swab_imputed))
  }
  
  if (drop_HC) {
    bool_keep <- bool_keep & df_comb$severity != "HC"
  } else {
    bool_keep <- bool_keep | df_comb$severity == "HC"
  }
  
  df_sub_comb <- df_comb[bool_keep, , drop = FALSE]
  
  if (!is.null(vec_var)) {
    df_sub <- df_sub_comb[, vec_var, drop = FALSE]
    
    na_df_sub <- apply(df_sub, 2, function(vv) all(is.na(vv)))
    if (any(na_df_sub)) {
      na_var <- names(na_df_sub)[na_df_sub]
      warning(paste0("Removed all NA variable(s): ", na_var, "\n"))
      vec_var <- vec_var[!(vec_var %in% na_var)]
      df_sub_comb <- df_sub_comb[, !(colnames(df_sub_comb) %in% na_var), drop = FALSE] 
      df_sub <- df_sub[, vec_var, drop = FALSE] 
    }
    
    
    var_df_sub <- apply(df_sub, 2, function(vv) var(vv, na.rm = TRUE))
    if (any(is.na(var_df_sub)) | any(var_df_sub == 0)) {
      cst_var <- names(var_df_sub)[is.na(var_df_sub) | var_df_sub == 0]
      warning(paste0("Removed constant variable(s): ", cst_var, "\n"))
      vec_var <- vec_var[!(vec_var %in% cst_var)]
      df_sub_comb <- df_sub_comb[, !(colnames(df_sub_comb) %in% cst_var), drop = FALSE] 
      df_sub <- df_sub[, vec_var, drop = FALSE] 
    }
    
    if (!is.null(nzf_thres)) {
      list_exclude <- exclude_low_nzf(df_sub_comb, vec_var, nzf_thres)
      vec_var <- list_exclude$vec_var
      df_sub_comb <- list_exclude$df_comb
      df_sub <- df_sub[, vec_var, drop = FALSE] 
    }
  } else {
    df_sub <- vec_var <- NULL
  }
  # rownames(df_sub_comb) <- paste0(df_sub_comb$sample_id_v2, "_", df_sub_comb$replicate_id)
  
  list("df_comb" = df_sub_comb, "df" = df_sub, "vec_var" = vec_var)
}



prepare_data <- function(df_comb, vec_var, min_nb_timepoints = 2, 
                         days_thres = NULL,
                         bool_rescale_time = FALSE,
                         bool_rm_na  = TRUE,
                         bool_rm_cst = TRUE) {
  
  common_subjects <- Reduce("intersect", lapply(vec_var, function(vv) {
    
    df_sub <- df_comb[!is.na(df_comb[,vv]),, drop = FALSE]
    tb <- table(df_sub$subject_id)
    names(tb)[tb >= min_nb_timepoints] # subjects with enough timepoints
    
  })) # all subjects having at least two measurements for each variable.
  
  if (length(common_subjects)==0) {
    stop("Stop. No subjects common to all variables.")
  }
  
  out <- lapply(common_subjects, function(subj) {
    
    time_obs <- df_comb[df_comb$subject_id == subj, "days_from_sx_or_swab_imputed"]
    reorder_time <- order(time_obs)
    time_obs <- time_obs[reorder_time]
    if (bool_rescale_time) {
      time_obs <- time_obs / (days_thres-1) # 49 days is max
    }
    Y <- df_comb[df_comb$subject_id == subj, vec_var, drop = FALSE]
    Y <- Y[reorder_time, ,drop = FALSE]
    
    if (bool_rm_na) {
      bool_na <- rowSums(is.na(Y))>0 # at least one variable is NA for that timepoint
      Y <- Y[!bool_na,, drop = FALSE]
      time_obs <- time_obs[!bool_na]
    }
    
    if (bool_rm_cst) {
      bool_cst <- all(apply(Y, 2, var)==0)
      
      if (bool_cst) {
        Y <- time_obs <- NULL
      }

    }
    create_named_list(Y, time_obs)
  })
  names(out) <- common_subjects
  
  bool_empty <- sapply(out, function(out_subj) is.null(out_subj$Y))
  out <- out[!bool_empty]
  common_subjects <- common_subjects[!bool_empty]
  
  Y <- lapply(out, "[[", "Y")
  time_obs <- lapply(out, "[[", "time_obs")
  
  p <- length(vec_var)
  N <- length(common_subjects)
 
  N_t <- sapply(time_obs, length)
  
  create_named_list(N, N_t, p, Y, time_obs)
  
}



prepare_data_irregular <- function(df_comb, vec_var, min_nb_timepoints = 2, 
                                   days_thres = NULL,
                                   bool_rescale_time = FALSE) {
  
  common_subjects <- Reduce("intersect", lapply(vec_var, function(vv) {
    
    df_sub <- df_comb[!is.na(df_comb[,vv]),, drop = FALSE]
    tb <- table(df_sub$subject_id)
    names(tb)[tb >= min_nb_timepoints] # subjects with enough timepoints
    
  })) # all subjects having at least two measurements for each variable.
  
  if (length(common_subjects)==0) {
    stop("Stop. No subjects common to all variables.")
  }
  
  out <- lapply(common_subjects, function(subj) {
    
    time_obs <- df_comb[df_comb$subject_id == subj, "days_from_sx_or_swab_imputed"]
    reorder_time <- order(time_obs)
    time_obs <- time_obs[reorder_time]
    if (bool_rescale_time) {
      time_obs <- time_obs / (days_thres-1) # 49 days is max
    }
    
    Y_mat <- df_comb[df_comb$subject_id == subj, vec_var, drop = FALSE]
    Y_mat <- Y_mat[reorder_time, ,drop = FALSE]
    
    list_obj_subj <- lapply(vec_var, function(vv) {
      Y_vv <- Y_mat[, vv]
      bool_na <- is.na(Y_vv)
      Y_vv <- Y_vv[!bool_na]
      time_obs_vv <- time_obs[!bool_na]
      create_named_list(time_obs_vv, Y_vv)
    })
    names(list_obj_subj) <- vec_var
    
    Y_subj <- lapply(list_obj_subj, "[[", "Y_vv")
    time_obs_subj <- lapply(list_obj_subj, "[[", "time_obs_vv")
    
    create_named_list(Y_subj, time_obs_subj)
  })
  names(out) <- common_subjects
  Y <- lapply(out, "[[", "Y_subj")
  time_obs <- lapply(out, "[[", "time_obs_subj")
  
  p <- length(vec_var)
  N <- length(common_subjects)
  
  create_named_list(N, p, Y, time_obs)
  
}



choose_severity_groups <- function(df_comb, severity_groups) {
  
  stopifnot(all(severity_groups %in% c("HC", "A", "B", "C", "D", "E")))
  
  df_comb <- df_comb[df_comb$severity %in% severity_groups,, drop = FALSE]
  df_comb$severity <- factor(df_comb$severity, levels = severity_groups)
  
  print(paste0("Retained levels: ", paste0(levels(df_comb$severity), collapse = ", ")))
  print(table(df_comb$severity))
  
  df_comb
  
}



plot_eigenfunctions_v3 <- function(time_g, Psi_hat, var_name, var_expl, days_thres, 
                                   bool_rescale_time = TRUE, vec_col = rep("black", 2), vec_lty = c(1, 2)) {
  
  disp_var_name <- gsub("log_", "", var_name)
  fac <- ifelse(bool_rescale_time, days_thres-1, 1) # scaling factor back to the original time window
  # xlim_max <- ifelse(bool_rescale_time, days_thres, 1)
  
  Psi_hat_var <- Psi_hat[[var_name]] 
  
  lim_eigen <- c(min(c(Psi_hat_var, 0)),
                 max(c(Psi_hat_var, 0)))
  plot(time_g*fac, Psi_hat_var[,1], 
       xlab = "Days from symptom onset",
       ylab = "Eigenfunctions",
       main = disp_var_name,
       ylim = lim_eigen, type = "l", col = vec_col[1], lwd = 2, lty = vec_lty[1])
  abline(h = 0, col = "grey", lwd = 2, lty = 3)
  lines(time_g*fac, Psi_hat_var[,2], 
        type = "l", lwd = 2, col = vec_col[2], lty = vec_lty[2])

  legend("bottomleft", 
         col = vec_col,
         lwd = 2,
         lty = vec_lty,
         bty = "n",
         legend = paste0(c("1st: ", "2nd: "), #, "3rd"), 
                         format(var_expl*100, digits = 2), "%")) 

  
}


plot_eigenfunctions_irregular_grid <- function(time_g, list_Psi_hat, var_name, var_expl, days_thres, L = 2, bool_rescale_time = TRUE,
                                               vec_col = rep("black", L), vec_lty = 1:L, bool_leg = T, ylab = NULL) {
  
  disp_var_name <- gsub("log_", "", var_name)
  fac <- ifelse(bool_rescale_time, days_thres-1, 1) # scaling factor back to the original time window
  # xlim_max <- ifelse(bool_rescale_time, days_thres, 1)
  
  Psi_hat_var <- sapply(list_Psi_hat[1:L], function(Psi_hat) Psi_hat[,var_name]) 
  
  if (is.null(ylab)) {
    ylab <- "Eigenfunctions"
  }
  lim_eigen <- c(min(c(Psi_hat_var, 0)),
                 max(c(Psi_hat_var, 0)))
  plot(time_g*fac, Psi_hat_var[,1], 
       xlab = "Days from symptom onset",
       ylab = ylab,
       main = disp_var_name,
       ylim = lim_eigen, type = "l", col = vec_col[1], lwd = 2, lty = vec_lty[1])
  abline(h = 0, col = "grey", lwd = 2, lty = 3)
  for (l in 2:L) {
    lines(time_g*fac, Psi_hat_var[,l], 
          type = "l", lwd = l, col = vec_col[l], lty = vec_lty[l])
  }

  if (bool_leg) {
    if (is.null(var_expl)) {
      leg <- paste0("Eigenfunction ", 1:L)
    } else {
      leg <- paste0(c("1st: ", "2nd: "), 
                    format(var_expl[1:L], digits = 3), "%")
    }
    
    legend("topleft", 
           col = vec_col,
           lwd = 2,
           lty = vec_lty,
           bty = "n",
           legend = leg,
           inset=c(-0.15,-0.21), xpd=TRUE) 
    
  }
 
  
}



plot_trajectories <- function(vec_subj, var_name, Y, Y_hat, Y_upp, Y_low, 
                              days_thres, time_obs, time_g, df_comb, scores, bool_rescale_time, 
                              ylim_offset = 1, mfrow = TRUE) {
  
  disp_var_name <- gsub("log_", "", var_name)
  fac <- ifelse(bool_rescale_time, days_thres-1, 1) # scaling factor back to the original time window
  xlim_max <- ifelse(bool_rescale_time, days_thres, 1)
  
  
  if (any(df_comb$severity == "HC")) {
    lo <- quantile(df_comb[df_comb$severity == "HC", var_name], probs = 0.25, na.rm = TRUE)
    up <- quantile(df_comb[df_comb$severity == "HC", var_name], probs = 0.75, na.rm = TRUE)
    
    mmax <- max(sapply(Y[vec_subj], function(y_subj) max(y_subj[, var_name], na.rm = T)), up, na.rm = T)
    mmin <- min(sapply(Y[vec_subj], function(y_subj) min(y_subj[, var_name], na.rm = T)), up, na.rm = T)
  } else {
    mmax <- max(sapply(Y[vec_subj], function(y_subj) max(y_subj[, var_name], na.rm = T)))
    mmin <- min(sapply(Y[vec_subj], function(y_subj) min(y_subj[, var_name], na.rm = T)))
  }
  
  # par(mfrow = c(1, length(vec_subj)), pty="s")
  if (mfrow) {
    par(mfrow = c(2, 2), pty="s")
  }

  for (subj_id in vec_subj) {
    main <- paste0(disp_var_name, ", subj. ", subj_id, "\n (sev.: ", 
                   unique(df_comb$severity[df_comb$subject_id %in% subj_id]), 
                   ", gender: ", unique(df_comb$gender[df_comb$subject_id %in% subj_id]), 
                   ", age: ", unique(df_comb$age[df_comb$subject_id %in% subj_id]), 
                   ")")
    
    sub <- paste0("Scores, FPC 1: ", format(scores$EF1[scores$subject_id %in% subj_id], digits =2),
                    ", FPC 2: ", format(scores$EF2[scores$subject_id %in% subj_id], digits =2))
    
    time_subj_id <- fac*time_obs[[subj_id]]
    plot(time_subj_id, Y[[subj_id]][,var_name], pch = 1, 
         xlim = c(0, xlim_max), 
         ylim = c(max(0.5, mmin) - ylim_offset, mmax + ylim_offset),
         main = main,
         xlab = "Days from symptom onset",  
         ylab = paste0(disp_var_name, " (log-scale)"),
         sub = sub)
    if (any(df_comb$severity == "HC")) {
      rect(-5, lo, days_thres+5, up, col = adjustcolor("gray", alpha.f=0.5), border = "gray")
    }
    
    df_comb$time <- df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab
    df_comb$time[!(df_comb$hospital_outcome %in% "dead")] <- df_comb$last_observation_days_from_start_sx_or_first_pos_swab[!(df_comb$hospital_outcome %in% "dead")]
    df_comb$hospital_outcome[df_comb$severity %in% c("HC", "A", "B")] <- "alive"
    df_comb$status <- ifelse(df_comb$hospital_outcome == "dead", 2, 1) # 1 = cencored, 2 = dead
    
    status <- unique(df_comb$status[df_comb$subject_id %in% subj_id])
    time <- unique(df_comb$time[df_comb$subject_id %in% subj_id])
    time <- time[!is.na(time)]
    
    stopifnot(length(status) == 1 & length(time) == 1)
    
    time_grid_subj_id <- time_g
    if (bool_rescale_time) {
      time_grid_subj_id <- time_grid_subj_id * fac
    }
    if (status == 2 & time < days_thres) { # cut prediction if patient died within the analysis window
      lines(time_grid_subj_id[time_grid_subj_id<=time],
            Y_hat[[subj_id]][time_grid_subj_id<=time, var_name],
            col="black",lwd=2)
      points(time_grid_subj_id[which.min(abs(time_grid_subj_id-time))], 
             Y_hat[[subj_id]][which.min(abs(time_grid_subj_id-time)), var_name], 
             pch = 4, cex = 1.25, lwd = 1.5)
      if (!is.null(Y_upp)) {
        lines(time_grid_subj_id[time_grid_subj_id<=time], 
              Y_upp[[subj_id]][time_grid_subj_id<=time, var_name],
              col="black",lwd=2,lty=2)
      }
      if (!is.null(Y_low)) {
        lines(time_grid_subj_id[time_grid_subj_id<=time], 
              Y_low[[subj_id]][time_grid_subj_id<=time, var_name],
              col="black",lwd=2,lty=2)
      }
    } else {
      lines(time_grid_subj_id,
            Y_hat[[subj_id]][, var_name],
            col="black",lwd=2)
      if (!is.null(Y_upp)) {  
        lines(time_grid_subj_id, 
            Y_upp[[subj_id]][, var_name],
            col="black",lwd=2,lty=2)
      }
      if (!is.null(Y_low)){
      lines(time_grid_subj_id, 
            Y_low[[subj_id]][, var_name],
            col="black",lwd=2,lty=2)
      }
    }
  }
  
}


plot_trajectories_irregular_grid <- function(vec_subj, var_name, Y, Y_hat, Y_upp, Y_low, 
                              days_thres, time_obs, time_g, df_comb, scores, bool_rescale_time, 
                              ylim_offset = 1, main = NULL, xlab = NULL, 
                              mfrow = TRUE, bool_sub = TRUE, mmin = NULL, mmax = NULL) {
  
  disp_var_name <- gsub("log_", "", var_name)
  fac <- ifelse(bool_rescale_time, days_thres-1, 1) # scaling factor back to the original time window
  xlim_max <- ifelse(bool_rescale_time, days_thres, 1)
  
  if (is.null(mmin) | is.null(mmax)) {
    if (any(df_comb$severity == "HC")) {
      lo <- quantile(df_comb[df_comb$severity == "HC", var_name], probs = 0.25, na.rm = TRUE)
      up <- quantile(df_comb[df_comb$severity == "HC", var_name], probs = 0.75, na.rm = TRUE)
      
      mmax <- max(sapply(Y[vec_subj], function(y_subj) max(y_subj[[var_name]], na.rm = T)), up, na.rm = T)
      mmin <- min(sapply(Y[vec_subj], function(y_subj) min(y_subj[[var_name]], na.rm = T)), up, na.rm = T)
    } else {
      mmax <- max(sapply(Y[vec_subj], function(y_subj) max(y_subj[[var_name]], na.rm = T)))
      mmin <- min(sapply(Y[vec_subj], function(y_subj) min(y_subj[[var_name]], na.rm = T)))
    }
  }
 
  if (mfrow) {
    par(mfrow = c(2, 2), pty="s") 
  }
  for (subj_id in vec_subj) {
    
    if (is.null(main)) {
      main <- paste0(disp_var_name, ", subj. ", subj_id, "\n (sev.: ", 
                     unique(df_comb$severity[df_comb$subject_id %in% subj_id]), 
                     ", gender: ", unique(df_comb$gender[df_comb$subject_id %in% subj_id]), 
                     ", age: ", unique(df_comb$age[df_comb$subject_id %in% subj_id]), 
                     ")")
    }
    
    if (is.null(xlab)) {
      xlab <- "Days from symptom onset"
    }
    
    if (bool_sub) {
      sub <- paste0("Scores, FPC 1: ", format(scores$EF1[scores$subject_id %in% subj_id], digits =2),
                    ", FPC 2: ", format(scores$EF2[scores$subject_id %in% subj_id], digits =2))
    } else {
      sub <- ""
    }
     
    time_subj_id_vv <- fac*time_obs[[subj_id]][[var_name]]
    plot(time_subj_id_vv, Y[[subj_id]][[var_name]], pch = 1, 
         xlim = c(0, xlim_max), 
         ylim = c(max(0.5, mmin) - ylim_offset, mmax + ylim_offset),
         main = main,
         xlab = xlab,  
         ylab = paste0(disp_var_name, " (log-scale)"),
         sub = sub)
    if (any(df_comb$severity == "HC")) {
      rect(-5, lo, days_thres+5, up, col = adjustcolor("gray", alpha.f=0.5), border = "gray")
    }
    main <- NULL
    
    df_comb$time <- df_comb$hospital_outcome_days_from_start_sx_or_first_pos_swab
    df_comb$time[!(df_comb$hospital_outcome %in% "dead")] <- df_comb$last_observation_days_from_start_sx_or_first_pos_swab[!(df_comb$hospital_outcome %in% "dead")]
    df_comb$hospital_outcome[df_comb$severity %in% c("HC", "A", "B")] <- "alive"
    df_comb$status <- ifelse(df_comb$hospital_outcome == "dead", 2, 1) # 1 = cencored, 2 = dead
    
    status <- unique(df_comb$status[df_comb$subject_id %in% subj_id])
    time <- unique(df_comb$time[df_comb$subject_id %in% subj_id])
    time <- time[!is.na(time)]
    
    stopifnot(length(status) == 1 & length(time) == 1)
    
    time_grid_subj_id <- time_g
    if (bool_rescale_time) {
      time_grid_subj_id <- time_grid_subj_id * fac
    }
    if (status == 2 & time < days_thres) {
      lines(time_grid_subj_id[time_grid_subj_id<=time],
            Y_hat[[subj_id]][[var_name]][time_grid_subj_id<=time],
            col="black",lwd=2)
      points(time_grid_subj_id[which.min(abs(time_grid_subj_id-time))], 
             Y_hat[[subj_id]][[var_name]][which.min(abs(time_grid_subj_id-time))], 
             pch = 4, cex = 1.25, lwd = 1.5)
      if (!is.null(Y_upp)){
      lines(time_grid_subj_id[time_grid_subj_id<=time], 
            Y_upp[[subj_id]][[var_name]][time_grid_subj_id<=time],
            col="black",lwd=2,lty=2)
      }
      if (!is.null(Y_low)){
      lines(time_grid_subj_id[time_grid_subj_id<=time], 
            Y_low[[subj_id]][[var_name]][time_grid_subj_id<=time],
            col="black",lwd=2,lty=2)
      }
    } else {
      lines(time_grid_subj_id,
            Y_hat[[subj_id]][[var_name]],
            col="black",lwd=2)
      if (!is.null(Y_upp)){
      lines(time_grid_subj_id, 
            Y_upp[[subj_id]][[var_name]],
            col="black",lwd=2,lty=2)
      }
      if (!is.null(Y_low)){
      lines(time_grid_subj_id, 
            Y_low[[subj_id]][[var_name]],
            col="black",lwd=2,lty=2)
      }
    }
  }
  
}


add_signif_label <- function(pval) {
  
  if (is.na(pval)) {
    pval_label <- ""
  } else if (pval < 0.001) {
    pval_label <- "****"
  } else if (pval < 0.01) {
    pval_label <- "***"
  } else  if (pval < 0.05) {
    pval_label <- "**"
  } else if (pval < 0.1) {
    pval_label <- "*"
  } else if (pval < 0.15) {
    pval_label <- "."
  } else {
    pval_label <- ""
  }
  
  pval_label
}

scatterplot_scores <- function(scores, var_expl, vec_col, zeta_ellipse = NULL, offset = 0.2) {
  
  par(mfrow = c(1, 1), mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
  plot(scores$EF1, scores$EF2, col = scores$col_severity_grey, pch = 20, main = "Scores",
       xlab = ifelse(is.null(var_expl), "Scores FPC 1", paste0("Scores FPC 1 (var. expl. ", format(var_expl[1]*100, digits = 3), "%)")), 
       ylab = ifelse(is.null(var_expl), "Scores FPC 2", paste0("Scores FPC 2 (var. expl. ", format(var_expl[2]*100, digits = 3), "%)")),
       xlim = c(min(scores$EF1)-offset, max(scores$EF1)+offset),
       ylim = c(min(scores$EF2)-offset, max(scores$EF2)+offset))
 
  points(scores$EF1[!is.na(scores$hospital_outcome) & scores$hospital_outcome == "dead"], 
         scores$EF2[!is.na(scores$hospital_outcome) & scores$hospital_outcome == "dead"], 
         pch = 13, col = adjustcolor("black", alpha.f = 0.4), cex = 1)
  
  if (!is.null(zeta_ellipse)) {
    
    for (subj_id in scores$subject_id) {
      points(zeta_ellipse[[subj_id]][,1], zeta_ellipse[[subj_id]][,2], 
             col = adjustcolor(scores$col_severity_grey[scores$subject_id == subj_id], alpha.f = 0.2),
             type = "l")
    }
  }
  
}


get_df_scores <- function(subject_ids, Zeta_hat, df_info, vec_col) {
  
  scores <- data.frame("subject_id" = subject_ids, 
                       "EF1" = Zeta_hat[,1], 
                       "EF2" = Zeta_hat[,2])
  rownames(scores) <- scores$subject_id
  
  scores$severity <- df_info$severity[match(scores$subject_id, df_info$subject_id)]
  scores$col_severity <- vec_col[scores$severity]
  scores$hospital_outcome <- df_info$hospital_outcome[match(scores$subject_id, df_info$subject_id)]
  
  vec_col_grey <- vec_col
  vec_col_grey["A"] <- "grey90"
  vec_col_grey["B"] <- "grey80"
  vec_col_grey["C"] <- "grey65"
  vec_col_grey["D"] <- "grey45"
  vec_col_grey["E"] <- "grey5"
  
  scores$col_severity_grey <- vec_col_grey[scores$severity]
  
  scores
}


plot_boxplots <- function(df_comb, vec_var, group_var,
                          vec_var_to_log = NULL,
                          method = "t.test", 
                          method_overall = "anova",
                          vec_col = c("orange2", "purple3", "grey"), 
                          save_path = NULL,
                          bool_log10_disp = FALSE,
                          bool_overall = TRUE,
                          width = 3.6,
                          height = 4.1,
                          xlab = "",
                          main = NULL,
                          ylab = NULL) {
  
  stopifnot(all(vec_var %in% names(df_comb)))
  stopifnot(all(group_var %in% names(df_comb)))
  
  df_comb_sub <- df_comb[, c(group_var, vec_var), drop = FALSE]
  
  for (var in vec_var_to_log) {
    
    df_comb_sub[,var] <- log2(df_comb_sub[,var] + 1) # + 1 as CRP can be close to zero
    
  }
  
  
  df_comb_sub[, group_var] <- as.factor(df_comb_sub[, group_var])
  
  vec_var <- vec_var[colSums(!is.na(df_comb_sub[, colnames(df_comb_sub) %in% vec_var, drop = FALSE])) > 0]
  
  res_tests <- lapply(vec_var, function(var) {
    
    if (!is.null(save_path)) {
      pdf(paste0(save_path, "boxplots_", var, "_conditional_on_", group_var, "_", method, 
                 ifelse(bool_log10_disp, "_log_disp", ""), ".pdf"), 
          width = width, height = height, paper='special')
    }
    
    keep_levels <- sapply(levels(df_comb_sub[, group_var]), 
                          function(ll) sum(!is.na(df_comb_sub[df_comb_sub[, group_var] == ll, var])) > 0)
    df_comb_sub_var <- df_comb_sub[df_comb_sub[, group_var] %in% names(keep_levels)[keep_levels], ]
    df_comb_sub_var[, group_var] <- droplevels(df_comb_sub_var[, group_var])
    
    if (bool_log10_disp) {
      eps <- .Machine$double.eps
      stopifnot(all(df_comb_sub_var[,var] > -eps, na.rm = TRUE))
      y_var <- df_comb_sub_var[,var]
      df_comb_sub_var[y_var<eps & !is.na(y_var),var] <- y_var[y_var<eps & !is.na(y_var)] + eps^0.05
      # print(summary( df_comb_sub_var[y_var<eps & !is.na(y_var),var]))
    }
    ymin <- min(df_comb_sub_var[, var], na.rm = TRUE)
    ymax <- max(df_comb_sub_var[, var], na.rm = TRUE)
    
    require(ggpubr)
    
    bool_violin <- FALSE
    if (bool_violin) {
      p <- ggviolin(df_comb_sub_var, 
                    x = group_var, 
                    y = var,
                    size = 0.5,
                    color = group_var, 
                    palette =  vec_col,
                    add = "jitter", 
                    ylim = c(ifelse(bool_log10_disp, round_power_down(ymin), ymin), 
                             ifelse(bool_log10_disp, round_power_up(ymax), ifelse(length(levels(df_comb_sub_var[, group_var]))>2, 1.12, 1.05)*ymax)),
                    title = ifelse(is.null(main), tools::toTitleCase(gsub("_", " ", var)), main), 
                    xlab = group_var, 
                    ylab = ifelse(is.null(ylab), paste0(ifelse(var %in% vec_var_to_log, "log2-", ""),
                                                        tools::toTitleCase(gsub("_", " ", var))), ylab))
    } else {
      p <- ggboxplot(df_comb_sub_var, 
                     x = group_var, 
                     y = var,
                     size = 0.5,
                     color = group_var, 
                     palette =  vec_col,
                     add = "jitter", 
                     ylim = c(ifelse(bool_log10_disp, round_power_down(ymin), ymin), 
                              ifelse(bool_log10_disp, round_power_up(ymax), ifelse(length(levels(df_comb_sub_var[, group_var]))>2, 1.12, 1.05)*ymax)),
                     title = ifelse(is.null(main), tools::toTitleCase(gsub("_", " ", var)), main),  
                     xlab = paste0("\n\n ", xlab),  #xlab
                     ylab = ifelse(is.null(ylab), paste0(ifelse(var %in% vec_var_to_log, "log2-", ""),
                                                         tools::toTitleCase(gsub("_", " ", var))), ylab))
    }
    
    p <- p + rremove("legend")
    p <- p + font("title", size = 15, face = "bold") +  theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.line.x = element_line(color = "black", size = 0.1),
      axis.line.y = element_line(color = "black", size = 0.1))
    
    
    if (length(levels(df_comb_sub_var[, group_var]))>2) {
      
      if (bool_overall) {
        # Add global anova/kruskal wallis p-value
        p <- p + stat_compare_means(method = method_overall, label.x = 1.1, 
                                    label.y = 1.12 *max(df_comb_sub_var[, var], na.rm = TRUE))
      }
      p <- p + stat_compare_means(
        label = "p.signif", 
        method = method, 
        ref.group = ".all.",
        hide.ns = TRUE)
    } else {
      p <- p + stat_compare_means(method = method, 
                                  label.x = 1.25)
    }
    
    if (bool_log10_disp) {
      
      p <- p + coord_trans(y = "log10")
          
      rymin <- round_power_down(ymin)
      rymax <- round_power_up(ymax)
     
      p <- p + scale_y_continuous(
        breaks = c(0, 0.1, 1, 10, 100, 1000, 10000, 100000),
        labels = format(c(0, 0.1, 1, 10, 100, 1000, 10000, 100000), 
                        scientific = FALSE),
        limits = c(rymin, rymax))
     
    }
    
    
    print(p)
    
    if (!is.null(save_path)) {
      dev.off()
    }
    
    var <-  gsub("\\,", ".", gsub("\\(", ".",  gsub("\\)", ".", gsub(" ", "-", var)))) 
    names(df_comb_sub_var) <- gsub("\\,", ".",gsub("\\(", ".",  gsub("\\)", ".", gsub(" ", "-", names(df_comb_sub_var)))))
    print(names(df_comb_sub_var))
    if (sum(colSums(table(df_comb_sub_var[,var], df_comb_sub_var[,group_var]))>0)>1) { # at least two groups have values for var
      
      res_test <-  tryCatch({
        as.vector(compare_means(as.formula(paste0("`", var, "` ~ `", group_var, "`")), 
                                data = df_comb_sub_var, method = method))
      }, error = function(e) { 
        NULL })
    } else {
      res_test <- NULL
    }
    
    print(res_test)
    res_test
    
  })
  
  names(res_tests) <- vec_var
  res_tests
  
}
