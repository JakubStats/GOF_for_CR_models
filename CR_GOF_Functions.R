prepare_recap_data_Mh <- function(cap_hist, covariates = NULL, 
                                  tau = ncol(cap_hist)) {
  
  if (any(rowSums(cap_hist) == 0)) {
    stop("Every row of cap_hist must contain at least one capture.")
  }
  if (!is.matrix(cap_hist) && !is.data.frame(cap_hist)) {
    stop("cap_hist must be a matrix or data frame.")
  }
  if (!is.null(covariates) && nrow(as.data.frame(covariates)) != nrow(cap_hist)) {
    stop("covariates must have one row/value per individual.")
  }
  cov_name <- deparse(substitute(covariates))
  
  # First capture (t1).
  t1 <- apply(cap_hist, 1, function(x) match(1, x)) 
  keep_idx <- which(t1 != tau)
  total_caps <- rowSums(cap_hist[keep_idx, , drop = FALSE])
  # Binomial response: Total captures minus the mandatory first capture.
  Y_recap <- total_caps - 1
  # Binomial size: Total occasions minus the occasion of first capture.
  n_opp <- tau - t1[keep_idx]

  df <- data.frame(
    Y_recap = Y_recap, 
    n_opp = n_opp, 
    ID = keep_idx
  )

  if (!is.null(covariates)) {
    if (is.matrix(covariates) || is.data.frame(covariates)) {
      X_recap <- covariates[keep_idx, , drop = FALSE]
      df <- cbind(df, X_recap)
    } else {
      df[[cov_name]] <- covariates[keep_idx]
    }
  }
  
  return(df)
}

prepare_recap_data_Mtbh <- function(cap_hist, covariates = NULL, 
                                    tau = ncol(cap_hist)) {
  if (any(rowSums(cap_hist) == 0)) {
    stop("Every row of cap_hist must contain at least one capture.")
  }
  if (!is.matrix(cap_hist) && !is.data.frame(cap_hist)) {
    stop("cap_hist must be a matrix or data frame.")
  }
  if (!is.null(covariates) && nrow(as.data.frame(covariates)) != nrow(cap_hist)) {
    stop("covariates must have one row/value per individual.")
  }
  cov_name <- deparse(substitute(covariates))
  
  # First capture occasion for each individual (f).
  f <- apply(cap_hist, 1, function(x) which(x == 1)[1])
  
  data_list <- lapply(which(f < tau), function(i) {
    idx_after <- (f[i] + 1):tau
    # Bernoulli response: The 0/1 outcomes after first capture.
    Y_response <- cap_hist[i, idx_after]
    # Previous capture status (used for behavioural effect).
    cap_prev <- cap_hist[i, f[i]:(tau - 1)]
    
    ind_df <- data.frame(
      ID = i, 
      Y_recap = Y_response, 
      cap_prev = cap_prev,
      occasion = factor(idx_after, levels = 2:tau)
    )
  
    if (!is.null(covariates)) {
      if (is.matrix(covariates) || is.data.frame(covariates)) {
        cov_rep <- covariates[i, , drop = FALSE]
        ind_df <- cbind(ind_df, cov_rep[rep(1, length(idx_after)), , drop = FALSE])
      } else {
        ind_df[[cov_name]] <- covariates[i]
      }
    }
    return(ind_df)
  })
  
  pl_data <- do.call(rbind, data_list)
  rownames(pl_data) <- NULL
  return(pl_data)
}

absorb <- function(m, k, D = nrow(m)/k) {
  m.dot <- matrix(0, D, ncol(m))
  m.a <- array(t(m), dim = c(ncol(m), k, D))
  for (j in 1:D) {
    if (is.vector(m.a[, , j])) m.dot[j, ] <- sum(m.a[, , j])
    else m.dot[j, ] <- apply(m.a[, , j], 1, sum)
  }
  m.dot
}