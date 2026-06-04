# ==============================================================================
# CR_GOF_Simulation_Scenario7.r
#
# Simulation scenario 7: Open population simulations.
#
# We evaluate model GOFs (for 3 incorrect closed population models) using DHARMa. 
# GOF residual plots are assessed visually and we use the KS and Dispersion 
# tests. 
#
# We repeat this 50 times.
#
# A full description of these data and their analysis appears in the 
# Supplementary Materials (Web Section C7) of the main text.
# ==============================================================================

library(DHARMa)

source("CR_GOF_Functions.R")

nsim <- 50   # Number of simulations.

# Set parameters.

N <- 250       # True population size.
tau <- 7       # Number of capture occasions.

phi <- 0.8  # Survival probabilities.

# IMPORTANT: Toggle below to change the distribution of X.

#dist.X <- "Normal"
dist.X <- "Binary"

beta <- c(0.5, 1)    # Intercept and covariate coef.

# Constant survival transition matrix.

M <- matrix(c(phi, 0, 1 - phi, 0, 0, 1, 0, 0, 1), 3, 3, byrow = TRUE)   

# Set up some functions.

fn.1 <- function(x) sum(x != 3)

# Assigns whether an individual survives based on M (time-invariant).

T.P <- function(x, M) {
  t.1 <- rmultinom(1, 1, M[x, ])
  which(t.1 > 0)
}

# Creates a (constant) survival history for one individual.

indiv.hist <- function(M, tau) {
  hist <- 1      
  for(k in 2:tau) {
    t.1 <- T.P(hist[k - 1], M)
    hist <- c(hist, t.1)
    if(t.1 == 3) break
  }
  
  t.2 <- length(hist)
  if(t.2 < tau) hist <- c(hist, rep(3, tau - t.2))
  hist
}

# Creates (constant) survival histories for all individuals.

pop.sim <- function(N.0, M, tau) {
  H <- NULL
  
  for(k in 1:N.0) {
    h.1 <- indiv.hist(M, tau)
    H <- rbind(H, h.1)
  }
  
  rownames(H) <- NULL
  H
}

# Creates a capture history matrix, then adds more individuals into the population by
# first calculating how many were missing on each occasion and then adding more.

pop.sim.all <- function(N.0, M, tau) {
  
  H <- pop.sim(N.0, M, tau)
  
  # This bit forces the population to be constant across each occasion. 
  
  for(j in 2:(tau - 1)) {
    D <- N.0 - fn.1(H[, j])     #  The number of missing individuals at time t.
    H.1 <- pop.sim(D, M, tau - j + 1)         # Constant.
    t.1 <- matrix(0, nrow(H.1), j - 1)
    t.2 <- cbind(t.1, H.1)
    H <- rbind(H, t.2)
  }
  
  D <- N.0 - fn.1(H[, tau])
  H.1 <- rep(1, D)
  t.1 <- matrix(0, length(H.1), (tau - 1))
  t.2 <- cbind(t.1, H.1)
  H <- rbind(H, t.2)
  colnames(H) <- NULL
  H
}

pop.sim.obs <- function(M, tau, N.0) {
  H <- pop.sim.all(N.0, M, tau + 1)
  H <- H[, -1]            # Remove first column.
  H1 <- H
  H1[H1 != 1] <- 0          # Change all 3s to 0s.
  H2 <- H
  H2[(H2 == 1) | (H2 == 2)] <- 1 
  H2[H2 == 3] <- 0
  N.all <- apply(H2, 2, sum)    # Population size.
  N.obs <- apply(H1, 2, sum)    # Observable population.
  hist <- H1
  
  list(hist = hist, H = H, N.all = N.all, N.obs = N.obs)
}

folder_name <- "Sim7_GOF"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

set.seed(123)

for (ii in 1:nsim) {

  sim.out <- pop.sim.obs(M, tau, N) # H contains status, 1 = obs and 3 = removed.
  hist <- sim.out$hist              # Present in population.
  
  NN <- nrow(hist)
  
  # IMPORTANT: Toggle here to change dist. of X.
  
  if (dist.X == "Normal") x_cov <- rnorm(NN, 0, 1)      # X ~ Normal covariate.
  if (dist.X == "Binary") x_cov <- rbinom(NN, 1, 0.5)   # X ~ Binary covariate.
  
  X <- cbind(rep(1, NN), x_cov)
  
  p <- plogis(drop(X%*%beta))
  
  # Capture functions.  
  
  cap.sim.0 <- function(x, p) {
    n <- sum(x)
    c <- rbinom(n, 1, p)
    x[x == 1] <- c
    x
  }  
  
  captures <- function(hist, p, X) {
    hist.out1 <- c()
    
    for(i in 1:nrow(hist)) {
      hist.out1 <- rbind(hist.out1, (cap.sim.0(hist[i, ], p = p[i])))
    }
    cap.hist.0 <- hist.out1
    
    
    y <- apply(cap.hist.0, 1, sum)
    cap.hist <- cap.hist.0[y > 0, ]
    X.obs <- X[y > 0, ]
    
    list(cap.hist = cap.hist, cap.hist.0 = cap.hist.0, X.obs = X.obs)
  }
  
  # Create capture histories here.
  
  cap.out <- captures(hist, p, X)
  cap.hist00 <- cap.out$cap.hist 
  cap.hist0 <- matrix(c(cap.hist00), nrow = nrow(cap.hist00), 
                     ncol = ncol(cap.hist00), byrow = FALSE)
  
  # Extract capture data.
  
  cap.hist <- cap.hist0[rowSums(cap.hist0) != 0, ]  # Remove individuals never captured.
 
  D <- nrow(cap.hist)               # Unique number of individuals captured.
  
  x.obs <- as.matrix(cap.out$X.obs[, 2])
  colnames(x.obs) <- c("x.obs")
  
  # Construct (Binomial) recapture-based quantities under the partial likelihood 
  # for fitting M_h-only models.
  
  pl_data <- prepare_recap_data_Mh(cap.hist, x.obs, tau)
  
  # Fit incorrectly specified (closed population) M_h-type (GLM) models here.
  
  # Incorrectly specified: Model with one covariate.
  
  model1 <- glm(cbind(Y_recap, n_opp - Y_recap) ~ x.obs,
                family = binomial(link = "logit"), data = pl_data)
  
  residuals_output1 <- simulateResiduals(model1, n = 1000)  
  
  # Construct (Bernoulli) recapture-based quantities under the partial likelihood
  # for fitting M_t- and M_b-type models.
  
  pl_data <- prepare_recap_data_Mtbh(as.matrix(cap.hist), x.obs, tau)
  
  # Fit incorrectly specified: closed population models: Models M_bh and M_th.
  
  # Incorrectly specified: Fit a M_bh-type GLM (accounts for behav. response and hete.).
  
  mod_mbh <- glm(Y_recap ~ x.obs + cap_prev, family = binomial(link = "logit"), 
                 data = pl_data)
  
  residuals_output2 <- simulateResiduals(mod_mbh, n = 1000)
  res_grouped2 <- recalculateResiduals(residuals_output2, group = pl_data$ID)
  
  # Incorrectly specified: Fit a M_th-type GLM (accounts for time effects and hete.).
  
  mod_mth <- glm(Y_recap ~ x.obs + occasion, family = binomial(link = "logit"), 
                 data = pl_data)
  
  residuals_output3 <- simulateResiduals(mod_mth, n = 1000)
  res_grouped3 <- recalculateResiduals(residuals_output3, group = pl_data$ID)
  
  # DHARMa plots: Side-by-side plots.
  
  figname <- file.path(folder_name, paste0("sim7_GOF_plots_dataset_", ii, ".png"))
  
  png(file = figname, width = 11, height = 9, units = "in", res = 300)
  
  par(mfrow = c(3, 1))
  
  plotQQunif(residuals_output1, main = paste("Dataset_", ii, ": Model 1 (M_h)"))
  plotQQunif(res_grouped2, main = paste("Dataset_", ii, ": Model 2 (M_th)"))
  plotQQunif(res_grouped3, main = paste("Dataset_", ii, ": Model 3 (M_bh)"))
  
  dev.off()
}

