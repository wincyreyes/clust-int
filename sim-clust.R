# Set library paths
.libPaths("/home/wincy.reyes@IUF.LAN/R/x86_64-pc-linux-gnu-library/4.3")
library("scrime")
library("MASS")
library("glmnet")
library("GMCM")
library("scales")
library("dplyr")
library("tidyr")
library("purrr")
library(mclust)
library("ks")
library("ggplot2")
library(kohonen)   # For SOM
library(PrEMIUM)  # For PrEMIUM
library(mvtnorm)  # Needed for PrEMIUM test assignment

# Helper: Train SOM and assign clusters for train and test
tain_test_som <- function(train_data, test_data, clus = 3, grid_dim = c(1,3), train_len = 100) {
  train_data <- scale(train_data)
  test_data <- scale(test_data)
  som_grid <- somgrid(xdim = grid_dim[1], ydim = grid_dim[2], topo = "rectangular")
  som_mod <- som(train_data, grid = som_grid, rlen = train_len)
  som_train <- as.numeric(factor(as.vector(som_mod$unit.classif)))
  som_test <- as.numeric(factor(predict(som_mod, newdata = test_data)$unit.classif))
  list(train = som_train, test = som_test)
}

# Helper: Train PrEMIUM and assign clusters for train and test
train_test_premium <- function(train_data, test_data, clus = 3) {
  train_data <- scale(as.matrix(train_data))
  test_data <- scale(as.matrix(test_data))
  fit <- premium(train_data, K = clus, verbose = FALSE)
  means <- fit$mu
  covs <- fit$sigma
  weights <- fit$pi
  logpdf <- function(x, mu, Sigma) {
    mvtnorm::dmvnorm(x, mean=mu, sigma=Sigma, log=TRUE)
  }
  zmat <- sapply(1:clus, function(k) {
    apply(test_data, 1, logpdf, mu=means[,k], Sigma=covs[,,k]) + log(weights[k])
  })
  premium_test <- max.col(zmat)
  train_assign <- apply(fit$Z, 1, which.max)
  list(train = train_assign, test = premium_test)
}

# ==== SIMULATION PARTIAL: Place inside your test/train split simulation loop, e.g. after ersmat_train, ersmat_test defined ====
#ersmat_train <- ENVdata[train_idx, ]
#ersmat_test <- ENVdata[test_idx, ]
# Add these lines:
out_som <- tain_test_som(ersmat_train, ersmat_test, clus=3)
ers_som_train <- out_som$train
ers_som_test <- out_som$test
out_premium <- train_test_premium(ersmat_train, ersmat_test, clus=3)
ers_premium_train <- out_premium$train
ers_premium_test <- out_premium$test
# In your GLM/results:
# Use ers_som_test, ers_premium_test as additional clusterings for test set analyses
# Use ers_som_train, ers_premium_train for training set evaluations if needed

# Update clus.inter, clus.inter2, and column/result bindings accordingly
# Add to cbind, result tables, and colnames ("som", "premium")
