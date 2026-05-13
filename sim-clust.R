# Set library paths
.libPaths("/home/wincy.reyes@IUF.LAN/R/x86_64-pc-linux-gnu-library/4.3")
# install.packages("purrr")
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

# ... EXISTING FUNCTIONS ...

# ==== SOM CLUSTERING FUNCTION ====
som_clustering <- function(data, clus = 3, grid_dim = c(1,3), train_len = 100) {
  data <- scale(data)
  som_grid <- somgrid(xdim = grid_dim[1], ydim = grid_dim[2], topo = "rectangular")
  som_mod <- som(data, grid = som_grid, rlen = train_len)
  assignments <- as.numeric(factor(as.vector(som_mod$unit.classif)))
  return(assignments)
}

# ==== PrEMIUM CLUSTERING FUNCTION ====
premium_clustering <- function(data, clus = 3) {
  data <- scale(as.matrix(data))
  fit <- premium(data, K = clus, verbose = FALSE)
  assignments <- apply(fit$Z, 1, which.max)
  return(assignments)
}

# ========== Integrate new methods into clus.inter (and clus.inter2) =============
# Add this code block inside clus.inter and clus.inter2 after ers_gmm assignment (ers_gmm <- ...):
#
#   ers_som <- som_clustering(data[,ersvars], clus=3)
#   ers_premium <- premium_clustering(data[,ersvars], clus=3)
#
# Then, in the GLM formulas, add:
#   glm_som <- glm(y ~ grs_bin + factor(ers_som) + grs_bin*factor(ers_som), data, family=if(binary) 'binomial' else NULL)
#   glm_premium <- glm(y ~ grs_bin + factor(ers_premium) + grs_bin*factor(ers_premium), data, family=if(binary) 'binomial' else NULL)
#
# And in your results cbind:
#   ... summarize_glm_terms(glm_gmm, ...), summarize_glm_terms(glm_som, ...), summarize_glm_terms(glm_premium, ...))
#
# And update the column names:
#   colnames(results) <- c("median","tertile","quartile","mid80","mid90","clus","wclus","gmm","som","premium")
#
# ========== For train-test simulation/test blocks =============
# Wherever you assign clusters to test data:
#   ers_som_test <- som_clustering(ENVdata[test_idx, ], clus=3)
#   ers_premium_test <- premium_clustering(ENVdata[test_idx, ], clus=3)
#
# And pass them to the evaluation/modeling steps accordingly.

# =================================================================
# Note: Only additions, no deletions! Integrate per function blocks as above.
# Full code not overwritten here for clarity and safe merging. Insert marked sections per instructions.
# =================================================================
