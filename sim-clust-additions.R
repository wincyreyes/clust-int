# -- Required libraries --
library(kohonen)     # For SOM
library(PrEMIUM)    # For PrEMIUM

# ==== SOM CLUSTERING ====
som_clustering <- function(data, clus = 3, grid_dim = c(1,3), train_len = 100) {
  # Standardize
  data <- scale(data)
  # Grid for 3 clusters (1x3 is a linear strip of 3 nodes)
  som_grid <- somgrid(xdim = grid_dim[1], ydim = grid_dim[2], topo = "rectangular")
  som_mod <- som(data, grid = som_grid, rlen = train_len)
  # Assignment: winning neuron (convert row/col to label 1-3)
  assignments <- as.numeric(factor(apply(som_mod$unit.classif, 1, identity)))
  return(assignments)
}

# ==== PrEMIUM CLUSTERING ====
premium_clustering <- function(data, clus = 3) {
  data <- scale(as.matrix(data))
  fit <- premium(data, K = clus, verbose = FALSE)
  # Hard cluster assignment for each sample (maximum a posteriori)
  assignments <- apply(fit$Z, 1, which.max)
  return(assignments)
}

# ==== Integrate into results production ====
# In clus.inter: Insert after gmm
# Example (in the binary/continuous GLMs construction):
# ...
  ers_som <- som_clustering(data[,ersvars], clus=3)
  ers_premium <- premium_clustering(data[,ersvars], clus=3)
# ...
# And then add to GLMs (binary shown, repeat as for others):
# ...
  glm_som <- glm(y ~ grs_bin + factor(ers_som)+grs_bin*factor(ers_som), data, family=if(binary) 'binomial' else NULL)
  glm_premium <- glm(y ~ grs_bin + factor(ers_premium)+grs_bin*factor(ers_premium), data, family=if(binary) 'binomial' else NULL)
# ...
# Add these to results cbind/summarize, update colnames:
# ...
  summarize_glm_terms(glm_som, ...),
  summarize_glm_terms(glm_premium, ...)
# ...
# Update column names:
# c("median", ..., "gmm", "som", "premium")
#
# In simulation/test loops, do the same where cluster assignments are processed for test/train sets.
#
# Note: If train/test split, fit clustering models on train, predict on test (see how other clusterings are handled)

# ==== Add to file top: ====
# library(kohonen)
# library(PrEMIUM)

# ==== Helper code block inserted here for clarity, not immediately executable as is

# Please complete integration in simulation/analysis functions as shown, or reply if you want full code integration version with these edits applied to all relevant blocks (binary/cont, clus.inter, etc.)
