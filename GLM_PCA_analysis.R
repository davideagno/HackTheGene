
# Response: matrix, dimension OTU x Samples (2771x2774)
Y_or <- read.csv("counts.csv", row.names = 1)
# Transpose it for Balsamico library
Y <- t(Y_or)

# Missing data
Y_oracle <- read.csv("counts_complete.csv", header = FALSE)

# Indices of NA values in orignial response matrix (i.e. test set)
ind_na <- which(is.na(Y_or), arr.ind = TRUE)
# Indices of observed values in original response matrix 
ind_obs <- which(!is.na(Y_or), arr.ind = TRUE)
# Split the dataset in train set and validation set
set.seed(123)
train <- sample(1:NROW(ind_obs), NROW(ind_obs)*0.9)
valid <- setdiff(1:NROW(ind_obs), train)
ind_train <- ind_obs[train,]
ind_valid <- ind_obs[valid,]
# Reconstruct Y_train matrix
Y_train_or <- matrix(NA, nrow = NROW(Y_or), ncol = NCOL(Y_or))
y_train_array <- Y_or[ind_train]
for (i in 1:NROW(ind_train)) {
  Y_train_or[ind_train[i,1], ind_train[i,2]] <- y_train_array[i]
}
Y_train <- t(Y_train_or)
mse <- function(pred, obs) sum((pred - obs)^2) / length(pred)

# BiocManager::install("pcaMethods", force = TRUE)
library(pcaMethods)
N <- NROW(Y_train)
P <- NCOL(Y_train)
PC_val <- 1:10
mse_glmpca <- rep(NA, length(PC_val))
out_glmpca_scores <- list()
out_glmpca_loadings <- list()
res <- matrix(NA, N, P)
mu <- colMeans(Y_train, na.rm = TRUE)
mu[mu == 0] <- 1e-9
for (i in 1:P) {
  res[,i] <- (Y_train[,i] - mu[i]) / sqrt(mu[i])
}
for (index in 10:length(PC_val)) {
  p <- PC_val[index]
  print(p)
  res_pca <- ppca(res, nPcs = p, seed = 123, scale = TRUE, center = TRUE)
  out_glmpca_scores[[index]] <- res_pca@scores
  out_glmpca_loadings[[index]] <- res_pca@loadings
  new_res <- res_pca@scores %*% t(res_pca@loadings)
  pred <- matrix(NA, nrow = N, ncol = P)
  for (i in 1:P) {
    pred[,i] <- new_res[,i]*sqrt(mu[i]) + mu[i]
  }
  pred_valid <- t(pred)[ind_valid]
  true_valid <- Y_oracle[ind_valid]
  mse_glmpca[index] <- mse(pred_valid, true_valid)
  print(mse_glmpca[index])
}
cbind(PC_val, mse_glmpca)
plot(PC_val, mse_glmpca, type = "l")

# Best
Y_obs <- matrix(NA, nrow = NROW(Y_or), ncol = NCOL(Y_or))
y_obs_array <- Y_or[ind_obs]
for (i in 1:NROW(ind_obs)) {
  Y_obs[ind_obs[i,1], ind_obs[i,2]] <- y_obs_array[i]
}
Y_obs <- t(Y_obs)
N <- NROW(Y_obs)
P <- NCOL(Y_obs)
res <- matrix(NA, N, P)
mu <- colMeans(Y_obs, na.rm = TRUE)
mu[mu == 0] <- 1e-9
for (i in 1:P) {
  res[,i] <- (Y_obs[,i] - mu[i]) / sqrt(mu[i])
}
res_pca <- ppca(res, nPcs = 2, seed = 123, scale = TRUE, center = TRUE)
new_res <- res_pca@scores %*% t(res_pca@loadings)
pred <- matrix(NA, nrow = N, ncol = P)
for (i in 1:P) {
  pred[,i] <- new_res[,i]*sqrt(mu[i]) + mu[i]
}
pred_na <- t(pred)[ind_na]
true_na <- Y_oracle[ind_na]
mse(pred_na, true_na)


