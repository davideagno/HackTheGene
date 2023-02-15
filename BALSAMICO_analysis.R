
# Covariates for samples: dataframe, 2774 obs, 7 variables
cc <- read.csv("covariates_col.csv")
# Discard ID_Sample, Body_Site, Sample_SRS
cc <- cc[, -c(1,5,7)]
X <- model.matrix(~ cc$N_Visita + cc$Sesso + cc$Centro_Seq + cc$Sottosito_Corporeo)

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
train <- sample(1:NROW(ind_obs), NROW(ind_obs)*0.7)
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

library(BALSAMICO)
N <- NROW(Y_train)
K <- NCOL(Y_train)
L_val <- 2:7
mse_balsamico <- rep(NA, length(L_val))
out_balsamico_U <- list()
out_balsamico_V <- list()
for (index in 1:length(L_val)) {
  l <- L_val[index]
  print(l)
  res <- VNMF(Y_train, X = X, L = l, maxit = 1e5)
  out_balsamico_U[[index]] <- res$W
  out_balsamico_V[[index]] <- res$H
  pred <- res$W %*% res$H
  pred_valid <- t(pred)[ind_valid]
  true_valid <- Y_oracle[ind_valid]
  mse_balsamico[index] <- mse(pred_valid, true_valid)
}
cbind(L_val, mse_balsamico)

# Best
Y_obs <- matrix(NA, nrow = NROW(Y_or), ncol = NCOL(Y_or))
y_obs_array <- Y_or[ind_obs]
for (i in 1:NROW(ind_obs)) {
  Y_obs[ind_obs[i,1], ind_obs[i,2]] <- y_obs_array[i]
}
Y_obs <- t(Y_obs)
res <- VNMF(Y_obs, X = X, L = 4, maxit = 1e5)
pred <- res$W %*% res$H
pred_na <- t(pred)[ind_na]
true_na <- Y_oracle[ind_na]
mse(pred_na, true_na)

