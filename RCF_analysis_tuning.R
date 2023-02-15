
library(lsa)

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
ind_known <- ind_obs
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
counts_train <- Y_train_or
counts <- Y_or
ind_val <- ind_valid
n <- dim(counts)[2]
rm(y_train_array, Y_or, Y_train_or)
answers_valid <- counts[ind_val]

na_train <- list()
for (j in 1:n) {
  na_train[[j]] <- which(is.na(counts_train[,j]))
}

# Column sd
csds_train <- as.numeric(matrixStats::colSds(as.matrix(counts_train), na.rm = TRUE))

# Rescaled train matrix
rescaled_counts_train <- counts_train
for (j in 1:n) {
  rescaled_counts_train[,j] <- rescaled_counts_train[,j] / csds_train[j]
}

# Cosine similarities (rescaled = original, normalized = centered)
csims_train <- matrix(0, nrow = n, ncol = n)
for(i in 1:n) {
  for (j in i:n) {
    nas_train <- unique(sort(c(na_train[[i]], na_train[[j]])))
    a <- counts_train[-nas_train,i]
    b <- counts_train[-nas_train,j]
    csims_train[i,j] <- cosine(a,b)
    csims_train[j,i] <- csims_train[i,j]
  }
}

# Error function
mse <- function(output, answers)
{
  err <- 0
  for (i in 1:length(answers)) {
    err <- err + (answers[i] - output[i])^2
    if (is.na(err))
    {
      print(i)
      break
    }
  }
  err <- err / length(answers) 
  return(err)
}

K_val <- 2:25
mse_counts <- rep(0, length(K_val))

for (index in 1:length(K_val)) {
  K <- K_val[index]
  print(K)
  top_csims_train <- apply(csims_train, 2, function(x) order(x, decreasing = TRUE)[1:(K+1)])
  final_counts_train <- counts_train
  final_rescaled_counts_train <- rescaled_counts_train
  for (j in 1:n) {
    nearest_users_train <- top_csims_train[-1,j]
    near_values_train <- final_counts_train[na_train[[j]], nearest_users_train]
    near_rescaled_values_train <- final_rescaled_counts_train[na_train[[j]], nearest_users_train]
    weight_csims_train <- csims_train[nearest_users_train,j]
    for (i in 1:length(na_train[[j]])) {
      final_counts_train[na_train[[j]][i],j] <- weighted.mean(near_values_train[i,], weight_csims_train, na.rm = TRUE)
    }
  }
  final_counts_train[,c(284, 1062, 1246, 2312)] <- 0
  error_nas <- which(is.na(final_counts_train), arr.ind = TRUE)
  if(dim(error_nas)[1] > 0)
  {
    for (i in 1:dim(error_nas)[1]) {
      final_counts_train[error_nas[i,1],error_nas[i,2]] <- mean(final_counts_train[error_nas[i,1],][!is.na(final_counts_train[error_nas[i,1],])])
    }
  }
  out_orig <- final_counts_train[ind_val]
  mse_counts[index] <- mse(out_orig, answers_valid)
}

cbind(K_val, mse_counts)
plot(K_val, mse_counts, type = "l")

