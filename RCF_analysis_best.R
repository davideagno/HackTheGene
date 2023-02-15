
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

Y_obs_or <- matrix(NA, nrow = NROW(Y_or), ncol = NCOL(Y_or))
y_obs_array <- Y_or[ind_obs]
for (i in 1:NROW(ind_obs)) {
  Y_obs_or[ind_obs[i,1], ind_obs[i,2]] <- y_obs_array[i]
}
Y_obs <- t(Y_obs_or)
rm(y_obs_array)
counts_obs <- Y_obs_or
n <- dim(counts_obs)[2]

na_obs <- list()
for (j in 1:n) {
  na_obs[[j]] <- which(is.na(counts_obs[,j]))
}

# column mean and sd
cmeans_obs <- as.numeric(colMeans(counts_obs, na.rm = TRUE))
csds_obs <- as.numeric(matrixStats::colSds(as.matrix(counts_obs), na.rm = TRUE))

# centered, rescaled and normalized matrices
centered_counts_obs <- counts_obs
rescaled_counts_obs <- counts_obs
normalized_counts_obs <- counts_obs
for (j in 1:n) {
  centered_counts_obs[,j] <- centered_counts_obs[,j] - cmeans_obs[j]
  rescaled_counts_obs[,j] <- rescaled_counts_obs[,j] / csds_obs[j]
  normalized_counts_obs[,j] <- (normalized_counts_obs[,j] - cmeans_obs[j]) / csds_obs[j]
} 

# Cosine similarities (rescaled = original, normalized = centered)
csims_obs <- matrix(0, nrow = n, ncol = n)
centered_csims_obs <- matrix(0, nrow = n, ncol = n)
for(i in 1:n) {
  for (j in i:n) {
    nas_obs <- unique(sort(c(na_obs[[i]], na_obs[[j]])))
    # original
    a <- counts_obs[-nas_obs,i]
    b <- counts_obs[-nas_obs,j]
    csims_obs[i,j] <- cosine(a,b)
    csims_obs[j,i] <- csims_obs[i,j]
    # centered
    nas_obs <- unique(sort(c(na_obs[[i]], na_obs[[j]])))
    a <- centered_counts_obs[-nas_obs,i]
    b <- centered_counts_obs[-nas_obs,j]
    centered_csims_obs[i,j] <- cosine(a,b)
    centered_csims_obs[j,i] <- centered_csims_obs[i,j]
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

K_orig <- 9
K_norm <- 5

top_csims_obs <- apply(csims_obs, 2, function(x) order(x, decreasing = TRUE)[1:(K_orig+1)])
top_centered_csims_obs <- apply(centered_csims_obs, 2, function(x) order(x, decreasing = TRUE)[1:(K_norm+1)])

final_counts_obs <- counts_obs
final_centered_counts_obs <- centered_counts_obs
final_rescaled_counts_obs <- rescaled_counts_obs
final_normalized_counts_obs <- normalized_counts_obs

for (j in 1:n) {
  nearest_users_obs <- top_csims_obs[-1,j]
  near_values_obs <- final_counts_obs[na_obs[[j]], nearest_users_obs]
  near_rescaled_values_obs <- final_rescaled_counts_obs[na_obs[[j]], nearest_users_obs]
  weight_csims_obs <- csims_obs[nearest_users_obs,j]
  
  nearest_centered_users_obs <- top_centered_csims_obs[-1,j]
  near_centered_values_obs <- final_centered_counts_obs[na_obs[[j]], nearest_centered_users_obs]
  near_normalized_values_obs <- final_normalized_counts_obs[na_obs[[j]], nearest_centered_users_obs]
  weight_centered_csims_obs <- centered_csims_obs[nearest_centered_users_obs, j]
  for (i in 1:length(na_obs[[j]])) {
    final_counts_obs[na_obs[[j]][i],j] <- weighted.mean(near_values_obs[i,], weight_csims_obs, na.rm = TRUE)
    final_normalized_counts_obs[na_obs[[j]][i],j] <- weighted.mean(near_normalized_values_obs[i,], weight_centered_csims_obs, na.rm = TRUE)
  }
}

for (j in 1:n) {
  final_normalized_counts_obs[,j] <- final_normalized_counts_obs[,j] * csds_obs[j] + cmeans_obs[j]
  final_normalized_counts_obs <- pmax(as.matrix(final_normalized_counts_obs), 0)
}

final_counts_obs[,c(284, 1062, 1246, 2312)] <- 0
error_nas <- which(is.na(final_counts_obs), arr.ind = TRUE)
if(dim(error_nas)[1] > 0)
{
  print('orig')
  print(dim(error_nas)[1])
  for (i in 1:dim(error_nas)[1]) {
    final_counts_obs[error_nas[i,1],error_nas[i,2]] <- mean(final_counts_obs[error_nas[i,1],][!is.na(final_counts_obs[error_nas[i,1],])])
  }
}

final_normalized_counts_obs[,c(284, 1062, 1246, 2312)] <- 0
error_nas <- which(is.na(final_normalized_counts_obs), arr.ind = TRUE)
if(dim(error_nas)[1] > 0)
{
  print('norm')
  print(dim(error_nas)[1])
  for (i in 1:dim(error_nas)[1]) {
    final_normalized_counts_obs[error_nas[i,1],error_nas[i,2]] <- mean(final_normalized_counts_obs[error_nas[i,1],][!is.na(final_normalized_counts_obs[error_nas[i,1],])])
  }
}

out_orig <- final_counts_obs[ind_na]
out_norm <- final_normalized_counts_obs[ind_na]
answers <- Y_oracle[ind_na]
mse(out_orig, answers)
mse(out_norm, answers)

