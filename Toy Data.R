n = 50
error_sd = 2

## Setting 1
Data_list = get(load("Setting_initial.rda"))
## Setting 2
# A = Data_list$Alpha_sparse
# A[11:13, 1] = Data_list$Alpha_sparse[1:3, 1]
# A[14:19, 2] = Data_list$Alpha_sparse[4:9, 2]
# A[1:10, ] = 0
# Data_list$Alpha_sparse = A
# rm(A)

beta_true =  Data_list$beta_sparse
Alpha_true = Data_list$Alpha_sparse
beta_true_pattern = beta_true != 0
Delta_true = Alpha_true != 0
p = nrow(Alpha_true)
k = ncol(Alpha_true)


X = mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(rep(1, p)))

Z = array(NA, dim = c(p, p, n))
Z_temp <- Z_index <- matrix(0, nrow = p, ncol = p)
Z_index[lower.tri(Z_index)] <- 1
Z_index <- which(Z_index == 1)
Z_index_n <- length(Z_index)
for (i in 1:n) {
  Z_temp[Z_index]  <- rnorm(Z_index_n, mean = 0, sd = 1)
  Z_temp[upper.tri(Z_temp)] = t(Z_temp)[upper.tri(Z_temp)]
  Z[, , i] = Z_temp
}


Y_true = X[, beta_true_pattern, drop = FALSE] %*% beta_true[beta_true_pattern, drop = FALSE]


for (i in 1:n) {
  for (r in 1:2) {
    subset <- Delta_true[, r]
    Y_true[i] = Y_true[i] + t(Alpha_true[subset, r, drop = FALSE]) %*% Z[subset, subset, i] %*% Alpha_true[subset, r, drop = FALSE]
    }
}

Y_true = Y_true + 1
error <- rnorm(n, 0, error_sd)
Y = Y_true + error


Data = list(Y = Y,
            Y_true = Y_true,
            Z = Z,
            X = X)

