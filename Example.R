library(mvtnorm)
library(invgamma)
library(truncnorm)
library(dplyr)
library(igraph)

# Data = get(load("ToyData.rda"))
source("Toy Data.R")
source("JNNTs.R")




X = Data$X
Y = Data$Y
Z = Data$Z
W = Data$W

n = ncol(X)
p = ncol(X)
intercept_index = TRUE

training_ratio = 0.9
training_n = round(n * 0.9)
training_ID = sample.int(n, training_n)
training_ID = sort(training_ID)

Y_training = Y[training_ID, ]
X_training = X[training_ID, ]
W_training = W[training_ID, ]
Z_training = Z[, , training_ID]

Y_test = Y[-training_ID, ]
X_test = X[-training_ID, ]
W_test = W[-training_ID, ]
Z_test = Z[, , -training_ID]

print("Fitting")
burnin = 2000
M = 5000
MPP = 0.5
a_rho = 0.2
b_rho = 0.8
a_lambda = 0
b_lambda = 1
sig2_i = 10
d_vec = 2:4
result = as.list(d_vec)
names(result) = d_vec
error_mat = matrix(NA, nrow = length(d_vec), ncol = 4)
colnames(error_mat) = c("In_MSE", "In_Rsq", "Out_MSE", "Out_Rsq")

for(d in 2:4) {
  ## model fitting
  result[[d-d_vec[1]+1]] = JNNTs(Y = Y_training,
                                  X = X_training,
                                  Z = Z_training,
                                  #  W = W_training,
                                  d = d, intercept = intercept_index,
                                  a_rho = a_rho, b_rho = b_rho,
                                  a_lambda = a_lambda, b_lambda = b_lambda,
                                  sig2_i = sig2_i,
                                  burnin = burnin, M = M,
                                  disp = TRUE)

  ## estimation
  if(nrow(result$ita_save) == 1) {
      intercept_esti = mean(result$ita_save)
  } else {
      intercept_esti = rowMeans(result$ita_save)
  }

  if(length(result$lambda_save) > 1) {
      lambda_seq = result$lambda_save
  } else {
      lambda_seq = rep(result$lambda_save[1], M)
  }

  beta_sel = apply(abs(result$gamma_save), 1, function(x, y) x - y, y = lambda_seq)
  beta_sel = colMeans(beta_sel > 0)
  beta_sel = which(beta_sel > MPP)
  beta_esti = rowMeans(result$beta_save)


  subnet_sel = array(NA, dim = c(p, p, d, M))
  subnet_esti = array(0, dim = c(p, p, d))
  for (iter in seq(M)) {
      a1 = abs(result$theta_save[, iter]) > lambda_seq[iter]
      for (r in seq(d)) {
      x <-  a1 * (abs(result$Theta_save[, r, iter]) > lambda_seq[iter])
      subnet_sel[, , r, iter] = matrix(x, ncol = 1) %*% x
      x = result$Alpha_save[, r, iter]
      subnet_esti[, , r] = subnet_esti[, , r] + matrix(x, ncol = 1) %*% x
      }
  }

  A_sel = rowMeans(subnet_sel, dims = 3)
  A_esti = subnet_esti / M


  ## Out-of-sample performance
  if(length(intercept_esti) == 1) {
      Out_Resi = rep(intercept_esti , n - training_n)
  } else {
      Out_Resi = cbind2(1, W_test) %*% intercept_esti
  }

  if(length(beta_sel) > 0) {
      Out_Resi = Out_Resi + X_test[, beta_sel, drop = FALSE] %*% beta_esti[beta_sel]
  }

  for (i in seq(length(Y_test))) {
    for(r in seq(d)) {
      A1 = which(A_sel[, , r]  > MPP, arr.ind = TRUE)
      A1 = A1[A1[,1] > A1[, 2], ]
      Out_Resi[i] = Out_Resi[i] + sum(A_esti[, , r][A1] * Z_test[, , i][A1]) * 2
    }
  }

  Out_Resi = Y_test - Out_Resi
  Out_MSE = sum(Out_Resi^2)/(n - training_n)
  Out_Rsq = 1 - sum(Out_Resi^2)/sum((Y_test-mean(Y_test))^2)

  ## In-sample performance
  if(length(intercept_esti) == 1) {
    In_Resi = rep(intercept_esti , training_n)
  } else {
    In_Resi = cbind2(1, W_training) %*% intercept_esti
  }

  if(length(beta_sel) > 0) {
    In_Resi = In_Resi + X_training[, beta_sel, drop = FALSE] %*% beta_esti[beta_sel]
  }

  for (i in seq(length(Y_training))) {
    for(r in seq(d)) {
    A1 = which(A_sel[, , r]  > MPP, arr.ind = TRUE)
    A1 = A1[A1[,1] > A1[, 2], ]
    In_Resi[i] = In_Resi[i] + sum(A_esti[, , r][A1] * Z_training[, , i][A1]) * 2
    }
  }

  In_Resi = Y_training - In_Resi
  In_MSE = sum(In_Resi^2)/(training_n)
  In_Rsq = 1 - sum(In_Resi^2)/sum((Y_training-mean(Y_training))^2)

  error_mat[d-d_vec[1]+1, ] = c(In_MSE = In_MSE, In_Rsq = In_Rsq,
                              Out_MSE = Out_MSE, Out_Rsq = Out_Rsq)

}

d_esti = which.max(error_mat[, "Out_Rsq"]) + d_vec[1] - 1

result_refit = JNNTs(Y = Y,
                     X = X,
                     Z = Z,
                     # W = W,
                     d = d_esti, intercept = intercept_index,
                     a_rho = a_rho, b_rho = b_rho,
                     a_lambda = a_lambda, b_lambda = b_lambda,
                     sig2_i = sig2_i,
                     burnin = burnin, M = M,
                     disp = TRUE)

if(nrow(result_refit$ita_save) == 1) {
  intercept_esti = mean(result_refit$ita_save)
} else {
  intercept_esti = rowMeans(result_refit$ita_save)
}

if(length(result_refit$lambda_save) > 1) {
  lambda_seq = result_refit$lambda_save
} else {
  lambda_seq = rep(result_refit$lambda_save[1], M)
}

cat("Estimation: Region;")
beta_sel = apply(abs(result_refit$gamma_save), 1, function(x, y) x - y > 0, y = lambda_seq)
beta_sel = colMeans(beta_sel)
beta_sel = which(beta_sel > MPP)
beta_esti = rowMeans(result_refit$beta_save)

cat("Network \r\n")
subnet_sel  = array(0, dim = c(p, p, d_esti))
subnet_esti = array(0, dim = c(p, p, d_esti))

for (iter in seq(dim(result_refit$Alpha_save)[3])) {
  theta_mask = abs(result_refit$theta_save[, iter]) > lambda_seq[iter]
  for (r in seq(d_esti)) {
    x <- theta_mask * (abs(result_refit$Theta_save[, r, iter]) > lambda_seq[iter])
    subnet_sel[, , r] = subnet_sel[, , r] + matrix(x, ncol = 1) %*% x
    x <- result_refit$Alpha_save[, r, iter]
    subnet_esti[, , r] = subnet_esti[, , r] + matrix(x, ncol = 1) %*% x
  }
}

subnet_sel =  subnet_sel / M
fit_plot = subnet_pp(subnet_sel, graph = TRUE, decision = MPP)
subnet_esti = subnet_esti/ M