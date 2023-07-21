# Working version as of 1/11/23

# See http://www.science.smith.edu/~jcrouser/SDS293/labs/lab11-r.html

install.packages("pls")
install.packages("dplyr")
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(tidyr)
library(pls)
library(plsVarSel)
library(vip)
library(plotly)
# library(systemsseRology)

set.seed(10097)

# Load in the data
filename <- 'hp1_pls_inputs_211215_withML.csv'
X <-  read.csv(filename,row.names=NULL,header=TRUE,check.names=FALSE)
colnames(X) <- make.unique(colnames(X))
X[is.na(X)] <- 0 # Replace NA with 0 (justified?)
x_full <- X

# Drop columns and write file containing them
X <- X %>% select(where(~n_distinct(.) > 1))
X <-  X[,colSums(X) > 8]

n_perms=100

X_init <- X

results <- as.data.frame(matrix(nrow=n_perms,ncol=4))
colnames(results)<-c("model_train_MSE","model_test_MSE","perm_train_MSE","perm_test_MSE")

# Begin loop
for (i in  1:n_perms){
  X <- X_init # Re-initialize X
  
  # Generate train and test sets
  train = X %>% sample_frac(0.8)
  test = X %>% setdiff(train)
  
  y <- log10(X$'Partition Ratio Probe')
  y_train <- log10(train$'Partition Ratio Probe')
  y_test <- log10(test$'Partition Ratio Probe')
  
  # Remove y from the X data matrices
  X <- X[,2:ncol(X)]
  train <- train[,2:ncol(train)]
  test <- test[,2:ncol(test)]
  
  x_train <-  model.matrix(y_train~., train)[,-1]
  x_test <-  model.matrix(y_test~., test)[,-1]
  
  
  # Train the PLSR model
  pls_fit <-  plsr(y_train~., data = train, scale = TRUE, validation = "CV")
  # summary(pls_fit)
  
  
  # Select dimensions for lowest cross-validation error
  MSEP_pls <-  MSEP(pls_fit)
  cvs_pls <- MSEP_pls$val[1,1, ]
  ind_pls <- which.min(cvs_pls)
  dims_pls <- ind_pls-1
  
  # Get the training set MSE
  train_mse_pls <- cvs_pls[ind_pls]
  
  # Compute the test set MSE
  pls_pred <-  predict(pls_fit, x_test, ncomp = dims_pls)
  test_mse_pls <- mean((pls_pred - y_test)^2)
  
  results[i,"model_train_MSE"] <- train_mse_pls
  results[i,"model_test_MSE"] <- test_mse_pls
  
  # Fit PLSR on the whole dataset (comment out for validation)
  # pls_fit2 <-  plsr(y~., data = X, scale = TRUE, ncomp = dims_pls)
  
  # Perform permutation test validation
  perm <- sample(1:length(y_train))
  y_perm <- y_train[perm]
  
  # Train the permuted (null) model
  perm_fit <-  plsr(y_perm~., data = train, scale = TRUE, validation = "CV")
  MSEP_perm <-  MSEP(perm_fit)
  cvs_perm <- MSEP_perm$val[1,1, ]
  
  # Get the training set MSE
  train_mse_perm <- cvs_perm[ind_pls] # Use number of dimensions previously computed
  
  # Compute the test MSE for permuted model
  perm_pred <-  predict(perm_fit, x_test, ncomp = dims_pls)
  test_mse_perm <- mean((perm_pred - y_test)^2)
  
  results[i,"perm_train_MSE"] <- train_mse_perm
  results[i,"perm_test_MSE"] <- test_mse_perm

} # End for loop

# Write results to file
write.csv(results,'validation_results_hp1.csv')

