#' Class probability estimation for giving data
#'
#' @param x.train Training data in matrix
#' @param y.train Training labels
#' @param x.test Test data in matrix
#' @param y.test Test labels
#' @param kernel Kernel to use for wSVMs
#' @param lambda Regularization params
#' @param eps small threshold
#'
#' @return List of probability estimation object

wsvm.prob <- function(x.train, y.train, x.test, y.test, kernel = list(type = 'linear', param1 = NULL, param2 = NULL), lambda , eps = 1e-10){

  if(!is.matrix(x.train)) x.train <- as.matrix(x.train)
  if(!is.vector(y.train)) y.train <- as.vector(y.train)
  if(!is.matrix(x.test)) x.test <- as.matrix(x.test)
  if(!is.vector(y.test)) y.test <- as.vector(y.test)

  n.train <- nrow(x.train)
  n.test <- nrow(x.test)

  m = floor(sqrt(n.train))
  PI <- rep(0, (m+1))
  PI_mid <- 1:(m-1)/m
  PI[2:m] <- PI_mid
  PI[1] <- 0
  PI[m+1]<- 1


  predicted_labels_train <- matrix(0, nrow = n.train, ncol =length(PI_mid))
  predicted_labels_test <- matrix(0, nrow = n.test, ncol = length(PI_mid))

  estimate_labels_train <- estimate_prob_train <- rep(0, n.train)
  estimate_labels_test <- estimate_prob_test <- rep(0, n.test)

  # store the series of wSVM classifier for future prediction
  c_matrix <- matrix(0, nrow = length(PI_mid), ncol = dim(x.train)[1])
  d_vector <- rep(0, length(PI_mid))
  for (i in 1:length(PI_mid)){
    tryCatch({
      wsvm.model <- wsvm(x.train, y.train, PI = PI_mid[i], kernel = kernel, lambda = lambda, eps = eps)
      c_matrix[i,] <- wsvm.model$c
      d_vector[i] <- wsvm.model$d
      predicted_labels_train[,i] <- wsvm.predict(x.train, x.train, wsvm.model$c, wsvm.model$d, wsvm.model$kernel)
      predicted_labels_test[,i] <- wsvm.predict(x.test, x.train, wsvm.model$c, wsvm.model$d, wsvm.model$kernel)
    }, error=function(e){})
  }

  # print("finish loop!!!!!")
  predicted_labels_train <- as.matrix(cbind(rep(1, n.train), predicted_labels_train, rep(-1, n.train)))
  predicted_labels_test <- as.matrix(cbind(rep(1, n.test), predicted_labels_test, rep(-1, n.test)))

  # print(predicted_labels_train)

  # print('check probability')
  estimate_prob_train <- as.vector(apply(predicted_labels_train,1,prob.estimate,PI))
  # print(estimate_prob_train)
  estimate_prob_test <- as.vector(apply(predicted_labels_test,1,prob.estimate,PI))
  # print(estimate_prob_test)

  estimate_labels_train <- estimate_prob_train
  estimate_labels_train[estimate_labels_train>=0.5] <- 1
  estimate_labels_train[estimate_labels_train<0.5] <- -1

  estimate_labels_test <- estimate_prob_test
  estimate_labels_test[estimate_labels_test>=0.5] <- 1
  estimate_labels_test[estimate_labels_test<0.5] <- -1

  train_error <- mean(estimate_labels_train!=y.train)
  test_error <- mean(estimate_labels_test!=y.test)

  wsvm.prob.model <- list(estimate_prob = estimate_prob_test, estimate_label = estimate_labels_test,
                          train.error = train_error, test.error = test_error, kernel = kernel, lambda = lambda, cmat = c_matrix, dvec =d_vector, PIseries = PI)

  return(wsvm.prob.model)

}
