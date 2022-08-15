#' Function to estimate the OVA multiclass probability
#'
#' @param train.data Training data as list object
#' @param tune.data Tuning data as list object
#' @param test.data Test data as list object
#' @param data_dim Data dimensionality
#' @param kernel Kernel to use in wSVMs
#' @param tuning.criteria Tuning by GKL or EGKL
#' @param seed Set random seed
#' @param eps Small threshold value
#'
#' @return OVA learning multiclass probability estimation as a list object

one2all.svm.probability <- function(train.data, tune.data, test.data, data_dim,
                                    kernel = list(type = 'linear', param1 = NULL, param2 = NULL), tuning.criteria = 'GKL', seed = 1, eps = 1e-10){

  # check if the provided train and test data is in a dataframe
  # data include the both tuning and testing set
  # tuning the parameters and training wsvm in tuning set, and calculate probability and evalutae in the test test
  if(!is.data.frame(train.data)) train.data <- as.data.frame(train.data)
  if(!is.data.frame(tune.data)) tune.data <- as.data.frame(tune.data)
  if(!is.data.frame(test.data)) test.data <- as.data.frame(test.data)

  n.train<- dim(train.data)[1]
  n.tune <- dim(tune.data)[1]
  n.test <- dim(test.data)[1]

  level <- levels(factor(test.data$ylabel))
  n_class <- length(level)

  train.one2all <- data_one2all(train.data, data.dim = data_dim, num_class = n_class)
  tune.one2all <- data_one2all(tune.data, data.dim = data_dim, num_class = n_class)


  test.all.x <-  test.data[,1:data_dim]
  test.all.y <-  as.vector(test.data[,data_dim+1])
  test.all.actual.prob <- as.matrix(test.data[, (data_dim+2):(data_dim+1+n_class)])
  colnames(test.all.actual.prob) <- NULL


  if(kernel$type == 'linear') param_size <-1
  else if(kernel$type == 'rbf') param_size <-2
  else if(kernel$type == 'poly' | kernel$type== 'nn') param_size <-3


  # params 5: class1, class 2, gkl, egkl and test error
  tuning_matrix <- matrix(0, nrow = n_class, ncol = param_size + 4)

  # matrix to hold all the binary class probability estimation
  # the column is the condintional probability for binary class
  # the first row is the class pair, for each pair column, left is class i and right is class j
  # the last (n_class + 2)th column indicates which base class to use use max voting class
  # the last (n_class + 1)-2th column is the n_class which includes the sum of probability difference
  # the last column the class with max sum of probability difference

  estimate_prob_multiclass <- matrix(0, nrow = n.test, ncol = n_class + 2)

  n_class_name_pred <- n_class_name_real <- vector()


  if (kernel$type == 'rbf'){
    # train the pairwise binary wsvm classifier and estimate the probability
    # i is plus class, j is minus class
    for (i in 1:n_class){

      print(i)

      class <- paste("class_", i, sep = "")

      train <- train.one2all[[class]]
      tune <- tune.one2all[[class]]

      x.train <- train[,1:data_dim]
      y.train <- train[,data_dim+1]

      x.tune <- tune[,1:data_dim]
      y.tune <- tune[,data_dim+1]

      # calculate the median distance between 2 classes
      # get sigmaM
      train.separate <- multiclass.separation(as.data.frame(train), data_dim)
      sigma_M <- median_class_distance(train.separate$class_list_X)$medianD

      # get theorical pairwise conditional probability
      px_test <- tune[, data_dim+(1+i)]


      #### For Gaussian (Radial) Kernel ####
      ######### paramter tunings ###########

      # set lambda grid
      lambda <- rep(0, 32)
      count <- 1
      for (jj in -8:7){
        low <- 10^jj
        up <- 10^(jj+1)
        lambda[count:(count+1)] <- seq(low+ (up-low)/2, up, length.out = 2)
        count <- count+2
      }

      sigma <- seq(1,6)*sigma_M/4

      Error_Radial <- matrix(0, nrow = length(lambda)*length(sigma), ncol = 5)

      print("start tuning")

      # count the row in matrix Error_Radial
      count <- 1

      # calculate the training and test error for all the tuning parameters
      for(s in 1:length(lambda)){
        for(t in 1:length(sigma)){
          # Record the parameters
          Error_Radial[count,1] <- lambda[s]
          Error_Radial[count,2] <- sigma[t]

          #fit linear svm
          prob_estimate_svm <- wsvm.prob (x.train, y.train, x.tune, y.tune, kernel = list(type = 'rbf', param1 = sigma[t], param2 = NULL), lambda = lambda[s])
          est_prob <- prob_estimate_svm$estimate_prob

          #calculate egkl and gkl
          egkl <- wsvm.egkl(y.tune, est_prob)

          # calculate the GKL
          gkl <- gkl.cal(px_test, est_prob)

          Error_Radial[count,3] <- egkl
          Error_Radial[count,4] <- gkl

          #testing error
          Error_Radial[count,5]  = prob_estimate_svm$test.error
          count <- count+1
        }
      }

      # set all NAs to 0
      Error_Radial[is.na(Error_Radial)] = 0

      # based on different tuning criteria
      if (tuning.criteria == 'EGKL'){
        Error_Radial[Error_Radial[,3] == 0,][,3] <- Inf
        best.tuning <- as.vector(Error_Radial[which.min(Error_Radial[,3]),])
        b.lambda <- Error_Radial[which.min(Error_Radial[,3]),1]
        b.sigma <- Error_Radial[which.min(Error_Radial[,3]),2]
      }

      else if(tuning.criteria == 'GKL'){
        # use GKL as criteria
        Error_Radial[Error_Radial[,4] == 0,][,4] <- Inf
        best.tuning <- as.vector(Error_Radial[which.min(Error_Radial[,4]),])
        b.lambda <- Error_Radial[which.min(Error_Radial[,4]),1]
        b.sigma <- Error_Radial[which.min(Error_Radial[,4]),2]
      }

      tuning_matrix[i, ]  <- c(i, best.tuning)

      print(best.tuning)
      print(c(b.lambda, b.sigma))

      colnames(Error_Radial) <- c("lambda","sigma","EGKL","GKL","test_error")
      # print(Error_Radial)

      # refit the wsvm with best lambda and sigma, and calculate the probability for the test.all data
      prob_estimate_svm <- wsvm.prob (x.train, y.train, test.all.x, test.all.y, kernel = list(type = 'rbf', param1 = b.sigma, param2 = NULL), lambda = b.lambda)


      # estimated binary class probability for the test.all data
      est_prob_class1 <- prob_estimate_svm$estimate_prob # probability estimation for class 1 (i)

      estimate_prob_multiclass [1:n.test, i] <- as.vector(est_prob_class1)

      n_class_name_pred <- c(n_class_name_pred, paste('cl_pred_p_', i, sep = ""))
      n_class_name_real <- c(n_class_name_real, paste('cl_real_p_', i, sep = ""))

    }

    # prepare the tuning matrix after tuning
    colnames(tuning_matrix) <- c("class1", "lambda","sigma","EGKL","GKL","tune_error")
    tuning_matrix <- as.data.frame(tuning_matrix)


    estimate_prob_multiclass[1:n.test,1:n_class] <- t(apply(estimate_prob_multiclass[1:n.test,1:n_class], 1, p_normalize))
    estimate_prob_multiclass[1:n.test,n_class+1] <- as.vector(apply(estimate_prob_multiclass[1:n.test,1:n_class],1,which.max))
    estimate_prob_multiclass[1:n.test,n_class + 2] <-  test.all.y
    estimate_prob_multiclass <- cbind(estimate_prob_multiclass, test.all.actual.prob)

    add_colnames <- c('pred_label_maxP', 'real_label')
    colnames(estimate_prob_multiclass) <- c(n_class_name_pred,add_colnames,n_class_name_real)
    # print(multiclass.prob)
    print(dim(estimate_prob_multiclass))
    # print(head(estimate_prob_multiclass))

    return(list(tuning_matrix = tuning_matrix, estimate_multiclass_prob_matrix = as.data.frame(estimate_prob_multiclass)))

  }

}
