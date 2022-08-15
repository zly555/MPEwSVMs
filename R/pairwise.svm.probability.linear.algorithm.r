#' Function to perform the baseline learning of multiclass probability estimation
#'
#' @param train.data Training data as list object
#' @param tune.data Tuning data as list object
#' @param test.data Test data as list object
#' @param data_dim Data dimensionality
#' @param kernel Kernel to use in wSVMs
#' @param linear.time.algorithm Set baseline learning method: 'largestClSize' for method 1; 'medianClDist' for method 2
#' @param tuning.criteria Tuning by GKL or EGKL
#' @param seed Set random seed
#' @param eps Small threshold value
#'
#' @return Baseline learning pairwise probability estimation as a list object

pairwise.svm.probability.linear.algorithm <- function(train.data, tune.data, test.data,  data_dim, kernel = list(type = 'linear', param1 = NULL, param2 = NULL),
                                                      linear.time.algorithm = list(type = 'largestClSize'), tuning.criteria = 'GKL', seed = 1, eps = 1e-10){

  # check if the provided train and test data is in a dataframe
  # data include the both tuning and testing set
  # tuning the parameters and training wsvm in tuning set, and calculate probability and evalutae in the test test
  if(!is.data.frame(train.data)) train.data <- as.data.frame(train.data)
  if(!is.data.frame(tune.data)) tune.data <- as.data.frame(tune.data)
  if(!is.data.frame(test.data)) test.data <- as.data.frame(test.data)

  n.train<- dim(train.data)[1]
  n.tune <- dim(tune.data)[1]
  n.test <- dim(test.data)[1]

  # separate the training and test data based on y labels
  # return the two lists includes XY data for each class and X data for each class
  data_train_list <- multiclass.separation(train.data, data_dim)
  data_tune_list <- multiclass.separation(tune.data, data_dim)

  data_train_XY <- data_train_list$class_list_XY
  data_train_X <- data_train_list$class_list_X

  data_tune_XY <- data_tune_list$class_list_XY
  data_tune_X <- data_tune_list$class_list_X

  if(kernel$type == 'linear') param_size <-1
  else if(kernel$type == 'rbf') param_size <-2
  else if(kernel$type == 'poly' | kernel$type== 'nn') param_size <-3

  # number of class
  n_class <- length(data_tune_XY)

  # split test data set to x and y
  test.all.x <- test.data[, 1:data_dim]
  test.all.y <- as.vector(test.data[, data_dim+1])

  test.all.actual.prob <- test.data[, (data_dim+2):(data_dim+1+n_class)]

  # # number of pair wise computation
  # n_pair <- n_class * (n_class-1)/2

  # get the base class
  base <- 0
  class_list <- seq(1:n_class)

  if(linear.time.algorithm$type == 'largestClSize'){
    base <- find.biggest.class (train.data, data.dim = data_dim)
  }

  else if(linear.time.algorithm$type == 'medianClDist'){
    base <- median_class_distance(data_train_X, calculate_between_class_dist = TRUE)$median_separated_class
  }

  else {stop("linear time algorithm not supoorted!")}

  other_classes <- class_list[class_list!=base]



  # matrix to hold all the binary class probability estimation
  # the column is the condintional probability for binary class
  # the first row is the class pair, for each pair column, left is class i and right is class j
  # the last (n_class + 2)th column indicates which base class to use use max voting class
  # the last (n_class + 1)-2th column is the n_class which includes the sum of probability difference
  # the last column the class with max sum of probability difference

  # number of pair wise computation
  n_pair <-  n_class - 1
  pair_indexes <- list()

  # params 5: class1, class 2, gkl, egkl and test error
  tuning_matrix <- matrix(0, nrow = n_pair, ncol = param_size + 5)

  estimate_prob_binary_classes <- matrix(0, nrow = n.test + 1, ncol = 2*n_pair)
  sum_prob_difference <- matrix(0, nrow = n.test, ncol = n_class)

  count_pair <- 1
  inc <- 0


  if (kernel$type == 'rbf'){
    # train the pairwise binary wsvm classifier and estimate the probability
    # i is plus class, j is minus class
    for (i in other_classes){
      print(i)

      data_XY1_train <- as.matrix(as.data.frame(data_train_XY[i]))
      data_XY2_train <- as.matrix(as.data.frame(data_train_XY[base]))
      dataXY.train <- covert_to_binary(data_XY1_train, data_XY2_train, data.dim = data_dim, seed = seed*(101*i+11))

      data_XY1_tune <- as.matrix(as.data.frame(data_tune_XY[i]))
      data_XY2_tune <- as.matrix(as.data.frame(data_tune_XY[base]))
      dataXY.tune <- covert_to_binary(data_XY1_tune, data_XY2_tune, data.dim = data_dim, seed = seed*(177*i+17))

      # calculate the median distance between 2 classes
      sigma_M <- median_class_distance(data_train_X[c(i,base)])$medianD

      x.train <- dataXY.train[,1:data_dim]
      y.train <- dataXY.train[,data_dim+1]
      x.tune <-  dataXY.tune[,1:data_dim]
      y.tune <-  dataXY.tune[,data_dim+1]

      # get theorical pairwise conditional probability
      p_ix <- dataXY.tune[, data_dim+(1+i)]
      p_jx <- dataXY.tune[, data_dim+(1+base)]
      px_test <- p_ix/(p_jx + p_ix)


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

      tuning_matrix[count_pair, ]  <- c(i, base, best.tuning)

      print(best.tuning)
      print(c(b.lambda, b.sigma))

      colnames(Error_Radial) <- c("lambda","sigma","EGKL","GKL","test_error")
      # print(Error_Radial)


      # refit the wsvm with best lambda and sigma, and calculate the probability for the test.all data
      prob_estimate_svm <- wsvm.prob (x.train, y.train, test.all.x, test.all.y, kernel = list(type = 'rbf', param1 = b.sigma, param2 = NULL), lambda = b.lambda)

      # estimated binary class probability for the test.all data
      est_prob_class1 <- prob_estimate_svm$estimate_prob # probability estimation for class 1 (i)
      est_prob_base <- 1 - est_prob_class1 # probability estimation for class 2 (j)

      estimate_prob_binary_classes [1, (count_pair+inc):(count_pair+inc+1)] <- c(i,base)
      estimate_prob_binary_classes [2:(n.test + 1), count_pair+inc] <- est_prob_class1
      estimate_prob_binary_classes [2:(n.test + 1), count_pair+inc+1] <- est_prob_base


      i_column_index <- count_pair+inc
      base_column_index <- count_pair+inc+1

      class.pair.forward <- paste('index_pair_', cantor_mapping(i,base), sep = "")
      class.pair.reverse <- paste('index_pair_', cantor_mapping(base,i), sep = "")

      # print("start debugging")
      # print(i)
      # print(j)
      # print(class.pair.forward)
      # print(class.pair.reverse)

      pair_indexes[[class.pair.forward]] <- c(i_column_index, base_column_index)
      pair_indexes[[class.pair.reverse]] <- c(base_column_index, i_column_index)

      # calculate the sum of probability difference for each class
      sum_prob_difference[,i] <- sum_prob_difference[,i] + (est_prob_class1 - est_prob_base)
      sum_prob_difference[,base] <- sum_prob_difference[,base] + (est_prob_base - est_prob_class1)

      count_pair <- count_pair + 1
      inc <- inc + 1

    }
  }


  # prepare the tuning matrix after tuning
  colnames(tuning_matrix) <- c("class1","class2", "lambda","sigma","EGKL","GKL","tune_error")
  tuning_matrix <- as.data.frame(tuning_matrix)

  # reduce the size of probability matrix
  estimate_prob_binary_classes <- as.data.frame(estimate_prob_binary_classes)

  return(list(tuning_matrix = tuning_matrix, base_class = base, estimate_prob_binary_classes_base = estimate_prob_binary_classes, k_class = n_class,
              actual_labels=test.all.y, actual_prob=test.all.actual.prob, pair_indexes = pair_indexes, sum_prob_difference = sum_prob_difference))

}
