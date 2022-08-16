#
# pk = [p(x1),p(x2) ... p(xn)]k
# p, phat is a n.data * nclass matrix
#' Evaluate the probability estimation performance
#'
#' @param p True probability matrix
#' @param phat Estimated probability matrix
#' @param method Evaluation methods
#'
#' @return Evaluation performance value
#' @export

prob.evaluation <- function(p, phat, method = list(type = 'egkl')) {

  if(!is.matrix(p)) p <- as.matrix(p)
  if(!is.matrix(phat)) phat <- as.matrix(phat)

  # ## Transfor p to log P
  # p <- log(p)
  # phat <- log(phat)

  n <- dim(p)[1]
  n_class <- dim(p)[2]

  if (method$type == 'l1_norm') {
    v <- 0
    for (i in 1:n_class){
      nv <- norm(as.matrix(p[,i]-phat[,i]), '1')
      v <- v + nv
    }
    K <- v/n
  }

  else if (method$type == 'l2_norm'){
    v <- 0
    for (i in 1:n_class){
      nv <- norm(as.matrix(p[,i]-phat[,i]), '2')^2
      v <- v + nv
    }
    K <- v/n
  }

  else if (method$type == 'egkl'){
    v <- 0
    for (i in 1:n_class){
      nv <- sum(p[,1]*log(p[,1]/phat[,1]))
      v <- v + nv
    }
    K <- v/n
  }

  else if (method$type == 'gkl'){
    v <- 0
    for (i in 1:n_class){
      nv <- sum(p[,1]*log(p[,1]/phat[,1]) + (1-p[,1])*log((1-p[,1])/(1-phat[,1])))
      v <- v + nv
    }
    K <- v/n
  }

  return(K)
}



#
# and output test error based on the evaluation set
#' Function to evaluate the multi-class SVM probability estimation performance
#'
#' @param multiclass.prob Multiclass probability estimation matrix
#' @param methods Multiclass probability estimation methods: e.g: "pairwise", "one2rest","liner_time"
#'
#' @return Performance evaluation matrix and classiifcation error matrix
#' @export

evaluate_performance <- function(multiclass.prob, methods = list(type = 'pairwise')){
  if (methods$type == 'pairwise'){
    if(!is.data.frame(multiclass.prob)) multiclass.prob <- as.data.frame(multiclass.prob)
    n_class <- (dim(multiclass.prob)[2]-4)/2

    phat <- as.matrix(multiclass.prob[,1:n_class])
    p <- as.matrix(multiclass.prob[, (n_class + 4 + 1):dim(multiclass.prob)[2]])

    colnames(phat) <- NULL
    rownames(phat) <- NULL
    colnames(p) <- NULL
    rownames(p) <- NULL

    evaluation_matrix <- matrix(0, nrow = 1, ncol = 4)
    evaluation_matrix[1,1] <- prob.evaluation(p, phat, method = list(type = 'l1_norm'))
    evaluation_matrix[1,2] <- prob.evaluation(p, phat, method = list(type = 'l2_norm'))
    evaluation_matrix[1,3] <- prob.evaluation(p, phat)
    evaluation_matrix[1,4] <- prob.evaluation(p, phat, method = list(type = 'gkl'))

    evaluation <- as.data.frame(evaluation_matrix)
    colnames(evaluation) <- c("L1_norm", 'L2_norm', 'EGKL', 'GKL')

    test_error_max <- mean(multiclass.prob[,n_class+1]!=multiclass.prob[,n_class+4])
    test_error_voting <- mean(multiclass.prob[,n_class+2]!=multiclass.prob[,n_class+4])
    test_error_maxPd <- mean(multiclass.prob[,n_class+3]!=multiclass.prob[,n_class+4])

    TestClassificationError = as.data.frame(matrix(c(test_error_max,test_error_voting, test_error_maxPd), nrow = 1))
    colnames(TestClassificationError) <- c("test_error_max", 'test_error_voting', 'test_error_maxPd')

    performance <- list(evaluationMat = evaluation, TestClassificationError = TestClassificationError)

    return(performance)
  }

  else if(methods$type == 'one2rest' | methods$type == 'liner_time' | methods$type == 'mllr' | methods$type == 'lda' |methods$type == 'tree' | methods$type == 'rf' | methods$type == 'xgb'){
    if(!is.data.frame(multiclass.prob)) multiclass.prob <- as.data.frame(multiclass.prob)

    n_class <- (dim(multiclass.prob)[2]-2)/2

    phat <- as.matrix(multiclass.prob[,1:n_class])
    p <- as.matrix(multiclass.prob[, (n_class + 2 + 1):dim(multiclass.prob)[2]])

    colnames(phat) <- NULL
    rownames(phat) <- NULL
    colnames(p) <- NULL
    rownames(p) <- NULL

    evaluation_matrix <- matrix(0, nrow = 1, ncol = 4)
    evaluation_matrix[1,1] <- prob.evaluation(p, phat, method = list(type = 'l1_norm'))
    evaluation_matrix[1,2] <- prob.evaluation(p, phat, method = list(type = 'l2_norm'))
    evaluation_matrix[1,3] <- prob.evaluation(p, phat)
    evaluation_matrix[1,4] <- prob.evaluation(p, phat, method = list(type = 'gkl'))

    evaluation <- as.data.frame(evaluation_matrix)
    colnames(evaluation) <- c("L1_norm", 'L2_norm', 'EGKL', 'GKL')

    test_error_maxP <- mean(multiclass.prob[,n_class+1]!=multiclass.prob[,n_class+2])

    TestClassificationError = as.data.frame(matrix(c(test_error_maxP), nrow = 1))
    colnames(TestClassificationError) <- c("test_error_maxP")

    performance <- list(evaluationMat = evaluation, TestClassificationError = TestClassificationError)

    return(performance)
  }

}


### generate the MP matrix for reliability diagram ####
#' Generate the MP matrix for reliability diagram
#'
#' @param mp Multiclass probability matrix
#' @param n_class Number of classes
#' @param n_data Number of test data
#' @param method Multilcass probability estimation method
#'
#' @return Probability matrix for relability diagram
#' @export

generate_MP <- function(mp, n_class, n_data, method){
  if(!is.data.frame(mp)) mp <- as.data.frame(mp)

  prob_est_df <- as.data.frame(matrix(0, nrow = n_data, ncol=5))
  colnames(prob_est_df) <- c('p_hat','pred_label','p_real','true_label','method')
  is_paiwise = method %in% MP_PAIRWISE

  if(is_paiwise){
    skip <- 2
    phat_label <- mp[, n_class+1]
    preal_label <- mp[, n_class+1+skip+1]
    prob_est_df[,1] <- unlist(lapply(1:length(phat_label), function(x) mp[x, phat_label[x]]))
    prob_est_df[,2] <- phat_label
    prob_est_df[,3] <- unlist(lapply(1:length(preal_label), function(x) mp[x, (n_class+2+skip)+preal_label[x]]))
    prob_est_df[,4] <- preal_label
    prob_est_df[,5] <- method
  }

  if(!is_paiwise){
    skip <- 0
    phat_label <- mp[, n_class+1]
    preal_label <- mp[, n_class+1+skip+1]
    prob_est_df[,1] <- unlist(lapply(1:length(phat_label), function(x) mp[x, phat_label[x]]))
    prob_est_df[,2] <- phat_label
    prob_est_df[,3] <- unlist(lapply(1:length(preal_label), function(x) mp[x, (n_class+2+skip)+preal_label[x]]))
    prob_est_df[,4] <- preal_label
    prob_est_df[,5] <- method
  }

  return(prob_est_df)

}
