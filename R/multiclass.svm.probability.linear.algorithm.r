#' Calculate the multiclass probability based on the baseline learning
#'
#' @param pairwise.prob Pairwise probability matrix from baseline learning
#' @param base_class Baseline class
#' @param k_class Number of classes
#' @param actual_labels True class labels
#' @param actual_prob True probability
#' @param pair_indexes Pairwise class indexes
#'
#' @return Multiclass class probaility estimation matrix

multiclass.svm.probability.linear.algorithm <- function(pairwise.prob, base_class, k_class, actual_labels, actual_prob, pair_indexes){
  #use the base class show on the last column of the pairwise.prob dataframe
  # calculate the denominator for all k base class
  ncol <- dim(pairwise.prob)[2]
  pairs <- pairwise.prob[1,]
  pairwise.prob.data <- pairwise.prob[-1, ]
  n.data <- dim(pairwise.prob.data)[1]
  class_list <- seq(1:k_class)

  multiclass.prob <- matrix(0, nrow = n.data, ncol = k_class + 2)

  k_class_name_pred <- vector()
  k_class_name_real <- vector()

  # base class is k

  base.class <- paste('base_class_', base_class, sep = "")
  pair <- class_list[class_list!=base_class]
  denominator_sum <- 1 #((k,k) = 1)

  for (t in pair){
    class.pair <- paste('index_pair_', cantor_mapping(base_class, t), sep = "")
    col_index <- pair_indexes[[class.pair]]
    denominator_sum <- as.vector(denominator_sum + pairwise.prob.data[,col_index[2]]/pairwise.prob.data[,col_index[1]])
    # print(denominator_sum)
  }


  for (i in 1:n.data){

    p_sum <- 0
    dem <- as.numeric(denominator_sum[i])

    for (t in pair){
      class.pair <- paste('index_pair_', cantor_mapping(t, base_class), sep = "")
      col_index <- pair_indexes[[class.pair]]
      num <- as.numeric(pairwise.prob.data[i,col_index[1]]/pairwise.prob.data[i,col_index[2]])
      p <- num/dem
      # print(p)
      multiclass.prob[i,t] <- p
      p_sum <- p_sum + p
    }

    multiclass.prob[i,base_class] <- 1-p_sum

  }


  for (k in 1:k_class){
    k_class_name_pred <- c(k_class_name_pred, paste('cl_pred_p_', k, sep = ""))
    k_class_name_real <- c(k_class_name_real, paste('cl_real_p_', k, sep = ""))
  }


  predict_label_max <- apply(multiclass.prob, 1, which.max)
  multiclass.prob[,k_class + 1] <- predict_label_max
  multiclass.prob[,k_class + 2] <- actual_labels
  multiclass.prob <-as.data.frame(multiclass.prob)
  multiclass.prob <- cbind(multiclass.prob, actual_prob)

  add_colnames <- c('pred_label_max', 'real_label')
  colnames(multiclass.prob) <- c(k_class_name_pred,add_colnames,k_class_name_real)
  # print(multiclass.prob)
  return(multiclass.prob)
}
