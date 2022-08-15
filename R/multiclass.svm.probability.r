#' Calculate the multiclass probability based on the pairwise probability matrix
#'
#' @param pairwise.prob Pairwise probability matrix
#' @param k_class Number of classes
#' @param actual_labels True class labels
#' @param actual_prob True probability
#' @param pair_indexes Pairwise class indexes
#'
#' @export
#' @return Multiclass class probaility estimation matrix

multiclass.svm.probability <- function(pairwise.prob, k_class, actual_labels, actual_prob, pair_indexes){
  #use the base class show on the last column of the pairwise.prob dataframe
  # calculate the denominator for all k base class
  ncol <- dim(pairwise.prob)[2]
  pairs <- pairwise.prob[1,]
  pairwise.prob.data <- pairwise.prob[-1, ]
  n.data <- dim(pairwise.prob.data)[1]
  class_list <- seq(1:k_class)

  base_dominator_sum <- list()

  multiclass.prob <- matrix(0, nrow = n.data, ncol = k_class + 4)

  k_class_name_pred <- vector()
  k_class_name_real <- vector()

  # base class is k
  for (k in 1:k_class){
    base.class <- paste('base_class_', k, sep = "")
    pair <- class_list[class_list!=k]
    denominator_sum <- 1 #((k,k) = 1)

    for (t in pair){
      class.pair <- paste('index_pair_', cantor_mapping(k, t), sep = "")
      col_index <- pair_indexes[[class.pair]]
      denominator_sum <- as.vector(denominator_sum + pairwise.prob.data[,col_index[2]]/pairwise.prob.data[,col_index[1]])
      # print(denominator_sum)
    }

    base_dominator_sum [[base.class]] <- denominator_sum

    k_class_name_pred <- c(k_class_name_pred, paste('cl_pred_p_', k, sep = ""))
    k_class_name_real <- c(k_class_name_real, paste('cl_real_p_', k, sep = ""))

  }

  # print(base_dominator_sum)

  for (i in 1:n.data){
    base <- as.numeric(pairwise.prob.data[i,][ncol-1])

    base.class <- paste('base_class_', base, sep = "")
    dem <- as.numeric(base_dominator_sum[[base.class]][i])

    pair <- class_list[class_list!=base]

    p_sum <- 0

    for (t in pair){
      class.pair <- paste('index_pair_', cantor_mapping(t, base), sep = "")
      col_index <- pair_indexes[[class.pair]]
      num <- as.numeric(pairwise.prob.data[i,col_index[1]]/pairwise.prob.data[i,col_index[2]])
      p <- num/dem
      # print(p)
      multiclass.prob[i,t] <- p
      p_sum <- p_sum + p
    }

    multiclass.prob[i,base] <- 1-p_sum

  }

  predict_label_max <- apply(multiclass.prob, 1, which.max)
  predict_label_voting <- as.vector(pairwise.prob.data[, ncol-1])
  predict_label_maxPd <- as.vector(pairwise.prob.data[, ncol])
  multiclass.prob[,k_class + 1] <- predict_label_max
  multiclass.prob[,k_class + 2] <- predict_label_voting
  multiclass.prob[,k_class + 3] <- predict_label_maxPd
  multiclass.prob[,k_class + 4] <- actual_labels
  multiclass.prob <-as.data.frame(multiclass.prob)
  multiclass.prob <- cbind(multiclass.prob, actual_prob)

  add_colnames <- c('pred_label_max', 'pred_label_voting','pred_label_maxPd','real_label')
  colnames(multiclass.prob) <- c(k_class_name_pred,add_colnames,k_class_name_real)
  # print(multiclass.prob)
  return(multiclass.prob)
}
