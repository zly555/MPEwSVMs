#
#' Reconstruct the pairwise conditional probability table based on baseline learning
#'
#' @param estimate_prob_binary_classes_base Pairwise probability with baseline learning
#' @param base_class Baseline class
#' @param num_class Number of classes
#' @param pair_indexes Index of pairwise classes
#' @param sum_prob_difference Sum of pairwise probability difference
#' @param actual_labels True labels
#' @param actual_prob True probability
#'
#' @return Reconstructed pairwise coupling probability estimation table from baseline learning as a list object

pairwise.svm.probability.linear.time.reconstruct <- function(estimate_prob_binary_classes_base, base_class, num_class, pair_indexes, sum_prob_difference, actual_labels, actual_prob){
  if(!is.data.frame(estimate_prob_binary_classes_base)) estimate_prob_binary_classes_base <- as.data.frame(estimate_prob_binary_classes_base)

  estimate_prob_binary_classes_base <- as.matrix(estimate_prob_binary_classes_base)
  estimate_prob_binary_classes_base_data <- estimate_prob_binary_classes_base[-1, ]

  colnames(estimate_prob_binary_classes_base) <- NULL
  rownames(estimate_prob_binary_classes_base) <- NULL

  class_list <- seq(1:num_class)
  n.data <- dim(estimate_prob_binary_classes_base_data)[1]
  ncol_pt <- dim(estimate_prob_binary_classes_base)[2]
  other_classes <- class_list[class_list!=base_class]
  num_other_classes <- length(other_classes)


  # number of pairwise computation
  n_pair <- num_other_classes * (num_other_classes-1)/2

  estimate_prob_binary_classes <- matrix(0, nrow = n.data + 1, ncol = 2*n_pair)

  count_pair <- 1
  inc <- 0

  #column index with pairs
  sum_prob_difference <- matrix(0, nrow = n.data, ncol = num_class)

  for (s in 1:(num_other_classes-1)){
    for (t in (s+1):num_other_classes){
      if(s<t & t<=num_other_classes){
        j <- other_classes[s]
        jp <- other_classes[t]

        class.pair.base.1 <- paste('index_pair_', cantor_mapping(j, base_class), sep = "")
        class.pair.base.2 <- paste('index_pair_', cantor_mapping(jp, base_class), sep = "")
        col_index.1 <- pair_indexes[[class.pair.base.1]]
        col_index.2 <- pair_indexes[[class.pair.base.2]]

        pro_jb <- as.vector(estimate_prob_binary_classes_base_data[,col_index.1[1]])
        pro_jpb <- as.vector(estimate_prob_binary_classes_base_data[,col_index.2[1]])

        est_prob_class1 <- pairwise_prob_infer(pro_jb, pro_jpb)
        est_prob_class2 <- 1 - est_prob_class1

        estimate_prob_binary_classes [1, (count_pair+inc):(count_pair+inc+1)] <- c(j,jp)
        estimate_prob_binary_classes [2:(n.data + 1), count_pair+inc] <- est_prob_class1
        estimate_prob_binary_classes [2:(n.data + 1), count_pair+inc+1] <- est_prob_class2

        j_column_index <- count_pair+inc + ncol_pt
        jp_column_index <- count_pair+inc+1 + ncol_pt

        class.pair.forward <- paste('index_pair_', cantor_mapping(j,jp), sep = "")
        class.pair.reverse <- paste('index_pair_', cantor_mapping(jp,j), sep = "")

        pair_indexes[[class.pair.forward]] <- c(j_column_index, jp_column_index)
        pair_indexes[[class.pair.reverse]] <- c(jp_column_index, j_column_index)

        # calculate the sum of probability difference for each class
        sum_prob_difference[,j] <- sum_prob_difference[,j] + (est_prob_class1 - est_prob_class2)
        sum_prob_difference[,jp] <- sum_prob_difference[,jp] + (est_prob_class2 - est_prob_class1)

        count_pair <- count_pair + 1
        inc <- inc + 1

      }
    }
  }

  # reconstructed pairwise table
  estimate_prob_binary_classes_full <- as.matrix(cbind(estimate_prob_binary_classes_base, estimate_prob_binary_classes))
  estimate_prob_binary_classes_data <- estimate_prob_binary_classes_full[-1, ]

  total_cols <- dim(estimate_prob_binary_classes_full)[2]
  odd_cols <- seq(1, total_cols-1, by = 2)
  compare_matrix = matrix(0, nrow=n.data, ncol =  total_cols)

  for (i in odd_cols){
    left <- as.numeric(estimate_prob_binary_classes_data[,i] >= estimate_prob_binary_classes_data[,i+1])
    right <- as.numeric(estimate_prob_binary_classes_data[,i] < estimate_prob_binary_classes_data[,i+1])
    compare_matrix[,i] <- left
    compare_matrix[,i+1] <- right
  }

  estimate_prob_binary_classes <- rbind(estimate_prob_binary_classes_full, compare_matrix)

  group <- as.vector(estimate_prob_binary_classes[1,])
  crossMat <- t(apply(estimate_prob_binary_classes[(n.data + 2):(n.data*2 + 1),], 1, function(x) x*group))


  # find the most abundent class with bigger pair wise probability
  mf <- as.vector(as.numeric(apply(crossMat, 1, function(x) {v <- as.vector(names(sort(table(x),decreasing=TRUE))); v[v!= 0][1]})))
  mf <- c(0, mf)

  estimate_prob_binary_classes <- estimate_prob_binary_classes [1:(n.data + 1), ]
  estimate_prob_binary_classes  <- cbind(estimate_prob_binary_classes, mf)

  base_class_max_prob_sum <- as.vector(apply(sum_prob_difference,1,which.max))
  base_class_max_prob_sum <- c(0, base_class_max_prob_sum)

  estimate_prob_binary_classes  <- as.data.frame(cbind(estimate_prob_binary_classes, base_class_max_prob_sum))

  result <- list(estimate_prob_binary_classes = estimate_prob_binary_classes, k_class = num_class,
                 actual_labels = actual_labels, actual_prob = actual_prob, pair_indexes = pair_indexes)

  return (result)

}
