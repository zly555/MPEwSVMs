#' Function to calculate the within class compactness
#'
#' @param classX Input class for calculate the Dcp
#'
#' @export
#' @return Within class compactness

within_class_compactness <- function(classX){
  X <- as.matrix(as.data.frame(classX))

  data_points <- seq(1:dim(X)[1])
  sum_distance <- rep(0, dim(X)[1])

  for (i in 1:length(data_points)){
    rest_points <- data_points[data_points!=i]
    for (j in rest_points){
      dist <- norm(as.matrix(X[i,]-X[j,]), type ='2')
      sum_distance[i] <- sum_distance[i] + dist
    }
  }

  median_point <- find_median_index(sum_distance)
  rest_points <- data_points[data_points!=median_point]
  distToMedian <- rep(0, dim(X)[1]-1)

  for (k in rest_points){
    dist <- norm(as.matrix(X[median_point,]-X[k,]), type ='2')
    distToMedian[k] <- dist
  }

  return(as.numeric(max(distToMedian)))

}


#
#' Calculate the mediam of between class distance
#'
#' @param list_classes_X List of X for all classes
#' @param calculate_between_class_dist Boolean to indicate whether to calculate the between class distance for baseline method 2
#'
#' @export
#' @return Theclass index of the median of between class distance or the SigmaM

median_class_distance <- function(list_classes_X, calculate_between_class_dist = FALSE){

  if(!calculate_between_class_dist){
    n_class <- length(list_classes_X)

    distance <- vector()
    count <- 0

    for (i in 1:n_class){
      for (j in (i+1):n_class){
        if(i<j & j<=n_class){
          X1 <- as.matrix(as.data.frame(list_classes_X[i]))
          X2 <- as.matrix(as.data.frame(list_classes_X[j]))

          for (s in 1: dim(X1)[1]){
            for (k in 1:dim(X2)[1]){
              dist <- norm(as.matrix(X1[s,]-X2[k,]), type ='2')
              distance <- c(distance, dist)
              count <- count + 1

            }
          }

        }

      }

    }

    medianD <- list(medianD =median(distance), num_pair = count)
    return(medianD)
  }


  else if(calculate_between_class_dist){
    n_class <- length(list_classes_X)

    distance <- vector()
    class_compactness <- rep(Inf, n_class)
    sum_class_distance <- rep(0, n_class)

    for(m in 1:n_class){
      class_compactness[m] <- within_class_compactness(list_classes_X[m])
    }

    count <- 0

    for (i in 1:n_class){
      for (j in (i+1):n_class){
        if(i<j & j<=n_class){
          X1 <- as.matrix(as.data.frame(list_classes_X[i]))
          X2 <- as.matrix(as.data.frame(list_classes_X[j]))

          between_class_dist <- vector()
          class_distance <- Inf

          for (s in 1: dim(X1)[1]){
            for (k in 1:dim(X2)[1]){
              dist <- norm(as.matrix(X1[s,]-X2[k,]), type ='2')
              between_class_dist <- c(between_class_dist, dist)
              distance <- c(distance, dist)
              count <- count + 1
            }
          }

          class_distance <- min(between_class_dist)

          sum_class_distance[i] <- sum_class_distance[i] + (class_distance/class_compactness[i])
          sum_class_distance[j] <- sum_class_distance[j] + (class_distance/class_compactness[j])

        }

      }

    }


    ave_sum_class_distance <- sum_class_distance/n_class
    median_separated_class <- find_median_index(ave_sum_class_distance)

    result <- list(medianD =median(distance), ave_sum_class_distance = ave_sum_class_distance, median_separated_class = median_separated_class, num_pair = count)
    return(result)
  }

}
