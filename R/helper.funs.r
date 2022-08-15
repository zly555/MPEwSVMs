#' Calculate the EGKL
#'
#' @param y Label vector
#' @param px Estimated probability
#'
#' @export
#' @return EGKL value

wsvm.egkl <- function(y, px) {
  if(!is.vector(y)) y <- as.vector(y)
  if(!is.vector(px)) px <- as.vector(px)

  return(-0.5*mean((y+1)*log(px) + (1-y)*log(1-px)))
}

##
#' Calculate the GKL
#'
#' @param p True probability
#' @param phat Estimated probability
#'
#' @export
#' @return GKL value

gkl.cal <- function(p, phat){
  if(!is.vector(p)) p <- as.vector(p)
  if(!is.vector(phat)) phat <- as.vector(phat)

  return(mean(p*log(p/phat) + (1-p)*log((1-p)/(1-phat))))
}


#' Calculate the density of multivariate normal
#'
#' @param x_mu Mean vector
#' @param Sigma Covariance matrix
#'
#' @export
#' @return Density of multivariate normal distribution

dmnorm <- function(x_mu, Sigma){
  return(exp(-0.5*t(x_mu)%*%solve(Sigma)%*%(x_mu))/sqrt(((2*pi)^dim(Sigma)[1])*det(Sigma)))
}


#' Calculate argmin function
#'
#' @param vec Input R vector
#'
#' @export
#' @return The minimum index with the minimum value of the given vector

argmin <- function(vec) {
  index <- 1:length(z)
  argmin <- min(index[z == min(z)])
  return(argmin)
}


#' Calculate argmax function
#'
#' @param vec Input R vector
#'
#' @export
#' @return The maximum index with the maximum value of the given vector

argmax <- function(vec){
  index <- 1:length(z)
  argmax <- max(index[z == max(z)])
  return(argmax)
}


#' Function to separate the data based on Y labels
#'
#' @param dataXY Full multiclass data set
#' @param data.dim Data dimesnionality
#'
#' @export
#' @return  Two lists, one includes y and other any includes X for each class

multiclass.separation <- function(dataXY, data.dim) {
  if(!is.data.frame(dataXY)) dataXY <- as.data.frame(dataXY)

  ncol <- dim(dataXY)[2]
  colnames(dataXY)[data.dim + 1] <- 'ylabel'

  level <- levels(factor(dataXY$ylabel))

  n_class <- length(level)

  classXY_list = list()
  classX_list = list()


  for (i in 1:n_class){
    class <- paste("class_", i, sep = "")
    classXY_list[[class]] <- dataXY[dataXY$ylabel==level[i],]
    classX_list[[class]] <- dataXY[dataXY$ylabel==level[i],][,1:(ncol-1-n_class)]
  }

  class_list <- list(class_list_XY = classXY_list, class_list_X = classX_list)

  return(class_list)
}



#' Function to find the class with most abundance data points
#'
#' @param dataXY Full multiclass data set
#' @param data.dim Data dimesnionality
#'
#' @export
#' @return Class label with the biggest class

find.biggest.class <- function(dataXY, data.dim){
  if(!is.data.frame(dataXY)) dataXY <- as.data.frame(dataXY)

  colnames(dataXY)[data.dim + 1] <- 'ylabel'
  dataXY$ylabel <- factor(dataXY$ylabel)
  count_per_class <- as.vector(summary(dataXY$ylabel))
  biggest_cl <- which.max(count_per_class)

  return(as.numeric(biggest_cl))
}



#' Function to format data label for one to all probability estimation
#'
#' @param dataXY Full multiclass data set
#' @param data.dim Data dimesnionality
#' @param sel_class Select class as {+1} class in binary wSVMs
#'
#' @export
#' @return Data set with two classes for OVA method

one2all_label_trans <- function(dataXY, data.dim, sel_class){
  if(!is.data.frame(dataXY)) dataXY <- as.data.frame(dataXY)

  colnames(dataXY)[data.dim + 1] <- 'ylabel'

  dataXY$ylabel[which(dataXY$ylabel != sel_class)] = -1
  dataXY$ylabel[which(dataXY$ylabel == sel_class)] = 1

  return(dataXY)
}


#' Function to generate K-size OVA data
#'
#' @param dataXY Full multiclass data set
#' @param data.dim Data dimesnionality
#' @param num_class Total number of classes, as K
#'
#' @export
#' @return A list include K OVA data

data_one2all <- function(dataXY, data.dim, num_class){
  data_one2all <- list()

  for (i in 1:num_class){
    data_temp <- dataXY
    class <- paste("class_", i, sep = "")
    data_one2all[[class]] <- one2all_label_trans(data_temp, data.dim = data.dim, sel_class = i)
  }

  return(data_one2all)
}


#' Change the class label to {+1} and {-1} on pairwise classes
#'
#' @param class1 First pairwise class
#' @param class2 Second pairwise class
#' @param data.dim Data dimensionality
#' @param seed Set random seed
#'
#' @export
#' @return Pairwise classes with correct class label as {+1} and {-1}

covert_to_binary <- function(class1, class2, data.dim, seed){
  if(!is.matrix(class1)) class1 <- as.matrix(class1)
  if(!is.matrix(class2)) class2 <- as.matrix(class2)

  set.seed(seed)

  class1[,data.dim+1] <- 1
  class2[,data.dim+1] <- -1

  binary_classes <- as.matrix(rbind(class1,class2))
  binary_classes <- binary_classes[sample(dim(binary_classes)[1]), ] # shuffle the data

  colnames(binary_classes) <- NULL
  rownames(binary_classes) <- NULL

  return(binary_classes)
}

#' Function to perform Softmax
#'
#' @param c Input vector
#'
#' @export
#' @return Softmax output vector

softmax <- function(c){
  return(exp(c)/sum(exp(c)))
}


#' Normalize the probability vector
#'
#' @param c Probability vector
#'
#' @export
#' @return Normalized probability vector

p_normalize <- function(c){
  return(c/sum(c))
}

#' Function to find the index of median value
#'
#' @param x Input vector
#'
#' @export
#' @return The index of median value

find_median_index <- function(x) {
  return (as.numeric(which.min(abs(x - median(x)))))
}


#' Function to implement cantor mapping for two ordered pair of numbers
#'
#' @param k1 First input
#' @param k2 Second input
#'
#' @export
#' @return Cantor mapping value

cantor_mapping <- function(k1, k2) {return(0.5*(k1+k2)*(k1+k2+1)+k2)}

#
#' Function to reconstruct the pairwise probability table from baseline learning
#'
#' @param pjb Pairwise conditional probability of pj given (j,base)
#' @param pjpb Pairwise conditional probability of pj' given (j',base)
#'
#' @export
#' @return Pairwise conditional probability of pj given (j,j')

pairwise_prob_infer <- function (pjb, pjpb) {return((pjb -pjb*pjpb)/(pjb+pjpb -2*pjb*pjpb))}
