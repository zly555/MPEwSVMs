#'  Predict the label with given wsvm classifer
#'
#' @param x New data point
#' @param x.train Training data matrix
#' @param c Fitted classiifier c
#' @param d Fitted classiifier d
#' @param kernel Kernel to use
#'
#' @export
#' @return Predicted label {+1,-1}

wsvm.predict <- function(x, x.train, c, d, kernel = list(type = 'linear', param1 = NULL, param2 = NULL)){

  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.matrix(c)) c <- as.matrix(c)
  if(!is.numeric(d)) d <- as.numeric(d)
  if(!is.matrix(x.train)) x.train <- as.matrix(x.train)
  # print('start predict')
  # print( dim(c))

  K <- wsvm.kernel(x, x.train, kernel)
  d <- as.numeric(d)

  return(as.vector(sign(K%*%c + d)))

}
