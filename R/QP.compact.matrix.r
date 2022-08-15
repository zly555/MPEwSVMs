#' Compute the compact matrix representation
#'
#' @param Amat Input Matrix
#'
#' @export
#' @return Compact matrix form for quodprog input

QP.compact.matrix <- function(Amat){
  nr <- nrow(Amat)
  nc <- ncol(Amat)
  Amat.compact <- matrix(0, nr, nc)
  Aind <- matrix(0, nr+1, nc)

  for (j in 1:nc){
    index <- (1:nr)[Amat[, j] != 0]
    number <- length(index)
    Amat.compact[1:number, j] <- Amat[index, j]
    Aind[1, j] <- number
    Aind[2:(number+1), j] <- index
  }

  max.number <- max(Aind[1, ])
  Amat.compact <- Amat.compact[1:max.number, ]

  Aind <- Aind[1:(max.number+1), ]
  compact <- list(Amat.compact = Amat.compact, Aind = Aind)

  return(compact)
}
