#' Estimate the binary class probability
#'
#' @param v Predicted Labels
#' @param p Vector of weights
#'
#' @return Estimated probability for {+1} class

prob.estimate<- function(v, p){
  df <- as.data.frame(cbind(v,p))
  colnames(df) <- c("labels","Pi")

  L_max<-aggregate(Pi ~ labels, data = df, max)
  L_min <-aggregate(Pi ~ labels, data = df, min)

  Pi_ls <- L_max$Pi[L_max$labels==1]
  Pi_us <- L_min$Pi[L_max$labels==-1]

  return((Pi_ls+Pi_us)/2)
}
