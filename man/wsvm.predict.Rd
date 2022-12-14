% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wsvm.predict.r
\name{wsvm.predict}
\alias{wsvm.predict}
\title{Predict the label with given wsvm classifer}
\usage{
wsvm.predict(
  x,
  x.train,
  c,
  d,
  kernel = list(type = "linear", param1 = NULL, param2 = NULL)
)
}
\arguments{
\item{x}{New data point}

\item{x.train}{Training data matrix}

\item{c}{Fitted classiifier c}

\item{d}{Fitted classiifier d}

\item{kernel}{Kernel to use}
}
\value{
Predicted label {+1,-1}
}
\description{
Predict the label with given wsvm classifer
}
