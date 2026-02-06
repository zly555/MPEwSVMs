
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MPEwSVMs

<!-- badges: start -->
<!-- badges: end -->

The goal of MPEwSVMs package is to perform the multiclass probability
estimation with weighted SVMs, using the pairwise coupling, One-vs-All
(OVA), and baseline learning algorithms from: “Zeng, L. and Zhang, H. H.
(2022). Linear Algorithms for Robust and Scalable Nonparametric Multiclass Probability Estimation”. [(arXiv Link)](https://arxiv.org/abs/2205.12460).

## Installation

You can install the development version of MPEwSVMs from
[GitHub](https://github.com/zly555/MPEwSVMs) with:

``` r
# install.packages("devtools")
devtools::install_github("zly555/MPEwSVMs")
```

## Example

Here we use an example to demonstrate the multiclass probability
estimation with “MPEwSVMs” package. We adopt a three-class non-linear
classification problem from [Zeng et
al.](https://arxiv.org/abs/2205.12460v1) simulation example 5 (page 18)
for demonstration.

### Install the Required R Library

``` r
library(quadprog)
library(data.table, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(tictoc, warn.conflicts = FALSE)
library(MPEwSVMs)
```

### Generate the Simulated Data

``` r
generator <- function(n, seed=1){
  
  set.seed (seed)

  x1 <- runif(n, min=-3, max=3)
  x2 <- runif(n, min=-6, max=6)

  x <- cbind(x1,x2)

  data.dim <- dim(x)[2]

  f1 <- function(x) {return(-x[1] +0.1*x[1]^2-0.05*x[2]^2+0.1)}
  f2 <- function(x) {return(-0.2*x[1]^2 +0.1*x[2]^2-0.2)}
  f3 <- function(x) {return(x[1]+0.1*x[1]^2-0.05*x[2]^2+0.1)}

  # calculate the pjx for each class
  pyx <- matrix(0, nrow=n, ncol =3)

  for (i in 1:3){
    pyx[,i] <- apply(x, 1, function(x) {exp(do.call(paste0("f",i), list(x)))/(exp(f1(x))+exp(f2(x))+exp(f3(x)))})
  }

  # generate moltinorminal random variable for y label, 1 indicate pick that class label
  y <- apply(pyx, 1, function(x) rmultinom(n=1, size =1, prob=x))

  # generate y label {1,2,3}
  y <- apply(y*c(1,2,3), 2, function(x) x[x > 0])

  data <- data.frame(cbind(x, y, pyx))
  colnames(data) <- c('x1','x2','ylabel','p1x','p2x','p3x')

  return(list(data = data, data_dim = data.dim))

}
```

### Running Setup

``` r
n <- 1000 # Training plus tuning data size
train.tune.ratio <- 0.5 # Training and tuning data ratio

set.seed (9036) # set the random seed
rand_seed <- sample.int(100000, 3)
```

### Generate the Training, Tuning and Test Data

``` r
# Training data = tuning data = 500, Test data = 10000
train.data <- generator(n=floor(n*train.tune.ratio), seed = rand_seed[1])
tune.data <- generator(n=floor((1-train.tune.ratio)*n), seed = rand_seed[2])
test.data <- generator(n=10*n, seed = rand_seed[3])
```

### Multiclass Probability Estimation with wSVMs

**1. Pairwise coupling with dynamic baseline choosing (P-SVM)**

``` r
  ##############
  ## Pairwise ##
  ##############

# Tuning with EGKL, with RBF kernel 
tic()
  result <- pairwise.svm.probability(train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, kernel = list(type = 'rbf', param1 = 1, param2 = NULL), tuning.criteria = 'EGKL')
  
  # Multiclass probability estimation matrix
  mp <- multiclass.svm.probability(result$estimate_prob_binary_classes, result$k_class, result$actual_labels, result$actual_prob, result$pair_indexes)
  # Evaluation performance matrix
  ep <- evaluate_performance(mp, methods = list(type = 'pairwise'))
    
T.diff <- toc()

elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

# The running time depends on the tuning size and the computation environment
print(elapsed_time) 
#[1] 2.32325

# get the probability evaulation performance
print(ep$evaluationMat)
#   L1_Norm    L2_Norm      EGKL       GKL
# 0.2072355 0.02634248 0.1321192 0.1109536

# get the test error
print(ep$TestClassificationError[1:2])
#  TestErr_MaxP TestErr_Voting
#     0.2398         0.2401
```

**2. One-vs-All (OVA) approach (A-SVM)**

``` r
  ################
  ## One-vs-All ##
  ################

# Tuning with EGKL, with RBF kernel 
tic()
  result <- one2all.svm.probability (train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, kernel = list(type = 'rbf', param1 = 1, param2 = NULL), tuning.criteria = 'EGKL')
  # Evaluation performance maxtrix
  ep <- evaluate_performance(result$estimate_multiclass_prob_matrix, methods = list(type ='one2rest'))

T.diff <- toc()
elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

# The running time depends on the tuning size and the computation environment
print(elapsed_time) 
#[1] 8.194167

# get the probability evaulation performance
print(ep$evaluationMat)
#   L1_Norm    L2_Norm      EGKL       GKL
# 0.2010323 0.02437899 0.1042973 0.1050824

# get the test error
print(ep$TestClassificationError[1:2])
#  TestErr_MaxP TestErr_Voting
#     0.2379             NA
```

**3. Baseline learning approach (B-SVM)**

``` r
  ##############
  ## Baseline ##
  ##############

# Tuning with EGKL, with RBF kernel, choose baseline with method 1 
tic()
  result <-pairwise.svm.probability.linear.algorithm(train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, kernel = list(type = 'rbf', param1 = 1, param2 = NULL), linear.time.algorithm = list(type = 'largestClSize'), tuning.criteria = 'EGKL')

  # Multiclass probability estimation matrix
  mp <- multiclass.svm.probability.linear.algorithm(pairwise.prob = result$estimate_prob_binary_classes_base, 
                                                    base_class = result$base_class, 
                                                    k_class = result$k_class, 
                                                    actual_labels = result$actual_labels, 
                                                    actual_prob = result$actual_prob, 
                                                    pair_indexes = result$pair_indexes)
  # Evaluation performance matrix
  ep <- evaluate_performance(mp, methods = list(type = 'liner_time'))
  
T.diff <- toc()
elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

# The running time depends on the tuning size and the computation environment
print(elapsed_time) 
#[1] 1.505583

# get the probability evaulation performance
print(ep$evaluationMat)
#   L1_Norm    L2_Norm     EGKL       GKL
# 0.2024622 0.02562124 0.133561 0.1041888

# get the test error
print(ep$TestClassificationError[1:2])
#  TestErr_MaxP TestErr_Voting
#     0.2403             NA

# get the select baseline class
print(result$base_class) 
# [1] 2
```

**4. Baseline learning with pairwise reconstruction (BP-SVM)**

``` r
  ######################################
  ## Baseline Pairwise Reconstruction ##
  ######################################

tic()
# Reconstruct the pairwise probability table based on baseline learning algorithm from 3
  result_pairwise <- pairwise.svm.probability.linear.time.reconstruct(estimate_prob_binary_classes_base =result$estimate_prob_binary_classes_base, 
                                                                        base_class = result$base_class, 
                                                                        num_class = result$k_class, 
                                                                        pair_indexes = result$pair_indexes, 
                                                                        sum_prob_difference = result$sum_prob_difference, 
                                                                        actual_labels = result$actual_labels, 
                                                                        actual_prob = result$actual_prob)


  # Multiclass probability estimation matrix
  mp <- multiclass.svm.probability(pairwise.prob = result_pairwise$estimate_prob_binary_classes, 
                                   k_class = result_pairwise$k_class, 
                                   actual_labels = result_pairwise$actual_labels, 
                                   actual_prob = result_pairwise$actual_prob, 
                                   pair_indexes = result_pairwise$pair_indexes)
  
  # Evaluation performance matrix
  ep <- evaluate_performance(mp, methods = list(type = 'pairwise'))
  
T.diff.reconst <- toc()
# The time includes the baseline learning from 3, and pairwise table reconstruction 
elapsed_time <- round(as.numeric(T.diff.reconst$toc - T.diff.reconst$tic)/60, 3) + elapsed_time

# The running time depends on the tuning size and the computation environment
print(elapsed_time) 
#[1] 1.508002

# get the probability evaulation performance
print(ep$evaluationMat)
#   L1_Norm    L2_Norm     EGKL       GKL
# 0.2024622 0.02562124 0.133561 0.1041888

# get the test error
print(ep$TestClassificationError[1:2])
#  TestErr_MaxP TestErr_Voting
#     0.2402         0.2388
```

### Conclusion

We developed the R package to perform the multiclass probability
estimation with wSVMs. The statistical performance needs to perform
Monte Carlo simulations. The baseline learning has close performance
with pairwise coupling method with fast computation time, and OVA
learning has the best performance. More examples and detail complexity
analysis check the paper at [Zeng et
al.](https://arxiv.org/abs/2205.12460).
