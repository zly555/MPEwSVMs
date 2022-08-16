
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MPEwSVMs

<!-- badges: start -->
<!-- badges: end -->

The goal of MPEwSVMs is to perform the multiclass probability estimation
with weighted SVMs, using pairwise coupling, OVA and linear time
algorithm from: “Zeng, L. and Zhang, H. H. (2022). Linear algorithms for
nonparametric multiclass probability estimation”.
[(arXiv)](https://arxiv.org/abs/2205.12460).

## Installation

You can install the development version of MPEwSVMs from
[GitHub](https://github.com/zly555/MPEwSVMs) with:

``` r
# install.packages("devtools")
devtools::install_github("zly555/MPEwSVMs")
```

## Example

Here we use a example to demonstrate the multiclass probability
estimation with “MPEwSVMs” package. We use a three class non-linear
simulation example from [Zeng et al](https://arxiv.org/abs/2205.12460v1)
of example 5 (page 18).

### Install the required library

``` r
library(quadprog)
library(data.table, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(tictoc, warn.conflicts = FALSE)
library(MPEwSVMs)
```

### The function to generate the simulation data

``` r
generator_ex5 <- function(n, seed=1){
  
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

### Running setup

``` r
n <- 1000 # Training plus tuning data size
train.tune.ratio <- 0.5 # Training and tuning data ratio
nsim <- 1 # number of MCMC simulations
n_methods <- 6 # Number of multiclass estimation methods: pairwise, one2all, baseline 1, baseline PR 1,  baseline 2, baseline PR 2  

set.seed (9036) # set the random seed
rand_seed <- sample.int(100000, 3*nsim)

exr_no <- 2
sim_no <- 1

# the matrix to hold all the results with all the siumulation results
Result_Matrix <- as.data.frame(matrix(0, nrow = (2*n_methods)*nsim, ncol = 12))
colnames(Result_Matrix) <- c('SimNum','Methods','ElapsedTime', 'TuningMethods','L1_Norm','L2_Norm','EGKL','GKL', 'TestErr_MaxP', 'TestErr_Voting', 'TestErr_MaxPD', 'BaseClass')

#### Matrix to hold probabililty for reliability diagram #######
MP_Matrix <- as.data.frame(matrix(0, nrow = (2*n_methods)*10*n*nsim, ncol = 5))
colnames(MP_Matrix) <- c('p_hat','pred_label','p_real','true_label','method')

## Pairwise and non-pairwise methods
MP_PAIRWISE <- c("pairwise_GKL","pairwise_EGKL","linearLCPR_GKL","linearLCPR_EGKL","linearMDPR_GKL","linearMDPR_EGKL")
MP_NON_PAIRWISE <- c("OVA_GKL","OVA_EGKL","linearLC_GKL","linearLC_EGKL", "linearMD_GKL","linearMD_EGKL")

seed.cnt <- 1
n_records <- 1
nrow_mp_cnt <- 1
```

### Running the simulations

``` r
##############################
### Start MCMC Simulations ###
##############################
tic()
for (sim in 1:nsim){

  train.data <- generator_ex5(n=floor(n*train.tune.ratio), seed = rand_seed[seed.cnt])
  tune.data <- generator_ex5(n=floor((1-train.tune.ratio)*n), seed = rand_seed[seed.cnt+1])
  test.data <- generator_ex5(n=10*n, seed = rand_seed[seed.cnt+2])
  seed.cnt <- seed.cnt + 3

  n_class <- length(levels(factor(train.data$data$ylabel)))
  n_data <- dim(test.data$data)[1]


  ##############
  ## pairwise ##
  ##############
  tic()
  result <- pairwise.svm.probability(train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, kernel = list(type = 'rbf', param1 = 1, param2 = NULL), tuning.criteria = 'GKL')
  # print(result$tuning_matrix)
  # print(result$estimate_prob_binary_classes)
  mp <- multiclass.svm.probability(result$estimate_prob_binary_classes, result$k_class, result$actual_labels, result$actual_prob, result$pair_indexes)
  ep <- evaluate_performance(mp, methods = list(type = 'pairwise'))
  
  T.diff <- toc()
  elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)
  
  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_PAIRWISE[1])

  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "pairwise"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "GKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9:11] <- test_errs
  Result_Matrix[n_records, 12] <- NA
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data

  
  tic()
  result <- pairwise.svm.probability(train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, kernel = list(type = 'rbf', param1 = 1, param2 = NULL), tuning.criteria = 'EGKL')
  # print(result$tuning_matrix)
  # print(result$estimate_prob_binary_classes)
  mp <- multiclass.svm.probability(result$estimate_prob_binary_classes, result$k_class, result$actual_labels, result$actual_prob, result$pair_indexes)
  ep <- evaluate_performance(mp, methods = list(type = 'pairwise'))

  T.diff <- toc()
  elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_PAIRWISE[2])

  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "pairwise"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "EGKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9:11] <- test_errs
  Result_Matrix[n_records, 12] <- NA
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data


  ###################
  ##  One to All   ##
  ###################
  
  tic()
  result <- one2all.svm.probability (train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, kernel = list(type = 'rbf', param1 = 1, param2 = NULL), tuning.criteria = 'GKL')
  # print(result$tuning_matrix)

  ep <- evaluate_performance(result$estimate_multiclass_prob_matrix, methods = list(type ='one2rest'))

  T.diff <- toc()
  elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=result$estimate_multiclass_prob_matrix, n_class=n_class, n_data=n_data, method=MP_NON_PAIRWISE[1])

  
  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "one2rest"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "GKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9] <- test_errs
  Result_Matrix[n_records, 10:12] <- NA
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data

  

  tic()
  result <- one2all.svm.probability (train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, kernel = list(type = 'rbf', param1 = 1, param2 = NULL), tuning.criteria = 'EGKL')
  # print(result$tuning_matrix)

  ep <- evaluate_performance(result$estimate_multiclass_prob_matrix, methods = list(type ='one2rest'))

  T.diff <- toc()
  elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=result$estimate_multiclass_prob_matrix, n_class=n_class, n_data=n_data, method=MP_NON_PAIRWISE[2])

  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "one2rest"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "EGKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9] <- test_errs
  Result_Matrix[n_records, 10:12] <- NA
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data



  ##########################
  ## Liner-Time Algorithm ##
  ##########################

  #########################
  ## largest class size  ##
  #########################
  
  tic()
  result <-pairwise.svm.probability.linear.algorithm(train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, 
    kernel = list(type = 'rbf', param1 = 1, param2 = NULL), linear.time.algorithm = list(type = 'largestClSize'), tuning.criteria = 'GKL')
  # print(result$tuning_matrix)
  # print(result$estimate_prob_binary_classes)

  mp <- multiclass.svm.probability.linear.algorithm(pairwise.prob = result$estimate_prob_binary_classes_base, 
                                                    base_class = result$base_class, 
                                                    k_class = result$k_class, 
                                                    actual_labels = result$actual_labels, 
                                                    actual_prob = result$actual_prob, 
                                                    pair_indexes = result$pair_indexes)
  
  ep <- evaluate_performance(mp, methods = list(type = 'liner_time'))
  
  T.diff <- toc()
  elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)


  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  BaseClass <- result$base_class

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_NON_PAIRWISE[3])
  
  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "largestClSize"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "GKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9] <- test_errs
  Result_Matrix[n_records, 10:11] <- NA
  Result_Matrix[n_records, 12] <- BaseClass
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data


  #### Reconstruct Pairwise table ####
  # recontruct the pairwise probability table based on linear time algorithm generated probability
  tic()
  result_pairwise <- pairwise.svm.probability.linear.time.reconstruct(estimate_prob_binary_classes_base =result$estimate_prob_binary_classes_base, 
                                                                      base_class = result$base_class, 
                                                                      num_class = result$k_class, 
                                                                      pair_indexes = result$pair_indexes, 
                                                                      sum_prob_difference = result$sum_prob_difference, 
                                                                      actual_labels = result$actual_labels, 
                                                                      actual_prob = result$actual_prob)



  mp <- multiclass.svm.probability(pairwise.prob = result_pairwise$estimate_prob_binary_classes, 
                                   k_class = result_pairwise$k_class, 
                                   actual_labels = result_pairwise$actual_labels, 
                                   actual_prob = result_pairwise$actual_prob, 
                                   pair_indexes = result_pairwise$pair_indexes)
  
  ep <- evaluate_performance(mp, methods = list(type = 'pairwise'))
  
  T.diff.reconst <- toc()
  elapsed_time <- round(as.numeric(T.diff.reconst$toc - T.diff.reconst$tic)/60, 3) + elapsed_time


  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_PAIRWISE[3])
  
  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "largestClSize_reconst"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "GKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9:11] <- test_errs
  Result_Matrix[n_records, 12] <- BaseClass
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data


  tic()
  result <-pairwise.svm.probability.linear.algorithm(train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, 
    kernel = list(type = 'rbf', param1 = 1, param2 = NULL), linear.time.algorithm = list(type = 'largestClSize'), tuning.criteria = 'EGKL')
  # print(result$tuning_matrix)
  # print(result$estimate_prob_binary_classes)

  mp <- multiclass.svm.probability.linear.algorithm(pairwise.prob = result$estimate_prob_binary_classes_base, 
                                                    base_class = result$base_class, 
                                                    k_class = result$k_class, 
                                                    actual_labels = result$actual_labels, 
                                                    actual_prob = result$actual_prob, 
                                                    pair_indexes = result$pair_indexes)
  
  ep <- evaluate_performance(mp, methods = list(type = 'liner_time'))

  T.diff <- toc()
  elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  BaseClass <- result$base_class

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_NON_PAIRWISE[4])
  
  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "largestClSize"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "EGKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9] <- test_errs
  Result_Matrix[n_records, 10:11] <- NA
  Result_Matrix[n_records, 12] <- BaseClass
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data


  ####################################
  ##   Reconstruct Pairwise table   ##
  ####################################
  
  # recontruct the pairwise probability table based on linear time algorithm generated probability
  tic()
  result_pairwise <- pairwise.svm.probability.linear.time.reconstruct(estimate_prob_binary_classes_base =result$estimate_prob_binary_classes_base, 
                                                                      base_class = result$base_class, 
                                                                      num_class = result$k_class, 
                                                                      pair_indexes = result$pair_indexes, 
                                                                      sum_prob_difference = result$sum_prob_difference, 
                                                                      actual_labels = result$actual_labels, 
                                                                      actual_prob = result$actual_prob)



  mp <- multiclass.svm.probability(pairwise.prob = result_pairwise$estimate_prob_binary_classes, 
                                   k_class = result_pairwise$k_class, 
                                   actual_labels = result_pairwise$actual_labels, 
                                   actual_prob = result_pairwise$actual_prob, 
                                   pair_indexes = result_pairwise$pair_indexes)
  
  ep <- evaluate_performance(mp, methods = list(type = 'pairwise'))

  T.diff.reconst <- toc()
  elapsed_time <- round(as.numeric(T.diff.reconst$toc - T.diff.reconst$tic)/60, 3) + elapsed_time

  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_PAIRWISE[4])

  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "largestClSize_reconst"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "EGKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9:11] <- test_errs
  Result_Matrix[n_records, 12] <- BaseClass
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data



  #############################
  ## Median class distance   ##
  #############################
  
  tic()
  result <-pairwise.svm.probability.linear.algorithm(train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, 
    kernel = list(type = 'rbf', param1 = 1, param2 = NULL), linear.time.algorithm = list(type = 'medianClDist'), tuning.criteria = 'GKL')
  # print(result$tuning_matrix)
  # print(result$estimate_prob_binary_classes)

  mp <- multiclass.svm.probability.linear.algorithm(pairwise.prob = result$estimate_prob_binary_classes_base, 
                                                    base_class = result$base_class, 
                                                    k_class = result$k_class, 
                                                    actual_labels = result$actual_labels, 
                                                    actual_prob = result$actual_prob, 
                                                    pair_indexes = result$pair_indexes)
  ep <- evaluate_performance(mp, methods = list(type = 'liner_time'))

  T.diff <- toc()
  elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

 # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  BaseClass <- result$base_class

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_NON_PAIRWISE[5])
  
  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "medianClDist"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "GKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9] <- test_errs
  Result_Matrix[n_records, 10:11] <- NA
  Result_Matrix[n_records, 12] <- BaseClass
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data



  ####################################
  ##   Reconstruct Pairwise table   ##
  ####################################
  
  # recontruct the pairwise probability table based on linear time algorithm generated probability
  tic()
  result_pairwise <- pairwise.svm.probability.linear.time.reconstruct(estimate_prob_binary_classes_base =result$estimate_prob_binary_classes_base, 
                                                                      base_class = result$base_class, 
                                                                      num_class = result$k_class, 
                                                                      pair_indexes = result$pair_indexes, 
                                                                      sum_prob_difference = result$sum_prob_difference, 
                                                                      actual_labels = result$actual_labels, 
                                                                      actual_prob = result$actual_prob)



  mp <- multiclass.svm.probability(pairwise.prob = result_pairwise$estimate_prob_binary_classes, 
                                   k_class = result_pairwise$k_class, 
                                   actual_labels = result_pairwise$actual_labels, 
                                   actual_prob = result_pairwise$actual_prob, 
                                   pair_indexes = result_pairwise$pair_indexes)
  ep <- evaluate_performance(mp, methods = list(type = 'pairwise'))

  T.diff.reconst <- toc()
  elapsed_time <- round(as.numeric(T.diff.reconst$toc - T.diff.reconst$tic)/60, 3) + elapsed_time

  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_PAIRWISE[5])

  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "medianClDist_reconst"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "GKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9:11] <- test_errs
  Result_Matrix[n_records, 12] <- BaseClass
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data


  tic()
  result <-pairwise.svm.probability.linear.algorithm(train.data$data, tune.data$data, test.data$data, data_dim =train.data$data_dim, 
    kernel = list(type = 'rbf', param1 = 1, param2 = NULL), linear.time.algorithm = list(type = 'medianClDist'), tuning.criteria = 'EGKL')
  # print(result$tuning_matrix)
  # print(result$estimate_prob_binary_classes)

  mp <- multiclass.svm.probability.linear.algorithm(pairwise.prob = result$estimate_prob_binary_classes_base, 
                                                    base_class = result$base_class, 
                                                    k_class = result$k_class, 
                                                    actual_labels = result$actual_labels, 
                                                    actual_prob = result$actual_prob, 
                                                    pair_indexes = result$pair_indexes)
  ep <- evaluate_performance(mp, methods = list(type = 'liner_time'))

  T.diff <- toc()
  elapsed_time <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

 # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  BaseClass <- result$base_class

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_NON_PAIRWISE[6])

  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "medianClDist"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "EGKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9] <- test_errs
  Result_Matrix[n_records, 10:11] <- NA
  Result_Matrix[n_records, 12] <- BaseClass
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data


  ####################################
  ##   Reconstruct Pairwise table   ##
  ####################################
  
  # recontruct the pairwise probability table based on linear time algorithm generated probability
  
  tic()
  result_pairwise <- pairwise.svm.probability.linear.time.reconstruct(estimate_prob_binary_classes_base =result$estimate_prob_binary_classes_base, 
                                                                      base_class = result$base_class, 
                                                                      num_class = result$k_class, 
                                                                      pair_indexes = result$pair_indexes, 
                                                                      sum_prob_difference = result$sum_prob_difference, 
                                                                      actual_labels = result$actual_labels, 
                                                                      actual_prob = result$actual_prob)



  mp <- multiclass.svm.probability(pairwise.prob = result_pairwise$estimate_prob_binary_classes, 
                                   k_class = result_pairwise$k_class, 
                                   actual_labels = result_pairwise$actual_labels, 
                                   actual_prob = result_pairwise$actual_prob, 
                                   pair_indexes = result_pairwise$pair_indexes)
  ep <- evaluate_performance(mp, methods = list(type = 'pairwise'))
  T.diff.reconst <- toc()
  elapsed_time <- round(as.numeric(T.diff.reconst$toc - T.diff.reconst$tic)/60, 3) + elapsed_time

  # get the probability evaulation performance
  pro_eva <- as.numeric(ep$evaluationMat)

  # get the test error
  test_errs <- as.numeric(ep$TestClassificationError)

  MP_Matrix[nrow_mp_cnt:(nrow_mp_cnt+n_data-1),] <- generate_MP(mp=mp, n_class=n_class, n_data=n_data, method=MP_PAIRWISE[6])

  Result_Matrix[n_records, 1] <- sim
  Result_Matrix[n_records, 2] <- "medianClDist_reconst"
  Result_Matrix[n_records, 3] <- elapsed_time
  Result_Matrix[n_records, 4] <- "EGKL"
  Result_Matrix[n_records, 5:8] <- pro_eva
  Result_Matrix[n_records, 9:11] <- test_errs
  Result_Matrix[n_records, 12] <- BaseClass
  n_records <- n_records + 1
  nrow_mp_cnt <- nrow_mp_cnt + n_data

}


T.diff <- toc()
time_elasped <- round(as.numeric(T.diff$toc - T.diff$tic)/60, 3)

sim_str <- paste(".ex", exr_no, ".sim", sim_no, ".csv", sep = "")

## Save the data output file
fwrite(MP_Matrix, file = paste("mp.rel", sim_str, sep = ""))
fwrite(Result_Matrix, file = paste("result", sim_str, sep = ""))
```
