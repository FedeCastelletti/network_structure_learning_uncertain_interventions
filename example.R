## Load required libraries

library(BCDAG)
library(mvtnorm)

## Generate a DAG (through its adjacency matrix A) and intervention targets

q  = 10   # number of nodes
pr = 2/q  # probability of edge inclusion

set.seed(123)
A = rDAG(q, w = pr)

## Generate DAG parameters B (regression coefficients) and Sigma (conditional variances)

B = A*matrix(runif(q*q, 0.1, 1), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE)

Sigma = diag(1, q, q)

## Number of datasets (intervention contexts) and corresponding sample sizes

K     = 4
n_all = rep(100, K)

## Generate intervention targets I.cal

I.cal.true = matrix(0, q, K)

set.seed(1234)

for(k in 1:K){
  
  I.cal.true[sample(1:q, 2),k] = 1
  
}

## Fix intervention conditional variances

sigma.int = rep(0.1, K)


#######################
## Generate the data ##
#######################

source("gen_data.R")

out_data = gen.dataset(A.true = A, B, Sigma, I.cal = I.cal.true, sigma.int, n = n_all)

X = out_data$X


###########################################################################
## Run the MCMC for posterior inference on DAGs and intervention targets ##
###########################################################################

## Fix hyperparameters of Beta(a_k,b_k) prior on the probabilities of intervention and prior probability of edge inclusion w

a_k = 1/q
b_k = 1
w   = 0.1

## Fix number of MCMC iterations and burn-in

S    = 6000
burn = 1000

source("mcmc_dags_targets.r")

out = mcmc_dag_targets(X = X, S = S, burn = burn, w = w, a = NULL, a_k = a_k, b_k = b_k, n_all = n_all)

dim(out$Graph_post)
dim(out$Targets_post)


#########################
## Posterior summaries ##
#########################

## Posterior probabilities of intervention

I.cal.hat = apply(X = out$Targets_post, FUN = mean, MARGIN = c(1,2))
round(I.cal.hat, 2)

## Posterior probabilities of edge inclusion

P_hat = apply(X = out$Graph_post, FUN = mean, MARGIN = c(1,2))
round(P_hat, 2)

## Median Probability DAG Model (MPM)

A_hat = round(P_hat)
A_hat

## Compare true DAG and posterior probabilities of edge inclusion

par(mfrow = c(1,2))

image(A)
image(P_hat)

## Compare true target and posterior probabilities of intervention

par(mfrow = c(1,2))

image(I.cal.true)
image(I.cal.hat)
