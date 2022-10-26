These R codes implement the Bayesian methodology of Castelletti and Peluso (2022, Journal of the American Statistical Association) for structure learning of DAGs from interventional data with unknown intervention targets

Specifically:

gen_data.R         : generates a dataset containing multivariate interventional data from an input DAG, intervention targets and DAG parameters (see "example.R" for its usage)
mcmc_dag_targets.R : contains the main MCMC algorithm for posterior inference on DAGs and intervention targets

move_dag.R      : performs one move from a DAG to an adjacent DAG (i.e. implements the proposal distribution over the space of DAGs)
marg_like_dag.R : computes the (log)marginal likelihood of a DAG model relative to a node in the DAG

example.R : implements mcmc_dag_targets.R on a simulated dataset generated using gen_data.R