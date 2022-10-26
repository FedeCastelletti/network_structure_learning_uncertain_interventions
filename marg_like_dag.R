#########################
## Auxiliary functions ##
#########################

library(gRbase)

pa = function (set, object){
  amat <- as(object,"matrix")
  rownames(amat) = colnames(amat) = 1:q
  if (is_ugMAT(amat)) 
    return(NULL)
  pa <- names(which(amat[, set] > 0))
  pa <- setdiff(pa, set)
  if (length(pa)) 
    as.numeric(pa)
  else NULL
}

fa = function(set, object){
  as.numeric(c(set, pa(set, object)))
}

###################
## Main function ##
###################

marg_like_j = function(j, dag, I, X, n, a){
  
  ## Compute the marginal likelihood for node j in dag
  
  g = 1/n
  
  dag.int = dag
  dag.int[,I] = 0
  
  pa_j = pa(j, dag.int)
  
  y = X[,j]
  XX = as.matrix(X[,pa_j])
  
  p_j   = length(pa_j)
  j_pos = q - p_j
  a_j   = (a + q - 2*j_pos + 3)/2 - p_j/2 - 1
  
  if(length(pa_j) == 0){
    
    m = 0.5*a_j*log(0.5*g) - 0.5*(a_j + n)*log(0.5*g + sum(y^2)/2) +
            lgamma(0.5*(a_j + n)) - lgamma(0.5*a_j) - 0.5*n*log(2*pi)
    
  } else{
    
    U_jj = (t(y)%*%XX)%*%chol2inv(chol(diag(g, p_j) + t(XX)%*%XX))%*%(t(XX)%*%y)
    
    m = - 0.5*n*log(2*pi) + 
            0.5*log(g^p_j) - 0.5*log(det(diag(g, p_j) + t(XX)%*%XX)) +
              0.5*a_j*log(0.5*g) - 0.5*(a_j + n)*log(g/2 + sum(y^2)/2 - U_jj/2) +
                lgamma(0.5*(a_j + n)) - lgamma(0.5*a_j)
    
  }
  
  return(m)
    
}

# k = 2
# 
# I = which(I.cal.true[,k] == 1)
# 
# marg_like_j(j = 2, dag = as(true_dag, "matrix"), I, X, n[1], a, g)
# marg_like_j(j = 3, dag = as(true_dag, "matrix"), I, X, n[1], a, g)
# marg_like_j(j = 4, dag = as(true_dag, "matrix"), I, X, n[1], a, g)
# marg_like_j(j = 5, dag = as(true_dag, "matrix"), I, X, n[1], a, g)
