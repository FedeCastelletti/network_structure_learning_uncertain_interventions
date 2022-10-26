gen.dataset = function(A.true, B, Sigma, I.cal, sigma.int, n_all){
  
  ###########
  ## INPUT ##
  ###########
  
  # A.true    : (q,q) adjacency matrix of the true DAG
  # B         : (q,q) matrix with regression coefficients in the observational data generating model X = BX + eps
  # Sigma     : (q,q) diagonal matrix with conditional variances of of X_1, ..., X_q
  # I.cal     : (q,K) 0-1 matrix indexing the intervened nodes
  # sigma.int : (K,1) vector with variances of the intervention densities of X_I
  # n         : (K,1) vector with sample sizes n_k
  
  ############
  ## OUTPUT ##
  ############
  
  # X       : a (n,q) dataset comprising K (n_k,q) interventional datasets
  # D       : a (n,q) design matrix (associate each observation to the target node) and interventional dataset
  
  
  gen.int.data = function(A.true, B, Sigma, I.k = NULL, sigma_I = NULL, n_k, k){
    
    A.int = A.true
    A.int[,I.k] = 0   # I remove all edges "pointing" to the intervened node (all j --> I)
    
    B_I = A.int*B # remove the regression coefficients associated to  j -> I
    
    Sigma_I = Sigma
    
    for(h in I.k){
      
      Sigma_I[h,h] = sigma_I
      
    }
    
    
    diag(B_I) = 1
    
    S_I = solve(t(B_I))%*%Sigma_I%*%solve(B_I)
    
    X = rmvnorm(n_k, rep(0, q), S_I)
    
    q = ncol(A.true)
    
    m = colMeans(X)
    
    X = t((t(X) - m))
    
    return(X = X)
    
  }
  
  K = length(n_all)
  
  out_X = lapply(1:K, function(k) gen.int.data(A.true, B, Sigma, I.k = which(I.cal[,k] == 1), sigma_I = sigma.int[k], n_k = n_all[k], k))
  out_D = lapply(1:K, function(k) matrix(0,n_all[k],q))
  
  for(k in 1:K){
    
    out_D[[k]][,which(I.cal[,k] == 1)] = k
    
  }
  
  return(list(X = do.call(rbind, out_X), D = do.call(rbind, out_D)))
  
}
