mcmc_dag_targets = function(X, S, burn, w, a = NULL, a_k = a_k, b_k = b_k, n_all = n_all){
  
  library(mvtnorm)
  
  ###########
  ## INPUT ##
  ###########
  
  # X : (n,q) data matrix
  
  # S    : number of MCMC iterations
  # burn : burn in period
  
  # w        : prior probability of edge inclusion in p(D)
  # a_k, b_k : hyper-parameters of the Beta(a_k, b_k) prior on the probability of intervention for any node
  # a        : common shape hyper-parameter for the DAG-Wishart prior
  
  # n_all : (K,1) vector with group sample sizes (n_1, ..., n_K)
  
  
  ############
  ## OUTPUT ##
  ############
  
  # X : (n,q) input data matrix
  
  # Graph_post  : (q,q,T) array collecting the T = (S - burn) adjacency matrices of the DAGs visited by the MCMC chain
  # Target_post : (q,K,T) array collecting the T = (S - burn) binary matrices representing the K intervention targets
  
  
  #########################
  ## Auxiliary functions ##
  #########################
  
  source("move_dag.R")
  source("marg_like_dag.R")
  
  q = dim(X)[2]
  
  if(is.null(a)){a = q}
  
  Graph_post   = array(NA, c(q, q, S))
  Targets_post = array(0, c(q, K, S))
  
  num_edges = c()
  
  # Initialize the chain
  
  Graph = matrix(0, q, q)
  I.cal = matrix(0, q, K)
  
  D = lapply(1:K, function(k) matrix(0,n_all[k],q))
  
  for(k in 1:K){
    
    D[[k]][,which(I.cal[,k]==1)] = k
    
  }
  
  D = do.call(rbind, D)
  
  
  Graph_post[,,1] = Graph
  
  out_D = lapply(1:K, function(k) matrix(0,n_all[k],q))
  
  for(k in 1:K){
    
    out_D[[k]][,which(I.cal[,k] == 1)] = k
    
  }
  
  cat("MCMC sampling")
  pb = utils::txtProgressBar(min = 2, max = S, style = 3)
  
  for(s in 1:S){
    
    ## Update the graph conditionally on the targets
    
    Graph_move = move(A = Graph)
    
    Graph_prop = Graph_move$A_new
    nodes_prop = Graph_move$nodes
    
    type.operator = Graph_move$type.operator
    
    # Distinguish 3 cases:
    
    if(type.operator == 1){
      
      # (1) Insert a directed edge
      
      logprior = log(w/(1-w))
      
      j_star = nodes_prop[2]
      
      marg_star = sum(sapply(unique(D[,j_star]), function(k) marg_like_j(j = j_star, dag = Graph_prop, I = which(I.cal[,k] == 1), X = X[D[,j_star] == k,], n = sum(D[,j_star] == k), a)))
      marg      = sum(sapply(unique(D[,j_star]), function(k) marg_like_j(j = j_star, dag = Graph, I = which(I.cal[,k] == 1), X = X[D[,j_star] == k,], n = sum(D[,j_star] == k), a)))
      
      
    }else{
      
      if(type.operator == 2){
        
        # (2) Delete a directed edge
        
        logprior = log((1-w)/w)
        
        j_star = nodes_prop[2]
        
        marg_star = sum(sapply(unique(D[,j_star]), function(k) marg_like_j(j = j_star, dag = Graph_prop, I = which(I.cal[,k] == 1), X = X[D[,j_star] == k,], n = sum(D[,j_star] == k), a)))
        marg      = sum(sapply(unique(D[,j_star]), function(k) marg_like_j(j = j_star, dag = Graph, I = which(I.cal[,k] == 1), X = X[D[,j_star] == k,], n = sum(D[,j_star] == k), a)))
        
      }else{
        
        # (3) Reverse a directed edge
        
        logprior = log(1)
        
        i_star = nodes_prop[1]
        j_star = nodes_prop[2]
        
        marg_star = sum(sapply(unique(D[,i_star]), function(k) marg_like_j(j = i_star, dag = Graph_prop, I = which(I.cal[,k] == 1), X = X[D[,i_star] == k,], n = sum(D[,i_star] == k), a))) +
          sum(sapply(unique(D[,j_star]), function(k) marg_like_j(j = j_star, dag = Graph_prop, I = which(I.cal[,k] == 1), X = X[D[,j_star] == k,], n = sum(D[,j_star] == k), a)))
        
        marg      = sum(sapply(unique(D[,i_star]), function(k) marg_like_j(j = i_star, dag = Graph, I = which(I.cal[,k] == 1), X = X[D[,i_star] == k,], n = sum(D[,i_star] == k), a))) +
          sum(sapply(unique(D[,j_star]), function(k) marg_like_j(j = j_star, dag = Graph, I = which(I.cal[,k] == 1), X = X[D[,j_star] == k,], n = sum(D[,j_star] == k), a)))
        
      }
      
    }
    
    # acceptance ratio
    
    ratio_D = min(0, marg_star - marg + logprior)
    
    # accept move
    
    if(log(runif(1)) < ratio_D){
      
      Graph = Graph_prop
      
    }
    
    Graph_post[,,s] = Graph
    
    
    # Update the targets given the DAG
    
    if((s%%5 == 0) & (s > burn)){
    
      # Permutation of intervention indicator for variable j and joint update of the targets
     
      I.cal.prop = I.cal
      
      node_int = sample(1:q, 1)
      
      ind_j = abs(I.cal.prop[node_int,] - 1)
      
      I.cal.prop[node_int,] = ind_j
      
      D.prop_list = out_D
      
      for(k in 1:K){
        
        D.prop_list[[k]][,which(I.cal.prop[,k] == 1)] = k
        D.prop_list[[k]][,which(I.cal.prop[,k] == 0)] = 0 # check here
        
      }
      
      D.prop = do.call(rbind, D.prop_list)
      
      j_star = node_int
      
      marg_I_prop = sum(sapply(unique(D.prop[,j_star]),
                               function(h) marg_like_j(j = j_star, dag = Graph, I = which(I.cal.prop[,h] == 1), X = X[D.prop[,j_star] == h,], n = sum(D.prop[,j_star] == h), a))) 
      
      marg_I = sum(sapply(unique(D[,j_star]),
                          function(h) marg_like_j(j = j_star, dag = Graph, I = which(I.cal[,h] == 1), X = X[D[,j_star] == h,], n = sum(D[,j_star] == h), a))) 
      
      log_prior = 0
      
      ratio = min(0, marg_I_prop - marg_I + logprior)
      
      # accept move
      
      if(log(runif(1)) < ratio){
        
        I.cal = I.cal.prop
        out_D = D.prop_list
        
      }
      
      D = do.call(rbind, out_D)
      
    
    Targets_post[,,s] = I.cal
    
      
    }else{
      
      
      # Proposal
      
      I.cal.prop = I.cal
      
      target_prop = rep(NA, K)
      
      
      for(k in sample(K)){
        
        # Propose a target I_k (add or remove a target node from the set I_k)
        
        target_prop[k] = sample(1:q, 1)
        
        if(I.cal[target_prop[k],k] == 0){
          
          I.cal.prop[target_prop[k],k] = 1
          
        }else{
          
          I.cal.prop[target_prop[k],k] = 0
          
        }
        
        D.k.prop = matrix(0, n_all[k], q)
        
        D.k.prop[,which(I.cal.prop[,k] == 1)] = k
        
        out_D_tmp = out_D
        
        out_D_tmp[[k]] = D.k.prop
        
        D.prop = do.call(rbind, out_D_tmp)
        
        
        j_star = target_prop[k]
        
        marg_I_k_prop = sum(sapply(unique(D.prop[,j_star]),
                                   function(h) marg_like_j(j = j_star, dag = Graph, I = which(I.cal.prop[,h] == 1), X = X[D.prop[,j_star] == h,], n = sum(D.prop[,j_star] == h), a))) 
        
        marg_I_k = sum(sapply(unique(D[,j_star]),
                              function(h) marg_like_j(j = j_star, dag = Graph, I = which(I.cal[,h] == 1), X = X[D[,j_star] == h,], n = sum(D[,j_star] == h), a))) 
        
        logprior_k = lgamma(a_k + sum(I.cal.prop[,k])) + lgamma(q - sum(I.cal.prop[,k]) + b_k) - 
          lgamma(a_k + sum(I.cal[,k])) - lgamma(q - sum(I.cal[,k]) + b_k)
        
        
        ratio_k = min(0, marg_I_k_prop - marg_I_k + logprior_k)
        
        # accept move
        
        if(log(runif(1)) < ratio_k){
          
          I.cal[,k] = I.cal.prop[,k]
          
        }
        
        out_D[[k]] = matrix(0,n_all[k],q)
        out_D[[k]][,which(I.cal[,k] == 1)] = k
        
        D = do.call(rbind, out_D)
        
      }
      
      Targets_post[,,s] = I.cal
      
    }
    
    utils::setTxtProgressBar(pb, s)
    close(pb)
    
  }
  
  return(out = list(X = X,
                    Graph_post  = Graph_post[,,(burn + 1):S],
                    Targets_post = Targets_post[,,(burn + 1):S]))
  
}

