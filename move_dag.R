#########################
## Auxiliary functions ##
#########################

n.edge = function(A){
  length(which(A[lower.tri(A)] == 1 | t(A)[lower.tri(A)] == 1))
}

names   = c("action","test","x","y")
actions = c("id","dd","rd")

# 3 types of possible moves indexed by {1,2,3}

# 1: insert a directed edge (id)
# 2: delete a directed edge (dd)
# 3: reverse a directed edge (rd)

# A     : adjacency matrix of the input DAG
# nodes : nodes (x,y) involved in the action


id = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 1
  return(A)
}

dd = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0
  return(A)
}

rd = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0 
  A[y,x] = 1
  return(A)
}

###################
## Main function ##
###################

move = function(A, q = q){
  
  ## Perform a random local move from the input DAG (with adjacency matrix A) to an adjacent DAG
  
  ###########
  ## INPUT ##
  ###########
  
  # A : adjacency matrix of the input DAG
  # q : number of vertices
  
  ############
  ## OUTPUT ##
  ############
  
  # A_new         : adjacency matrix of the sampled (new) DAG
  # type.operator : type of operator applied to the input DAG
  # nodes         : nodes involved in the operation
  
  A_na = A
  diag(A_na) = NA
  
  id_set = c()
  dd_set = c()
  rd_set = c()
  
  # Set of possible nodes for id
  
  set_id = which(A_na == 0, TRUE)
  
  if(length(set_id) != 0){
    id_set = cbind(1, rbind(set_id))
  }
  
  # Set of possible nodes for dd
  
  set_dd = which(A_na == 1, TRUE)
  
  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }
  
  # Set of possible nodes for rd
  
  set_rd = which(A_na == 1, TRUE)
  
  if(length(set_rd != 0)){
    rd_set = cbind(3, set_rd)
  }
  
  # Set of all possible operators
  
  O = rbind(id_set, dd_set, rd_set)
  
  # Sample one of the possible graphs, each obtained by applying an operator in O
  # Check that the proposed graph is a DAG
  
  repeat {
    
    i = sample(dim(O)[1],1)
    
      act_to_exe  = paste0(actions[O[i,1]],"(A=A,c(",as.vector(O[i,2]),",",as.vector(O[i,3]),"))")
      A_succ      = eval(parse(text = act_to_exe))
      act_to_eval = paste0("is.DAG(A_succ)")
      val = eval(parse(text = act_to_eval))
    
    if (val != 0){
      break
    }
  }
  
  A_new = A_succ
  
  return(list(A_new = A_new, type.operator = O[i,1], nodes = O[i,2:3]))
  
}

