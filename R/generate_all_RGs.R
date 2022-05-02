# generate_all_models_3Ts() copied from the following link
# https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R
# The function generate_all_models_3Ts() was written by James Watson and edited by Aimee Taylor
generate_all_models_3Ts = function(M1=2, M2=2, M3=2){
  if(M1>4 | M2>4 | M3>4 | sum(c(M1,M2,M3))>6) stop('Whooh there, too complex')
  
  # # Set up all the within timepoint models
  
  # AT: I think this is a only slightly more elegant way for within time points
  T1models = if(M1 == 1){matrix(0,1,1)}else{permutations(n = 2, r = choose(M1,2), v = c(0,.5), repeats.allowed = TRUE)}
  T2models = if(M2 == 1){matrix(0,1,1)}else{permutations(n = 2, r = choose(M2,2), v = c(0,.5), repeats.allowed = TRUE)}
  T3models = if(M3 == 1){matrix(0,1,1)}else{permutations(n = 2, r = choose(M3,2), v = c(0,.5), repeats.allowed = TRUE)}
  
  # and the across timepoint models
  T12models = as.matrix(expand.grid( rep(list(c(0,.5, 1)), M1*M2) ))
  T23models = as.matrix(expand.grid( rep(list(c(0,.5, 1)), M2*M3) ))
  T13models = as.matrix(expand.grid( rep(list(c(0,.5, 1)), M1*M3) ))
  
  # compute size of problem
  KK1 = nrow(T1models)
  KK2 = nrow(T2models)
  KK3 = nrow(T3models)
  
  KK12= nrow(T12models)
  KK23= nrow(T23models)
  KK13= nrow(T13models)
  
  KK = (KK1 * KK2 * KK3) * (KK12 * KK23 * KK13)
  
  print(paste('complexity of problem is: ', KK))
  
  
  Pbar = txtProgressBar(min=1, max = KK)
  ind = 1
  correct_ind = 1; correct_models = list()
  
  # Now for the iteration with brute force...
  for(i1 in 1:KK1){
    for(i2 in 1:KK2){
      for(i3 in 1:KK3){
        for(i12 in 1:KK12){
          for(i23 in 1:KK23){
            for(i13 in 1:KK13){
              # plotting progress...
              setTxtProgressBar(Pbar, ind); ind = ind+1
              # Lets make the sub adjacency matrices
              Adj_T1 = matrix(rep(0, M1^2), ncol = M1)
              Adj_T2 = matrix(rep(0, M2^2), ncol = M2)
              Adj_T3 = matrix(rep(0, M3^2), ncol = M3)
              
              Adj_T1[lower.tri(Adj_T1)] = as.vector(T1models[i1,])
              Adj_T2[lower.tri(Adj_T2)] = as.vector(T2models[i2,])
              Adj_T3[lower.tri(Adj_T3)] = as.vector(T3models[i3,])
              
              Adj_T1_T2 = matrix(T12models[i12,], ncol = M1) # this has M1 cols and M2 rows
              Adj_T2_T3 = matrix(T23models[i23,], ncol = M2) # this has M2 cols and M3 rows
              Adj_T1_T3 = matrix(T13models[i13,], ncol = M1) # this has M1 cols and M3 rows
              
              Adj_all = matrix(rep(0,(M1+M2+M3)^2), ncol = M1+M2+M3)
              Adj_all[1:M1, 1:M1] = Adj_T1
              Adj_all[(M1+1):(M1+M2), (M1+1):(M1+M2)] = Adj_T2
              Adj_all[(M1+M2+1):(M1+M2+M3), (M1+M2+1):(M1+M2+M3)] = Adj_T3
              
              Adj_all[(M1+1):(M1+M2), 1:M1] = Adj_T1_T2
              Adj_all[(M1+M2+1):(M1+M2+M3), (M1+1):(M1+M2)] = Adj_T2_T3
              Adj_all[(M1+M2+1):(M1+M2+M3), 1:M1] = Adj_T1_T3
              
              Adj_all = Adj_all + t(Adj_all)
              colnames(Adj_all) = 1:(M1+M2+M3)
              
              # turn this into an igraph item
              G1 = graph_from_adjacency_matrix(adjmatrix = Adj_all, mode='undirected', 
                                               diag=F, weighted = 'w', add.colnames = T)
              # set the group attributes for G1
              G1 = set_vertex_attr(G1, "group", value = c(rep(1,M1), rep(2,M2), rep(3,M3)))
              G1 = set_vertex_attr(G1, "label", value = 1:(M1+M2+M3))
              
              # test correctness and store if correct
              if(test_correct_graph(G = G1)){
                correct_models[[correct_ind]] = G1
                correct_ind=correct_ind+1
              }
            }
          }
        }
      }
    }
  }
  return(correct_models)
}
