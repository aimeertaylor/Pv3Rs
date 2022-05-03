# test_correct_graph() copied from the following link
# https://github.com/jwatowatson/RecurrentVivax/blob/master/Genetic_Model/iGraph_functions.R
# This function was written by James Watson and edited by Aimee Taylor

# This implements the Viable graph brute-force search algorithm described in the appendix
test_correct_graph = function(G){
  
  INCORRECT = FALSE # Presume correct until proven incorrect
  
  # First check that all connected components are cliques
  S = components(graph = G)
  for(k in 1:S$no){
    G_sub = induced_subgraph(graph = G, vids = vertex_attr(G)$label[S$membership==k])
    MC = max_cliques(graph = G_sub, min = vcount(G_sub)) 
    if(length(MC) == 0) INCORRECT = TRUE
  }
  
  # If all connected components are cliques, then iterate through all cliques of size 3
  if(!INCORRECT){ 
    all_cliques = cliques(graph = G, min = 3, max=3)
    K = length(all_cliques)
    if(K > 0){
      i = 1
      while(!INCORRECT & i < (K+1) ){ # stop searching as soon as we find it's incorrect
        G_sub = induced_subgraph(graph = G, vids = all_cliques[[i]])
        Score_sub = sum(edge_attr(G_sub)$w) # the sum of the edges tells us if it's correct
        if(Score_sub == 2.5){ INCORRECT = TRUE } # 2.5 is the only possible incorrect score so we reject
        i = i+1
      }
    }
  }
  if(INCORRECT) return(FALSE) else return(TRUE)
}