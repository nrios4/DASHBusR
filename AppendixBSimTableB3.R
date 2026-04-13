# AppendixBSimTableB3.R
# This R code runs the simulations used to produce Table B3 in Appendix B.

rm(list = ls()) # clear environment

# Set working directory to the root directory of the archive, e.g.,
# setwd("C:/Users/riosn/OneDrive/Desktop/DASHBusRCodeR1")

source("source_code.R") # load custom functions
library(crossdes) # for find.BIB

# Simulation Settings
tau2 = 0.01 # 0.01
m = 8 
p = 2 

num_edges = choose(m,2)/2  
# num_edges = ceiling(2*choose(m,2)/3)

n = num_edges + p + 1 
# n = 2*(num_edges + p) + 1
# n = 3*(num_edges + p) + 1

N = 2^p*num_edges*(num_edges-1)
L = 6 # 4,5,6 cap on path length
n_sim = 100 # number of simulations
R = diag(num_edges + p)*tau2


set.seed(123) 

SA_bayesD = numeric(n_sim)
fedorov_bayesD = numeric(n_sim)
rand_bayesD = numeric(n_sim)
coord_bayesD = numeric(n_sim)
SA_CPU = numeric(n_sim)
fedorov_CPU = numeric(n_sim)
rand_CPU = numeric(n_sim)
coord_CPU = numeric(n_sim)

for(iter in 1:n_sim){
  
  print(paste("Beginning Simulation",iter, sep = " "))
  
  # Step 1: Generate the graph
  graph_found = FALSE
  while(!graph_found){
    G = generate_random_connected_graph(m = m, n = num_edges)
    V(G)$name = LETTERS[1:m]
    edge_names <- apply(as_edgelist(G), 1, function(e) paste(sort(e), collapse = "_"))
    
    # Step 2: Find all paths of length L and convert them to X_edge design matrices
    all_paths = get_all_paths_with_L_edges(G,L)
    graph_found = length(all_paths) > 0
    
  }
  
  
  num_paths = length(all_paths)
  X_paths = matrix(0, nrow = num_paths, ncol = num_edges)
  colnames(X_paths) = edge_names
  
  for(i in 1:num_paths){
    
    cols_in_path = which(edge_names %in% all_paths[[i]])
    X_paths[i,cols_in_path] = 1
  }
  
  # Step 3: enumerate all possible design points under the constraints
  factorial_design = as.matrix(generate_2p_factorial_design(p))
  X_all = matrix(NA, nrow = num_paths*(2^p), ncol = num_edges + p)
  counter = 1
  for(i in 1:num_paths){
    x_edges = X_paths[i,]
    for(j in 1:(2^p)){
      X_all[counter,1:num_edges] = x_edges
      X_all[counter,-(1:num_edges)] = factorial_design[j,]
      counter = counter + 1
    }
    
  }
  
  # Step 4: Find optimal design using Theorem 1
  blockdesign = find.BIB(trt = num_edges, b = num_edges*(num_edges-1), k = L)
  X_optpaths = matrix(0, nrow = nrow(blockdesign), ncol = num_edges)
  colnames(X_optpaths) = edge_names
  
  for(i in 1:nrow(blockdesign)){
    
    X_optpaths[i,blockdesign[i,]] = 1
  }
  
  
  
  Xopt = matrix(NA, nrow = nrow(blockdesign)*(2^p), ncol = num_edges + p)
  counter = 1
  for(i in 1:nrow(blockdesign)){
    x_edges = X_optpaths[i,]
    for(j in 1:(2^p)){
      Xopt[counter,1:num_edges] = x_edges
      Xopt[counter,-(1:num_edges)] = factorial_design[j,]
      counter = counter + 1
    }
    
  }
  
  optD = get_bayes_D(Xopt,R)
  
  # Step 4: Generate an initial design with n runs
  X_init = X_all[sample(1:nrow(X_all),n),]
  
  # Step 5: Use SA to find a design of size n
  SA_start = Sys.time()
  # SA_out <- SA_design(X_init = X_init, p=p, L=L, R=R, X_paths=X_paths, n_iter = 50000)
  if(num_edges == choose(m,2)/2){
    SA_out <- SA_design(X_init = X_init, p=p, L=L, R=R, X_paths=X_paths, n_iter = 50000)
  }
  if(num_edges == ceiling(2*choose(m,2)/3)){
    if(p == 2){
      SA_out <- SA_design(X_init = X_init, p=p, L=L, R=R, X_paths=X_paths, n_iter = 150000)  
    } else{
      SA_out <- SA_design(X_init = X_init, p=p, L=L, R=R, X_paths=X_paths, n_iter = 250000)
    }
    
  }
  SA_bayesD[iter] = SA_out$phi/optD
  SA_end = Sys.time()
  SA_CPU[iter] = as.numeric(SA_end - SA_start, units = "secs")
  
  # Step 6: Use optFederov to find a design of size n 
  fedorov_start = Sys.time()
  fedorov_out = fedorov_bayes(X_init = X_init, X_all = X_all, R=R, n_repeat = 2)
  fedorov_bayesD[iter] = fedorov_out$phi/optD
  fedorov_end = Sys.time()
  fedorov_CPU[iter] = as.numeric(fedorov_end - fedorov_start, units = "secs")
  
  
  # Step 7: Find a completely random design for comparison
  rand_start = Sys.time()
  X_rand = X_all[sample(1:nrow(X_all),n),]
  randD = get_bayes_D(X_rand, R)
  rand_bayesD[iter] = randD/optD
  rand_end = Sys.time()
  rand_CPU[iter] = as.numeric(rand_end - rand_start, units = "secs")  
  
  # Step 8: Use coordinate exchange algorithm for comparison
  coord_start = Sys.time()
  coord_out = coordinate_exchange_bayes(X_init = X_init, X_all = X_all, R = R,
                                        n_repeat = 2)
  coord_bayesD[iter] = coord_out$phi/optD
  coord_end = Sys.time()
  coord_CPU[iter] = as.numeric(coord_end - coord_start, units = "secs")
  
}



res_mat = cbind(SA_bayesD, fedorov_bayesD, rand_bayesD, coord_bayesD,
                SA_CPU, fedorov_CPU, rand_CPU, coord_CPU)
round(colMeans(res_mat),3)


# change directory to the SimResults folder and save results here
setwd("./SimResults")
namestr = paste("simresAppendixBTableB3_m",m,"_n",n,"_p",p,"num_edges",num_edges,"L",L,".csv",sep="")
write.csv(res_mat,namestr)
setwd("..") # return to root directory to make it so other files can be run more easily

