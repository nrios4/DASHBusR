# Section4.2Sim.R
# This R code runs the simulations described in Section 4.2.

rm(list = ls()) # clear environment

# Set working directory to the root directory of the archive, e.g.,
# setwd("C:/Users/riosn/OneDrive/Desktop/DASHBusRCodeR1")

source("source_code.R") # load custom functions

# Simulation settings - change these to get a different simulation
scenario = "C"  # one of "A","B","C"
m = 27 # one of 25,27
n = 100 # one of 50,100 
n_sim = 100 # number of simulated datasets


namestr = paste("Sim4.2v2scenario",scenario,"m",m,"n",n,".csv",sep = "")

tau2 = 100
mu = 10 # changed mu from 0 -> 10
num_edges = 30  # |E|
p = 2
Rinv = diag(num_edges + p)*tau2
R = diag(num_edges + p)*(1/tau2)
L = 2^p + 1
sigma2 = 0.01

mean_CRPS_prior1 = numeric(n_sim)  
mean_CRPS_prior2 = numeric(n_sim)  

mean_coverage_prior1 = numeric(n_sim)
mean_coverage_prior2 = numeric(n_sim)

mean_CRPS_prior1_rand = numeric(n_sim)
mean_CRPS_prior2_rand = numeric(n_sim)
mean_coverage_prior1_rand = numeric(n_sim)
mean_coverage_prior2_rand = numeric(n_sim)


lambda1 = 0.05
lambda2 = 0.01
lambda3 = 0.01
changepoint = floor(num_edges/2)
R2 = diag(c(rep(lambda1, changepoint), rep(lambda2, num_edges - changepoint), rep(lambda3, p)))



set.seed(123)


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
  
  print("Generated Graph")
  
  num_paths = length(all_paths)
  X_paths = matrix(0, nrow = num_paths, ncol = num_edges)
  colnames(X_paths) = edge_names
  
  for(i in 1:num_paths){
    
    cols_in_path = which(edge_names %in% all_paths[[i]])
    X_paths[i,cols_in_path] = 1
  }
  
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
  
  # Step 3: Generate an initial design with n runs
  X_init = X_all[sample(1:nrow(X_all),n),]
  
  print("Found initial design")
  
  # Step 4: Use SA to find a design of size n
  SA_out <- SA_design(X_init = X_init, p=p, L=L, R=R, X_paths=X_paths, n_iter = 100000)  
  
  print("Find optimal design")
  
  # Step 5: Generate responses according to the design and the scenarios
  X = SA_out$model.matrix
  y = numeric(n)
  beta1 = rnorm(1,2,sqrt(0.01))
  beta2 = rnorm(1,-2,sqrt(0.01))
  edge_costs = numeric(num_edges)
  if(scenario == "A"){
    true_mus = numeric(num_edges)
    true_sds = rep(sqrt(0.01),num_edges)
    counter = 1
    for(j in 1:(m-1)){
      for(k in (j+1):m){
        if(G[j,k] == 1){
          # print(paste("(",j,",",k,")",sep=""))
          true_mus[counter] = 10*log(j*k)
          edge_costs[counter] = rnorm(1,true_mus[counter],sqrt(0.01))
          counter = counter + 1
        }
        
      }
    }
    
  }
  if(scenario == "B"){
    
    true_mus = numeric(num_edges)
    true_sds = numeric(num_edges)
    counter = 1
    for(j in 1:(m-1)){
      for(k in (j+1):m){
        if(G[j,k] == 1){
          true_mus[counter] = 10*log(j*k)
          true_sds[counter] = sqrt(j*k/25)
          edge_costs[counter] = rnorm(1,true_mus[counter],true_sds[counter])
          counter = counter + 1
        }
        
      }
    }
    
    
  }
  if(scenario == "C"){
    
    true_mus = numeric(num_edges)
    true_sds = numeric(num_edges)
    counter = 1
    for(j in 1:(m-1)){
      for(k in (j+1):m){
        if(G[j,k] == 1){
          true_mus[counter] = 10*log(j*k)
          edge_costs[counter] = abs(rt(1,df=4, ncp = 1/true_mus[counter])) 
          counter = counter + 1
        }
        
      }
    }
    
  }
  # print(paste("edge costs:", edge_costs, sep = ""))
  true_theta = matrix(c(edge_costs, beta1, beta2),ncol = 1)
  y = X%*%true_theta + rnorm(n, mean = 0, sd = sqrt(sigma2))
  
  # fit the proposed model
  Y = matrix(y, ncol = 1)
  
  
  theta_hat = solve(t(X)%*%X + R)%*%(t(X)%*%Y + t(R)%*%rep(mu,num_edges+p)) 
  var_hats = diag(solve(t(X)%*%X + R))
  
  theta_hat2 = solve(t(X)%*%X + R2)%*%(t(X)%*%Y + t(R2)%*%rep(mu,num_edges+p)) 
  var_hats2 = diag(solve(t(X)%*%X + R2))
  
  # find credible intervals
  z_crit = qnorm(0.025, lower.tail = FALSE)
  lower_credible_interval = theta_hat - z_crit*sqrt(var_hats)
  upper_credible_interval = theta_hat + z_crit*sqrt(var_hats)
  mean_coverage_prior1[iter] = mean( (lower_credible_interval[1:num_edges] <= true_theta[1:num_edges]) & 
                                       (true_theta[1:num_edges] <= upper_credible_interval[1:num_edges])   )
  
  lower_credible_interval2 = theta_hat2 - z_crit*sqrt(var_hats2)
  upper_credible_interval2 = theta_hat2 + z_crit*sqrt(var_hats2)
  mean_coverage_prior2[iter] = mean( (lower_credible_interval2[1:num_edges] <= true_theta[1:num_edges]) & 
                                       (true_theta[1:num_edges] <= upper_credible_interval2[1:num_edges])   )
  
  # find mean CRPS over all edges
  CRPS_vec_prior1 = numeric(num_edges)
  CRPS_vec_prior2 = numeric(num_edges)
  for(j in 1:num_edges){
    
    CRPS_vec_prior1[j] = CRPS_analytic_normal(Y=edge_costs[j],means = theta_hat[j], 
                                              variances =  var_hats[j])
    CRPS_vec_prior2[j] = CRPS_analytic_normal(Y=edge_costs[j],means = theta_hat2[j], 
                                              variances =  var_hats2[j])
    
  }
  
  mean_CRPS_prior1[iter] = mean(CRPS_vec_prior1)
  mean_CRPS_prior2[iter] = mean(CRPS_vec_prior2)
  
  # repeat for random design
  # find random sample
  X_rand = X_all[sample(1:nrow(X_all),n),]
  # re-generate the resposnes y_rand
  y_rand = X_rand%*%true_theta + rnorm(n, mean = 0, sd = sqrt(sigma2))
  Y_rand = matrix(y_rand, ncol = 1)
  # re-fit the proposed model for the random sample
  
  theta_hat_rand = solve(t(X_rand)%*%X_rand + R)%*%(t(X_rand)%*%Y_rand + t(R)%*%rep(mu,num_edges+p)) 
  var_hats_rand = diag(solve(t(X_rand)%*%X_rand + R))
  
  theta_hat2_rand = solve(t(X_rand)%*%X_rand + R2)%*%(t(X_rand)%*%Y_rand + t(R2)%*%rep(mu,num_edges+p)) 
  var_hats2_rand = diag(solve(t(X_rand)%*%X_rand + R2))
  
  # find credible intervals
  z_crit = qnorm(0.025, lower.tail = FALSE)
  lower_credible_interval_rand = theta_hat_rand - z_crit*sqrt(var_hats_rand)
  upper_credible_interval_rand = theta_hat_rand + z_crit*sqrt(var_hats_rand)
  mean_coverage_prior1_rand[iter] = mean( (lower_credible_interval_rand[1:num_edges] <= true_theta[1:num_edges]) & 
                                            (true_theta[1:num_edges] <= upper_credible_interval_rand[1:num_edges])   )
  
  lower_credible_interval2_rand = theta_hat2_rand - z_crit*sqrt(var_hats2_rand)
  upper_credible_interval2_rand = theta_hat2_rand + z_crit*sqrt(var_hats2_rand)
  mean_coverage_prior2_rand[iter] = mean( (lower_credible_interval2_rand[1:num_edges] <= true_theta[1:num_edges]) & 
                                            (true_theta[1:num_edges] <= upper_credible_interval2_rand[1:num_edges])   )
  
  # find mean CRPS over all edges
  CRPS_vec_prior1_rand = numeric(num_edges)
  CRPS_vec_prior2_rand = numeric(num_edges)
  for(j in 1:num_edges){
    
    CRPS_vec_prior1_rand[j] = CRPS_analytic_normal(Y=edge_costs[j],means = theta_hat_rand[j], 
                                                   variances =  var_hats_rand[j])
    CRPS_vec_prior2_rand[j] = CRPS_analytic_normal(Y=edge_costs[j],means = theta_hat2_rand[j], 
                                                   variances =  var_hats2_rand[j])
    
  }
  
  mean_CRPS_prior1_rand[iter] = mean(CRPS_vec_prior1_rand)
  mean_CRPS_prior2_rand[iter] = mean(CRPS_vec_prior2_rand)
  
  
}


res_out = cbind(mean_CRPS_prior1, mean_CRPS_prior2, mean_coverage_prior1, mean_coverage_prior2,
                mean_CRPS_prior1_rand, mean_CRPS_prior2_rand, mean_coverage_prior1_rand, mean_coverage_prior2_rand)
colnames(res_out) = c("CRPS Prior 1", "CRPS Prior 2", "Coverage Prior 1", "Coverage Prior 2",
                      "CRPS Prior 1 RAND", "CRPS Prior 2 RAND", "Coverage Prior 1 RAND", "Coverage Prior 2 RAND")

# to save data, set working directory to SimResults subfolder relative to the root
# directory of the archive
setwd("./SimResults")
write.csv(res_out, namestr)

# Return to the root directory (to make it easier to run other files)
setwd("..")



