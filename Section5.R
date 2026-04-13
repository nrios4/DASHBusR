# Section5.R

# This code reproduces the analysis in Section 5.
# It will reproduce Table 3, Table 4, and also Tables C1 and C2 in Appendix C. 

rm(list = ls()) # Clear environemnt
# Set working directory to the root directory of the archive, e.g.,
# setwd("C:/Users/riosn/OneDrive/Desktop/DASHBusRCodeR1")

# load source code
source("source_code.R")

# load libraries
library(scoringRules) # for crps_t() function

# load the data
busdata = readRDS("bus_data_X_Y_format.rds")

# load the ggmaps prior
load("Edge_prior.Rdata")

X = as.matrix(busdata[,1:(ncol(busdata)-1)])
Y = busdata$Y
n = nrow(X)

# Prior 1:
mu1 = c(rep(10,ncol(X)-2),0,0) 
R = 10*diag(ncol(X))

# Prior 2:
mu2 = c(edge_prior$km,0,0)


###### Table 3: Fit model to all n = 1000 data points under both prior
model1 = fit_model(X,Y,mu1,R)
model2 = fit_model(X,Y,mu2,R)
table3 = cbind(model1$summary, model2$summary)
table3

##### Repeated Train-Test splits for predictive analysis
n_split = 100 # number of splits
normal_crps_vec_1 = numeric(n_split)
normal_crps_vec_2 = numeric(n_split)
normal_coverage_vec_1 = numeric(n_split)
normal_coverage_vec_2 = numeric(n_split)
t1_crps_vec_1 = numeric(n_split)
t1_crps_vec_2 = numeric(n_split)
t1_coverage_vec_1 = numeric(n_split)
t1_coverage_vec_2 = numeric(n_split)
t2_crps_vec_1 = numeric(n_split)
t2_crps_vec_2 = numeric(n_split)
t2_coverage_vec_1 = numeric(n_split)
t2_coverage_vec_2 = numeric(n_split)
t3_crps_vec_1 = numeric(n_split)
t3_crps_vec_2 = numeric(n_split)
t3_coverage_vec_1 = numeric(n_split)
t3_coverage_vec_2 = numeric(n_split)


beta0 = 1
# alpha0 will be 50, 500, 5000

set.seed(123)

for(i in 1:n_split){
  
  # Training-Test Split
  # 90% Training, 10% Testing
  train_inds = sample(1:n, floor(0.9*n))
  X_train = X[train_inds,]
  X_test = X[-train_inds,]
  Y_train = Y[train_inds]
  Y_test = Y[-train_inds]
  
  # Fit model to training data with prior 1
  model1 = fit_model(X = X_train, Y = Y_train, mu = mu1, R = R)
  # Fit model to training data with prior 2
  model2 = fit_model(X = X_train, Y = Y_train, mu = mu2, R = R)
  
  # Use posterior to get predictions on test set (for plug-in estimator of sigma2)
  normal_pred_out1 = predict_model(model = model1, X = X_test)
  normal_pred_out2 = predict_model(model = model2, X = X_test)
  
  # Use posterior to get predictions on test set (assuming inverse-gamma prior)
  t1_pred_out1 = predict_model_t(model = model1, X = X_test, alpha0 = 50, beta0 = beta0)
  t1_pred_out2 = predict_model_t(model = model2, X = X_test, alpha0 = 50, beta0 = beta0)
  t2_pred_out1 = predict_model_t(model = model1, X = X_test, alpha0 = 500, beta0 = beta0)
  t2_pred_out2 = predict_model_t(model = model2, X = X_test, alpha0 = 500, beta0 = beta0)
  t3_pred_out1 = predict_model_t(model = model1, X = X_test, alpha0 = 5000, beta0 = beta0)
  t3_pred_out2 = predict_model_t(model = model2, X = X_test, alpha0 = 5000, beta0 = beta0)
  
  # find mean CRPS
  normal_crps_vec_1[i] = mean(CRPS_analytic_normal(Y = Y_test, means = normal_pred_out1[,1], variances = normal_pred_out1[,2]))
  normal_crps_vec_2[i] = mean(CRPS_analytic_normal(Y = Y_test, means = normal_pred_out2[,1], variances = normal_pred_out2[,2]))
  alpha_star1 = 50 + 0.5*nrow(X_train)
  alpha_star2 = 500 + 0.5*nrow(X_train)
  alpha_star3 = 5000 + 0.5*nrow(X_train)
  t1_crps_vec_1[i] = mean(crps_t(Y_test, df = 2*alpha_star1, location = t1_pred_out1[,1], scale = sqrt(t1_pred_out1[,2])))
  t1_crps_vec_2[i] = mean(crps_t(Y_test, df = 2*alpha_star1, location = t1_pred_out2[,1], scale = sqrt(t1_pred_out2[,2])))
  t2_crps_vec_1[i] = mean(crps_t(Y_test, df = 2*alpha_star2, location = t2_pred_out1[,1], scale = sqrt(t2_pred_out1[,2])))
  t2_crps_vec_2[i] = mean(crps_t(Y_test, df = 2*alpha_star2, location = t2_pred_out2[,1], scale = sqrt(t2_pred_out2[,2])))
  t3_crps_vec_1[i] = mean(crps_t(Y_test, df = 2*alpha_star3, location = t3_pred_out1[,1], scale = sqrt(t3_pred_out1[,2])))
  t3_crps_vec_2[i] = mean(crps_t(Y_test, df = 2*alpha_star3, location = t3_pred_out2[,1], scale = sqrt(t3_pred_out2[,2])))
  
  # find mean coverage
  normal_coverage_vec_1[i] = mean( (normal_pred_out1[,3] <= Y_test) & (Y_test <= normal_pred_out1[,4])  )
  normal_coverage_vec_2[i] = mean( (normal_pred_out2[,3] <= Y_test) & (Y_test <= normal_pred_out2[,4])  )
  t1_coverage_vec_1[i] = mean( (t1_pred_out1[,3] <= Y_test) & (Y_test <= t1_pred_out1[,4])   )
  t1_coverage_vec_2[i] = mean(  (t1_pred_out2[,3] <= Y_test) & (Y_test <= t1_pred_out2[,4])  )
  t2_coverage_vec_1[i] = mean( (t2_pred_out1[,3] <= Y_test) & (Y_test <= t2_pred_out1[,4])   )
  t2_coverage_vec_2[i] = mean(  (t2_pred_out2[,3] <= Y_test) & (Y_test <= t2_pred_out2[,4])  )
  t3_coverage_vec_1[i] = mean( (t3_pred_out1[,3] <= Y_test) & (Y_test <= t3_pred_out1[,4])   )
  t3_coverage_vec_2[i] = mean(  (t3_pred_out2[,3] <= Y_test) & (Y_test <= t3_pred_out2[,4])  )
  
}

resmat_crps = cbind(t1_crps_vec_1, t2_crps_vec_1, t3_crps_vec_1, normal_crps_vec_1,
                    t1_crps_vec_2, t2_crps_vec_2, t3_crps_vec_2, normal_crps_vec_2)
resmat_coverage = cbind(t1_coverage_vec_1,t2_coverage_vec_1, t3_coverage_vec_1, normal_coverage_vec_1,
                        t1_coverage_vec_2,t2_coverage_vec_2, t3_coverage_vec_2, normal_coverage_vec_2)



crps_means = colMeans(resmat_crps)
crps_SEs = apply(resmat_crps, 2, function(x) sd(x)/sqrt(n_split))
crps_CI_lower = crps_means - qt(0.025, df = n_split-1, lower.tail = FALSE)*crps_SEs
crps_CI_upper = crps_means + qt(0.025, df = n_split-1, lower.tail = FALSE)*crps_SEs
coverage_means = colMeans(resmat_coverage)

table4 = round(cbind( crps_means, crps_CI_lower, crps_CI_upper, coverage_means  ),3)
cases_str = c("alpha0 = 50, mu = mu1", 
           "alpha0 = 500, mu = mu1",
           "alpha0 = 5000, mu = mu1",
           "plug-in, mu = mu1",
           "alpha0 = 50, mu = mu2",
           "alpha0 = 500, mu = mu2",
           "alpha0 = 5000, mu = mu2",
           "plug-in, mu = mu2")
rownames(table4) = cases_str
table4

#################################################
# Finding an efficient follow-up experiment (L = 5,6)
#################################################

n = 50 # target sample size
L = 5 # path length constraint (L = 5,6)


R_updated = t(X)%*%X + 0.1*diag(ncol(X))

# Load the graph
bus_graph = readRDS("bus_graph.rds")
bus_graph = simplify(bus_graph)
plot(bus_graph)

p = 2
num_edges = length(E(bus_graph))
all_paths = get_all_paths_with_L_edges(bus_graph, L)


num_paths = length(all_paths)
X_paths = matrix(0, nrow = num_paths, ncol = num_edges)
edge_names <- apply(as_edgelist(bus_graph), 1, function(e) paste(sort(e), collapse = "_"))
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

set.seed(123)
# Generate an initial design with n runs
X_init = X_all[sample(1:nrow(X_all),n),]

# find 30 different designs, look at the relative phi-efficiency compared to fedorov
n_design = 30 # number of designs to search for with SA
SA_designs = list()
SA_elapsed_times = numeric(n_design)

for(i in 1:n_design){
  
  print(paste("Finding design",i,sep=" "))
  # Use SA to find a design of size n
  t1 <- Sys.time()
  SA_designs[[i]] <- SA_design(X_init = X_init, p=p, L=L, R=R_updated, X_paths=X_paths, n_iter = 25000)  
  t2 <- Sys.time()
  SA_elapsed_times[i] = as.numeric(t2 - t1, units = "secs")
  
  
}

t3 <- Sys.time()
fedorov_out = fedorov_bayes(X_init = X_init, X_all = X_all, R = R_updated, n_repeat = 2)
t4 <- Sys.time()
fedorov_elapsed = as.numeric(t4 - t3, units = "secs")

SA_Deffs = unlist(lapply(SA_designs, function(x) x$phi))

namestr = paste("SA_rel_Deffs_n",n,"_L",L,".csv", sep = "")
relDeffs = SA_Deffs/fedorov_out$phi
# write.csv(relDeffs, namestr)

# CPU times
mean(SA_elapsed_times)
fedorov_elapsed

# Find and store the best design in a readable format
X_best = SA_designs[[which.max(SA_Deffs)]]$model.matrix
namestr2 = paste("bestSAdesign_n",n,"L",L,".csv", sep = "")

X_best_edges = X_best[,1:num_edges]
colnames(X_best_edges) = edge_names
best_paths = NULL
# best_paths = c()
for(i in 1:n){
  tmp1 = edge_names[X_best_edges[i,] == 1] 
  tmp2 = get_path_from_strings(tmp1)
  best_paths = rbind(best_paths, tmp2$vertex_order)
}
colnames(best_paths) = paste("Bus Stop", 1:(L+1))
best_design_df = data.frame(best_paths, bus_size = X_best[,num_edges+1],temperature = X_best[,num_edges+2])

# write.csv(best_design_df, namestr2)

library(xtable)
xtable(best_design_df)
