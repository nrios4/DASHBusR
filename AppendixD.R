# Section Appendix D

# This code reproduces the analysis in Appendix D.
# It will reproduce Table D1 and D2

rm(list = ls()) # Clear environemnt
# Set working directory to the location of "bus_data_X_Y_format.rds" and "source_code.R"

# load source code
source("source_code.R")

# load libraries
library(scoringRules) # for crps_t() function
library(tidyverse)
library(knitr)
# load the data
busdata = readRDS('bus_data_X_Y_format.rds')
busdata_ext = readRDS("bus_data_intersection_X_Y_format.rds")
intxn_col_index<-which(!colnames(busdata_ext)%in%colnames(busdata)) #find intersection term col index: which p in busdata_ext not in busdata
# load the ggmaps prior
load("Edge_prior.Rdata")
bus_graph = readRDS("bus_graph.rds")
X = as.matrix(busdata_ext[,1:(ncol(busdata_ext)-1)])
#X_intxn =get_inter_paths_matrix(bus_graph, L = 3)
Y = busdata_ext$Y
n = nrow(X)

# Prior 1:
#mu1 = rep(0,ncol(X)) #Archived
mu1 = c(rep(10,intxn_col_index[1]-1),rep(0,length(intxn_col_index)),0,0) #New
R = 10*diag(ncol(X))


# Prior 2:
mu2 = c(edge_prior$km,rep(0,length(intxn_col_index)),0,0)


###### Table 3: Fit model to all n = 1000 data points under both prior
model1 = fit_model(X,Y,mu1,R)
model2 = fit_model(X,Y,mu2,R)
tableD1 = cbind(model1$summary, model2$summary)


tableD1_df<-tableD1 %>% 
  round(2) %>% 
  as.data.frame
colnames(tableD1_df)<-c('pos_mean1','pos_var1','pos_low1','pos_up1','pos_mean2','pos_var2','pos_low2','pos_up2')
tableD1_df %>% 
  transmute(
    pos_mean1=pos_mean1,
    pos_var1=pos_var1,
    ci95_1 = sprintf("(%.2f, %.2f)", pos_low1, pos_up1),
    pos_mean2=pos_mean2,
    pos_var2=pos_var2,
    ci95_2 = sprintf("(%.2f, %.2f)", pos_low2, pos_up2))%>%
  rownames_to_column(var = "Route") %>%
  mutate(
    Route = Route %>%
      str_replace_all("\\.", " ") %>%
      str_replace_all("_", " to ")
  ) %>%
  filter(pos_mean1!=0,pos_mean2!=0) %>% 
  kable(format = 'latex',booktabs=T,linesep = "")

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

tableD2 = cbind(paste0(round(crps_means,3),
                      " (",round(crps_CI_lower,3),
                      ", ",round(crps_CI_upper,3),
                      ")"),
               round(coverage_means,3))
cases_str = c("alpha0 = 50, mu = mu1", 
              "alpha0 = 500, mu = mu1",
              "alpha0 = 5000, mu = mu1",
              "plug-in, mu = mu1",
              "alpha0 = 50, mu = mu2",
              "alpha0 = 500, mu = mu2",
              "alpha0 = 5000, mu = mu2",
              "plug-in, mu = mu2")
rownames(tableD2) = cases_str
tableD2 %>% 
  kable(format = 'latex',booktabs=T,linesep = "")
