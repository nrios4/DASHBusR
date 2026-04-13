# Table2.R
# Produces Table 2

# clear environment
rm(list = ls())

# set working directory to SimResults subfolder relative to the root
# directory of the archive
setwd("./SimResults")

# Create Table 2
scenarios = c(rep("A",8),rep("B",8),rep("C",8))
mvec = rep(c(rep(25,4),rep(27,4)),3)
nvec = rep(c(50,50,100,100),6)
design_vec = rep(c("SA","Rand"),12)
n_sim = 100

Table2 = data.frame(Scenario = scenarios, m = mvec, n = nvec, Design = design_vec,
                       CRPS_R1 = rep("",length(scenarios)),
                       CRPS_R2 = rep("",length(scenarios)),
                       Coverage_R1 = numeric(length(scenarios)),
                       Coverage_R2 = numeric(length(scenarios)))
metrics_mat = matrix(NA, nrow = 12, ncol = 8)
metrics_mat_rand = metrics_mat
for(i in 1:length(scenarios)){
  m = mvec[i]
  n = nvec[i]
  scenario = scenarios[i]
  design = design_vec[i]
  namestr = paste("Sim4.2v2scenario",scenario,"m",m,"n",n,".csv",sep = "")
  print(namestr)
  tmp_file = read.csv(namestr)[,-1]
  if(design == "SA"){
    metrics_vec = round(colMeans(tmp_file[,1:4]),3)
    crps_SEs = apply(tmp_file[,1:2], 2, function(x) sd(x)/sqrt(n_sim))
    crps_CI_lower = round(colMeans(tmp_file[,1:2]) - qt(0.025, df = n_sim-1, lower.tail = FALSE)*crps_SEs,3)
    crps_CI_upper = round(colMeans(tmp_file[,1:2]) + qt(0.025, df = n_sim-1, lower.tail = FALSE)*crps_SEs,3)
    Table2[i,5] = paste(metrics_vec[1], " (", crps_CI_lower[1], ",", crps_CI_upper[1], ")",
                                sep = "")
    Table2[i,6] = paste(metrics_vec[2], " (", crps_CI_lower[2], ",", crps_CI_upper[2], ")",
                        sep = "")
    Table2[i,c(7,8)] = metrics_vec[3:4]
  } else{
    metrics_vec = round(colMeans(tmp_file[,-(1:4)]),3)
    crps_SEs = apply(tmp_file[,5:6], 2, function(x) sd(x)/sqrt(n_sim))
    crps_CI_lower = round(colMeans(tmp_file[,5:6]) - qt(0.025, df = n_sim-1, lower.tail = FALSE)*crps_SEs,3)
    crps_CI_upper = round(colMeans(tmp_file[,5:6]) + qt(0.025, df = n_sim-1, lower.tail = FALSE)*crps_SEs,3)
    Table2[i,5] = paste(metrics_vec[1], " (", crps_CI_lower[1], ",", crps_CI_upper[1], ")",
                        sep = "")
    Table2[i,6] = paste(metrics_vec[2], " (", crps_CI_lower[2], ",", crps_CI_upper[2], ")",
                        sep = "")
    Table2[i,c(7,8)] = metrics_vec[3:4]
  }

  
}

Table2

# library(kableExtra)
# kable(Table2, format = "latex" )

# Return to the root directory (to make it easier to run other files)
setwd("..")
